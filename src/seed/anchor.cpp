// ==================================================================
// chain.cpp - minimap2 风格的锚点收集（seeding）与预处理
// ==================================================================
//
// 功能说明：
// 本文件实现了序列比对中的锚点（anchor）收集功能，参考 minimap2 的设计。
// 锚点是 ref 和 query 序列之间共享的种子匹配位置，用于后续的链化（chaining）
// 和精确比对。
//
// 参考来源：
// - minimap2/seed.c: mm_collect_matches() - 收集种子匹配
// - minimap2/hit.c: 锚点处理和过滤
// - minimap2/lchain.c: 链化算法（DP 动态规划）
//
// 核心数据流程（minimap2 风格，参考 README Algorithm overview）：
// 1) 对参考序列 minimizers 建索引（hash -> occurrences）
// 2) 对每个 query minimizer 查索引；若该参考 minimizer 不在 top -f 最频繁集合内，
//    则收集其在参考序列的所有出现位置作为 seeds/anchors
// 3) seeds/anchors 随后按参考坐标排序并做 DP chaining
//
// 重要：
// - minimap2 的频次过滤（-f/-U/--q-occ-frac/-e）发生在“展开 occurrences 之前”，
//   否则重复区域会导致 anchors 数量爆炸。
//
// 设计说明：
// - 使用 C++ 标准库（std::sort, std::unordered_map）替代 minimap2 的 ksort/khash
// - 保持与 minimap2 相同的语义：hash 匹配表示潜在的同源区域
// - 模板化设计：支持任意继承 SeedHitBase 的种子类型
//
// 性能优化：
// - 预先对 ref_hits 按 hash 排序
// - 避免重复的 hash 计算（使用缓存的 hash 值）
// - 预分配输出容器，减少内存重新分配
// ==================================================================

#include "seed.h"
#include <algorithm>
#include <unordered_map>
#include <cstdint>
#include <vector>
#include <cmath>
#include <limits>

namespace anchor {

// ==================================================================

// ------------------------------------------------------------------
// 函数：compute_occ_cutoff_top_frac
// 功能：根据 distinct minimizer 的 occurrence 列表，计算“最高频 top f_top_frac”对应的频次阈值。
//
// 输入：
// - occs: 每个 distinct minimizer 在参考中出现的次数（例如：{1,1,2,10,3,...}）
// - f_top_frac: top fraction（例如 2e-4 表示忽略最频繁的 0.02% minimizers）
//
// 输出：
// - 返回 occ_cutoff：当 a minimizer 的 occurrence >= occ_cutoff 时，可认为其属于 top 高频区域
//   （上层可据此决定是否过滤/稀疏采样）。
//
// 关键语义（对齐 minimap2 的直觉）：
// - f_top_frac==0：不做基于 top fraction 的过滤 => 返回 +inf
// - occs 很少或算出来 n_skip==0：也相当于“不过滤” => 返回 +inf
// - f_top_frac>=1：极端情况，几乎全部都要过滤 => 返回 1
//
// 实现细节：
// - 我们需要“第 n_skip 大的 occurrence 值”，用 nth_element 做 O(N) 期望时间选择。
// - 使用 std::greater<std::size_t>() 让 nth_element 得到降序意义下的第 n_skip 大。
//
// 复杂度：
// - 时间：O(N) 期望（nth_element），最坏 O(N log N)
// - 空间：O(N)（拷贝一份 tmp，避免修改输入）
// ------------------------------------------------------------------
std::size_t compute_occ_cutoff_top_frac(const std::vector<std::size_t>& occs,
                                        double f_top_frac)
{
    if (occs.empty()) return std::numeric_limits<std::size_t>::max();
    if (f_top_frac <= 0.0) return std::numeric_limits<std::size_t>::max();
    if (f_top_frac >= 1.0) return 1; // 极端情况：几乎全丢

    // top f fraction of DISTINCT minimizers
    const std::size_t n = occs.size();

    // n_skip 表示需要“忽略掉的 distinct minimizer 数量”
    // 例：n=10000, f=2e-4 => n_skip=f*n=2 => 忽略出现次数排名前 2 的 minimizer
    const std::size_t n_skip = static_cast<std::size_t>(std::floor(f_top_frac * static_cast<double>(n)));
    if (n_skip == 0) return std::numeric_limits<std::size_t>::max();

    // nth_element 会重排 tmp，使得 tmp[n_skip-1] 是“第 n_skip 大”的元素
    std::vector<std::size_t> tmp = occs;
    std::nth_element(tmp.begin(), tmp.begin() + static_cast<std::ptrdiff_t>(n_skip - 1), tmp.end(),
                     std::greater<std::size_t>());
    return tmp[n_skip - 1];
}

// ------------------------------------------------------------------
// 函数：compute_ref_occ_threshold
// 功能：计算最终“参考端 occurrence 阈值”。
//
// minimap2 里常见逻辑可理解为：
// - 先用 -f 估计一个阈值 f_cutoff（忽略 top fraction 的高频 minimizer）
// - 再用 -U 设置上下界，把阈值限制到 [u_floor, u_ceil] 范围内
// - 最终阈值 = max(u_floor, min(u_ceil, f_cutoff))
//
// 这样做的好处：
// - 当参考序列很大且重复很多，-f 能自适应地忽略最频繁的 minimizer
// - 同时 -U 能确保阈值不会过小（误杀正常 seeds）也不会过大（让 repeats 爆炸）
//
// 复杂度：
// - 主要来自 compute_occ_cutoff_top_frac 的 nth_element
// ------------------------------------------------------------------
std::size_t compute_ref_occ_threshold(const std::vector<std::size_t>& occs,
                                      const SeedFilterParams& p)
{
    const std::size_t f_cutoff = compute_occ_cutoff_top_frac(occs, p.f_top_frac);
    const std::size_t capped = std::min(p.u_ceil, f_cutoff);
    return std::max(p.u_floor, capped);
}

// ==================================================================
// 辅助函数：按对角线排序锚点（minimap2 chaining 预处理常用）
// ==================================================================
//
// 为什么要按对角线排序？
// - chaining 本质是把“在相同对角线附近”的 anchors 串起来形成一条链。
// - 对角线（diagonal）可理解为 ref_pos - qry_pos（正向）
//   对于同源区域，锚点通常会聚集在相近的 diagonal 上。
//
// 注意反向链（is_rev==true）的对角线定义不同：
// - minimap2 内部会把反向链 query 坐标映射到反向互补坐标系。
// - 在不知道 query 总长度 qlen 的情况下，我们无法直接算 (ref - q_rc)。
// - 但排序只需要“单调等价”的 key：
//   q_rc = qlen - (q + span)  => ref - q_rc = ref + q + span - qlen
//   qlen 是常数，可以忽略，因此用 (ref + q + span) 作为反向链的 diagonal key。
//
// 排序键：
// 1) rid_ref（不同参考序列分开）
// 2) is_rev（正向/反向分开，避免混链）
// 3) diagonal key（见上）
// 4) pos_ref / pos_qry（稳定化）
//
// 复杂度：O(A log A)，A=anchors.size()
// ------------------------------------------------------------------
void sortAnchorsByDiagonal(Anchors& anchors)
{
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.rid_ref != b.rid_ref) return a.rid_ref < b.rid_ref;
            if (a.is_rev != b.is_rev)   return a.is_rev < b.is_rev;

            // 正向：diag = r - q
            // 反向：等价于 r - q_rc，其中 q_rc = qlen - (q + span)
            //       qlen 为常数，排序可用 (r + q + span) 代替
            const int64_t diag_a = a.is_rev
                ? (static_cast<int64_t>(a.pos_ref) + static_cast<int64_t>(a.pos_qry) + static_cast<int64_t>(a.span))
                : (static_cast<int64_t>(a.pos_ref) - static_cast<int64_t>(a.pos_qry));
            const int64_t diag_b = b.is_rev
                ? (static_cast<int64_t>(b.pos_ref) + static_cast<int64_t>(b.pos_qry) + static_cast<int64_t>(b.span))
                : (static_cast<int64_t>(b.pos_ref) - static_cast<int64_t>(b.pos_qry));
            if (diag_a != diag_b) return diag_a < diag_b;

            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        });
}

// ==================================================================
// 辅助函数：按位置排序锚点
// ==================================================================
void sortAnchorsByPosition(Anchors& anchors)
{
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.rid_ref != b.rid_ref) return a.rid_ref < b.rid_ref;
            if (a.is_rev != b.is_rev)   return a.is_rev < b.is_rev;
            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        });
}

// ==================================================================
// 辅助函数：过滤高频锚点（后过滤；不等价于 minimap2 的 -f/-U）
// ==================================================================
void filterHighFrequencyAnchors(Anchors& anchors, std::size_t max_occ)
{
    if (anchors.empty() || max_occ == 0) return;

    std::unordered_map<hash_t, std::size_t> hash_count;
    for (const auto& anchor : anchors) {
        hash_count[anchor.hash]++;
    }

    // 注意：这属于"后过滤"，语义上不同于 minimap2 的 -f/-U（后者在展开 occurrences 前过滤）。
    auto new_end = std::remove_if(anchors.begin(), anchors.end(),
        [&hash_count, max_occ](const Anchor& anchor) {
            return hash_count[anchor.hash] > max_occ;
        });

    anchors.erase(new_end, anchors.end());
}

// ==================================================================
// 链化（Chaining）相关函数实现
// ==================================================================
//
// 参考 minimap2/lchain.c 的 mg_lchain_dp 函数实现。
//
// minimap2 的链化算法核心思想：
// 1. 锚点按 (target_pos, query_pos) 排序
// 2. 使用 DP 计算到达每个锚点的最优累积得分 f[i]
// 3. 对于每个锚点 i，考虑所有可能的前驱锚点 j (j < i)
// 4. 两锚点可链接的条件：
//    - 在同一条参考序列（rid_ref 相同）
//    - 链方向相同（is_rev 相同）
//    - 参考距离和查询距离都在允许范围内
//    - 对角线偏移在带宽内
// 5. 链化得分 = 基础得分 - 惩罚项
// 6. 使用 max_skip 和 max_iter 进行剪枝优化
// 7. 回溯找出所有满足条件的链
// ==================================================================

// ------------------------------------------------------------------
// 辅助函数：计算 log2（用于惩罚项，参考 minimap2 的 mg_log2）
// ------------------------------------------------------------------
static inline float mg_log2(float x) {
    // 使用标准库的 log2，对于小数值返回 0
    return x >= 1.0f ? std::log2(x) : 0.0f;
}

// ------------------------------------------------------------------
// chainScoreSimple - 计算两个锚点之间的链化得分
// ------------------------------------------------------------------
// 实现细节（对应 minimap2/lchain.c 的 comput_sc_simple）：
//
// 输入约定：
// - ai 是"后面"的锚点（位置更大）
// - aj 是"前面"的锚点（位置更小）
//
// 可链接性检查：
// - 必须在同一参考序列：ai.rid_ref == aj.rid_ref
// - 链方向必须相同：ai.is_rev == aj.is_rev
// - 查询方向 gap (dq) 必须 > 0 且 <= max_dist_x
// - 参考方向 gap (dr) 必须 > 0（除非是反向链的特殊情况）
// - 对角线偏移 |dr - dq| 必须 <= bw
//
// 得分计算：
// - 基础得分 = min(span_j, min(dq, dr))
// - 对角线惩罚 = gap_penalty * |dr - dq| + skip_penalty * min(dr, dq)
// - 对数惩罚 = 0.5 * log2(|dr - dq| + 1)
// - 最终得分 = 基础得分 - 对角线惩罚 - 对数惩罚
// ------------------------------------------------------------------
std::int32_t chainScoreSimple(const Anchor& ai, const Anchor& aj, const ChainParams& params)
{
    // 检查是否在同一参考序列、同一链方向
    if (ai.rid_ref != aj.rid_ref) return INT32_MIN;
    if (ai.is_rev != aj.is_rev)   return INT32_MIN;

    // 计算 query 方向的 gap
    const std::int32_t dq = static_cast<std::int32_t>(ai.pos_qry) - static_cast<std::int32_t>(aj.pos_qry);

    // dq 必须 > 0（ai 必须在 aj 之后）且在允许范围内
    if (dq <= 0 || dq > params.max_dist_x) return INT32_MIN;

    // 计算 reference 方向的 gap
    const std::int32_t dr = static_cast<std::int32_t>(ai.pos_ref) - static_cast<std::int32_t>(aj.pos_ref);

    // 对于正向链，dr 也必须 > 0
    // 对于反向链，由于坐标系不同，规则可能不同
    // 这里简化处理：要求 dr 在合理范围内
    if (dr <= 0 || dr > params.max_dist_y) return INT32_MIN;

    // 计算对角线偏移（diagonal deviation）
    const std::int32_t dd = (dr > dq) ? (dr - dq) : (dq - dr);

    // 对角线偏移必须在带宽内
    if (dd > params.bw) return INT32_MIN;

    // 计算 gap 的较小值
    const std::int32_t dg = (dr < dq) ? dr : dq;

    // 基础得分：min(span_j, dg)
    // span 代表前一个锚点覆盖的长度
    const std::int32_t q_span = static_cast<std::int32_t>(aj.span);
    std::int32_t sc = (q_span < dg) ? q_span : dg;

    // 惩罚项（只有当存在 gap 偏差或 gap 大于 span 时才惩罚）
    if (dd > 0 || dg > q_span) {
        // 线性惩罚
        const float lin_pen = params.gap_penalty * static_cast<float>(dd)
                            + params.skip_penalty * static_cast<float>(dg);
        // 对数惩罚（对 gap 偏差取对数）
        const float log_pen = (dd >= 1) ? mg_log2(static_cast<float>(dd + 1)) : 0.0f;

        // 总惩罚 = 线性惩罚 + 0.5 * 对数惩罚
        sc -= static_cast<std::int32_t>(lin_pen + 0.5f * log_pen);
    }

    return sc;
}

// ------------------------------------------------------------------
// chainAnchors - 使用 DP 对锚点进行链化并返回最佳链
// ------------------------------------------------------------------
// 实现细节（对应 minimap2/lchain.c 的 mg_lchain_dp）：
//
// 算法流程：
// 1. 对锚点按位置排序：(rid_ref, is_rev, pos_ref, pos_qry)
// 2. 初始化 DP 数组：
//    - f[i] : 以锚点 i 结尾的最大链得分
//    - p[i] : 锚点 i 的最优前驱（-1 表示无前驱）
//    - v[i] : 到达锚点 i 之前的峰值得分（用于回溯）
// 3. DP 转移：
//    对于每个锚点 i，遍历其前面的锚点 j，计算链化得分
//    f[i] = max{f[j] + score(j, i)} for all valid j
// 4. 回溯：
//    找出得分最高的锚点，回溯提取链中的所有锚点
// 5. 返回最佳链的锚点列表（如果满足 min_cnt 和 min_score）
//
// 优化策略：
// - max_iter: 限制每个锚点考虑的前驱数量
// - max_skip: 连续跳过的锚点数达到阈值后提前终止
//
// 复杂度：
// - 最坏 O(N^2)，但通过 max_iter/max_skip 剪枝后通常接近 O(N * max_iter)
//
// 修改说明：
// - 不再返回 Chains 结构体，直接返回最佳链的 Anchors
// - 简化了逻辑，避免构建和管理多个链
// ------------------------------------------------------------------
Anchors chainAnchors(Anchors& anchors, const ChainParams& params)
{
    const std::int64_t n = static_cast<std::int64_t>(anchors.size());

    if (n == 0) return Anchors{};

    // ------------------------------------------------------------------
    // 步骤 1：按位置排序锚点
    // 排序键：(rid_ref, is_rev, pos_ref, pos_qry)
    // 这样同一参考区域的锚点会聚在一起，便于 DP 处理
    // ------------------------------------------------------------------
    std::sort(anchors.begin(), anchors.end(),
        [](const Anchor& a, const Anchor& b) {
            if (a.rid_ref != b.rid_ref) return a.rid_ref < b.rid_ref;
            if (a.is_rev != b.is_rev)   return a.is_rev < b.is_rev;
            if (a.pos_ref != b.pos_ref) return a.pos_ref < b.pos_ref;
            return a.pos_qry < b.pos_qry;
        });

    // ------------------------------------------------------------------
    // 步骤 2：分配 DP 数组
    // ------------------------------------------------------------------
    std::vector<std::int32_t> f(n);      // f[i] = 以锚点 i 结尾的最大链得分
    std::vector<std::int64_t> p(n);      // p[i] = 锚点 i 的最优前驱索引（-1 表示无前驱）
    std::vector<std::int32_t> t(n, 0);   // t[i] = 标记数组（用于 max_skip 优化）

    // ------------------------------------------------------------------
    // 步骤 3：DP 填表
    // ------------------------------------------------------------------
    // st : 滑动窗口的起始位置（排除距离过远的锚点）
    std::int64_t st = 0;
    std::int32_t max_f_global = 0;       // 全局最大得分
    std::int64_t best_end_idx = -1;      // 得分最高的链的结束锚点索引

    for (std::int64_t i = 0; i < n; ++i) {
        const Anchor& ai = anchors[static_cast<std::size_t>(i)];

        // 初始化：f[i] = span（单锚点自身的得分）
        std::int32_t max_f = static_cast<std::int32_t>(ai.span);
        std::int64_t max_j = -1;
        std::int32_t n_skip = 0;

        // 移动窗口起点：排除距离过远或不同参考序列的锚点
        while (st < i) {
            const Anchor& ast = anchors[static_cast<std::size_t>(st)];
            // 不同参考序列或参考距离超限，则移动窗口
            if (ast.rid_ref != ai.rid_ref ||
                ast.is_rev != ai.is_rev ||
                static_cast<std::int32_t>(ai.pos_ref - ast.pos_ref) > params.max_dist_x) {
                ++st;
            } else {
                break;
            }
        }

        // 限制迭代次数
        std::int64_t iter_start = (i - st > params.max_iter) ? (i - params.max_iter) : st;

        // 从后往前遍历候选前驱
        for (std::int64_t j = i - 1; j >= iter_start; --j) {
            const Anchor& aj = anchors[static_cast<std::size_t>(j)];

            // 计算链化得分
            const std::int32_t sc = chainScoreSimple(ai, aj, params);
            if (sc == INT32_MIN) continue;

            // 累加前驱的得分
            const std::int32_t total_sc = f[static_cast<std::size_t>(j)] + sc;

            if (total_sc > max_f) {
                max_f = total_sc;
                max_j = j;
                // 找到更好的前驱，重置 skip 计数
                if (n_skip > 0) --n_skip;
            } else if (t[static_cast<std::size_t>(j)] == static_cast<std::int32_t>(i)) {
                // 这个锚点已经在当前 i 的遍历中被访问过
                // 如果连续跳过太多次，提前终止
                if (++n_skip > params.max_skip) break;
            }

            // 标记 j 的前驱已被访问（用于 max_skip 优化）
            if (p[static_cast<std::size_t>(j)] >= 0) {
                t[static_cast<std::size_t>(p[static_cast<std::size_t>(j)])] = static_cast<std::int32_t>(i);
            }
        }

        // 记录 DP 结果
        f[static_cast<std::size_t>(i)] = max_f;
        p[static_cast<std::size_t>(i)] = max_j;

        // 更新全局最大得分和对应的结束索引
        if (max_f > max_f_global) {
            max_f_global = max_f;
            best_end_idx = i;
        }
    }

    // ------------------------------------------------------------------
    // 步骤 4：回溯提取最佳链
    // ------------------------------------------------------------------
    // 如果全局最大得分低于阈值，返回空
    if (max_f_global < params.min_score || best_end_idx < 0) {
        return Anchors{};
    }

    // 回溯收集链中的锚点索引
    std::vector<std::int64_t> chain_indices;
    std::int64_t cur = best_end_idx;
    while (cur >= 0) {
        chain_indices.push_back(cur);
        cur = p[static_cast<std::size_t>(cur)];
    }

    // 检查链长度是否满足要求
    if (static_cast<std::int32_t>(chain_indices.size()) < params.min_cnt) {
        return Anchors{};
    }

    // 反转使得锚点按位置顺序排列（从前到后）
    std::reverse(chain_indices.begin(), chain_indices.end());

    // 提取锚点构建结果
    Anchors result;
    result.reserve(chain_indices.size());
    for (std::int64_t idx : chain_indices) {
        result.push_back(anchors[static_cast<std::size_t>(idx)]);
    }

    return result;
}


} // namespace anchor

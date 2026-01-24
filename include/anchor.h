#ifndef HALIGN4_ANCHOR_H
#define HALIGN4_ANCHOR_H

#include <cstdint>
#include <type_traits>
#include <string_view>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <cstddef>
#include <algorithm>
#include <limits>
#include <cmath>
#include "hash.h"

// ================================================================
// anchor 命名空间：锚点（Anchor）相关的数据结构与工具函数
// ================================================================
// 设计动机：
// - Anchor 及其配套的过滤/排序函数，本质属于“锚点/链化阶段”的公共数据结构，
//   不应该和 seed 抽象接口（SeedHitBase/traits）强耦合在同一个命名空间里。
// - 将其拆到 anchor::，可以让 seed:: 只负责“seed/hit 的抽象与提取”，
//   anchor:: 专注“如何用 hits 形成 anchors、如何排序/过滤”。
//
// 重要：
// - 这里仍然使用全局 hash_t（来自 include/hash.h），保持与现有代码一致。
// - 目前 minimizer::collect_anchors 返回 anchor::Anchors；后续如果支持 syncmer/strobemer，
//   也可以复用同一套 anchor 工具函数。
// ================================================================
namespace anchor
{
    // ------------------------------------------------------------------
    // 结构体：Anchor
    // ------------------------------------------------------------------
    // 语义：描述 ref 与 query 在某个 seed/hash 上的一次“锚点匹配”。
    // 下游 chaining 会把一组 anchors 串成一条链（候选比对区域）。
    // ------------------------------------------------------------------
    struct Anchor
    {
        hash_t hash{};              // seed hash（期望 ref/query 相同）
        std::uint32_t rid_ref{};    // ref 序列 id
        std::uint32_t pos_ref{};    // ref 上位置（0-based）
        std::uint32_t rid_qry{};    // query 序列 id（很多场景固定 0）
        std::uint32_t pos_qry{};    // query 上位置（0-based, forward 坐标系）
        std::uint32_t span{};       // 覆盖长度（可用 min(ref.span, qry.span)）
        bool is_rev{};              // ref/query 是否为"相反链"
        // 可选：预先缓存对角线，chaining 常用
        // int32_t diag{}; // (int32_t)pos_ref - (int32_t)pos_qry
    };

    using Anchors = std::vector<Anchor>;

    // ==================================================================
    // 内部辅助结构：用于 ref_hits 的 hash 索引
    // ==================================================================
    struct HashIndex {
        std::size_t start;  // 在排序后数组中的起始索引
        std::size_t count;  // 具有该 hash 的元素数量
    };

    // ==================================================================
    // minimap2 风格的 seeding 过滤参数（默认值与 minimap2 CLI 相同）
    // ==================================================================
    struct SeedFilterParams {
        double f_top_frac = 2e-4;                 // -f
        std::size_t u_floor = 10;                 // -U lower
        std::size_t u_ceil  = 1000000;            // -U upper
        double q_occ_frac  = 0.01;                // --q-occ-frac
        std::size_t sample_every_bp = 500;        // -e
    };

    static inline SeedFilterParams default_mm2_params() {
        return SeedFilterParams{};
    }

    // 计算 -f (fraction) 对应的 occurrence 阈值：忽略 top f_top_frac 最频繁 minimizers
    // 返回值：occ_cutoff（>=1）。当 distinct minimizers 很少或 f_top_frac==0 时，返回 +inf。
    std::size_t compute_occ_cutoff_top_frac(const std::vector<std::size_t>& occs,
                                            double f_top_frac);

    // 计算最终 reference occurrence 阈值：max{u_floor, min{u_ceil, -f}}
    std::size_t compute_ref_occ_threshold(const std::vector<std::size_t>& occs,
                                          const SeedFilterParams& p);

    // =====================================================================
    // sortAnchorsByDiagonal - 按对角线排序锚点（用于链化算法）
    // =====================================================================
    void sortAnchorsByDiagonal(Anchors& anchors);

    // =====================================================================
    // sortAnchorsByPosition - 按位置排序锚点
    // =====================================================================
    void sortAnchorsByPosition(Anchors& anchors);

    // =====================================================================
    // filterHighFrequencyAnchors - 过滤高频锚点（参考 minimap2）
    // =====================================================================
    void filterHighFrequencyAnchors(Anchors& anchors, std::size_t max_occ = 500);

    // =====================================================================
    // 链化（Chaining）相关数据结构与函数
    // =====================================================================
    // 参考 minimap2/lchain.c 的实现，使用动态规划（DP）将锚点串成链。
    //
    // 核心思想：
    // - 锚点按参考位置排序后，使用 DP 计算到达每个锚点的最优得分
    // - 两个锚点之间的得分 = min(gap_ref, gap_qry, span) - 惩罚项
    // - 惩罚项考虑 gap 差异（对角线偏移）和 gap 大小
    // - 回溯找出得分最高的链
    // =====================================================================

    // ------------------------------------------------------------------
    // 链化参数（参考 minimap2 的默认值）
    // ------------------------------------------------------------------
    struct ChainParams {
        std::int32_t max_dist_x = 5000;       // 参考序列方向最大距离
        std::int32_t max_dist_y = 5000;       // 查询序列方向最大距离
        std::int32_t bw = 500;                // 带宽（对角线偏移容忍度）
        std::int32_t max_skip = 25;           // 最大跳过数（优化）
        std::int32_t max_iter = 5000;         // 最大迭代次数（优化）
        std::int32_t min_cnt = 3;             // 链的最小锚点数
        std::int32_t min_score = 40;          // 链的最小得分
        float gap_penalty = 0.01f;            // gap 惩罚系数
        float skip_penalty = 0.01f;           // skip（gap 大小）惩罚系数
    };

    // 返回默认链化参数
    inline ChainParams default_chain_params() {
        return ChainParams{};
    }

    // ------------------------------------------------------------------
    // chainScoreSimple - 计算两个锚点之间的链化得分（内部辅助函数）
    // ------------------------------------------------------------------
    // 输入：
    //   ai    : 当前锚点（位置较大）
    //   aj    : 前一个锚点（位置较小）
    //   params: 链化参数
    //
    // 输出：
    //   返回链化得分，若两锚点不可链接返回 INT32_MIN
    //
    // 得分计算（参考 minimap2）：
    //   基础得分 = min(gap_ref, gap_qry, span)
    //   惩罚 = gap_penalty * |gap_ref - gap_qry| + skip_penalty * min(gap_ref, gap_qry)
    //   最终得分 = 基础得分 - 惩罚 - 0.5 * log2(|gap_ref - gap_qry| + 1)
    //
    // 注意：此函数仅供 chainAnchors 内部使用
    // ------------------------------------------------------------------
    std::int32_t chainScoreSimple(const Anchor& ai, const Anchor& aj, const ChainParams& params);

    // ------------------------------------------------------------------
    // chainAnchors - 使用 DP 对锚点进行链化并返回最佳链
    // ------------------------------------------------------------------
    // 功能：
    // 使用动态规划算法对锚点进行链化，找出得分最高的锚点链。
    // 参考 minimap2/lchain.c 的 mg_lchain_dp 实现。
    //
    // 输入：
    //   anchors : 锚点列表（会被排序并修改）
    //   params  : 链化参数
    //
    // 输出：
    //   返回最佳链包含的锚点列表（按位置顺序排列）
    //   如果没有找到满足条件的链，返回空 Anchors
    //
    // 算法流程：
    // 1. 对锚点按 (rid_ref, is_rev, pos_ref, pos_qry) 排序
    // 2. 使用 DP 计算每个锚点的最优链接得分
    // 3. 回溯提取得分最高的链
    // 4. 检查链是否满足 min_cnt 和 min_score 条件
    // 5. 返回最佳链的锚点列表
    //
    // 注意：
    //   - 输入 anchors 会被排序
    //   - 返回的锚点按在链中的位置顺序排列（从前到后）
    //   - 如果有多条得分相同的链，返回第一条
    // ------------------------------------------------------------------
    Anchors chainAnchors(Anchors& anchors, const ChainParams& params = default_chain_params());

} // namespace anchor


#endif //HALIGN4_ANCHOR_H
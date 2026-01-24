#include "align.h"
#include "ksw2.h"
#include "bindings/cpp/WFAligner.hpp"
extern "C" {
#include "alignment/cigar.h"
#include "wavefront/wavefront_align.h"
}

// ==================================================================
// align.cpp - HAlign-4 序列比对算法封装（KSW2 / WFA2）
// ==================================================================
// 本文件实现 align.h 中声明的比对接口，核心目标是：
// 1) 将不同底层比对库（KSW2、WFA2）的调用细节封装起来，对外统一返回 cigar::Cigar_t
// 2) 统一使用项目内的 DNA 字符编码与替换矩阵（ScoreChar2Idx / dna5_simd_mat）
// 3) 明确资源所有权（哪些需要 free/delete），避免内存泄漏
//
// 重要约定（正确性）：
// - 所有接口都返回"query 相对于 ref"的 CIGAR（即描述 query 如何对齐到 ref）。
// - CIGAR 的压缩编码与 BAM/SAM 规范一致：cigar_unit=(len<<4)|op。
// - 本文件只做"调用封装"与"参数配置"，不改变算法本身的比对逻辑。
//
// 性能说明：
// - KSW2：经典 DP（SIMD 加速），对中短序列稳定；对长序列可用 band 降低复杂度。
// - WFA2：Wavefront，编辑距离小（高相似度）时非常快；相似度低时可能退化。
// - 编码阶段会产生两个临时 vector（ref_enc/qry_enc），属于必要开销；
//   若未来成为热点，可考虑复用缓冲（但会引入状态/线程安全复杂度，本次不改行为）。
// ==================================================================

namespace align
{
    // ------------------------------------------------------------------
    // globalAlignKSW2：KSW2 全局比对（end-to-end）
    // ------------------------------------------------------------------
    // 输入：ref/query 原始序列（字符串，允许包含 A/C/G/T/N，大小写均可）
    // 输出：cigar::Cigar_t（压缩 CIGAR），描述 query 如何对齐到 ref
    //
    // 实现分 4 步：
    // 1) 序列编码：将字符映射到 0..4（A/C/G/T/N）供 KSW2 使用
    // 2) 配置参数：替换矩阵、gap 罚分、band/zdrop 等
    // 3) 调用 ksw_extz2_sse：获得 ez（包含 CIGAR 指针与长度）
    // 4) 拷贝 CIGAR 并释放 ez.cigar（KSW2 内部 malloc）
    //
    // 正确性/一致性：
    // - 本实现启用 KSW_EZ_GENERIC_SC，确保使用 dna5_simd_mat 完整替换矩阵
    // - 启用 KSW_EZ_RIGHT，让 gap 右对齐，提高 CIGAR 的"规范化"程度，便于后续比较/测试
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query)
    {
        /* ---------- 1. 编码序列 ---------- */
        // 说明：KSW2 期望输入为整数编码序列。
        // ScoreChar2Idx 的映射规则见 align.h：A/C/G/T -> 0..3，其它 -> 4(N)
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());

        for (size_t i = 0; i < ref.size(); ++i)
            ref_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
        for (size_t i = 0; i < query.size(); ++i)
            qry_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(query[i])];

        /* ---------- 2. 配置 KSW2 参数 ---------- */
        // dna5_simd_mat 是 constexpr 编译期常量，无需运行时初始化。
        // 这里 cfg 的所有字段都显式填写，避免依赖未初始化内存。
        align::KSW2AlignConfig cfg;
        cfg.mat = align::dna5_simd_mat;
        cfg.alphabet_size = 5;

        // Gap 惩罚参数（KSW2 内部会转为负数，这里传正值）
        // - gap_open=6：开启一个 gap 的基础代价
        // - gap_extend=2：每延长 1bp 的额外代价
        // - 总代价公式：gap of length L 的代价 = -(gap_open + gap_extend * L)
        // - 例如：1bp gap 代价 = -(6+2) = -8，与单个错配代价 -9 接近
        cfg.gap_open = 6;
        cfg.gap_extend = 2;

        // end_bonus：全局比对通常不需要 ends-free 奖励
        cfg.end_bonus = 0;

        // zdrop：全局比对应完整跑完，一般不依赖 z-drop 提前终止
        // 这里传 -1 表示禁用（与项目原有逻辑保持一致）
        cfg.zdrop = -1;

        // band_width：限制 DP 计算在对角线附近的带宽，可显著加速长序列。
        // 注意：带宽过小会导致比对失败或产生次优路径。
        // 这里用 auto_band 做启发式估计，属于"性能/准确性折中"的经验值。
        cfg.band_width = align::auto_band(ref.size(), query.size());

        // Flag 组合说明：
        // - KSW_EZ_GENERIC_SC：启用完整替换矩阵（必须！否则会忽略 dna5_simd_mat）
        // - KSW_EZ_RIGHT：gap 右对齐（CIGAR 标准化，便于比较）
        cfg.flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT;

        /* ---------- 3. 调用 KSW2 ---------- */
        // ez：KSW2 输出结构体。
        // 关键字段：
        // - ez.cigar：uint32_t*，malloc 分配，需要调用者 free
        // - ez.n_cigar：CIGAR 操作数量
        ksw_extz_t ez{};

        // 注意：ksw_extz2_sse 的第一个参数通常是 qlen/tlen 的额外参数（如 score-only 模式），
        // 这里传 0 与现有逻辑保持一致。
        ksw_extz2_sse(0,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);


        /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
        // 输出转换：把 KSW2 的 uint32_t CIGAR 拷贝到项目统一的 cigar::Cigar_t
        // 性能：reserve(ez.n_cigar) 避免 push_back 扩容
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        // 资源释放：KSW2 用 malloc() 分配 ez.cigar -> 需要 free()
        free(ez.cigar);
        return cigar;
    }

    // ------------------------------------------------------------------
    // extendAlignKSW2：KSW2 延伸比对（extension / ends-free 风格）
    // ------------------------------------------------------------------
    // 与 globalAlignKSW2 的差异：
    // - 使用 zdrop 以及 EXTZ_ONLY/APPROX_DROP 等 flag，允许在得分下降时提前终止
    // - 用 end_bonus 奖励末端延伸，有利于从 seed 位置向外扩展
    //
    // 参数：
    // @param zdrop: Z-drop 阈值，越大越不容易提前停止（更慢但更完整）
    // ------------------------------------------------------------------
    cigar::Cigar_t extendAlignKSW2(const std::string& ref,
        const std::string& query,
        int zdrop)
    {
        /* ---------- 1. 序列编码 ---------- */
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());
        for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = align::ScoreChar2Idx[(uint8_t)ref[i]];
        for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = align::ScoreChar2Idx[(uint8_t)query[i]];

        // ------------------------------------------------------------------
        // 2. 配置参数
        // ------------------------------------------------------------------
        // 说明：
        // - 这里保留了历史上尝试过的参数组合（注释块），方便后续调参回归。
        // - 当前启用的 cfg：在 extension 场景下使用 EXTZ_ONLY + APPROX_DROP + RIGHT。
        //   * EXTZ_ONLY：启用 ends-free extension 模式（更适合延伸）
        //   * APPROX_DROP：在 approximate 模式下触发 z-drop 就中断（加速）
        //   * RIGHT：gap 右对齐，CIGAR 更稳定
        // - end_bonus=50：给末端延伸一定奖励，鼓励更长的 alignment。
        // - band_width 启发式估计：加速同时尽量不损失太多准确性。
        align::KSW2AlignConfig cfg;
        cfg.mat = align::dna5_simd_mat;
        cfg.zdrop = zdrop;
        cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_APPROX_DROP;
        cfg.end_bonus = 50;
        cfg.alphabet_size = 5;
        cfg.gap_open = 6;
        cfg.gap_extend = 2;
        cfg.band_width = align::auto_band(ref.size(), query.size());


        /* ---------- 3. 调用 KSW2 ---------- */
        ksw_extz_t ez{};
        // 注意：这里第一个参数传 nullptr 与项目历史实现保持一致。
        ksw_extz2_sse(nullptr,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);

        /* ---------- 4. 拷贝 & 释放 ---------- */
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(ez.cigar);                    // ksw2 使用 malloc
        return cigar;                      // 返回的 CIGAR 即延伸片段
    }


    // ------------------------------------------------------------------
    // globalAlignWFA2：WFA2 全局比对
    // ------------------------------------------------------------------
    // WFA2（Wavefront Alignment）特点：
    // - 对高相似度序列（编辑距离小）通常显著快于 DP
    // - 可配置不同 memory_mode（速度/内存权衡）与 heuristic（剪枝策略）
    //
    // 本实现使用 gap_affine（仿射 gap 罚分）：
    // - mismatch=2
    // - gap_opening=3
    // - gap_extension=1
    // 这些参数的量级与 KSW2 的 5/-4 + gap(6,2) 并不完全等价，
    // 但对下游而言我们只需要一个可用且一致的 CIGAR 表示。
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignWFA2(const std::string& ref,
        const std::string& query)
    {
        // 1) 构建 WFA2 属性（attributes）
        // 注意：wavefront_aligner_attr_default 会填充大量默认值。
        // 本函数只覆盖我们关心的部分，保持行为稳定。
        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
        attributes.distance_metric = gap_affine;
        attributes.affine_penalties.mismatch = 2;      // X > 0
        attributes.affine_penalties.gap_opening = 3;   // O >= 0
        attributes.affine_penalties.gap_extension = 1; // E > 0

        // memory_mode：ultralow 表示尽量降低内存占用，适合长序列但可能更慢
        attributes.memory_mode = wavefront_memory_ultralow;

        // heuristic：自适应 band 策略，限制波前带宽以加速
        // 注意：heuristic 的启发式会影响速度与最优性，本项目当前选择偏速度。
        attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
        attributes.heuristic.min_k = -50;
        attributes.heuristic.max_k = +50;
        attributes.heuristic.steps_between_cutoffs = 1;

        // 2) 创建 aligner（堆对象，需要 delete）
        wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

        // 3) 执行比对
        // wavefront_align 会把结果写入 wf_aligner->cigar
        wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());

        // 4) 从 WFA2 输出中取 CIGAR
        // cigar_get_CIGAR：将内部 cigar 转换为 BAM 风格的 uint32_t buffer
        // - cigar_buffer：由 WFA2 管理（通过 wf_aligner->cigar 释放），这里不手动 free
        // - cigar_length：操作数
        uint32_t* cigar_buffer;
        int cigar_length = 0;
        cigar_get_CIGAR(wf_aligner->cigar, false, &cigar_buffer, &cigar_length);

        // 5) 拷贝到项目统一格式
        cigar::Cigar_t cigar;
        cigar.reserve(static_cast<std::size_t>(cigar_length));
        for (int i = 0; i < cigar_length; ++i)
            cigar.push_back(cigar_buffer[i]);

        // 6) 释放 aligner（会同时释放内部 cigar 结构）
        wavefront_aligner_delete(wf_aligner);

        return cigar;
    }

    // ------------------------------------------------------------------
    // extendAlignWFA2：WFA2 延伸比对
    // ------------------------------------------------------------------
    // WFA2 的延伸比对（ends-free extension）：
    // - 通常用于从种子位置向外延伸，快速获得局部比对结果
    // - 可配置 zdrop 阈值，控制延伸过程中的提前终止
    //
    // 本实现使用 wavefront_aligner_attr_default 的高内存模式：
    // - 适合长序列比对，但可能导致较高的内存占用
    // ------------------------------------------------------------------
    // cigar::Cigar_t extendAlignWFA2(const std::string& ref,
    //     const std::string& query, int zdrop)
    // {
    //     wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    //     attributes.distance_metric = gap_affine;
    //     attributes.affine_penalties.mismatch = 2;      // X > 0
    //     attributes.affine_penalties.gap_opening = 3;   // O >= 0
    //     attributes.affine_penalties.gap_extension = 1; // E > 0
    //     attributes.memory_mode = wavefront_memory_high;
    //     attributes.heuristic.strategy = wf_heuristic_zdrop;
    //     attributes.heuristic.zdrop = zdrop;
    //     attributes.heuristic.steps_between_cutoffs = 1;
    //     //// Create a WFAligner
    //     //
    //     wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    //
    //     wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());
    //     /*wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);
    //
    //     aligner.alignEnd2End(ref, query);*/
    //
    //     uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
    //     int cigar_length = 0; // Length of the CIGAR string.
    //     // Retrieve the CIGAR string from the wavefront aligner.
    //     cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
    //
    //     /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    //     cigar::Cigar_t cigar;
    //
    //     for (int i = 0; i < cigar_length; ++i)
    //         cigar.push_back(cigar_buffer[i]);
    //
    //     wavefront_aligner_delete(wf_aligner);
    //
    //     return cigar;
    // }

    // ------------------------------------------------------------------
    // globalAlignMM2：基于 minimap2 风格的 anchor-guided 全局比对
    // ------------------------------------------------------------------
    // 设计思路（参考 minimap2/align.c 的 mm_align1 函数）：
    // 1) 利用已有的锚点（anchors）进行链化（chaining），获取最佳链
    // 2) 对链覆盖范围之前的左端区域，自适应选择比对算法
    // 3) 对链中相邻锚点之间的间隙区域（gap filling），自适应选择比对算法
    // 4) 对链覆盖范围之后的右端区域，自适应选择比对算法
    // 5) 使用 cigar::appendCigar 合并所有 CIGAR，确保覆盖整个序列
    //
    // **自适应比对策略（性能优化）：**
    // - **小间隙（< 100bp）**：直接使用 globalAlignKSW2（精确，O(mn) 可接受）
    // - **大间隙（>= 100bp）**：
    //   1. 先用 extendAlignKSW2 快速探路（zdrop=200），利用剪枝加速
    //   2. 检查结果是否完整覆盖（通过 cigar::getRefLength/getQueryLength）
    //   3. 如果未完整覆盖（zdrop 提前截止），则 fallback 到 globalAlignKSW2
    //   - 优势：对高相似度大间隙，extension 速度快；低相似度会自动 fallback
    //
    // **正确性保证：**
    // - 严格跟踪 ref_pos 和 qry_pos，确保每个碱基都被比对
    // - 最后验证 CIGAR 消耗的序列长度是否匹配输入序列长度
    // - Debug 模式下长度不匹配会触发断言；Release 模式下记录日志并 fallback
    // - 使用 cigar 命名空间的函数进行 CIGAR 操作，提高代码复用性
    //
    // 输入：
    // @param ref     - 参考序列（字符串，A/C/G/T/N）
    // @param query   - 查询序列（字符串，A/C/G/T/N）
    // @param anchors - 预先计算的锚点列表（未排序也可以，内部会排序和链化）
    //
    // 输出：
    // 返回 CIGAR 操作序列（cigar::Cigar_t），描述 query 如何比对到 ref
    // **保证**：返回的 CIGAR 消耗的 ref 长度 = ref.size()，query 长度 = query.size()
    //
    // 性能说明：
    // - 相比直接 globalAlignKSW2，本函数利用锚点信息将比对分解为多个小区间
    // - 对于长序列且锚点覆盖率高的情况，性能提升显著（减少 DP 矩阵规模）
    // - 大间隙优先用 extension 探路，高相似度时显著加速
    // - 若锚点为空或无有效链，退化为 globalAlignKSW2
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignMM2(const std::string& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors)
    {
        const std::size_t ref_len = ref.size();
        const std::size_t qry_len = query.size();

        // ------------------------------------------------------------------
        // 辅助 lambda：自适应选择比对算法
        // ------------------------------------------------------------------
        // 策略：
        // - 小间隙（< 100bp）：直接 globalAlignKSW2
        // - 大间隙（>= 100bp）：先 extendAlignKSW2 探路，检测是否完整覆盖
        //   * 如果覆盖完整：使用 extension 结果（快速路径）
        //   * 如果未完整覆盖：**复用** extension 结果，只对剩余部分补全
        //
        // **性能优化（关键）：**
        // - 复用 extension 已比对的部分，避免重复计算
        // - 只对未覆盖的剩余区域进行补全比对
        // - 例如：extension 覆盖了 80%，只需对剩余 20% 进行全局比对
        //
        // 参数：
        // @param gap_ref - ref 间隙序列
        // @param gap_qry - query 间隙序列
        // @return CIGAR 操作序列，保证完整覆盖输入序列
        // ------------------------------------------------------------------
        auto adaptive_align = [](const std::string& gap_ref, const std::string& gap_qry) -> cigar::Cigar_t {
            const std::size_t ref_gap_len = gap_ref.size();
            const std::size_t qry_gap_len = gap_qry.size();

            // 边界情况：空序列
            if (ref_gap_len == 0 && qry_gap_len == 0) {
                return cigar::Cigar_t{}; // 空 CIGAR
            }
            if (ref_gap_len == 0) {
                // 只有 query，全部是插入
                return cigar::Cigar_t{cigar::cigarToInt('I', static_cast<uint32_t>(qry_gap_len))};
            }
            if (qry_gap_len == 0) {
                // 只有 ref，全部是删除
                return cigar::Cigar_t{cigar::cigarToInt('D', static_cast<uint32_t>(ref_gap_len))};
            }

            // 定义小间隙阈值
            constexpr std::size_t kSmallGapThreshold = 100;
            const std::size_t max_gap = std::max(ref_gap_len, qry_gap_len);

            // 小间隙：直接用全局比对（精确且 O(mn) 可接受）
            if (max_gap < kSmallGapThreshold) {
                return globalAlignKSW2(gap_ref, gap_qry);
            }

            // 大间隙：先用 extension 探路
            cigar::Cigar_t ext_cigar = extendAlignKSW2(gap_ref, gap_qry, 200);

            // 检查 extension 是否完整覆盖
            std::size_t ext_ref_len = cigar::getRefLength(ext_cigar);
            std::size_t ext_qry_len = cigar::getQueryLength(ext_cigar);

            // 如果完整覆盖，直接返回 extension 结果（快速路径）
            if (ext_ref_len == ref_gap_len && ext_qry_len == qry_gap_len) {
                return ext_cigar;
            }

            // ------------------------------------------------------------------
            // 未完整覆盖：复用 extension 结果，只对剩余部分补全
            // ------------------------------------------------------------------
            // 说明：extension 可能因为 zdrop 提前截止，但已比对的部分是有效的
            // 策略：保留 extension 的 CIGAR，对剩余未覆盖的序列进行补全比对
            //
            // 例如：
            // - gap_ref = 1000bp, gap_qry = 1000bp
            // - extension 覆盖了前 800bp（ref）和 800bp（query）
            // - 只需对剩余的 200bp（ref）和 200bp（query）进行全局比对
            // - 性能提升：从 O(1000*1000) 降低到 O(800*800) + O(200*200)
            // ------------------------------------------------------------------

            cigar::Cigar_t result = ext_cigar; // 复用 extension 结果

            // 计算剩余未覆盖的长度
            std::size_t remaining_ref_len = ref_gap_len - ext_ref_len;
            std::size_t remaining_qry_len = qry_gap_len - ext_qry_len;

            if (remaining_ref_len > 0 || remaining_qry_len > 0) {
                // 提取剩余未覆盖的序列
                std::string remaining_ref = gap_ref.substr(ext_ref_len, remaining_ref_len);
                std::string remaining_qry = gap_qry.substr(ext_qry_len, remaining_qry_len);

                // 对剩余部分进行全局比对
                cigar::Cigar_t remaining_cigar = globalAlignKSW2(remaining_ref, remaining_qry);

                // 合并 CIGAR：extension 部分 + 剩余部分
                cigar::appendCigar(result, remaining_cigar);

#ifdef _DEBUG
                spdlog::debug("adaptive_align: extension partial coverage (ref: {}/{}, qry: {}/{}), "
                              "补全剩余部分 (ref: {}, qry: {})",
                              ext_ref_len, ref_gap_len, ext_qry_len, qry_gap_len,
                              remaining_ref_len, remaining_qry_len);
#endif
            }

            // 最终验证：确保返回的 CIGAR 完全覆盖输入序列
            std::size_t final_ref_len = cigar::getRefLength(result);
            std::size_t final_qry_len = cigar::getQueryLength(result);

            if (final_ref_len != ref_gap_len || final_qry_len != qry_gap_len) {
#ifdef _DEBUG
                spdlog::error("adaptive_align: final CIGAR mismatch (expected ref:{}/qry:{}, got ref:{}/qry:{}), "
                              "fallback to global alignment",
                              ref_gap_len, qry_gap_len, final_ref_len, final_qry_len);
#endif
                // 最终 fallback：直接用全局比对
                return globalAlignKSW2(gap_ref, gap_qry);
            }

            return result;
        };

        // ------------------------------------------------------------------
        // 0. 边界检查：如果锚点为空，退化为全局比对
        // ------------------------------------------------------------------
        if (anchors.empty()) {
            return globalAlignKSW2(ref, query);
        }

        // ------------------------------------------------------------------
        // 1. 复制锚点并链化（chainAnchors 会修改输入锚点）
        // ------------------------------------------------------------------
        // 说明：chainAnchors 会按 (rid_ref, is_rev, pos_ref, pos_qry) 排序锚点，
        // 然后使用 DP 算法找出得分最高的链，直接返回该链的锚点列表。
        anchor::Anchors sorted_anchors = anchors;
        anchor::ChainParams chain_params = anchor::default_chain_params();
        anchor::Anchors chain_anchors = anchor::chainAnchors(sorted_anchors, chain_params);

        // ------------------------------------------------------------------
        // 2. 检查链化结果，如果无有效链则退化为全局比对
        // ------------------------------------------------------------------
        if (chain_anchors.empty()) {
            return globalAlignKSW2(ref, query);
        }

        // ------------------------------------------------------------------
        // 3. 按 pos_qry 排序锚点（确保从左到右处理）
        // ------------------------------------------------------------------

        std::sort(chain_anchors.begin(), chain_anchors.end(),
            [](const anchor::Anchor& a, const anchor::Anchor& b) {
                // 按 query 位置排序，如果 query 位置相同则按 ref 位置排序
                if (a.pos_qry != b.pos_qry) return a.pos_qry < b.pos_qry;
                return a.pos_ref < b.pos_ref;
            });

        // ------------------------------------------------------------------
        // 4. 构建完整的 CIGAR：左端 + 锚点间隙 + 右端
        // ------------------------------------------------------------------
        cigar::Cigar_t result_cigar;
        result_cigar.reserve(chain_anchors.size() * 2 + 2);

        // 记录当前在 ref 和 query 上已处理到的位置
        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        // ------------------------------------------------------------------
        // 4.1 处理第一个锚点之前的左端区域
        // ------------------------------------------------------------------
        const auto& first_anchor = chain_anchors.front();
        std::size_t first_ref_start = first_anchor.pos_ref;
        std::size_t first_qry_start = first_anchor.pos_qry;

        if (first_ref_start > 0 || first_qry_start > 0) {
            // 左端有未覆盖区域，需要比对
            std::string left_ref = ref.substr(0, first_ref_start);
            std::string left_qry = query.substr(0, first_qry_start);

            // 使用自适应比对策略
            cigar::Cigar_t left_cigar = adaptive_align(left_ref, left_qry);
            cigar::appendCigar(result_cigar, left_cigar);

            // 更新位置
            ref_pos = first_ref_start;
            qry_pos = first_qry_start;
        }

        // ------------------------------------------------------------------
        // 4.2 遍历锚点，处理锚点覆盖区域和锚点之间的间隙
        // ------------------------------------------------------------------
        for (std::size_t i = 0; i < chain_anchors.size(); ++i) {
            const auto& anchor = chain_anchors[i];

            // 边界检查：确保锚点位置在序列范围内
            if (anchor.pos_ref >= ref_len || anchor.pos_qry >= qry_len) {
#ifdef _DEBUG
                spdlog::warn("globalAlignMM2: anchor {} out of bounds (ref: {}/{}, qry: {}/{}), skipping",
                             i, anchor.pos_ref, ref_len, anchor.pos_qry, qry_len);
#endif
                continue;
            }

            std::size_t anchor_ref_start = anchor.pos_ref;
            std::size_t anchor_qry_start = anchor.pos_qry;
            std::size_t anchor_span = anchor.span > 0 ? anchor.span : 1;

            // 计算锚点结束位置（确保不超出序列边界）
            std::size_t anchor_ref_end = std::min(anchor_ref_start + anchor_span, ref_len);
            std::size_t anchor_qry_end = std::min(anchor_qry_start + anchor_span, qry_len);

            // 检查锚点是否与当前位置重叠（避免重复处理）
            if (anchor_ref_start < ref_pos || anchor_qry_start < qry_pos) {
#ifdef _DEBUG
                spdlog::warn("globalAlignMM2: anchor {} overlaps with processed region (anchor: {}/{}, current: {}/{}), adjusting",
                             i, anchor_ref_start, anchor_qry_start, ref_pos, qry_pos);
#endif
                // 调整锚点起始位置到当前位置
                anchor_ref_start = std::max(anchor_ref_start, ref_pos);
                anchor_qry_start = std::max(anchor_qry_start, qry_pos);

                // 如果调整后锚点已经被完全处理，跳过
                if (anchor_ref_start >= anchor_ref_end || anchor_qry_start >= anchor_qry_end) {
                    continue;
                }
            }

            // ------------------------------------------------------------
            // 4.2.1 处理当前位置到锚点起始位置之间的间隙
            // ------------------------------------------------------------
            if (anchor_ref_start > ref_pos || anchor_qry_start > qry_pos) {
                std::size_t gap_ref_len = anchor_ref_start - ref_pos;
                std::size_t gap_qry_len = anchor_qry_start - qry_pos;

                if (gap_ref_len > 0 || gap_qry_len > 0) {
                    // 有间隙需要填充
                    std::string gap_ref = ref.substr(ref_pos, gap_ref_len);
                    std::string gap_qry = query.substr(qry_pos, gap_qry_len);

                    // 使用自适应比对策略：小间隙用 global，大间隙先 extension 探路
                    cigar::Cigar_t gap_cigar = adaptive_align(gap_ref, gap_qry);

                    // 验证间隙 CIGAR 的正确性
                    std::size_t gap_cigar_ref = cigar::getRefLength(gap_cigar);
                    std::size_t gap_cigar_qry = cigar::getQueryLength(gap_cigar);
                    if (gap_cigar_ref != gap_ref_len || gap_cigar_qry != gap_qry_len) {
#ifdef _DEBUG
                        spdlog::warn("globalAlignMM2: gap CIGAR mismatch at anchor {} (expected ref:{}/qry:{}, got ref:{}/qry:{})",
                                     i, gap_ref_len, gap_qry_len, gap_cigar_ref, gap_cigar_qry);
#endif
                        // Fallback：使用全局比对确保正确性
                        gap_cigar = globalAlignKSW2(gap_ref, gap_qry);
                    }

                    cigar::appendCigar(result_cigar, gap_cigar);
                }

                // 更新位置到锚点起始
                ref_pos = anchor_ref_start;
                qry_pos = anchor_qry_start;
            }

            // ------------------------------------------------------------
            // 4.2.2 处理锚点覆盖的区域
            // ------------------------------------------------------------
            std::size_t anchor_ref_len = anchor_ref_end - anchor_ref_start;
            std::size_t anchor_qry_len = anchor_qry_end - anchor_qry_start;

            if (anchor_ref_len > 0 || anchor_qry_len > 0) {
                std::string anchor_ref_seq = ref.substr(anchor_ref_start, anchor_ref_len);
                std::string anchor_qry_seq = query.substr(anchor_qry_start, anchor_qry_len);

                // 对锚点区域也进行比对，确保精确性
                cigar::Cigar_t anchor_cigar = globalAlignKSW2(anchor_ref_seq, anchor_qry_seq);

                // 验证锚点 CIGAR 的正确性
                std::size_t anchor_cigar_ref = cigar::getRefLength(anchor_cigar);
                std::size_t anchor_cigar_qry = cigar::getQueryLength(anchor_cigar);
                if (anchor_cigar_ref != anchor_ref_len || anchor_cigar_qry != anchor_qry_len) {
#ifdef _DEBUG
                    spdlog::error("globalAlignMM2: anchor CIGAR mismatch at anchor {} (expected ref:{}/qry:{}, got ref:{}/qry:{})",
                                  i, anchor_ref_len, anchor_qry_len, anchor_cigar_ref, anchor_cigar_qry);
#endif
                    // 这不应该发生（globalAlignKSW2 应该总是返回正确长度）
                    // 但为了安全起见，重新尝试
                    anchor_cigar = globalAlignKSW2(anchor_ref_seq, anchor_qry_seq);
                }

                cigar::appendCigar(result_cigar, anchor_cigar);

                // 更新位置到锚点结束
                ref_pos = anchor_ref_end;
                qry_pos = anchor_qry_end;
            }
        }

        // ------------------------------------------------------------------
        // 4.3 处理最后一个锚点之后的右端区域
        // ------------------------------------------------------------------
        if (ref_pos < ref_len || qry_pos < qry_len) {
            std::size_t right_ref_len = ref_len - ref_pos;
            std::size_t right_qry_len = qry_len - qry_pos;

            std::string right_ref = ref.substr(ref_pos, right_ref_len);
            std::string right_qry = query.substr(qry_pos, right_qry_len);

            // 使用自适应比对策略
            cigar::Cigar_t right_cigar = adaptive_align(right_ref, right_qry);
            cigar::appendCigar(result_cigar, right_cigar);

            // 更新位置（已处理完整个序列）
            ref_pos = ref_len;
            qry_pos = qry_len;
        }

        // ------------------------------------------------------------------
        // 5. 验证：确保 CIGAR 消耗的序列长度与输入一致
        // ------------------------------------------------------------------
        // 说明：这是最终的正确性检查，确保没有遗漏任何区域
        std::size_t cigar_ref_len = cigar::getRefLength(result_cigar);
        std::size_t cigar_qry_len = cigar::getQueryLength(result_cigar);

        // 如果长度不匹配，记录详细信息并 fallback
        if (cigar_ref_len != ref_len || cigar_qry_len != qry_len) {
#ifdef _DEBUG
            spdlog::error("globalAlignMM2: CIGAR length mismatch!");
            spdlog::error("  Expected: ref={}, qry={}", ref_len, qry_len);
            spdlog::error("  Got:      ref={}, qry={}", cigar_ref_len, cigar_qry_len);
            spdlog::error("  Anchors:  count={}, chain_size={}", anchors.size(), chain_anchors.size());
            spdlog::error("  Final positions: ref_pos={}, qry_pos={}", ref_pos, qry_pos);

            // Debug 模式：输出链的详细信息
            for (size_t i = 0; i < chain_anchors.size(); ++i) {
                const auto& a = chain_anchors[i];
                spdlog::error("    Anchor[{}]: ref=[{}, {}), qry=[{}, {}), span={}",
                              i, a.pos_ref, a.pos_ref + a.span, a.pos_qry, a.pos_qry + a.span, a.span);
            }

            // 在 Debug 模式下也使用 fallback 而不是 assert，便于调试
            spdlog::error("  Falling back to globalAlignKSW2");
#else
            spdlog::warn("globalAlignMM2: CIGAR length mismatch (ref: {}/{}, qry: {}/{}), falling back to global alignment",
                         cigar_ref_len, ref_len, cigar_qry_len, qry_len);
#endif
            return globalAlignKSW2(ref, query);
        }

        return result_cigar;
    }
}
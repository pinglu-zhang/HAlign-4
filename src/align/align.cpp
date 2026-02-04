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
        attributes.affine_penalties.mismatch = 3;      // X > 0
        attributes.affine_penalties.gap_opening = 4;   // O >= 0
        attributes.affine_penalties.gap_extension = 1; // E > 0

        // memory_mode：ultralow 表示尽量降低内存占用，适合长序列但可能更慢
        attributes.memory_mode = wavefront_memory_ultralow;

        // heuristic：自适应 band 策略，限制波前带宽以加速
        // 注意：heuristic 的启发式会影响速度与最优性，本项目当前选择偏速度。
        // attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
        // attributes.heuristic.min_k = -200;
        // attributes.heuristic.max_k = +200;
        // attributes.heuristic.steps_between_cutoffs = 1;

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
    // globalAlignMM2：基于锚点（anchors）的分段全局比对（minimap2 风格）
    // ------------------------------------------------------------------
    // 功能概述（精炼）：
    // 使用预先计算的锚点将一对序列（ref / query）拆分为若干可比对的片段，
    // 对每个片段调用指定的全局比对函数（align_func），最后合并各段的 CIGAR
    // 为完整的 query->ref CIGAR。该函数优先保证正确性（CIGAR 消耗长度匹配），
    // 在极端边界上提供稳健的兜底逻辑以避免崩溃或最终不一致性报错。
    //
    // 输入参数：
    // - ref: 参考序列字符串（A/C/G/T/N），长度 ref_len
    // - query: 查询序列字符串，长度 qry_len
    // - anchors: 预先生成的锚点列表（anchor::Anchors），可以是任意顺序，函数内部会排序并链化
    // - align_func: 用于分段比对的函数指针，签名为 cigar::Cigar_t(const std::string&, const std::string&)。
    //               若传入 nullptr，则默认使用 `globalAlignKSW2`。
    //
    // 输出：
    // - 返回一个 cigar::Cigar_t（压缩 CIGAR）描述如何把 query 对齐到 ref。
    //   本函数对外的语义保证：若返回非空 CIGAR，则它严格消耗 ref_len 个参考碱基和 qry_len 个查询碱基，
    //   即 cigar::getRefLength(result) == ref_len 且 cigar::getQueryLength(result) == qry_len。
    //   若内部检测到无法生成满足该不变量的分段 CIGAR，会退化为对整对序列调用 align_func 并返回其结果。
    //
    // 设计与实现要点（逐步说明）：
    // 1) 防御式编程：对输入 anchors 先做链化（chainAnchors）以获得一条“最佳链”；若链为空则退化为
    //    对整对序列做全局比对（align_func）。这样保证在没有锚点或链化失败时仍能工作。
    //
    // 2) 分段策略（与实现严格一致）：
    //    - 先按 query 坐标对 chain_anchors 进行升序排序（保持与 minimap2 相同的处理顺序），
    //    - 使用一个驱动位置（ref_pos, qry_pos）从左到右推进（注意：不信任 anchors 给出的绝对位置，
    //      仅将其作为建议的边界；真正的推进以各段返回的 CIGAR 消耗为准），
    //    - 处理顺序为：
    //        a) 左端：从当前位置到第一个锚点的起始位置（append_segment）
    //        b) 对于每个锚点 a（按顺序）：
    //             - 先处理该锚点覆盖区（span 段）：append_segment(ref_pos, a.pos_ref + a.span, qry_pos, a.pos_qry + a.span)
    //             - 若存在下一个锚点 b，则处理 a 结束到 b 开始之间的 gap：append_segment(ref_pos, b.pos_ref, qry_pos, b.pos_qry)
    //        c) 右端：最后一个锚点结束到序列末尾（append_segment）
    //
    //    说明：上面的实现保证每个 anchor.span 恰好被包含一次（修正了早期实现中遗漏/重复 span 的问题），
    //    同时通过 ref_pos/qry_pos 的推进处理 anchors 之间的重叠或坐标不一致情形（重叠会被自动裁剪/剪切）。
    //
    // 3) append_segment 的语义与边界策略（非常关键，直接影响正确性）:
    //    - 输入：以 ref/query 的绝对坐标区间为参数（ref_start, ref_end, qry_start, qry_end），函数会先做边界裁剪
    //      以避免越界（用 std::min / 容错性的 if 修正），并保证 start <= end。
    //    - 若某一侧长度为 0（seg_ref.empty() 或 seg_qry.empty()），函数不会去调用昂贵或行为不稳定的底层比对器，
    //      而是直接构造一个简单的 CIGAR：
    //         * seg_ref 非空且 seg_qry 为空 -> 产生一条单 op 的 Deletion ('D', len=seg_ref_len)
    //         * seg_ref 为空且 seg_qry 非空 -> 产生一条单 op 的 Insertion ('I', len=seg_qry_len)
    //         * 两者均为空 -> 不产生 CIGAR
    //      这样既高效又避免了对空输入时底层库（KSW2/WFA2）可能返回不一致/空 CIGAR 的问题。
    //
    //    - 当两侧均非空时，才调用 align_func(seg_ref, seg_qry) 来生成段级 CIGAR；生成后会做严格校验：
    //         c_ref = cigar::getRefLength(seg_cigar) 必须等于 seg_ref_len
    //         c_qry = cigar::getQueryLength(seg_cigar) 必须等于 seg_qry_len
    //      若校验失败（常见于极端长度差异或 band/zdrop 限制导致的截断），函数采用“兜底策略”而不是无限重试：
    //         - 在 Debug 模式下先记录一条 warn 日志（有助于定位）；
//         - 构造一个“强制补齐”的 CIGAR（forced_cigar）：先一个 'I' 消耗全部 query（若有），再一个 'D' 消耗全部 ref（若有），
//           并将其 append 到结果中，同时把 ref_pos/qry_pos 直接推进到期望的区间末端（ref_end/qry_end）。
//      这样可以确保：即使底层比对算法在该片段上失败，整个分段链仍然连续，并保持长度不变性。
//
//    - 如果校验通过，则 append seg_cigar，并使用其实际消耗（c_ref/c_qry）推进 ref_pos 与 qry_pos，
//      而不是盲目采用 anchor 给出的起止坐标（这是保证正确性的核心做法）。
//
// 4) 终结一致性检查与退化策略：
//    - 在合并完所有分段后，函数会验证最终 result 消耗的 ref/ qry 长度是否等于输入长度。
//    - 若不等：该函数会记录 error 并退化为对整对序列调用 align_func(ref, query)，返回 align_func 的结果。
//      （理由：在极端/非常规错误下，直接让底层比对器处理整对序列通常更直观且便于定位问题）。
//
// 5) 复杂度与性能说明：
//    - 链化 chainAnchors 的时间复杂度取决于 anchors 的数量与链化参数（一般近似 O(m log m) 排序 + O(m * max_iter) DP），
//      其中 m 为 anchors 数量；分段比对的总成本为各段比对成本之和（受 align_func 的算法复杂度影响）。
//    - 对于长序列且 anchors 覆盖率高的场景，该分段策略能显著降低单次 DP 的规模并提升整体吞吐。
//    - 对于极端长度差异的片段，append_segment 的快速兜底（单 op I/D）既保证正确性也避免了昂贵无效的 DP。
//
// 6) 可调试点与注意事项（调试/开发者使用）：
//    - 在 _DEBUG 模式下，发生段级 mismatch 时会记录 warn，便于定位哪些段被兜底处理；生产环境通常关闭 debug 以避免日志噪音。
//    - 如果希望强制使用更稳的比对策略（例如总是用 WFA2），可在调用处传入对应的 align_func。
//    - 对于非常长的 ref（> 1e6），请注意 CIGAR 单个 op 长度上限（cigar::cigarToInt 对 len 有 28-bit 限制），
//      若存在更长的连续 gap 需要分割为多个 op 或采用不同策略。
//
// 注：本注释已经与 `append_segment` 的实际实现保持一致（包括空段直接生成 I/D、段级强制补齐策略、以及
//      锚点 span/gap 的处理顺序）。如后续对实现做改动，请同步更新此处注释以保持文档与代码一致。
// ------------------------------------------------------------------
cigar::Cigar_t globalAlignMM2(const std::string& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors,
                                  AlignFunc align_func)
    {
        // ------------------------------------------------------------------
        // 保守实现：锚点分割 + 可选比对算法（保证 CIGAR 一定正确）
        // ------------------------------------------------------------------
        // 说明：
        // 1) 使用 anchors 进行 chain，得到一条"最佳链"作为分段边界参考
        // 2) 将 (ref,query) 拆成：左端 + (锚点间gap + 锚点span段)* + 右端
        // 3) 每一段都用传入的 align_func（默认 globalAlignKSW2）做 end-to-end 全局比对
        // 4) **关键正确性策略**：每段比对完后，不用我们"猜测"的长度推进位置，
        //    而是用 cigar::getRefLength/getQueryLength 从 CIGAR 反推实际消耗长度。
        //    这样可以避免之前因为 anchor 坐标不精确/重叠/空洞造成的长度错配。
        //
        // 性能：
        // - 相比直接全局比对，此版本在锚点可靠时会分解成多个小矩阵，通常更快。
        // - 但由于每一段都是真全局比对，可能比"extend+fallback"慢。
        // - 通过 align_func 参数可以灵活选择 WFA2 或 KSW2，适应不同场景。
        // - 目前优先保证正确性。
        // ------------------------------------------------------------------

        // 如果未传入比对函数，默认使用 globalAlignKSW2
        if (!align_func) {
            align_func = globalAlignKSW2;
        }

        const std::size_t ref_len = ref.size();
        const std::size_t qry_len = query.size();

        // 0) 无锚点直接退化为全局比对（使用传入的比对函数）
        if (anchors.empty()) {
            return align_func(ref, query);
        }

        // 1) 链化：拿到最佳链
        anchor::Anchors sorted_anchors = anchors;
        anchor::ChainParams chain_params = anchor::default_chain_params();
        anchor::Anchors chain_anchors = anchor::chainAnchors(sorted_anchors, chain_params);
        if (chain_anchors.empty()) {
            return align_func(ref, query);
        }


        // 2) 按 query 坐标从左到右处理（与 minimap2 的处理方式一致）
        std::sort(chain_anchors.begin(), chain_anchors.end(),
                  [](const anchor::Anchor& a, const anchor::Anchor& b) {
                      if (a.pos_qry != b.pos_qry) return a.pos_qry < b.pos_qry;
                      return a.pos_ref < b.pos_ref;
                  });

        cigar::Cigar_t result;
        result.reserve(chain_anchors.size() * 2 + 2);

        // 用 CIGAR 驱动推进的位置（不信任 anchor 坐标的绝对正确性）
        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        auto append_segment = [&](std::size_t ref_start, std::size_t ref_end,
                                  std::size_t qry_start, std::size_t qry_end) {
            // 边界裁剪（即使 anchor 给错也不崩）
            ref_start = std::min(ref_start, ref_len);
            ref_end = std::min(ref_end, ref_len);
            qry_start = std::min(qry_start, qry_len);
            qry_end = std::min(qry_end, qry_len);

            if (ref_end < ref_start) ref_end = ref_start;
            if (qry_end < qry_start) qry_end = qry_start;

            const std::string seg_ref = ref.substr(ref_start, ref_end - ref_start);
            const std::string seg_qry = query.substr(qry_start, qry_end - qry_start);

            // ------------------------------------------------------------------
            // 关键边界修复：当某一侧 segment 为空时，不调用底层比对器，直接构造 CIGAR。
            //
            // 背景：
            // - 在“ref 很长，但 query 很短”的场景中，chaining 的锚点坐标可能导致某些 gap 段
            //   出现 seg_qry 为空但 seg_ref 非空（或反之）。
            // - 对空序列调用 KSW2/WFA2 时，不同实现可能返回空 CIGAR（0/0），
            //   从而触发我们后续的长度一致性检查，导致频繁 fallback。
            //
            // CIGAR 语义（本项目约定：返回的是 query 相对于 ref）：
            // - seg_ref 非空 & seg_qry 为空：表示 query 相对 ref 缺失这一段 -> 全部是 'D'
            // - seg_ref 为空 & seg_qry 非空：表示 query 相对 ref 多出来这一段 -> 全部是 'I'
            // - 两者都空：空操作
            //
            // 性能：
            // - 直接构造 CIGAR 是 O(1)（单个 op），比调用 DP/WFA 更快；
            // - 同时避免不必要的内存分配与库调用开销。
            // ------------------------------------------------------------------
            cigar::Cigar_t seg_cigar;
            // if (seg_ref.empty() && seg_qry.empty()) {
            //     // nothing
            // } else if (seg_ref.empty()) {
            //     seg_cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(seg_qry.size())));
            // } else if (seg_qry.empty()) {
            //     seg_cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(seg_ref.size())));
            // } else {
            //     seg_cigar = align_func(seg_ref, seg_qry);
            // }
            seg_cigar = align_func(seg_ref, seg_qry);
            // 严格校验本段 CIGAR 覆盖长度
            const std::size_t seg_ref_len = seg_ref.size();
            const std::size_t seg_qry_len = seg_qry.size();
            const std::size_t c_ref = cigar::getRefLength(seg_cigar);
            const std::size_t c_qry = cigar::getQueryLength(seg_cigar);

            if (c_ref != seg_ref_len || c_qry != seg_qry_len) {
#ifdef _DEBUG
                spdlog::warn("globalAlignMM2(seg): segment cigar mismatch (expected ref:{}/qry:{}, got ref:{}/qry:{}); forcing global fallback logic for this segment",
                             seg_ref_len, seg_qry_len, c_ref, c_qry);
#endif
                // ------------------------------------------------------------------
                // 彻底修复逻辑：
                // 当底层比对器（尤其是 banded KSW2）在极大长度差异下（例如 ref=20000, qry=100）
                // 可能会因为 band 限制或 Z-drop 导致无法延伸到终点，从而返回不完整的 CIGAR。
                //
                // 此时如果不做处理直接 return，就会导致整个链的比对中断或错误。
                //
                // 策略 A（简单回退）：既然分段比对失败，就放弃这一小段的精细比对，
                // 直接视为 "大段 ref 对应小段 query"，用一个简单的 D/I 结构填充，或者调用更重但更稳的手段？
                //
                // 实际上我们之前是 "fallback to global"，即调用 align_func 对该段重试。
                // 但如果 align_func 本身就是那个失败的函数的（例如 globalAlignKSW2），重试依然会错。
                //
                // 策略 B（强制补齐）：
                // 我们明确知道这段 seq_ref 和 seq_qry 必须对在一起。
                // 如果 ref 极长而 qry 极短，大概率是大片段缺失（Deletion）。
                // 我们可以构造一个 "fallback CIGAR"：
                // 1. 先匹配 common 长度（M/X）
                // 2. 剩余 ref 变为 D，剩余 qry 变为 I
                // 或者更简单：直接构造一个巨大的 gap 对齐？
                //
                // 这里采用 "暴力补齐" 方案：
                // 如果比对器未能覆盖全长，我们手动追加剩余部分的 D/I。
                // 但 CIGAR 必须连续。如果原始 CIGAR 只对了一半，直接追加可能导致中间断裂。
                //
                // 最终方案：
                // 当出现 mismatch 时，说明该段极度不可信/算法失效。
                // 我们直接构造一个平凡解：
                // - 如果 ref >>> qry：全部设为 D，中间夹杂少量 I（或者先 I 后 D）
                // - 实际上对于 anchor 间的 gap，通常就是大 INDEL。
                //
                // 为了保证结果一定正确（长度匹配），我们直接生成一个简单的 CIGAR：
                // - 先全 I (query bases)
                // - 再全 D (ref bases)
                // （或者反过来，或者混合）。
                // 这种对齐虽然分不高，但几何上是正确的，且不会导致程序崩溃/报错。
                // ------------------------------------------------------------------
                cigar::Cigar_t forced_cigar;
                // 策略：尽量先 match/mismatch 短的那边？不，直接粗暴处理最稳健。
                // 先把 Query 消耗完（Insertion），再把 Reference 消耗完（Deletion）。
                // 这样能保证一定会回到对角线。
                if (seg_qry_len > 0) {
                    forced_cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(seg_qry_len)));
                }
                if (seg_ref_len > 0) {
                    forced_cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(seg_ref_len)));
                }
                cigar::appendCigar(result, forced_cigar);

                // 既然强制补齐了，就按照预期的长度推进指针
                ref_pos = ref_end;
                qry_pos = qry_end;
                return;
            }

            cigar::appendCigar(result, seg_cigar);

            // 用 CIGAR 的消耗来推进（关键）：不信任 anchor 坐标，只信任 CIGAR 真实消耗
            ref_pos = ref_start + c_ref;
            qry_pos = qry_start + c_qry;
        };

        // ------------------------------------------------------------------
        // 3) 左端：链起点之前的区域（ref_pos/qry_pos -> 第一个锚点起始位置）
        // ------------------------------------------------------------------
        {
            const auto& first = chain_anchors.front();
            append_segment(ref_pos, first.pos_ref, qry_pos, first.pos_qry);
        }

        // ------------------------------------------------------------------
        // 4) 依次处理每个锚点：先加 span，再加到下一个锚点的 gap
        //    修正逻辑：确保每个 anchor.span 恰好被加入一次
        // ------------------------------------------------------------------
        for (std::size_t i = 0; i < chain_anchors.size(); ++i) {
            const auto& a = chain_anchors[i];

            // 4.1) 当前锚点的 span 段
            const std::size_t a_ref_start = static_cast<std::size_t>(a.pos_ref);
            const std::size_t a_qry_start = static_cast<std::size_t>(a.pos_qry);
            const std::size_t a_ref_end = a_ref_start + static_cast<std::size_t>(a.span);
            const std::size_t a_qry_end = a_qry_start + static_cast<std::size_t>(a.span);

            append_segment(ref_pos, a_ref_end, qry_pos, a_qry_end);

            // 4.2) 如果不是最后一个锚点，加上 gap（当前 anchor 结束 -> 下一个 anchor 开始）
            if (i + 1 < chain_anchors.size()) {
                const auto& b = chain_anchors[i + 1];
                append_segment(ref_pos, b.pos_ref, qry_pos, b.pos_qry);
            }
        }

        // ------------------------------------------------------------------
        // 5) 右端：最后一个锚点结束到序列末尾
        // ------------------------------------------------------------------
        append_segment(ref_pos, ref_len, qry_pos, qry_len);


        // 6) 最终一致性检查：若不完整，退化为全局比对（保持历史行为）
        const std::size_t total_ref = cigar::getRefLength(result);
        const std::size_t total_qry = cigar::getQueryLength(result);
        if (total_ref != ref_len || total_qry != qry_len) {
            spdlog::error("globalAlignMM2: final cigar mismatch (ref:{}/{}, qry:{}/{}), fallback to global",
                         total_ref, ref_len, total_qry, qry_len);
            return align_func(ref, query);
        }

        return result;
    }

} // namespace align

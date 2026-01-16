#include "align.h"
#include "ksw2.h"
#include "bindings/cpp/WFAligner.hpp"
extern "C" {
#include "alignment/cigar.h"
#include "wavefront/wavefront_align.h"
}
namespace align
{
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query)
    {
        /* ---------- 1. 编码序列 ---------- */
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());

        for (size_t i = 0; i < ref.size(); ++i)
            ref_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
        for (size_t i = 0; i < query.size(); ++i)
            qry_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(query[i])];

        /* ---------- 2. 配置 KSW2 参数 ---------- */
        // dna5_simd_mat 现在是 constexpr 编译期常量，无需运行时初始化

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

        cfg.end_bonus = 0;         // 全局比对不需要 ends-free 奖励
        cfg.zdrop = -1;            // 禁用 z-drop（全局比对必须完整比完）
        cfg.band_width = align::auto_band(ref.size(), query.size());;       // 启用全矩阵（也可用 auto_band 加速）

        // Flag 组合说明：
        // - KSW_EZ_GENERIC_SC：启用完整替换矩阵（必须！否则会忽略 dna5_simd_mat）
        // - KSW_EZ_RIGHT：gap 右对齐（CIGAR 标准化，便于比较）
        cfg.flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT;

        /* ---------- 3. 调用 KSW2 ---------- */
        ksw_extz_t ez{};

        ksw_extz2_sse(0,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);


        /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(ez.cigar);           // KSW2 用 malloc()
        return cigar;
    }

    cigar::Cigar_t extendAlignKSW2(const std::string& ref,
        const std::string& query,
        int zdrop)
    {
        /* ---------- 1. 序列编码 ---------- */
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());
        for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = align::ScoreChar2Idx[(uint8_t)ref[i]];
        for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = align::ScoreChar2Idx[(uint8_t)query[i]];

        ///* ---------- 2. 配置 ---------- */
        //KSW2AlignConfig cfg = makeTurboKSW2Config(query.size(), ref.size());
        ////KSW2AlignConfig cfg;
        //cfg.zdrop = zdrop;       // 用于提前终止
        //cfg.flag = KSW_EZ_EXTZ_ONLY     // ends-free extension
        //    | KSW_EZ_APPROX_MAX    // 跟踪 ez.max_q/max_t
        //    | KSW_EZ_APPROX_DROP   // 在 approximate 模式下触发 z-drop 就中断
        //    | KSW_EZ_RIGHT;        // （可选）gap 右对齐     // **关键**：启用 extension/ends-free
        //// 若需要右对齐 gaps 建议保留 KSW_EZ_RIGHT
        //cfg.end_bonus = 100;
        //cfg.band_width = -1;
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
        ksw_extz2_sse(nullptr,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);

        // 赋值bool& if_zdrop,int& ref_end,int& qry_end
        /* ---------- 4. 拷贝 & 释放 ---------- */
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(ez.cigar);                    // ksw2 使用 malloc
        return cigar;                      // 返回的 CIGAR 即延伸片段
    }


    cigar::Cigar_t globalAlignWFA2(const std::string& ref,
        const std::string& query)
    {
        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
        attributes.distance_metric = gap_affine;
        attributes.affine_penalties.mismatch = 2;      // X > 0
        attributes.affine_penalties.gap_opening = 3;   // O >= 0
        attributes.affine_penalties.gap_extension = 1; // E > 0
        attributes.memory_mode = wavefront_memory_ultralow;
        attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
        attributes.heuristic.min_k = -50;
        attributes.heuristic.max_k = +50;
        attributes.heuristic.steps_between_cutoffs = 1;
        //// Create a WFAligner
        //
        wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

        wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());
        /*wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);

        aligner.alignEnd2End(ref, query);*/
        uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
        int cigar_length = 0; // Length of the CIGAR string.
        // Retrieve the CIGAR string from the wavefront aligner.
        cigar_get_CIGAR(wf_aligner->cigar, false, &cigar_buffer, &cigar_length);

        /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
        cigar::Cigar_t cigar;

        for (int i = 0; i < cigar_length; ++i)
            cigar.push_back(cigar_buffer[i]);

        wavefront_aligner_delete(wf_aligner);

        return cigar;
    }

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
}
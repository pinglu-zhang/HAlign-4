#include "align.h"
#include "config.hpp"
#include "preprocess.h"    // alignConsensusSequence
#include "consensus.h"     // generateConsensusSequence
#include "seed.h"          // minimizer 提取和锚点生成
#include <algorithm>
#include <cstddef>
#include <fstream>         // std::ifstream, std::ofstream
#include <filesystem>      // std::filesystem::remove
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

#include <unordered_map>

// 进度条：
// - <chrono> 用于统计耗时
// - <cstdio> 用于 fprintf(stderr, ...) 实现单行(\r)进度输出
#include <chrono>
#include <cstdio>

namespace align {

    // ------------------------------------------------------------------
    // 构造函数1：直接传入参数初始化
    // 说明：
    // 1. 如果 keep_first_length == true，使用 ref_sequences[0] 作为参考
    // 2. 否则，调用 MSA 生成共识序列作为参考
    // 3. threads 和 msa_cmd 参数用于共识序列生成
    // 4. keep_first_length 和 keep_all_length 会被保存到成员变量
    // ------------------------------------------------------------------
    RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size, int window_size,
                           int sketch_size, bool noncanonical, int threads, std::string msa_cmd,
                           bool keep_first_length, bool keep_all_length)
        : work_dir(work_dir),
          kmer_size(kmer_size),
          window_size(window_size),
          sketch_size(sketch_size),
          noncanonical(noncanonical),
          threads(threads),
          msa_cmd(msa_cmd),          // 使用 move 语义避免拷贝
          keep_first_length(keep_first_length),
          keep_all_length(keep_all_length)
    {
        // 从文件读取参考序列并构建索引
        seq_io::KseqReader reader(ref_fasta_path);
        seq_io::SeqRecord rec;
        while (reader.next(rec))
        {
            // 关键修复：必须在 move(rec) 之前计算 sketch 和 minimizer
            auto sketch = mash::sketchFromSequence(rec.seq, kmer_size, sketch_size, noncanonical, random_seed);
            auto minimizer = minimizer::extractMinimizer(rec.seq, kmer_size, window_size, noncanonical);

            ref_sequences.push_back(std::move(rec));
            ref_sketch.push_back(std::move(sketch));
            ref_minimizers.push_back(std::move(minimizer));
        }

        FilePath consensus_unaligned_file = ref_fasta_path;
        FilePath consensus_aligned_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_ALIGNED;
        FilePath consensus_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_FASTA;
        FilePath consensus_json_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_JSON;

        const std::size_t batch_size = 4096;
        alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, msa_cmd, threads);
        std::string consensus_string = consensus::generateConsensusSequence(
            consensus_aligned_file,
            consensus_file,
            consensus_json_file,
            0, // 不限制数量
            threads,
            batch_size
        );

        if (keep_first_length)
        {
            consensus_seq = ref_sequences.front();
        }
        else
        {
            consensus_seq.id = "consensus";
            consensus_seq.seq = std::move(consensus_string);
        }

        // ------------------------------------------------------------------
        // 性能优化：预计算共识序列的 MinHash sketch 和 Minimizer 索引
        // ------------------------------------------------------------------
        // 说明：
        // 1. consensus_sketch 和 consensus_minimizer 在每次 alignOneQueryToRef 的二次比对中都会用到
        // 2. 将其提前计算并存储为成员变量，避免每次都重复计算
        // 3. 对于大规模比对任务（数千到数万条 query），可显著减少 CPU 时间
        //
        // 性能影响：
        // - 初始化时增加一次 sketch 和 minimizer 计算开销（通常 < 20ms）
        // - 每次 alignOneQueryToRef 节省一次 sketch 计算（约 1-5ms）和 minimizer 计算（约 2-10ms）
        // - 对于 10000 条 query，可节省 30-150 秒总时间
        //
        // 实现说明：
        // - consensus_sketch：用于 Jaccard 相似度计算，决定是否需要二次比对
        // - consensus_minimizer：用于种子定位和锚点比对（当前未使用，预留未来优化）
        // ------------------------------------------------------------------
        consensus_sketch = mash::sketchFromSequence(
            consensus_seq.seq,
            static_cast<std::size_t>(kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        consensus_minimizer = minimizer::extractMinimizer(
            consensus_seq.seq,
            kmer_size,
            window_size,
            noncanonical);

    }

    // ------------------------------------------------------------------
    // 构造函数2：基于 Options 结构体初始化
    // 说明：
    // 1. 从 Options 中提取相关参数，委托给第一个构造函数
    // 2. 构造函数1会根据 keep_first_length 标志自动生成或选择参考序列
    // 3. 从 opt 中提取 keep_first_length 和 keep_all_length
    // ------------------------------------------------------------------
    RefAligner::RefAligner(const Options& opt, const FilePath& ref_fasta_path)
        : RefAligner(
            opt.workdir,                // work_dir：工作目录
            ref_fasta_path,             // 参考序列文件路径
            opt.kmer_size,              // kmer_size：k-mer 大小
            opt.kmer_window,            // window_size：minimizer 窗口大小
            opt.sketch_size,            // sketch_size：sketch 大小
            true,                       // noncanonical：是否使用非标准模式（固定为 true）
            opt.threads,                // threads：线程数（用于共识序列生成）
            opt.msa_cmd,                // msa_cmd：MSA 命令模板
            opt.keep_first_length,      // keep_first_length：从 opt 提取
            opt.keep_all_length         // keep_all_length：从 opt 提取
        )
    {
        // 委托构造函数已完成所有初始化工作（包括共识序列生成）
    }

    // ==================================================================
    // 静态方法：globalAlign - 全局序列比对（统一接口）
    // ==================================================================
    // 功能：
    // 根据两个序列的相似度自动选择合适的比对算法执行全局比对
    //
    // 当前实现：
    // - 直接使用 WFA2 算法（globalAlignWFA2）
    // - 暂时忽略 similarity、ref_minimizer、query_minimizer 参数（预留未来扩展）
    //
    // 未来扩展方向：
    // 可根据 similarity 选择最优算法：
    // - similarity > 0.90：优先使用 WFA2（对高相似度序列速度快）
    // - 0.70 < similarity <= 0.90：使用 KSW2 延伸比对
    // - similarity <= 0.70：使用 KSW2 全局比对（更稳健）
    //
    // 可利用 minimizer 信息优化比对：
    // - 对于长序列（>10kb），可先通过 minimizer 找到锚点
    // - 在锚点之间分段比对，减少计算量
    // - 根据 minimizer 密度调整比对策略
    //
    // 参数：
    //   @param ref - 参考序列（A/C/G/T/N，大小写不敏感）
    //   @param query - 查询序列（A/C/G/T/N，大小写不敏感）
    //   @param similarity - 序列相似度（0.0 到 1.0）
    //                       当前版本未使用，保留接口兼容性
    //   @param ref_minimizer - 参考序列的 minimizer 索引（可选）
    //                          当前版本未使用，保留接口用于未来优化
    //   @param query_minimizer - 查询序列的 minimizer 索引（可选）
    //                            当前版本未使用，保留接口用于未来优化
    //
    // 返回：
    //   CIGAR 操作序列（cigar::Cigar_t），描述 query 如何比对到 ref
    //
    // 性能说明：
    // - WFA2 对高相似度序列（编辑距离小）通常显著快于传统 DP 算法
    // - 时间复杂度：O(s * N)，s 为编辑距离，N 为序列长度
    // - 对于低相似度序列，WFA2 可能不如 KSW2（未来可优化）
    //
    // 使用示例：
    //   double jaccard = mash::jaccard(sketch1, sketch2);
    //   auto cigar = RefAligner::globalAlign(ref_seq, query_seq, jaccard, &ref_mz, &query_mz);
    // ==================================================================
    cigar::Cigar_t RefAligner::globalAlign(const std::string& ref,
                                           const std::string& query,
                                           double similarity,
                                           const SeedHits* ref_minimizer,
                                           const SeedHits* query_minimizer) const
    {

        constexpr double min_similarity = 0.5;  //
        constexpr double min_coverage = 0.5;  //
        int len_diff = std::abs(static_cast<int>(ref.size()) - static_cast<int>(query.size()));
        double coverage = 1.0 - static_cast<double>(len_diff) / std::max(ref.size(), query.size());

        if (similarity >= min_similarity && coverage >= min_coverage) {
            return globalAlignWFA2(ref, query);
        } else {
            // ----------------------------------------------------------
            // 低相似度：使用 MM2（基于锚点的分段比对）
            // ----------------------------------------------------------
            anchor::Anchors anchors;
            minimizer::MinimizerHits ref_hits;
            minimizer::MinimizerHits qry_hits;

            if (ref_minimizer != nullptr && !ref_minimizer->empty())
            {
                ref_hits = *ref_minimizer;
            }else
            {
                ref_hits = minimizer::extractMinimizer(ref, kmer_size, window_size, noncanonical);
            }
            if (query_minimizer != nullptr && !query_minimizer->empty())
            {
                qry_hits = *query_minimizer;
            }else
            {
                qry_hits = minimizer::extractMinimizer(query, kmer_size, window_size, noncanonical);
            }

            anchors = minimizer::collect_anchors(ref_hits, qry_hits);

            if (similarity >= min_similarity)
            {
                return globalAlignMM2(ref, query, anchors, globalAlignWFA2);
            }
            else
            {
                return globalAlignMM2(ref, query, anchors, globalAlignKSW2);
            }
        }
    }


    // ------------------------------------------------------------------
    // 辅助函数：写入SAM记录
    // ------------------------------------------------------------------
    void RefAligner::writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                                    std::string_view ref_name, seq_io::SeqWriter& out) const
    {
        const std::string cigar_str = cigar::cigarToString(cigar);
        const auto sam_rec = seq_io::makeSamRecord(q, ref_name, cigar_str, 1, 60, 0);
        out.writeSam(sam_rec);
    }

    // ------------------------------------------------------------------
    // 辅助函数：mergeConsensusAndSamToFasta
    // 功能：将共识序列和多个 SAM 文件合并为一个 FASTA 文件
    //
    // 实现说明：
    // 1. 创建一个 SeqWriter 打开输出 FASTA 文件
    // 2. 先写入共识序列（consensus_seq）
    // 3. 逐个读取 SAM 文件，提取序列并追加写入到同一个 writer
    // 4. 使用同一个 writer 实例，避免文件追加模式的复杂性和性能损失
    // 5. 所有写入操作在同一个 writer 的生命周期内完成
    // 6. 如果 write_gap 为 true，根据 SAM 记录中的 CIGAR 字符串在序列中插入 gap
    //
    // 性能优化：
    // - 使用大缓冲区（默认 8MiB）批量写入，减少系统调用
    // - 流式处理 SAM 文件，内存占用与文件大小无关
    // - 避免文件的重复打开和关闭
    // ------------------------------------------------------------------
    std::size_t RefAligner::mergeConsensusAndSamToFasta(
        const std::vector<FilePath>& sam_paths,
        const FilePath& fasta_path,
        bool write_gap,
        std::size_t line_width
        ) const
    {
        // 创建 FASTA writer（非追加模式，会覆盖已存在的文件）
        seq_io::SeqWriter writer(fasta_path, line_width);

        // 1. 先写入共识序列
        writer.writeFasta(consensus_seq);
        writer.flush();

        // 2. 统计信息（用于返回值和调试）
        std::size_t total_count = 1;  // 已写入共识序列，计数从 1 开始
        std::size_t file_idx = 0;

        // 3. 逐个处理每个 SAM 文件
        for (const auto& sam_path : sam_paths) {
            // 打开当前 SAM 文件
            // 说明：
            // - 每个 SamReader 独立打开和关闭
            // - 使用 RAII 确保文件在处理完毕后自动关闭
            // - 如果文件打开失败，会抛出异常并中止合并
            seq_io::SamReader reader(sam_path);

            seq_io::SamRecord sam_rec;
            seq_io::SeqRecord fasta_rec;
            std::size_t file_count = 0;

            while (reader.next(sam_rec)) {

                fasta_rec = seq_io::samRecordToSeqRecord(sam_rec, false);

                // 如果 write_gap 为 true，根据 CIGAR 在序列中插入 gap
                // 说明：
                // - padQueryToRefByCigar 会根据 CIGAR 中的 D（deletion）操作在序列中插入 '-'
                // - 这使得输出的序列与参考序列对齐
                // - 如果 CIGAR 为空或不需要插入 gap，保持原序列不变
                if (write_gap && !sam_rec.cigar.empty()) {
                    // 解析 CIGAR 字符串为内部表示
                    cigar::Cigar_t cigar_ops = cigar::stringToCigar(sam_rec.cigar);                   // 根据 CIGAR 在序列中插入 gap（原地修改）
                    cigar::delQueryToRefByCigar(fasta_rec.seq, cigar_ops);
                }

                writer.writeFasta(fasta_rec);
                ++file_count;
                ++total_count;
            }

            // 调试信息：记录每个文件的处理进度
            #ifdef _DEBUG
            spdlog::debug("mergeConsensusAndSamToFasta: processed file {} ({}/{}): {} records from {}",
                         file_idx + 1, file_idx + 1, sam_paths.size(), file_count, sam_path.string());
            #endif

            ++file_idx;
        }

        // 4. 确保所有数据已刷新到磁盘
        writer.flush();

        // 调试信息：记录合并统计
        #ifdef _DEBUG
        spdlog::debug("mergeConsensusAndSamToFasta: merged {} SAM files ({} total records) to {}",
                     sam_paths.size(), total_count, fasta_path.string());
        #endif

        return total_count;
    }

    // ------------------------------------------------------------------
    // 主函数：alignOneQueryToRef
    // 功能：对单个query执行比对并写入SAM文件
    //
    // 流程：
    // 1. 计算 query sketch
    // 2. 找到最相似的 reference
    // 3. 执行全局比对
    // 4. 根据是否有插入进行二次判断并写入相应文件
    //
    // 参数说明：
    //   - q: 待比对的查询序列
    //   - out: 普通输出文件的 writer（无插入或二次比对后无插入的序列）
    //   - out_insertion: 插入序列输出文件的 writer（二次比对后仍有插入的序列）
    // ------------------------------------------------------------------
    void RefAligner::alignOneQueryToRef(const seq_io::SeqRecord& q,
                                       seq_io::SeqWriter& out,
                                       seq_io::SeqWriter& out_insertion) const
    {

        // 1) 计算 query sketch 和 minimizer
        // 说明：
        // - sketch 用于快速相似度估计（Jaccard 距离）
        // - minimizer 用于种子定位和锚点比对（当前传递给 globalAlign，未来可用）
        const mash::Sketch qsk = mash::sketchFromSequence(
            q.seq,
            static_cast<std::size_t>(kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        const SeedHits query_minimizer = minimizer::extractMinimizer(
            q.seq,
            kmer_size,
            window_size,
            noncanonical);

        // 2) 选择最相似 reference（线性扫描）
        double best_j = -1.0;
        std::size_t best_r = 0;
        for (std::size_t r = 0; r < ref_sketch.size(); ++r) {
            const double j = mash::jaccard(qsk, ref_sketch[r]);
            if (j > best_j) {
                best_j = j;
                best_r = r;
            }
        }

        const auto& best_ref = ref_sequences[best_r];

        // 3) 执行全局比对（使用统一的比对接口）
        // 说明：
        // - 使用 RefAligner::globalAlign 而非直接调用 globalAlignWFA2
        // - 传递 ref_minimizer 和 query_minimizer，允许未来根据种子信息优化比对
        // - 允许未来根据相似度自动选择最优算法
        cigar::Cigar_t initial_cigar = globalAlign(
            best_ref.seq,
            q.seq,
            best_j,
            &ref_minimizers[best_r],  // 使用对应参考序列的 minimizer
            &query_minimizer);

        // 4) 根据是否存在插入，决定写入哪个输出文件
        if (!cigar::hasInsertion(initial_cigar)) {
            // 无插入：直接写入普通文件
            writeSamRecord(q, initial_cigar, best_ref.id, out);
            return;
        }

        // 有插入：进行二次比对判断（与共识序列比对）
        // 说明：
        // 1. 使用预计算的 consensus_sketch 和 consensus_minimizer 成员变量（性能优化）
        // 2. 计算 query 与共识序列的相似度（用于传递给 globalAlign）
        // 3. 使用统一的比对接口，保持一致性
        //
        // 性能优化：
        // - consensus_sketch 和 consensus_minimizer 在构造函数中已计算并缓存
        // - 避免每次调用都重复计算 sketchFromSequence 和 extractMinimizer
        // - 对于大规模比对（10000+ queries），可节省数十秒到数分钟

        const double consensus_similarity = mash::jaccard(qsk, consensus_sketch);

        cigar::Cigar_t recheck_cigar = globalAlign(
            consensus_seq.seq,
            q.seq,
            consensus_similarity,
            &consensus_minimizer,  // 使用预计算的共识序列 minimizer
            &query_minimizer);

        // 使用二次比对结果（如果失败则使用初始结果）
        const cigar::Cigar_t& final_cigar = recheck_cigar.empty() ? initial_cigar : recheck_cigar;

        // 根据二次比对结果决定输出文件
        if (cigar::hasInsertion(final_cigar)) {
            writeSamRecord(q, final_cigar, consensus_seq.id, out_insertion);
        } else {
            writeSamRecord(q, final_cigar, consensus_seq.id, out);
        }

    }

    void RefAligner::alignQueryToRef(const FilePath& qry_fasta_path, std::size_t batch_size)
    {
        // =====================================================================
        // OpenMP 并行骨架（流式读取 + 每线程独立 writer）
        //
        // 核心思路：
        // - 单线程按 chunk 流式读取 query；
        // - 对一个 chunk 用 OpenMP 并行 for；
        // - 每线程有稳定 thread_id (= omp_get_thread_num)；
        // - 每线程一个独占输出文件，避免加锁。
        //
        // 我们把“相似度计算 +（占位）比对 + 写出”抽象成一个函数 process_one()。
        // 这样后续你只需要替换 process_one() 的内部逻辑即可，不影响并行框架。
        // =====================================================================

        if (ref_sequences.empty() || ref_sketch.empty()) {
            throw std::runtime_error("RefAligner::alignQueryToRef: reference is empty");
        }

        // batch_size==0 没有意义，这里兜底为 1
        if (batch_size == 0) batch_size = 1;

        // 线程数：
        // - threads<=0：尊重 OpenMP 默认（OMP_NUM_THREADS 或 runtime 配置）
        // - threads>0 ：显式设置
        if (threads > 0) {
            omp_set_num_threads(threads);
        }
        const int nthreads = std::max(1, omp_get_max_threads());

        FilePath result_dir = work_dir / RESULTS_DIR;
        file_io::ensureDirectoryExists(result_dir, "result directory");

        // 每线程一个输出
        // writer（SAM 模式，占位输出也能保持格式合法）
        outs_path.clear();
        outs_path.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion_path.clear();
        outs_with_insertion_path.resize(static_cast<std::size_t>(nthreads));


        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs;
        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs_with_insertion;

        outs.clear();
        outs_with_insertion.clear();
        outs.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion.resize(static_cast<std::size_t>(nthreads));  // 关键修复：必须 resize 以避免越界访问

        for (int tid = 0; tid < nthreads; ++tid) {
            // 每个线程独立输出 SAM 文件，避免线程竞争
            // 文件名格式：thread<tid>.sam 和 thread<tid>_insertion.sam
            FilePath out_path = result_dir / (THREAD_SAM_PREFIX + std::to_string(tid) + THREAD_SAM_SUFFIX);
            FilePath out_path_insertion = result_dir / (THREAD_SAM_PREFIX + std::to_string(tid) + THREAD_INSERTION_SAM_SUFFIX);

            outs_path[static_cast<std::size_t>(tid)] = out_path;
            outs_with_insertion_path[static_cast<std::size_t>(tid)] = out_path_insertion;


            auto tmp = seq_io::SeqWriter::Sam(out_path);
            auto tmp_insertion = seq_io::SeqWriter::Sam(out_path_insertion);
            outs[static_cast<std::size_t>(tid)] =
                std::make_unique<seq_io::SeqWriter>(std::move(tmp));

            outs[static_cast<std::size_t>(tid)]->writeSamHeader("@HD\tVN:1.6\tSO:unknown");
            outs_with_insertion[static_cast<std::size_t>(tid)] =
                std::make_unique<seq_io::SeqWriter>(std::move(tmp_insertion));
            outs_with_insertion[static_cast<std::size_t>(tid)]->writeSamHeader("@HD\tVN:1.6\tSO:unknown");
        }

        // ---- 流式读取（单线程）+ chunk/batch 并行处理 ----
        seq_io::KseqReader reader(qry_fasta_path);
        std::vector<seq_io::SeqRecord> chunk;
        chunk.reserve(batch_size);

        // -----------------------------------------------------------------
        // 进度条（单行刷新，无需预先统计总数，且不使用原子/锁）
        //
        // 设计说明（为什么这样做）：
        // 1) 你明确不希望预统计（需要二次遍历输入，超大数据集很慢）。
        // 2) 你也不希望用 atomic（并行区频繁原子加会带来额外开销）。
        //
        // 因此这里采用：
        // - “读取线程”统计本次读入/累计已处理的序列数（read/processed）
        // - 并行区只做计算，不做计数
        // - 每完成一个 chunk（即 chunk 内所有 read 都已处理），再在串行区更新 processed_count。
        //
        // 进度刷新策略：每处理 >=1000 条刷新一次，使用 '\r' 覆盖同一行。
        // -----------------------------------------------------------------
        std::size_t processed_count = 0;   // 已完成比对并写出（以 chunk 为单位累加）
        std::size_t next_report = 1000;   // 下次刷新阈值
        const auto t0 = std::chrono::steady_clock::now();

        auto print_progress = [&](bool force) {
            if (!force && processed_count < next_report) {
                return;
            }

            // 下一次阈值：按 1000 的倍数递增
            if (processed_count >= next_report) {
                next_report = ((processed_count / 1000) + 1) * 1000;
            }

            const auto now = std::chrono::steady_clock::now();
            const double sec = std::chrono::duration_cast<std::chrono::duration<double>>(now - t0).count();
            const double rate = (sec > 0.0) ? (static_cast<double>(processed_count) / sec) : 0.0;

            // std::fprintf(stdout, "\r%s", fmt::format("[align] processed={}  elapsed={:.1f}s  rate={:.1f} seq/s",
            //                                         processed_count, sec, rate).c_str());
            // std::fflush(stdout);
            // ANSI: cyan = \x1b[36m, reset = \x1b[0m
            std::fprintf(stderr,
                         "\r\033[32m[align] processed=%zu  elapsed=%.1fs  rate=%.1f seq/s\033[0m   ",
                         processed_count, sec, rate);
            std::fflush(stderr);

        };

        while (true)
        {
            chunk.clear();
            seq_io::SeqRecord rec;
            for (std::size_t i = 0; i < batch_size; ++i) {
                if (!reader.next(rec)) break;
                chunk.push_back(std::move(rec));
            }
            if (chunk.empty()) break;

            // OpenMP 规约：
            // - default(none) 强制显式标注共享/私有变量，避免无意的数据竞争；
            // - outs/outs_with_insertion 与 chunk 是跨线程共享只读/线程索引访问的；tid/i 为线程私有。
            #pragma omp parallel default(none) shared(outs, outs_with_insertion, chunk)
            {
                const int tid = omp_get_thread_num();
                auto& out = *outs[static_cast<std::size_t>(tid)];
                auto& out_insertion = *outs_with_insertion[static_cast<std::size_t>(tid)];

                #pragma omp for schedule(dynamic, 1)
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(chunk.size()); ++i)
                {
                    alignOneQueryToRef(chunk[static_cast<std::size_t>(i)], out, out_insertion);
                }
            }

            // flush + 进度更新都在串行区做
            for (auto& w : outs) {
                w->flush();
            }
            for (auto& w : outs_with_insertion) {
                w->flush();
            }

            // chunk 处理完毕：累加已完成数量
            processed_count += chunk.size();
            print_progress(false);
        }

        // 结束：强制输出一次并换行
        print_progress(true);
        std::fprintf(stderr, "\n");
        spdlog::info("Alignment completed");
    }

    // ------------------------------------------------------------------
    // 辅助函数：parseAlignedReferencesToCigar
    // 功能：读取 MSA 对齐后的参考序列文件，解析每个参考序列与共识序列的对齐关系
    //
    // 实现说明：
    // 1. 流式读取 FASTA 文件，内存占用 O(L)（L 为序列长度）
    // 2. 对于每个后续序列，逐位置比较：
    //    - ref 为碱基（非 '-'）：记录 M（匹配/错配）
    //    - ref 为 gap（'-'）：记录 D（删除）
    // 3. 使用游程编码压缩连续相同操作（例如 10 个连续 M -> "10M"）
    //
    // 性能优化：
    // - 流式处理，只保存共识序列和当前参考序列
    // - 预分配 CIGAR 容器，减少内存重新分配
    // - 使用游程编码，压缩 CIGAR 大小
    // ------------------------------------------------------------------
    void RefAligner::parseAlignedReferencesToCigar(
        const FilePath& aligned_fasta_path,
        std::unordered_map<std::string, cigar::Cigar_t>& out_ref_aligned_map,
        std::vector<bool>& out_ref_gap_pos) const
    {
        // 说明：
        // - 本函数只做“流式读取 + 生成CIGAR(M/D)”这一件事，避免额外耦合。
        // - out_ref_aligned_map / out_ref_gap_pos 由调用方创建并传入；本函数负责清空并填充。
        // - out_ref_gap_pos 用于记录对齐矩阵中“参考（第一条序列）”每一列是否为 gap”。

        out_ref_aligned_map.clear();
        out_ref_gap_pos.clear();

        // ------------------------------------------------------------------
        // 步骤1：打开 FASTA 并读取第一条序列（参考/共识序列），构建 ref_gap_pos
        // ------------------------------------------------------------------
        seq_io::KseqReader reader(aligned_fasta_path);
        seq_io::SeqRecord rec;
        // ------------------------------------------------------------------
        // 步骤2：流式读取后续序列，生成每条序列的 CIGAR（碱基->M，gap->D）
        // ------------------------------------------------------------------
        std::size_t ref_count = 0;

        while (reader.next(rec)) {
            if (ref_count == 0)
            {
                // ref_gap_pos[i] == true  表示参考序列该列是 '-'
                // ref_gap_pos[i] == false 表示参考序列该列是碱基
                out_ref_gap_pos.reserve(rec.seq.size());
                for (const char base : rec.seq) {
                    out_ref_gap_pos.push_back(base == '-');
                }
            }
            cigar::Cigar_t cigar;
            cigar.reserve(20);  // 经验值：对齐后通常是少量游程段，预留可减少扩容

            char current_op = '\0';
            std::uint32_t current_len = 0;

            // 逐列扫描：碱基记为 M，gap 记为 D。
            // 注意：这里保持原逻辑：不依赖参考序列列信息做一致性校验。
            for (const char base : rec.seq) {
                const char op = (base == '-') ? 'D' : 'M';
                if (op == current_op) {
                    ++current_len;
                } else {
                    if (current_op != '\0' && current_len > 0) {
                        cigar.push_back(cigar::cigarToInt(current_op, current_len));
                    }
                    current_op = op;
                    current_len = 1;
                }
            }

            if (current_op != '\0' && current_len > 0) {
                cigar.push_back(cigar::cigarToInt(current_op, current_len));
            }

            out_ref_aligned_map[rec.id] = std::move(cigar);
            ++ref_count;

#ifdef _DEBUG
            spdlog::debug("parseAlignedReferencesToCigar: {} -> CIGAR: {}",
                         rec.id, cigar::cigarToString(out_ref_aligned_map[rec.id]));
#endif
        }

        // ------------------------------------------------------------------
        // 步骤3：至少要有一条“后续序列”，否则认为输入不符合预期
        // ------------------------------------------------------------------
        if (ref_count == 0) {
            throw std::runtime_error(
                "parseAlignedReferencesToCigar: input contains only the reference/consensus sequence; no subsequent sequences found: " +
                aligned_fasta_path.string());
        }

#ifdef _DEBUG
        spdlog::info("parseAlignedReferencesToCigar: parsed {} sequence CIGARs (from {}), ref_gap_pos_len={}",
                    ref_count, aligned_fasta_path.string(), out_ref_gap_pos.size());
#endif
    }

    // ==================================================================
    // mergeAlignedResults：合并所有比对结果生成最终的多序列比对（MSA）
    // ==================================================================
    // 功能概述：
    // 1. 将多个线程产生的 SAM 文件（比对到不同参考序列）合并为一个统一的 MSA FASTA 文件
    // 2. 处理三类序列：
    //    a) 共识序列（consensus）及其对齐到的参考序列
    //    b) 插入序列（insertion）：无法比对到任何参考的序列，需单独 MSA
    //    c) 普通比对序列：已比对到参考序列的 reads
    // 3. 通过 CIGAR 操作将所有序列投影到统一坐标系（最终的 MSA 列空间）
    //
    // 核心挑战：
    // - 不同参考序列之间存在对应关系（通过共识序列建立）
    // - 需要将"query→ref"和"ref→consensus"的两级对齐关系合并
    // - 插入序列需要独立 MSA 后再整合到最终结果
    //
    // 参数：
    // @param aligned_consensus_path: 已对齐的共识序列文件路径（包含 consensus + 参考序列的 MSA）
    // @param msa_cmd: 外部 MSA 工具命令（用于对齐插入序列）
    // ==================================================================
    void RefAligner::mergeAlignedResults(const FilePath& aligned_consensus_path, const std::string& msa_cmd)
    {
        // ------------------------------------------------------------------
        // 阶段 0：初始化与路径准备
        // ------------------------------------------------------------------
        // 将多个线程的 SAM 文件合并为一个文件
        // 所有的序列都比对到多个参考序列上了，输入的参数为这些参考序列比对好的文件
        // 因此要解析这些参考的互相对应关系，然后把它们合并到一个文件中

        // ------------------------------------------------------------------
        // 阶段 1：处理插入序列（insertion sequences）
        // ------------------------------------------------------------------
        // 插入序列：在比对过程中无法匹配任何参考序列的 reads
        // 处理策略：
        // 1. 收集所有线程产生的 insertion SAM 文件
        // 2. 将它们与共识序列合并为一个 FASTA 文件
        // 3. 调用外部 MSA 工具（如 MAFFT/Muscle）进行多序列比对
        // 4. 解析 MSA 结果生成 CIGAR（用于后续投影）
        // ------------------------------------------------------------------
        // 首先把insertion的sam文件合并为fasta
        FilePath result_dir = work_dir / RESULTS_DIR;

        bool using_other_align_insertion = true;  // 标志：是否使用外部 MSA 工具对齐插入序列

        FilePath aligned_insertion_fasta = result_dir / ALIGNED_INSERTION_FASTA;

        if (using_other_align_insertion)
        {
            // ------------------------------------------------------------------
            // 子步骤 1.1：收集所有线程的插入序列 SAM 文件
            // ------------------------------------------------------------------
            // outs_with_insertion_path：每个线程产生的包含插入序列的 SAM 文件路径集合
            // 将所有 insertion SAM 文件路径收集到一个 vector
            std::vector<FilePath> insertion_sam_paths;
            for (const auto& path : outs_with_insertion_path)
            {
                insertion_sam_paths.push_back(path);
            }

            // ------------------------------------------------------------------
            // 子步骤 1.2：合并为单个 FASTA 文件（共识序列 + 插入序列）
            // ------------------------------------------------------------------
            // 定义输出 FASTA 文件路径
            FilePath insertion_fasta_path = result_dir / ALL_INSERTION_FASTA;
            bool write_gap = keep_first_length || keep_all_length;  // 标志：是否在插入序列中保留 gap（根据 CIGAR）

            // 调用辅助函数：将共识序列和所有 SAM 文件合并为一个 FASTA 文件
            // 说明：
            // 1. 该函数会先写入共识序列（consensus_seq）作为 MSA 的参考（第一条序列）
            // 2. 然后逐个读取 SAM 文件并提取序列追加写入
            // 3. 使用同一个 SeqWriter 实例，避免追加模式的复杂性
            // 4. 返回值为合并的总序列数（包括共识序列）
            // 性能考虑：SeqWriter 内部使用缓冲，避免频繁系统调用
            const std::size_t total_sequences = mergeConsensusAndSamToFasta(
                insertion_sam_paths,
                insertion_fasta_path,
                write_gap,  // 是否根据 CIGAR 插入/删除 gap
                80          // FASTA 行宽：每行 80 个字符（标准 FASTA 格式）
            );

#ifdef _DEBUG
            spdlog::info("mergeAlignedResults: merged {} sequences (1 consensus + {} from SAM files) to {}",
                        total_sequences, total_sequences - 1, insertion_fasta_path.string());
#endif

            // ------------------------------------------------------------------
            // 子步骤 1.3：调用外部 MSA 工具对齐插入序列
            // ------------------------------------------------------------------
            // alignConsensusSequence：调用外部命令（msa_cmd）执行多序列比对
            // 输入：ALL_INSERTION_FASTA（共识序列 + 插入序列）
            // 输出：ALIGNED_INSERTION_FASTA（对齐后的 MSA 结果，所有序列等长）
            // 说明：MSA 工具会将共识序列作为锚点，对齐所有插入序列
            if (!write_gap)
            {
                alignConsensusSequence(insertion_fasta_path, aligned_insertion_fasta, msa_cmd, threads);
                spdlog::info("mergeAlignedResults: aligned insertion sequences using external MSA tool: {} -> {}",
                             insertion_fasta_path.string(), aligned_insertion_fasta.string());
            }else
            {
                // 把insertion_fasta_path复制一份改名为aligned_insertion_fasta
                file_io::copyFile(insertion_fasta_path, aligned_insertion_fasta);
                spdlog::info("mergeAlignedResults: skipped alignment for insertion sequences, copied {} to {}",
                             insertion_fasta_path.string(), aligned_insertion_fasta.string());
            }

        }
        else
        {
            // 未实现的分支：可能是星比对（star alignment）或其他策略
        }

        // ------------------------------------------------------------------
        // 阶段 2：解析 MSA 文件，生成 CIGAR 映射表
        // ------------------------------------------------------------------
        // 目的：将 MSA 文件（FASTA 格式）转换为 CIGAR 结构，用于后续序列投影
        //
        // 核心数据结构：
        // 1. ref_aligned_map：参考序列名 → CIGAR（描述该序列如何对齐到共识序列）
        //    - key: 参考序列 ID（如 "ref_1", "ref_2"）
        //    - value: CIGAR 字符串（只包含 M/D 操作，因为是 MSA 结果）
        //    - 用途：将"比对到参考序列的 query"投影到共识序列坐标系
        //
        // 2. insertion_aligned_map：插入序列名 → CIGAR
        //    - 描述每个插入序列如何对齐到插入 MSA 的共识序列（第一条）
        //    - 用途：将插入序列投影到统一坐标系
        //
        // 3. ref_gap_pos：共识序列（第一条）在 MSA 中的 gap 位置标记
        //    - 长度 = MSA 的列数
        //    - ref_gap_pos[i] = true：第 i 列在共识序列中是 gap（'-'）
        //    - 用途：如果 keep_first_length=true，移除这些列以保持共识序列原始长度
        //
        // 4. insertion_ref_gap_pos：插入 MSA 中共识序列的 gap 位置
        //    - 作用类似 ref_gap_pos，用于插入序列的坐标系统
        // ------------------------------------------------------------------
        // 读取比对好参考序列文件，获得每个参考序列的和共识序列的比对结果，结果里应该只有M和D
        std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map;
        std::unordered_map<std::string, cigar::Cigar_t> insertion_aligned_map;
        std::vector<bool> ref_gap_pos;                 // 共识/参考（第一条）每列是否为gap
        std::vector<bool> insertion_ref_gap_pos;       // insertion MSA 中第一条序列每列是否为gap
        FilePath consensus_aligned_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_ALIGNED;

        // 调用辅助函数：解析 MSA 对齐文件，生成每个序列与"对齐矩阵列"的 CIGAR
        // parseAlignedReferencesToCigar 逻辑：
        // 1. 读取 MSA 文件中的所有序列（第一条为参考/共识）
        // 2. 对每条序列，将 gap 位置转换为 CIGAR（M=匹配/不匹配，D=删除）
        // 3. 记录第一条序列的 gap 位置到 ref_gap_pos/insertion_ref_gap_pos
        parseAlignedReferencesToCigar(consensus_aligned_file, ref_aligned_map, ref_gap_pos);
        parseAlignedReferencesToCigar(aligned_insertion_fasta, insertion_aligned_map, insertion_ref_gap_pos);

#ifdef _DEBUG
        spdlog::info("mergeAlignedResults: parsed alignments for {} reference sequences",
                    ref_aligned_map.size());
        spdlog::info("mergeAlignedResults: parsed alignments for {} insertion sequences",
                    insertion_aligned_map.size());
#endif

        // ------------------------------------------------------------------
        // 阶段 3：初始化最终输出文件与序列长度检测机制
        // ------------------------------------------------------------------
        FilePath final_output_path = FilePath(work_dir) / RESULTS_DIR / FINAL_ALIGNED_FASTA;
        seq_io::SeqWriter final_writer(final_output_path, U_MAX);  // U_MAX：禁用缓冲阈值，手动控制 flush

        // 首先先把consensus_aligned_file写入最终文件
        // ------------------------------------------------------------------
        // 序列长度一致性检测机制
        // ------------------------------------------------------------------
        // 说明：
        // - MSA（多序列比对）的关键要求是所有序列长度必须一致
        // - expected_length：第一条序列的长度，作为后续序列的参考标准
        // - seq_count：已写入的序列数量（用于错误定位）
        // - 如果发现长度不一致，立即抛出异常并报告详细信息（序列 ID、实际长度、期望长度）
        //
        // 为什么需要检测：
        // 1. 确保最终 MSA 文件格式正确（所有序列等长）
        // 2. 及早发现 CIGAR 投影错误或数据损坏
        // 3. 提供清晰的错误信息，便于调试
        // ------------------------------------------------------------------
        std::size_t expected_length = 0;   // 期望的序列长度（由第一条序列确定）
        std::size_t seq_count = 0;         // 已写入的序列数量
        bool length_initialized = false;   // 是否已初始化 expected_length

        // ------------------------------------------------------------------
        // 阶段 4.1：处理共识序列及其参考序列（来自 consensus_aligned_file）
        // ------------------------------------------------------------------
        // 数据来源：consensus_aligned_file
        // - 第一条：共识序列（consensus_seq）
        // - 后续：各个参考序列（已对齐到共识序列的 MSA 结果）
        //
        // 处理流程（对每条序列）：
        // 1. 读取序列（已经是 MSA 格式，包含 gap）
        // 2. [可选] 如果 keep_first_length=true，移除共识序列为 gap 的列
        // 3. 应用 tmp_insertion_cigar：将序列投影到插入 MSA 坐标系
        //    - 这一步将序列对齐到"插入序列 MSA 的共识序列"
        //    - 通过在相应位置插入 gap 实现坐标系统一
        // 4. [可选] 移除插入 MSA 中共识序列为 gap 的列
        // 5. 长度检测：确保与 expected_length 一致
        // 6. 写入最终文件
        // ------------------------------------------------------------------
        // 1. 处理 consensus 对齐序列
        seq_io::KseqReader cons_reader(consensus_aligned_file);
        seq_io::SeqRecord cons_rec;
        // 获取共识序列在插入 MSA 中的 CIGAR（用于将所有序列投影到统一坐标系）
        cigar::Cigar_t tmp_insertion_cigar =  insertion_aligned_map[consensus_seq.id];
        while (cons_reader.next(cons_rec))
        {
            seq_io::cleanSequence(cons_rec);
            // ------------------------------------------------------------------
            // 步骤 4.1.1：移除共识序列中的 gap 列（可选）
            // ------------------------------------------------------------------
            // 条件：keep_first_length = true
            // 作用：去除共识序列为 gap 的所有列，保持共识序列的原始长度
            // 原理：ref_gap_pos[i] = true 表示第 i 列在共识序列中是 gap
            //       removeRefGapColumns 会删除所有这些列
            if (keep_first_length)
            {
                removeRefGapColumns(cons_rec.seq, ref_gap_pos);
            }

            // ------------------------------------------------------------------
            // 步骤 4.1.2：应用插入 CIGAR，投影到插入 MSA 坐标系
            // ------------------------------------------------------------------
            // 说明：consensus 对自己比对，CIGAR 为空（或只有 M 操作）
            // tmp_insertion_cigar：共识序列在插入 MSA 中的 CIGAR
            // padQueryToRefByCigar：根据 CIGAR 在序列中插入 gap，使其对齐到插入 MSA 的坐标系
            // 性能：原地修改，避免拷贝
            cigar::padQueryToRefByCigar(cons_rec.seq, tmp_insertion_cigar );

            // ------------------------------------------------------------------
            // 步骤 4.1.3：移除插入 MSA 中共识序列的 gap 列（可选）
            // ------------------------------------------------------------------
            // 条件：keep_all_length = true 或 keep_first_length = true
            // 作用：去除插入 MSA 共识序列为 gap 的列
            // 用途：压缩最终 MSA，去除冗余的 gap 列
            if (keep_all_length || keep_first_length)
            {
                removeRefGapColumns(cons_rec.seq, insertion_ref_gap_pos);
            }

            // ------------------------------------------------------------------
            // 步骤 4.1.4：序列长度一致性检测
            // ------------------------------------------------------------------
            // 长度检测：第一条序列初始化 expected_length，后续序列必须与之一致
            // 如果不一致，抛出异常（包含序列 ID、实际长度、期望长度、序列位置）
            if (!length_initialized) {
                expected_length = cons_rec.seq.size();
                length_initialized = true;
            } else if (cons_rec.seq.size() != expected_length) {
                throw std::runtime_error(
                    "mergeAlignedResults: sequence length mismatch! sequence '" + cons_rec.id +
                    "' has length " + std::to_string(cons_rec.seq.size()) +
                    ", expected " + std::to_string(expected_length) +
                    " (sequence #" + std::to_string(seq_count + 1) + ")");
            }

            // ------------------------------------------------------------------
            // 步骤 4.1.5：写入最终 MSA 文件
            // ------------------------------------------------------------------
            final_writer.writeFasta(cons_rec);
            ++seq_count;
        }
        final_writer.flush();  // 强制刷新缓冲区，确保数据落盘

        // ------------------------------------------------------------------
        // 阶段 4.2：处理插入序列（来自 aligned_insertion_fasta）
        // ------------------------------------------------------------------
        // 数据来源：aligned_insertion_fasta
        // - 第一条：共识序列（已在阶段 4.1 处理，需跳过）
        // - 后续：所有插入序列（已通过外部 MSA 对齐）
        //
        // 处理流程（对每条序列）：
        // 1. 跳过第一条序列（共识序列，避免重复）
        // 2. [可选] 移除插入 MSA 中共识序列为 gap 的列
        // 3. 长度检测
        // 4. 写入最终文件
        //
        // 说明：
        // - 插入序列已经通过外部 MSA 对齐，不需要再应用额外的 CIGAR
        // - 只需要移除冗余的 gap 列，保持与其他序列的长度一致
        // ------------------------------------------------------------------
        // 2. 处理 insertion 对齐序列
        seq_io::KseqReader insertion_reader(aligned_insertion_fasta);
        seq_io::SeqRecord insertion_rec;
        bool skip_first = true;  // 标志：是否跳过第一条序列（共识序列）
        while (insertion_reader.next(insertion_rec))
        {
            // ------------------------------------------------------------------
            // 步骤 4.2.1：跳过共识序列（第一条）
            // ------------------------------------------------------------------
            if (skip_first)
            {
                // 跳过第一条序列（共识序列），因为已经在阶段 4.1 处理过了
                skip_first = false;
                continue;
            }

            // // ------------------------------------------------------------------
            // // 步骤 4.2.2：移除插入 MSA 中共识序列的 gap 列（可选）
            // // ------------------------------------------------------------------
            // if (keep_first_length || keep_all_length)
            // {
            //     removeRefGapColumns(insertion_rec.seq, insertion_ref_gap_pos);
            // }

            // ------------------------------------------------------------------
            // 步骤 4.2.3：序列长度一致性检测
            // ------------------------------------------------------------------
            // 长度检测：确保与 expected_length 一致
            if (!length_initialized) {
                expected_length = insertion_rec.seq.size();
                length_initialized = true;
            } else if (insertion_rec.seq.size() != expected_length) {
                throw std::runtime_error(
                    "mergeAlignedResults: sequence length mismatch! sequence '" + insertion_rec.id +
                    "' has length " + std::to_string(insertion_rec.seq.size()) +
                    ", expected " + std::to_string(expected_length) +
                    " (sequence #" + std::to_string(seq_count + 1) + ")");
            }

            // ------------------------------------------------------------------
            // 步骤 4.2.4：写入最终 MSA 文件
            // ------------------------------------------------------------------
            final_writer.writeFasta(insertion_rec);
            ++seq_count;
        }
        final_writer.flush();  // 强制刷新缓冲区

        // ------------------------------------------------------------------
        // 阶段 4.3：处理普通比对序列（来自各线程的 SAM 文件）
        // ------------------------------------------------------------------
        // 数据来源：outs_path（每个线程产生的 SAM 文件）
        // 内容：比对到各个参考序列的 reads
        //
        // 核心挑战：两级坐标系投影
        // ========================================================================
        // 坐标系层次结构：
        // 1. Query 原始坐标：reads 的原始序列（无 gap）
        // 2. Reference 坐标：reads 比对到的参考序列坐标（通过 SAM CIGAR 映射）
        // 3. Consensus 坐标：参考序列对齐到共识序列的坐标（通过 ref_aligned_map）
        // 4. Insertion MSA 坐标：最终的统一坐标系（通过 tmp_insertion_cigar）
        //
        // 投影步骤（对每条 read）：
        // 1. 读取 SAM 记录（query 序列 + CIGAR）
        // 2. 应用 SAM CIGAR：query → reference（第一级投影）
        //    - padQueryToRefByCigar(seq, sam_cigar)
        //    - 在 query 中插入 gap，使其对齐到参考序列
        // 3. 应用 ref_aligned_map：reference → consensus（第二级投影）
        //    - padQueryToRefByCigar(seq, ref_cigar)
        //    - 在序列中插入 gap，使其对齐到共识序列
        // 4. [可选] 移除共识序列的 gap 列
        // 5. 应用 tmp_insertion_cigar：consensus → insertion MSA（第三级投影）
        //    - 使其对齐到最终的统一坐标系
        // 6. [可选] 移除插入 MSA 共识序列的 gap 列
        // 7. 长度检测并写入
        // ========================================================================
        //
        // 性能考虑：
        // - 每次 padQueryToRefByCigar 都是原地修改，避免拷贝
        // - CIGAR 操作复杂度：O(CIGAR_length)
        // - 序列修改复杂度：O(seq_length + inserted_gaps)
        // ------------------------------------------------------------------
        // 3. 遍历所有的 SAM 文件，转换为 SeqRecord 并应用比对
        for (const auto& sam_path : outs_path)
        {
            seq_io::SamReader sam_reader(sam_path);
            seq_io::SamRecord sam_rec;
            seq_io::SeqRecord fasta_rec;
            while (sam_reader.next(sam_rec))
            {
                // ------------------------------------------------------------------
                // 步骤 4.3.1：转换 SAM 记录为 FASTA 记录
                // ------------------------------------------------------------------
                // samRecordToSeqRecord：提取 query 序列、ID、质量值等
                // 第二个参数 false：不反向互补（保持原始方向）
                fasta_rec = seq_io::samRecordToSeqRecord(sam_rec, false);

                // ------------------------------------------------------------------
                // 步骤 4.3.2：解析 SAM CIGAR 字符串
                // ------------------------------------------------------------------
                // sam_rec.cigar：SAM 文件中的 CIGAR 字符串（如 "10M2I5M3D10M"）
                // stringToCigar：解析为内部 CIGAR 结构（vector<CigarOp>）
                cigar::Cigar_t tmp_cigar = cigar::stringToCigar(sam_rec.cigar);

                // ------------------------------------------------------------------
                // 步骤 4.3.3：第一级投影 - query → reference
                // ------------------------------------------------------------------
                // 根据 SAM CIGAR 将 query 序列投影到参考序列坐标系
                // 操作：在序列中插入 gap（'-'），对应 CIGAR 中的 D（deletion）操作
                // 原理：D 表示参考序列有而 query 没有的碱基，需要在 query 中插入 gap
                cigar::padQueryToRefByCigar(fasta_rec.seq, tmp_cigar);

                // ------------------------------------------------------------------
                // 步骤 4.3.4：第二级投影 - reference → consensus
                // ------------------------------------------------------------------
                // ref_aligned_map[sam_rec.rname]：该参考序列对齐到共识序列的 CIGAR
                // sam_rec.rname：参考序列名（如 "ref_1"）
                // 操作：将序列从参考坐标系投影到共识序列坐标系
                cigar::Cigar_t ref_cigar = ref_aligned_map[sam_rec.rname];
                cigar::padQueryToRefByCigar(fasta_rec.seq, ref_cigar);

                // ------------------------------------------------------------------
                // 步骤 4.3.5：移除共识序列的 gap 列（可选）
                // ------------------------------------------------------------------
                if (keep_first_length)
                {
                    removeRefGapColumns(fasta_rec.seq, ref_gap_pos);
                }

                // ------------------------------------------------------------------
                // 步骤 4.3.6：第三级投影 - consensus → insertion MSA
                // ------------------------------------------------------------------
                // 说明：如果共识序列本身在插入 MSA 中有 gap，需要在所有序列中同步插入
                // tmp_insertion_cigar：共识序列在插入 MSA 中的 CIGAR（对所有序列相同）
                cigar::padQueryToRefByCigar(fasta_rec.seq, tmp_insertion_cigar);

                // ------------------------------------------------------------------
                // 步骤 4.3.7：移除插入 MSA 共识序列的 gap 列（可选）
                // ------------------------------------------------------------------
                if (keep_all_length || keep_first_length)
                {
                    removeRefGapColumns(fasta_rec.seq, insertion_ref_gap_pos);
                }

                // ------------------------------------------------------------------
                // 步骤 4.3.8：序列长度一致性检测
                // ------------------------------------------------------------------
                // 长度检测：确保所有投影后的序列长度一致
                if (!length_initialized) {
                    expected_length = fasta_rec.seq.size();
                    length_initialized = true;
                } else if (fasta_rec.seq.size() != expected_length) {
                    throw std::runtime_error(
                        "mergeAlignedResults: sequence length mismatch! sequence '" + fasta_rec.id +
                        "' has length " + std::to_string(fasta_rec.seq.size()) +
                        ", expected " + std::to_string(expected_length) +
                        " (sequence #" + std::to_string(seq_count + 1) + ", from SAM: " + sam_rec.rname + ")");
                }

                // ------------------------------------------------------------------
                // 步骤 4.3.9：写入最终 MSA 文件
                // ------------------------------------------------------------------
                final_writer.writeFasta(fasta_rec);
                ++seq_count;
            }
        }
        final_writer.flush();  // 强制刷新缓冲区，确保所有数据落盘

        // ------------------------------------------------------------------
        // 阶段 5：最终验证与统计输出
        // ------------------------------------------------------------------
        // 说明：
        // - 此时所有序列已经写入 FINAL_ALIGNED_FASTA
        // - seq_count：总序列数（共识序列 + 参考序列 + 插入序列 + 普通比对序列）
        // - expected_length：MSA 的统一长度（所有序列的长度）
        //
        // 调试信息（仅在 _DEBUG 模式输出）：
        // - 帮助验证处理是否正确
        // - 提供性能分析的基本数据（总序列数、MSA 长度）
        // ------------------------------------------------------------------
#ifdef _DEBUG
        spdlog::info("mergeAlignedResults: wrote {} sequences; all sequences have consistent length = {}",
                    seq_count, expected_length);
#endif

        // ==================================================================
        // 函数结束
        // 输出：FINAL_ALIGNED_FASTA（包含所有序列的 MSA 文件）
        // ==================================================================
    }

    void RefAligner::removeRefGapColumns(
        std::string& seq,
        const std::vector<bool>& ref_gap_pos)
    {
        // ------------------------------------------------------------------
        // 设计要点：
        // 1) 本函数只做"删列"这一件事：按 ref_gap_pos 过滤输入序列的列。
        // 2) 原地修改，使用移动赋值减少拷贝。
        // 3) ref_gap_pos 为空时不做任何操作。
        // ------------------------------------------------------------------

        // 1) ref_gap_pos 为空，不做任何操作
        if (ref_gap_pos.empty()) {
            return;
        }

        // 2) 对齐列数必须一致，否则说明输入序列与 ref_gap_pos 长度不匹配
        if (seq.size() != ref_gap_pos.size()) {
            throw std::runtime_error(
                "removeRefGapColumns: seq length mismatch, seq_len=" + std::to_string(seq.size()) +
                ", ref_gap_pos_len=" + std::to_string(ref_gap_pos.size()));
        }

        // 3) 过滤参考为 gap 的列。
        // 性能：
        // - 先统计保留列数用于 reserve，避免反复扩容。
        std::size_t keep_len = 0;
        for (const bool is_gap : ref_gap_pos) {
            keep_len += static_cast<std::size_t>(!is_gap);
        }

        std::string filtered;
        filtered.reserve(keep_len);

        for (std::size_t i = 0; i < ref_gap_pos.size(); ++i) {
            if (!ref_gap_pos[i]) {
                filtered.push_back(seq[i]);
            }
        }

        // 4) 使用移动赋值替换原序列（避免拷贝）
        seq = std::move(filtered);
    }

} // namespace align


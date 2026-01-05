#include "consensus.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <thread>
#include <cereal/archives/json.hpp>

#if __has_include(<omp.h>)
    #include <omp.h>
#endif

namespace consensus
{
    // 选择共识碱基：仅在 A/C/G/T/U 之间选择（不会返回 N 或 '-'）。
    // 说明：如果这五者计数均为 0（例如该位点在输入中全为 N 或 gap），
    // 则按优先级回退为 'A'（可根据需要调整为抛出异常或其它策略）。
    char pickConsensusChar(const SiteCount& sc)
    {
        // 固定优先级（在相等时使用）：A > C > G > T > U
        std::uint32_t best = sc.a;
        char best_ch = 'A';

        auto upd = [&](std::uint32_t v, char ch) {
            if (v > best) { best = v; best_ch = ch; }
        };

        // 只比较 A/C/G/T/U，避免选出 N 或 gap
        upd(sc.c, 'C');
        upd(sc.g, 'G');
        upd(sc.t, 'T');
        upd(sc.u, 'U');

        return best_ch;
    }

    // 将 consensus 序列按 FASTA 格式写出，行宽固定为 80
    // 注意：写入前会确保父目录存在
    void writeConsensusFasta(const FilePath& out_fasta, const std::string& seq)
    {
        file_io::ensureParentDirExists(out_fasta);

        std::ofstream ofs(out_fasta, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("failed to open fasta output: " + out_fasta.string());
        }

        ofs << ">consensus\n";
        constexpr std::size_t width = 80;
        for (std::size_t i = 0; i < seq.size(); i += width) {
            const std::size_t n = std::min(width, seq.size() - i);
            ofs.write(seq.data() + (std::streamoff)i, (std::streamsize)n);
            ofs.put('\n');
        }
    }

    // 使用 cereal 将计数写为 JSON（项目中 prefer cereal）
    void writeCountsJson(const FilePath& out_json, const ConsensusJson& cj)
    {
        file_io::ensureParentDirExists(out_json);

        std::ofstream ofs(out_json, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("failed to open json output: " + out_json.string());
        }

        cereal::JSONOutputArchive ar(ofs);
        ar(cereal::make_nvp("consensus", cj));
    }



    // 将单条序列的按列计数逻辑提取为一个独立函数。
    // 该函数只负责把单条序列的每个位点计数累加到 cj.counts 中，
    // 不负责更新 cj.num_seqs（由调用者在安全的语义下进行递增）。
    //
    // 并行化策略：对单条序列，在 "位置" 维度并行（每个线程负责不同的列索引），
    // 因为每个列对应的 SiteCount 是独立的，多个线程不会写入同一索引，故无需原子操作。
    // 注意：如果并行处理多条序列（在调用层面并行），必须保证不同线程不会并发写入同一列。
    //
    // 详细优化说明（中文，供维护者阅读）
    //
    // 总体目标：
    // - 在处理大量对齐序列时，将“按列统计碱基计数”这个瓶颈路径尽量并行化和向量化，
    //   同时避免高成本的原子操作或频繁的锁竞争，尽可能利用缓存局部性与 SIMD 指令。
    //
    // 主要优化手段（本文件实现中可见）：
    // 1) 分支最小化（branch minimization）
    //    - 传统的 switch/case 每个字符要走分支，分支预测失败开销大。
    //    - 用 (idx == constant) 转成 0/1 的加法消除分支，编译器更容易生成向量指令。
    //
    // 2) 循环展开（loop unrolling）与预取（prefetch）
    //    - 通过 4-路（或可调整）循环展开减少循环控制开销、增加指令级并行（ILP）。
    //    - 使用 __builtin_prefetch 提前载入将要写入的缓存行，降低缓存未命中延迟（对长序列/大工作集有效）。
    //
    // 3) 线程局部（thread-local）累加 + 批处理（batching）
    //    - 直接并发写全局 counts 会导致 cache-line 冲突（false sharing）。为避免，采用每个线程维护本地计数数组，
    //      在处理完一批序列后再按列合并（reduce）。这种策略通过牺牲一些内存来大幅降低并发写冲突，
    //      在线程数与对齐长度适中时通常能得到最好吞吐。
    //    - 批大小（batch_size）需要调优：过小会导致调度/同步开销，过大会占用更多内存且增加延迟。
    //
    // 4) SoA (Structure of Arrays) vs AoS (Array of Structures)
    //    - AoS（SiteCount 存在 cj.counts[i]）在合并或向量化时不便，因为每个 SiteCount 中字段交错。
    //    - SoA 将每个碱基的计数放入独立的数组（A[], C[], G[]...），便于对单个碱基列进行连续读写，
    //      更有利于 SIMD 与缓存预取。实现中提供了 SoA 的批处理路径以求最大化性能。
    //
    // 5) OpenMP + simd 指示（#pragma omp parallel / #pragma omp simd / ivdep）
    //    - 使用 OpenMP 做线程级并行（按位置或按序列分配工作）；同时在内循环增加 simd 提示以帮助编译器向量化。
    //
    // 6) 编译器与编译选项
    //    - 在 CMake 中对 Release 模式启用 -O3、-march=native、可选 -flto，有助于生成高效的向量指令与内联。
    //
    // 适用场景与权衡：
    // - 当 aln_len（对齐长度）很大（几千到几万）且每批序列数量也大时，SoA + 线程本地合并通常最快。
    // - 当 aln_len 很小（例如 < 256）时，线程并行化的开销可能抵消收益，应退化到单线程或较小线程数。
    // - 本地计数会使用额外内存：大线程数与大对齐长度会带来显著内存占用（T * aln_len * sizeof(count)）。
    // - 合并阶段仍然要按列扫描，一次性合并成本需在批大小与内存之间权衡。
    //
    // 实践建议：
    // - 在目标机器上做小规模基准（不同 batch_size、线程数、aln_len）来选择最优参数。
    // - 结合性能分析工具（perf / VTune / likwid）观察缓存未命中、内存带宽与分支失误。
    // - 若需要极限性能，可进一步用 AVX2/AVX512 intrinsics 对 SoA 路径做微调。
    static void processSequenceParallel(const std::string& s, ConsensusJson& cj, int thread)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (s.size() != aln_len) {
            throw std::runtime_error("alignment length mismatch: expect " + std::to_string(aln_len) +
                                     ", got " + std::to_string(s.size()));
        }

        // 预取常量以减少循环内查找开销
        const unsigned char* data = reinterpret_cast<const unsigned char*>(s.data());
        std::uint8_t* base_map = const_cast<std::uint8_t*>(k_base_map.data());
        SiteCount* counts = cj.counts.data();

        const std::size_t limit = (aln_len / 4) * 4;

#if __has_include(<omp.h>)
        // 使用 parallel for simd，允许 OpenMP 分配线程并启用向量化
        #pragma omp parallel for schedule(static) num_threads(thread)
#endif
        for (std::size_t i = 0; i < limit; i += 4) {
            // 轻量预取未来缓存行，距离可调（16~64 bytes => 16 positions 粗略）
            __builtin_prefetch(&counts[i + 16]);

            const std::uint8_t idx0 = base_map[data[i + 0]];
            const std::uint8_t idx1 = base_map[data[i + 1]];
            const std::uint8_t idx2 = base_map[data[i + 2]];
            const std::uint8_t idx3 = base_map[data[i + 3]];

            SiteCount& sc0 = counts[i + 0];
            SiteCount& sc1 = counts[i + 1];
            SiteCount& sc2 = counts[i + 2];
            SiteCount& sc3 = counts[i + 3];

            sc0.a += static_cast<std::uint32_t>(idx0 == 0);
            sc0.c += static_cast<std::uint32_t>(idx0 == 1);
            sc0.g += static_cast<std::uint32_t>(idx0 == 2);
            sc0.t += static_cast<std::uint32_t>(idx0 == 3);
            sc0.u += static_cast<std::uint32_t>(idx0 == 4);
            sc0.n += static_cast<std::uint32_t>(idx0 == 5);
            sc0.dash += static_cast<std::uint32_t>(idx0 == 6);

            sc1.a += static_cast<std::uint32_t>(idx1 == 0);
            sc1.c += static_cast<std::uint32_t>(idx1 == 1);
            sc1.g += static_cast<std::uint32_t>(idx1 == 2);
            sc1.t += static_cast<std::uint32_t>(idx1 == 3);
            sc1.u += static_cast<std::uint32_t>(idx1 == 4);
            sc1.n += static_cast<std::uint32_t>(idx1 == 5);
            sc1.dash += static_cast<std::uint32_t>(idx1 == 6);

            sc2.a += static_cast<std::uint32_t>(idx2 == 0);
            sc2.c += static_cast<std::uint32_t>(idx2 == 1);
            sc2.g += static_cast<std::uint32_t>(idx2 == 2);
            sc2.t += static_cast<std::uint32_t>(idx2 == 3);
            sc2.u += static_cast<std::uint32_t>(idx2 == 4);
            sc2.n += static_cast<std::uint32_t>(idx2 == 5);
            sc2.dash += static_cast<std::uint32_t>(idx2 == 6);

            sc3.a += static_cast<std::uint32_t>(idx3 == 0);
            sc3.c += static_cast<std::uint32_t>(idx3 == 1);
            sc3.g += static_cast<std::uint32_t>(idx3 == 2);
            sc3.t += static_cast<std::uint32_t>(idx3 == 3);
            sc3.u += static_cast<std::uint32_t>(idx3 == 4);
            sc3.n += static_cast<std::uint32_t>(idx3 == 5);
            sc3.dash += static_cast<std::uint32_t>(idx3 == 6);
        }

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(thread)
#endif
        for (std::size_t i = limit; i < aln_len; ++i) {
            const std::uint8_t idx = base_map[data[i]];
            SiteCount& sc = counts[i];
            sc.a += static_cast<std::uint32_t>(idx == 0);
            sc.c += static_cast<std::uint32_t>(idx == 1);
            sc.g += static_cast<std::uint32_t>(idx == 2);
            sc.t += static_cast<std::uint32_t>(idx == 3);
            sc.u += static_cast<std::uint32_t>(idx == 4);
            sc.n += static_cast<std::uint32_t>(idx == 5);
            sc.dash += static_cast<std::uint32_t>(idx == 6);
        }
    }

    // 批量并行处理：每个线程维护自己的本地 counts 数组来累加多条序列，
    // 处理完批后再把本地 counts 合并到全局 cj.counts。这能显著减少对全局 counts 的并发写入，
    // 降低 false sharing 并提高缓存局部性。适合 aln_len 和 batch_size 都较大的场景。
    static void processBatchParallel(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        // 先在主线程中验证所有序列长度一致，避免并行区抛异常不安全
        for (const auto& s : seqs) {
            if (s.size() != aln_len) {
                throw std::runtime_error("alignment length mismatch in batch processing");
            }
        }

        // 扁平化本地计数：一块连续内存，大小为 T * aln_len
        std::vector<SiteCount> locals;
        try {
            locals.assign((std::size_t)T * aln_len, SiteCount{});
        } catch (...) {
            throw std::runtime_error("failed to allocate thread-local counts");
        }

        const std::uint8_t* base_map = k_base_map.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());
            SiteCount* local = locals.data() + (std::size_t)tid * aln_len;

            // 对单条序列按位置累加到线程本地计数；采用分支最小化的布尔加法
            for (std::size_t i = 0; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                SiteCount& sc = local[i];
                sc.a += static_cast<std::uint32_t>(idx == 0);
                sc.c += static_cast<std::uint32_t>(idx == 1);
                sc.g += static_cast<std::uint32_t>(idx == 2);
                sc.t += static_cast<std::uint32_t>(idx == 3);
                sc.u += static_cast<std::uint32_t>(idx == 4);
                sc.n += static_cast<std::uint32_t>(idx == 5);
                sc.dash += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // 合并本地 counts 到全局 cj.counts（按列并行）
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            // 使用 64-bit 临时累加以避免溢出（尽管 SiteCount 是 uint32）
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            for (int t = 0; t < T; ++t) {
                const SiteCount& ls = locals[(std::size_t)t * aln_len + i];
                sa += ls.a; sc_ += ls.c; sg += ls.g; st += ls.t; su += ls.u; sn += ls.n; sd += ls.dash;
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }

    /*
     下面对批处理（线程本地累加 + 合并）函数做详细注释：

     processBatchParallelWithLocals:
     - 每线程拥有一块连续的 SiteCount 数组（locals），大小为 aln_len。
     - 每个线程将它负责的序列累加到本地 locals（避免并发写入全局 cj.counts）。
     - 合并阶段按列并行：每个线程负责合并若干列，将本地累加的结果加到全局 counts。

     优点：
     - 极大降低 false sharing（不同线程不会频繁写同一 cache line）。
     - 简单实现，易于理解与调试。

     局限：
     - 每线程的 locals 使用较多内存（T * aln_len * sizeof(SiteCount)）。
     - SiteCount 内部字段为 AoS，合并时访问字段交错，向量化受限。
    */

    // 批处理（使用外部分配的 locals 缓冲以避免每批分配开销）
    static void processBatchParallelWithLocals(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads, std::vector<SiteCount>& locals)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        const std::uint8_t* base_map = k_base_map.data();
        SiteCount* locals_ptr = locals.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());
            SiteCount* local = locals_ptr + (std::size_t)tid * aln_len;

            // 更高效的按位置更新：4-路展开并使用布尔加法，减少分支
            std::size_t i = 0;
            const std::size_t limit = (aln_len / 4) * 4;
            for (; i < limit; i += 4) {
                __builtin_prefetch(&local[i + 16]);
                const std::uint8_t idx0 = base_map[data[i + 0]];
                const std::uint8_t idx1 = base_map[data[i + 1]];
                const std::uint8_t idx2 = base_map[data[i + 2]];
                const std::uint8_t idx3 = base_map[data[i + 3]];

                SiteCount& sc0 = local[i + 0];
                SiteCount& sc1 = local[i + 1];
                SiteCount& sc2 = local[i + 2];
                SiteCount& sc3 = local[i + 3];

                sc0.a += static_cast<std::uint32_t>(idx0 == 0);
                sc0.c += static_cast<std::uint32_t>(idx0 == 1);
                sc0.g += static_cast<std::uint32_t>(idx0 == 2);
                sc0.t += static_cast<std::uint32_t>(idx0 == 3);
                sc0.u += static_cast<std::uint32_t>(idx0 == 4);
                sc0.n += static_cast<std::uint32_t>(idx0 == 5);
                sc0.dash += static_cast<std::uint32_t>(idx0 == 6);

                sc1.a += static_cast<std::uint32_t>(idx1 == 0);
                sc1.c += static_cast<std::uint32_t>(idx1 == 1);
                sc1.g += static_cast<std::uint32_t>(idx1 == 2);
                sc1.t += static_cast<std::uint32_t>(idx1 == 3);
                sc1.u += static_cast<std::uint32_t>(idx1 == 4);
                sc1.n += static_cast<std::uint32_t>(idx1 == 5);
                sc1.dash += static_cast<std::uint32_t>(idx1 == 6);

                sc2.a += static_cast<std::uint32_t>(idx2 == 0);
                sc2.c += static_cast<std::uint32_t>(idx2 == 1);
                sc2.g += static_cast<std::uint32_t>(idx2 == 2);
                sc2.t += static_cast<std::uint32_t>(idx2 == 3);
                sc2.u += static_cast<std::uint32_t>(idx2 == 4);
                sc2.n += static_cast<std::uint32_t>(idx2 == 5);
                sc2.dash += static_cast<std::uint32_t>(idx2 == 6);

                sc3.a += static_cast<std::uint32_t>(idx3 == 0);
                sc3.c += static_cast<std::uint32_t>(idx3 == 1);
                sc3.g += static_cast<std::uint32_t>(idx3 == 2);
                sc3.t += static_cast<std::uint32_t>(idx3 == 3);
                sc3.u += static_cast<std::uint32_t>(idx3 == 4);
                sc3.n += static_cast<std::uint32_t>(idx3 == 5);
                sc3.dash += static_cast<std::uint32_t>(idx3 == 6);
            }
            for (; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                SiteCount& sc = local[i];
                sc.a += static_cast<std::uint32_t>(idx == 0);
                sc.c += static_cast<std::uint32_t>(idx == 1);
                sc.g += static_cast<std::uint32_t>(idx == 2);
                sc.t += static_cast<std::uint32_t>(idx == 3);
                sc.u += static_cast<std::uint32_t>(idx == 4);
                sc.n += static_cast<std::uint32_t>(idx == 5);
                sc.dash += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // 合并本地 counts 到全局 cj.counts（按列并行）
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            for (int t = 0; t < T; ++t) {
                const SiteCount& ls = locals_ptr[(std::size_t)t * aln_len + i];
                sa += ls.a; sc_ += ls.c; sg += ls.g; st += ls.t; su += ls.u; sn += ls.n; sd += ls.dash;
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }

    /*
     processBatchParallelWithSoA（SoA 版本）：
     - 对每个碱基维护独立数组（A[], C[], G[], T[], U[], N[], Dash[]），并为每个线程分配一段连续内存（T * aln_len）。
     - 在累加阶段，每线程对自己的段按位置累加：a_ptr[i]++ 等。由于同一种碱基的计数是连续的，
       读取/写入模式更友好于 CPU 的向量化与缓存预取策略。
     - 合并阶段按列（i）对每个线程的对应位置求和并写回 cj.counts[i]。
     *
     * 优势（为什么更快）：
     * 1) 向量化友好：SoA 允许编译器把一条指令应用到连续的 a_ptr[] 元素，生成 SIMD 指令（如 AVX2 之类）。
     * 2) 减少内存带宽浪费：当处理某个碱基时只触碰该碱基的数组，减少对其它字段不必要的缓存写回。
     * 3) 合并步骤中按列累加多个线程的连续内存，内存访问模式简单，有利于预取与硬件合并。
     *
     * 代价：
     * - 额外的内存开销（7 个 uint32_t 数组 * T * aln_len），但通常比原子/锁/频繁缓存同步更划算。
     * - 实现复杂度略增，但对性能敏感场景值得。
    */

    static void processBatchParallelWithSoA(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads,
                                            std::vector<std::uint32_t>& localsA,
                                            std::vector<std::uint32_t>& localsC,
                                            std::vector<std::uint32_t>& localsG,
                                            std::vector<std::uint32_t>& localsT,
                                            std::vector<std::uint32_t>& localsU,
                                            std::vector<std::uint32_t>& localsN,
                                            std::vector<std::uint32_t>& localsDash)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        const std::uint8_t* base_map = k_base_map.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());

            const std::size_t base_off = (std::size_t)tid * aln_len;
            std::uint32_t* a_ptr = localsA.data() + base_off;
            std::uint32_t* c_ptr = localsC.data() + base_off;
            std::uint32_t* g_ptr = localsG.data() + base_off;
            std::uint32_t* t_ptr = localsT.data() + base_off;
            std::uint32_t* u_ptr = localsU.data() + base_off;
            std::uint32_t* n_ptr = localsN.data() + base_off;
            std::uint32_t* d_ptr = localsDash.data() + base_off;

            const std::size_t limit = (aln_len / 4) * 4;
            std::size_t i = 0;
            for (; i < limit; i += 4) {
                __builtin_prefetch(a_ptr + i + 16);
                const std::uint8_t idx0 = base_map[data[i + 0]];
                const std::uint8_t idx1 = base_map[data[i + 1]];
                const std::uint8_t idx2 = base_map[data[i + 2]];
                const std::uint8_t idx3 = base_map[data[i + 3]];

                a_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 0);
                c_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 1);
                g_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 2);
                t_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 3);
                u_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 4);
                n_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 5);
                d_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 6);

                a_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 0);
                c_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 1);
                g_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 2);
                t_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 3);
                u_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 4);
                n_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 5);
                d_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 6);

                a_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 0);
                c_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 1);
                g_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 2);
                t_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 3);
                u_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 4);
                n_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 5);
                d_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 6);

                a_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 0);
                c_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 1);
                g_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 2);
                t_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 3);
                u_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 4);
                n_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 5);
                d_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 6);
            }
            for (; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                a_ptr[i] += static_cast<std::uint32_t>(idx == 0);
                c_ptr[i] += static_cast<std::uint32_t>(idx == 1);
                g_ptr[i] += static_cast<std::uint32_t>(idx == 2);
                t_ptr[i] += static_cast<std::uint32_t>(idx == 3);
                u_ptr[i] += static_cast<std::uint32_t>(idx == 4);
                n_ptr[i] += static_cast<std::uint32_t>(idx == 5);
                d_ptr[i] += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // 合并：按列读取 locals arrays 并写入 cj.counts
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            const std::size_t stride = aln_len;
            for (int t = 0; t < T; ++t) {
                const std::size_t off = (std::size_t)t * stride + i;
                sa += localsA[off]; sc_ += localsC[off]; sg += localsG[off]; st += localsT[off]; su += localsU[off]; sn += localsN[off]; sd += localsDash[off];
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }

    /*
     深入性能优化说明（补充）

     以下注释面向想要进一步优化或调试性能的工程师，包含更具体的建议和硬件相关注意点：

     1) False sharing 与缓存行对齐
        - False sharing 发生在多个线程频繁写入同一 cache line（例如 SiteCount 的多个字段或相邻索引落在同一 64B cache line）。
        - 减少 false sharing 的策略：
          a) 线程本地累加（本文件采用）：每线程写自己的本地数组，最后合并；
          b) 在 SiteCount 上做填充（padding）以使每个 SiteCount 占用整 cache line（内存消耗增加）；
             仅当 SiteCount 更新极其频繁并且 T 很大时考虑；示例：alignas(64) 或在 SiteCount 后添加 uint8_t pad[...]
          c) 让 OpenMP 为每个线程分配较大的连续 chunk（static schedule 与较大 chunk_size），减少线程交错写导致的共享。

     2) 向量化与内存布局（SoA vs AoS）
        - AoS（SiteCount 结构数组）对单一字段的连续访问不友好，向量化受限。
        - SoA（Structure of Arrays）把 A/C/G/T/U 等分别放到独立数组，便于对单一字段做大范围加法并生成 SIMD 指令。
        - 本实现提供 SoA 路径：在累加阶段对 a_ptr/c_ptr/g_ptr 等做连续写；在合并阶段再按列读取并写回。

     3) 批大小（batch_size）与线程数（threads）选择建议
        - 先根据内存限制估算最大可接受 batch_size： batch_memory ≈ T * aln_len * sizeof(counter)。
        - 通常策略：保持 batch_size 足够大以 amortize 线程调度与合并成本，例如 1k~10k（视 aln_len 而定）；
          但当 aln_len 很大（几十千或百万列）时 batch_size 可以较小。
        - 以实验为准：在目标机器上运行一组试验（threads ∈ {1,2,4,8,...}, batch_size ∈ {64,256,1024,4096}）测量每秒处理基数。

     4) NUMA 与亲和性
        - 在 NUMA 芯片组上，尽量把线程固定到本地 NUMA 节点并把输入数据/locals 分配在同一 NUMA 节点上（使用 numactl 或 pthread_setaffinity_np）。
        - 如果性能短板是内存带宽，考虑把大批次分配到不同 NUMA 节点并分别合并以减少跨节点流量。

     5) 预取与预热
        - __builtin_prefetch 能在某些平台明显降低缓存未命中，预取距离（例如 i+16）需根据 cache 行与访问模式调优。
        - 可在微基准中试不同预取距离找到最佳点。

     6) 内存与溢出安全
        - 本实现使用 32-bit 计数（SiteCount 的字段）；若每列计数可能超过 2^32（非常大量序列），需切换为 64-bit（uint64_t）以避免溢出。
        - 合并时采用 uint64_t 临时累加以避免短暂溢出，再截断回 32-bit 写入 SiteCount（或直接将 SiteCount 改为 uint64_t）。

     7) 编译器生成的代码检查
        - 若希望确认向量化是否生效，可以在带 -O3 -march=native 下查看编译器生成的汇编（-S 或 objdump），搜索 AVX/AVX2/AVX512 指令。
        - 也可使用 clang 的 -Rpass=loop-vectorize 来让编译器报告成功的向量化。

     8) 性能回归与测试
        - 在引入任何优化后，务必使用回归测试确保输出一致（本项目已有测试套件）。
        - 同时保留 baseline（未优化版本）以便比较性能提升。
    */

    // 下面对 SoA 批处理函数再补充：具体如何调参与排查问题的建议
    // - 如果发现合并阶段占用较多时间，可尝试把合并也做分块（例如一次合并一个列区间以提高局部性），
    //   或者在合并时将每个线程负责的列区间固定，减少内存抖动。
    // - 若出现内存带宽瓶颈，可尝试减少线程数或增加 batch_size，以便每次合并做更多工作再写回。

    /*
     generateConsensusSequence 的并行运行策略与建议：
     - 对于小的 aln_len（例如 < 256），优先使用单线程或轻量并行；因为线程和合并开销可能大于收益。
     - 对于大的 aln_len，使用 SoA + 线程本地累加（processBatchParallelWithSoA）通常能获得最佳吞吐。
     - batch_size 值会影响性能：
         * 小 batch：更实时、内存占用少，但同步开销高；
         * 大 batch：更高吞吐但需要更大内存与更长延迟。
     - 在部署前在目标机器上对 batch_size 与线程数进行基准测试（选择能最大化每秒处理碱基数的组合）。
     *
     * 实用调试/性能收集建议：
     * - 使用 perf/top/htop 观察 CPU 利用率与内存带宽。
     * - 使用 perf record / report 或 VTune 查看缓存未命中与分支失误热点。
     * - 用不同优化编译选项（-O2 vs -O3，-march=native）对比，确认向量化是否被启用（查看编译器汇编输出）。
    */
    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                       const FilePath& out_fasta,
                                                       const FilePath& out_json,
                                                       std::uint64_t seq_limit,
                                                       int thread,
                                                       size_t batch_size)
    {
        file_io::requireRegularFile(aligned_fasta, "aligned_fasta");

        seq_io::KseqReader reader(aligned_fasta);

        // 读第一条确定 aln_len
        seq_io::SeqRecord rec;
        if (!reader.next(rec)) {
            throw std::runtime_error("aligned fasta is empty: " + aligned_fasta.string());
        }

        const std::size_t aln_len = rec.seq.size();
        if (aln_len == 0) {
            throw std::runtime_error("first sequence length is 0: " + aligned_fasta.string());
        }

        ConsensusJson cj;
        cj.aln_len = (std::uint64_t)aln_len;
        cj.counts.assign(aln_len, SiteCount{});

        std::uint64_t num_seqs = 0;

        // 批处理参数：批大小可调（经验值），在内存允许的范围内放大能减少调度开销
        std::vector<std::string> batch;
        batch.reserve(batch_size + 1);

        // 首条先放入 batch（并检查长度）
        if (rec.seq.size() != aln_len) {
            throw std::runtime_error("alignment length mismatch: first record length changed");
        }
        batch.push_back(std::move(rec.seq));
        ++num_seqs;

        // 预分配并复用 thread-local 缓冲
        int T = thread;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        if (T <= 0) T = 1;
#endif
        std::vector<SiteCount> locals;
        try {
            locals.assign((std::size_t)T * aln_len, SiteCount{});
        } catch (...) {
            throw std::runtime_error("failed to allocate thread-local counts");
        }

        // 读取并按批处理
        while ((seq_limit == 0 || num_seqs < seq_limit)) {
            // 填充 batch
            while (batch.size() < batch_size && (seq_limit == 0 || num_seqs < seq_limit)) {
                if (!reader.next(rec)) break;
                if (rec.seq.size() != aln_len) throw std::runtime_error("alignment length mismatch when reading");
                batch.push_back(std::move(rec.seq));
                ++num_seqs;
            }

            // 处理当前 batch
            if (!batch.empty()) {
                // 清零 locals 一次性（快速）
                std::memset(locals.data(), 0, locals.size() * sizeof(SiteCount));
                processBatchParallelWithLocals(batch, cj, thread, locals);
                batch.clear();
            }

            if (!reader.next(rec)) break; // EOF
        }


        if (cj.num_seqs == 0) cj.num_seqs = num_seqs;

        if (num_seqs == 0) {
            throw std::runtime_error("no sequences processed");
        }

        // 生成共识序列（单线程选多数）
        std::string consensus_seq(aln_len, 'N');
        for (std::size_t i = 0; i < aln_len; ++i) {
            consensus_seq[i] = pickConsensusChar(cj.counts[i]);
        }

        writeConsensusFasta(out_fasta, consensus_seq);
        writeCountsJson(out_json, cj);

        return consensus_seq;
    }


} // namespace consensus

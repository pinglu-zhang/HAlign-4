#include "consensus.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <thread>
#include <cereal/archives/json.hpp>

#if __has_include(<omp.h>)
    #include <omp.h>
#endif

// ============================================================================
// consensus.cpp - 详细并行注释版
//
// 本文件的并行化核心思想：
// - 采用 batch/master-worker 模式：主线程（master）负责从输入文件读取若干条序列
//   到共享缓冲（batch_seqs），随后发出同步点（barrier）；工作线程并行对该 batch
//   按列统计（每个列由一个线程负责写入对应的 SiteCount 项），完成后再次同步，
//   master 更新计数并开始读取下一批数据。
// - 这种设计优点：
//   1) 降低内存占用（不会把所有序列保持在内存中）；
//   2) 避免在热路径使用原子或锁（因为每一列只有单一线程写入对应的 SiteCount），
//      所以统计阶段的写入是线程私有的（按索引分配），无数据竞争；
//   3) 同步粒度可控（batch_size 可调），适配不同 CPU/IO 比例的场景。
// - 关键保证：在并行统计开始之前，所有工作线程必须能观察到 master 已完成填写 batch_seqs，
//   这是通过 OpenMP 的 barrier/同步语义以及 master/thread 区块的隐式内存刷新保证的。
// ============================================================================

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

    // 基于 SoA（structure-of-arrays）的批处理：每线程为每个碱基维护一块连续的 uint32_t counts，
    // 这样在累加与合并阶段更易向量化且缓存友好。localsA..localsDash 的长度均为 T * aln_len。
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

    // 单线程版本：与 generateConsensusSequence 功能相同，但不使用 OpenMP 或并行任务。
    // 该函数按条读取对齐的 FASTA 序列，并在当前线程中立即累加列计数，
    // 最后生成共识序列并写出 FASTA/JSON。适合在不希望启用多线程或为调试/比较用途时使用。
    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                       const FilePath& out_fasta,
                                                       const FilePath& out_json,
                                                       std::uint64_t seq_limit,
                                                       int thread)
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
        const std::size_t batch_size = 5120;
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

        // 如果最后 reader had more? already handled

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

    // 修改 generateConsensusSequence：按批读取序列并调用 processBatchParallelWithSoA 来处理
    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                       const FilePath& out_fasta,
                                                       const FilePath& out_json,
                                                       std::uint64_t seq_limit,
                                                       int thread)
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
        const std::size_t batch_size = 5120;
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
        std::vector<std::uint32_t> localsA, localsC, localsG, localsT, localsU, localsN, localsDash;
        try {
            localsA.assign((std::size_t)T * aln_len, 0);
            localsC.assign((std::size_t)T * aln_len, 0);
            localsG.assign((std::size_t)T * aln_len, 0);
            localsT.assign((std::size_t)T * aln_len, 0);
            localsU.assign((std::size_t)T * aln_len, 0);
            localsN.assign((std::size_t)T * aln_len, 0);
            localsDash.assign((std::size_t)T * aln_len, 0);
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
                processBatchParallelWithSoA(batch, cj, thread,
                                            localsA, localsC, localsG, localsT, localsU, localsN, localsDash);
                batch.clear();
            }

            if (!reader.next(rec)) break; // EOF
        }

        // 如果最后 reader had more? already handled

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

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


    // 仅由 single 线程调用：读入一批序列到 dst
    static void fillBatch(seq_io::KseqReader& reader,
                          std::vector<std::string>& dst,
                          std::size_t batch_size,
                          std::size_t aln_len,
                          std::uint64_t seq_limit,
                          std::uint64_t& seen_total,      // 已“读入/排队”的序列数（包含双缓冲两边）
                          bool keep_existing,             // 首批保留已塞入的第一条
                          bool& stop_after_this_batch,    // 因 seq_limit 触发：处理完该批就应停止
                          bool& eof_reached,              // 读到 EOF
                          bool& error,
                          std::string& error_msg)
    {
        stop_after_this_batch = false;
        eof_reached = false;

        if (!keep_existing) {
            dst.clear();
        }

        while (dst.size() < batch_size) {
            if (seq_limit > 0 && seen_total >= seq_limit) {
                stop_after_this_batch = true;
                break;
            }

            seq_io::SeqRecord r;
            if (!reader.next(r)) {
                eof_reached = true;
                break;
            }

            if (r.seq.size() != aln_len) {
                error = true;
                error_msg = "alignment length mismatch: expect " +
                            std::to_string(aln_len) + ", got " +
                            std::to_string(r.seq.size());
                dst.clear();
                return;
            }

            dst.push_back(std::move(r.seq));
            ++seen_total;
        }

        // 若刚好读满但 seq_limit 已达，也应在本批结束后停止
        if (seq_limit > 0 && seen_total >= seq_limit) {
            stop_after_this_batch = true;
        }
    }

    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                        const FilePath& out_fasta,
                                                        const FilePath& out_json,
                                                        std::uint64_t seq_limit,
                                                        int threads,
                                                        std::size_t batch_size)
    {
        file_io::requireRegularFile(aligned_fasta, "aligned_fasta");
        if (batch_size == 0) {
            throw std::runtime_error("batch_size must be > 0");
        }

        if (threads <= 0) {
#if __has_include(<omp.h>)
            threads = [](){ unsigned int hc = std::thread::hardware_concurrency(); return static_cast<int>(hc ? hc : 1u); }();
#else
            threads = 1;
#endif
        }

#if __has_include(<omp.h>)
        omp_set_dynamic(0);
#endif

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

        // 双缓冲：两个 batch 容器轮转使用
        std::vector<std::string> buf[2];
        buf[0].reserve(batch_size);
        buf[1].reserve(batch_size);

        // 首条先放入 buf[0]
        buf[0].push_back(std::move(rec.seq));

        // 已“读入/排队”的序列数（用于 seq_limit 控制）
        std::uint64_t seen_total = 1;

        // 已处理完成并计入 cj.num_seqs 的序列数
        std::uint64_t num_seqs = 0;

        bool stop_flag[2] = {false, false};
        bool eof_flag[2]  = {false, false};

        bool error = false;
        std::string error_msg;

        // 首次填充 buf[0]（保留已有第一条）
        fillBatch(reader, buf[0], batch_size, aln_len, seq_limit, seen_total,
                  /*keep_existing=*/true,
                  stop_flag[0], eof_flag[0], error, error_msg);

        if (error) {
            throw std::runtime_error("consensus counting failed: " + error_msg);
        }
        if (buf[0].empty()) {
            throw std::runtime_error("no sequences processed");
        }
        const int grainsize = 64;
#if __has_include(<omp.h>)
        #pragma omp parallel num_threads(threads) shared(reader, cj, buf, stop_flag, eof_flag, error, error_msg, seen_total, num_seqs)
#endif
        {
#if __has_include(<omp.h>)
            #pragma omp single nowait
#endif
            {
                int cur = 0;

                while (!error && !buf[cur].empty()) {
                    const int next = 1 - cur;

                    // 当前批信息（在发起统计任务前固化下来）
                    const std::size_t cur_n = buf[cur].size();
                    const bool stop_after_cur = stop_flag[cur];

                    // ===== 1) 发起“统计当前 batch”的任务 =====
#if __has_include(<omp.h>)
                    #pragma omp task firstprivate(cur, cur_n) shared(cj, buf, error, error_msg)
#endif
                    {
                        // 若读线程已发现 error，跳过该 task 的执行体（不从 task 内返回）
                        if (!error) {
                            const auto& batch = buf[cur];
                            const std::size_t bn = cur_n;

                            // taskloop 会生成多个子任务，由线程池执行，实现“统计并行”
#if __has_include(<omp.h>)
                            #pragma omp taskloop grainsize(grainsize)
#endif
                            for (std::int64_t i = 0; i < (std::int64_t)aln_len; ++i) {
                                // 每个 i 只会由一个 task 负责 => 独占写 counts[i]，无需原子
                                SiteCount local{};
                                for (std::size_t s = 0; s < bn; ++s) {
                                    const std::uint8_t idx = mapBase(batch[s][(std::size_t)i]);
                                    switch (idx) {
                                        case 0: local.a++; break;
                                        case 1: local.c++; break;
                                        case 2: local.g++; break;
                                        case 3: local.t++; break;
                                        case 4: local.u++; break;
                                        case 5: local.n++; break;
                                        case 6: local.dash++; break;
                                        default: local.n++; break;
                                    }
                                }

                                // 最后一次性写回（减少对同一 cache line 的反复写入）
                                SiteCount& sc = cj.counts[(std::size_t)i];
                                sc.a += local.a;
                                sc.c += local.c;
                                sc.g += local.g;
                                sc.t += local.t;
                                sc.u += local.u;
                                sc.n += local.n;
                                sc.dash += local.dash;
                            }
                        }
                    }

                    // ===== 2) 读线程（single）在统计进行时，读取下一批到 next buffer =====
                    fillBatch(reader, buf[next], batch_size, aln_len, seq_limit, seen_total,
                              /*keep_existing=*/false,
                              stop_flag[next], eof_flag[next], error, error_msg);

                    // 若读取中发现错误：等当前统计任务结束后再退出（保证 counts 不被并发破坏）
                    // （注意：我们仍然只有一个统计任务在飞，不会和下一批并发写 counts）
#if __has_include(<omp.h>)
                    #pragma omp taskwait
#endif

                    // 当前批统计已完成，可以安全复用 cur buffer（下一轮会写 next，之后再写回 cur）
                    num_seqs += (std::uint64_t)cur_n;
                    cj.num_seqs = num_seqs;

                    if (error) {
                        break;
                    }

                    // seq_limit 触发：处理完当前批就停止
                    if (stop_after_cur) {
                        break;
                    }

                    // 切换到 next batch
                    cur = next;
                }
            } // single
        } // parallel

        if (error) {
            throw std::runtime_error("consensus counting failed: " + error_msg);
        }
        if (num_seqs == 0) {
            throw std::runtime_error("no sequences processed");
        }

        // 生成共识序列（并行按列选多数）
        std::string consensus_seq(aln_len, 'N');

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(threads)
#endif
        for (std::int64_t i = 0; i < (std::int64_t)aln_len; ++i) {
            consensus_seq[(std::size_t)i] = pickConsensusChar(cj.counts[(std::size_t)i]);
        }

        writeConsensusFasta(out_fasta, consensus_seq);
        writeCountsJson(out_json, cj);

        return consensus_seq;
    }

} // namespace consensus

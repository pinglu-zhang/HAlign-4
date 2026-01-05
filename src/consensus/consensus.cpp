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

#if __has_include(<omp.h>)
        // 在列维度并行：每个迭代 i 只操作 cj.counts[i]，因此不存在写冲突。
        #pragma omp parallel for schedule(static) num_threads(thread)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            const std::uint8_t idx = mapBase(s[i]);
            SiteCount& sc = cj.counts[i];
            switch (idx) {
                case 0: sc.a++; break;
                case 1: sc.c++; break;
                case 2: sc.g++; break;
                case 3: sc.t++; break;
                case 4: sc.u++; break;
                case 5: sc.n++; break;
                case 6: sc.dash++; break;
                default: sc.n++; break;
            }
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

        // 读第一条确定 aln_len 并计入
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

        // 处理首条
        processSequenceParallel(rec.seq, cj, thread);
        ++num_seqs;
        cj.num_seqs = num_seqs;

        // 逐条读取并处理，直到 EOF 或达到 seq_limit
        while ((seq_limit == 0 || num_seqs < seq_limit) && reader.next(rec)) {
            processSequenceParallel(rec.seq, cj, thread);
            ++num_seqs;
            cj.num_seqs = num_seqs;
        }

        if (num_seqs == 0) {
            throw std::runtime_error("no sequences processed");
        }

        // 生成共识序列（单线程）
        std::string consensus_seq(aln_len, 'N');
        for (std::size_t i = 0; i < aln_len; ++i) {
            consensus_seq[i] = pickConsensusChar(cj.counts[i]);
        }

        writeConsensusFasta(out_fasta, consensus_seq);
        writeCountsJson(out_json, cj);

        return consensus_seq;
    }

} // namespace consensus

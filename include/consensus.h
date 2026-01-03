#ifndef HALIGN4_CONSENSUS_H
#define HALIGN4_CONSENSUS_H

#include <cstddef>
#include "utils.h"
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>

// 维护“最长的 k 条序列”，用于共识构建的候选集合
class TopKLongestSelector
{
public:
    explicit TopKLongestSelector(std::size_t k = 0);

    void reset(std::size_t k);

    // 传值：调用方可选择复制传入或 std::move 传入，内部会 move 到堆里
    void consider(seq_io::SeqRecord rec);

    std::size_t size() const;
    std::size_t capacity() const;
    bool empty() const;

    // 取出结果（按长度降序；同长度按输入出现顺序升序），并清空内部状态
    std::vector<seq_io::SeqRecord> takeSortedDesc();

private:
    struct Item
    {
        std::size_t len{0};
        std::uint64_t order{0};   // 输入顺序（用于稳定 tie-break）
        seq_io::SeqRecord rec;
    };

    // “更差/更小”判定：长度更短更差；长度相同则 order 更大（更晚出现）更差
    static bool worseThan(const Item& a, const Item& b);

    // 候选是否比当前最差（堆顶）更好
    static bool betterThan(const Item& cand, const Item& worst);

    void siftUp(std::size_t idx);
    void siftDown(std::size_t idx);

private:
    std::size_t k_{0};
    std::uint64_t order_counter_{0};

    // 自实现 min-heap：heap_[0] 永远是“最差”的那条（最短/同长最晚）
    std::vector<Item> heap_;
};
namespace consensus
{
    // 计数：AGCTUN-
    // 若你担心极端情况下每列计数可能超过 2^32，可把 uint32_t 改为 uint64_t
    struct SiteCount
    {
        std::uint32_t a = 0;
        std::uint32_t c = 0;
        std::uint32_t g = 0;
        std::uint32_t t = 0;
        std::uint32_t u = 0;
        std::uint32_t n = 0;
        std::uint32_t dash = 0;

        template <class Archive>
        void serialize(Archive& ar)
        {
            ar(cereal::make_nvp("A", a),
               cereal::make_nvp("C", c),
               cereal::make_nvp("G", g),
               cereal::make_nvp("T", t),
               cereal::make_nvp("U", u),
               cereal::make_nvp("N", n),
               cereal::make_nvp("-", dash));
        }
    };

    struct ConsensusJson
    {
        std::uint64_t num_seqs = 0;
        std::uint64_t aln_len = 0;
        std::vector<SiteCount> counts;

        template <class Archive>
        void serialize(Archive& ar)
        {
            ar(cereal::make_nvp("num_seqs", num_seqs),
               cereal::make_nvp("aln_len", aln_len),
               cereal::make_nvp("counts", counts));
        }
    };

    // 你要求“makeBaseMap 固定在头文件中”：这里用 inline constexpr 表
    // idx: 0 A, 1 C, 2 G, 3 T, 4 U, 5 N, 6 '-'
    inline constexpr std::array<std::uint8_t, 256> k_base_map = []() consteval {
        std::array<std::uint8_t, 256> m{};

        for (std::size_t i = 0; i < m.size(); ++i) m[i] = 5; // default N

        auto set = [&](unsigned char ch, std::uint8_t idx) { m[ch] = idx; };

        set((unsigned char)'A', 0); set((unsigned char)'a', 0);
        set((unsigned char)'C', 1); set((unsigned char)'c', 1);
        set((unsigned char)'G', 2); set((unsigned char)'g', 2);
        set((unsigned char)'T', 3); set((unsigned char)'t', 3);
        set((unsigned char)'U', 4); set((unsigned char)'u', 4);
        set((unsigned char)'N', 5); set((unsigned char)'n', 5);

        // gap：'-' 和 '.' 都当作 gap
        set((unsigned char)'-', 6);
        set((unsigned char)'.', 6);

        return m;
    }();

    inline std::uint8_t mapBase(char ch)
    {
        return k_base_map[(unsigned char)ch];
    }

    // 共识选择：固定 tie-break（A C G T U N -）
    char pickConsensusChar(const SiteCount& sc);

    // 输出 fasta：只写 1 条 >consensus
    void writeConsensusFasta(const FilePath& out_fasta, const std::string& seq);

    // 输出 json：用 cereal JSON archive
    void writeCountsJson(const FilePath& out_json, const ConsensusJson& cj);

    // batch 折中方案：
    // - 1 个线程（master）读取 batch_size 条序列到内存
    // - 其余线程并行按列统计 counts
    //
    // seq_limit: 0 表示全部；>0 表示最多处理前 seq_limit 条
    // threads: 总线程数（<=0 自动取 omp_get_max_threads）
    // batch_size: 每批序列数量（建议 64~1024，按 aln_len 与内存调参）
    //
    // 返回共识序列（长度=对齐长度，可能包含 '-'）

    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                    const FilePath& out_fasta,
                                                    const FilePath& out_json,
                                                    std::uint64_t seq_limit,
                                                    int threads,
                                                    std::size_t batch_size = 512);

} // namespace consensus

#endif // HALIGN4_CONSENSUS_H

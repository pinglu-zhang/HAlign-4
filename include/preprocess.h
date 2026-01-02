#ifndef HALIGN4_PREPROCESS_H
#define HALIGN4_PREPROCESS_H

#include <cstddef>
#include "config.hpp"
#include "utils.h"

// 预处理输入 FASTA，并返回处理的序列数量（total records processed）
uint_t preprocessInputFasta(const std::string input_path, const std::string workdir, const int cons_n = 1000);

// 对外：对输入 fasta 文件进行 MSA（consensus 对齐），将结果写入 output 文件
void alignConsensusSequence(const FilePath& input_file, const FilePath& output_file,
                            const std::string& msa_cmd, const std::string& workdir, int threads);

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

#endif //HALIGN4_PREPROCESS_H

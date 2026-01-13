#ifndef HALIGN4_ALIGN_H
#define HALIGN4_ALIGN_H
#include "utils.h"
#include "mash.h"
#include "seed.h"
#include <filesystem>
#include <string>
#include <vector>
#include <functional>

namespace cigar
{

    // ------------------------------------------------------------------
    // CIGAR 表示与转换
    // ------------------------------------------------------------------

    // 单个 CIGAR 操作的压缩编码（uint32_t）
    // 高位表示操作长度，低 4 位为操作类型编码
    using CigarUnit = uint32_t;

    // 整个 CIGAR 操作序列（压缩形式）
    using Cigar_t = std::vector<CigarUnit>;

    // ------------------------------------------------------------------
    // 函数：cigarToInt
    // 功能：将 CIGAR 操作字符（如 'M'）与其长度编码成一个整数
    // 编码方式：高 28 位表示长度，低 4 位为操作类型（0=M, 1=I, 等）
    // ------------------------------------------------------------------
    CigarUnit cigarToInt(char operation, uint32_t len);

    // ------------------------------------------------------------------
    // 函数：intToCigar
    // 功能：将一个压缩整数还原为操作字符与其长度
    // 示例：0x50 -> ('M', 5)
    // ------------------------------------------------------------------
    void intToCigar(CigarUnit cigar, char& operation, uint32_t& len);
}

namespace align {
    using SeedHit = minimizer::MinimizerHit;
    using SeedHits = std::vector<SeedHit>;
    static constexpr seed::SeedKind kSeedKind = seed::SeedKind::minimizer;

    class RefAligner
    {
        public:
        // 初始化函数
        RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size = 21, int window_size = 10,
                    int sketch_size = 2000, bool noncanonical = true);

        // 说明：
        // - threads <= 0 表示使用 OpenMP 运行时默认线程数（例如由 OMP_NUM_THREADS 控制）
        // - batch_size 用于控制“流式读取”的批次大小，越大吞吐越高但占用内存更多
        void alignQueryToRef(const FilePath& qry_fasta_path, int threads = 0, std::size_t batch_size = 256);


        private:
        // 把“相似度计算 +（占位）比对 + 写出”抽象成成员函数，便于后续替换实现而不影响并行框架。
        // 约束：该函数不应修改共享 reference 数据结构（除非自行加锁）。
        void alignOneQueryToRef(const seq_io::SeqRecord& q, seq_io::SeqWriter& out) const;

        FilePath work_dir;
        seq_io::SeqRecords ref_sequences;
        mash::Sketches ref_sketch;
        std::vector<SeedHits> ref_minimizers;

        int kmer_size = 21;
        int window_size = 10;
        int sketch_size = 2000;
        int random_seed = 42;

        bool keep_first_length = false;
        bool keep_all_length = false;
        bool noncanonical = true;


    };

} // namespace align

#endif //HALIGN4_ALIGN_H

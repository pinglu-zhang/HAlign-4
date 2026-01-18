#ifndef HALIGN4_ALIGN_H
#define HALIGN4_ALIGN_H
#include "utils.h"
#include "mash.h"
#include "seed.h"
#include <filesystem>
#include <string>
#include <vector>
#include <functional>
#include "config.hpp"  // 包含 Options 结构体的完整定义

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

    // ------------------------------------------------------------------
    // 函数：hasInsertion
    // 功能：检测 CIGAR 序列中是否存在插入操作（'I'）
    // 参数：cigar - CIGAR 操作序列（压缩形式）
    // 返回：true 表示存在至少一个插入操作，false 表示不存在
    // 性能：O(N)，N 为 CIGAR 操作数量；短路优化，找到第一个 'I' 即返回
    // ------------------------------------------------------------------
    bool hasInsertion(const Cigar_t& cigar);

    // ------------------------------------------------------------------
    // 函数：cigarToString
    // 功能：将 CIGAR 从压缩格式转换为 SAM 标准字符串格式
    // 参数：cigar - CIGAR 操作序列（压缩形式）
    // 返回：SAM 格式的 CIGAR 字符串，例如 "100M5I95M"
    // 性能优化：
    //   1. 预分配字符串空间（cigar.size() * 5），减少内存重新分配
    //   2. 复杂度 O(N)，N 为 CIGAR 操作数量
    //   3. 使用 std::to_string 进行整数到字符串转换（编译器优化）
    // 示例：
    //   输入：[cigarToInt('M', 100), cigarToInt('I', 5), cigarToInt('M', 95)]
    //   输出："100M5I95M"
    // ------------------------------------------------------------------
    std::string cigarToString(const Cigar_t& cigar);
}

namespace align {
    using SeedHit = minimizer::MinimizerHit;
    using SeedHits = std::vector<SeedHit>;
    static constexpr seed::SeedKind kSeedKind = seed::SeedKind::minimizer;

    // ------------------------------------------------------------------
    // DNA 字符到索引的映射表（编译期常量）
    // ------------------------------------------------------------------
    // 说明：将 DNA 字符（A/C/G/T/N，大小写不敏感）映射到 0-4 的索引
    // - 'A'/'a' -> 0
    // - 'C'/'c' -> 1
    // - 'G'/'g' -> 2
    // - 'T'/'t' -> 3
    // - 'N'/'n' 或其他 -> 4（未知碱基）
    // 用途：KSW2 需要将序列编码为整数数组才能进行比对
    // ------------------------------------------------------------------
    static constexpr uint8_t ScoreChar2Idx[256] = {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 0-15
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 16-31
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 32-47 (空格等)
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 48-63 (数字)
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,  // 64-79  (@,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O)
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,  // 80-95  (P,Q,R,S,T,U,V,W,X,Y,Z,...)
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,  // 96-111 (`,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,  // 112-127(p,q,r,s,t,u,v,w,x,y,z,...)
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 128-143
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 144-159
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 160-175
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 176-191
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 192-207
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 208-223
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 224-239
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4   // 240-255
    };

    struct KSW2AlignConfig {
        const int8_t* mat;                  // 一维替换矩阵 (flattened 5x5)
        int alphabet_size;                 // 通常为 5
        int gap_open;                      // gap open penalty (positive)
        int gap_extend;                    // gap extend penalty (positive)
        int end_bonus;                     // 末端奖励分
        int zdrop = 100;                   // Z-drop 剪枝参数
        int band_width = -1;               // -1 表示全矩阵
        int flag = 0;      // 默认使用全替换矩阵
    };;

    // ------------------------------------------------------------------
    // DNA5 替换矩阵（编译期常量，用于 KSW2 全局/延伸比对）
    // ------------------------------------------------------------------
    // 说明：
    // 1. 矩阵维度：5x5，索引 0-3 为 A/C/G/T，索引 4 为 N（未知碱基/通配符）
    // 2. 一维展开：mat[i*5+j] 对应二维 mat[i][j]
    // 3. 评分规则：
    //    - 精确匹配（A-A, C-C, G-G, T-T）：+5 分（奖励正确配对）
    //    - 错配（A-C, A-G 等）：-4 分（惩罚碱基不匹配）
    //    - 涉及 N（未知碱基）：0 分（既不奖励也不惩罚）
    // 4. 参数平衡：
    //    - match/mismatch 比例约 5:4，与常见 Blast 参数（+1/-1 或 +5/-4）一致
    //    - 配合 gap_open=6, gap_extend=2 时，indel 代价约为 6+2k
    //    - 单个错配代价 9（从 +5 降到 -4），单个 1bp gap 代价 8，略偏好 gap
    // 5. 使用要求：
    //    - **必须**配合 KSW_EZ_GENERIC_SC flag 使用（启用完整替换矩阵）
    //    - 不带该 flag 时，KSW2 只支持简单 match/mismatch 模式，会忽略此矩阵
    // ------------------------------------------------------------------
    static constexpr int8_t dna5_simd_mat[25] = {
        // A   C   G   T   N
           5, -4, -4, -4,  0,  // A (i=0)
          -4,  5, -4, -4,  0,  // C (i=1)
          -4, -4,  5, -4,  0,  // G (i=2)
          -4, -4, -4,  5,  0,  // T (i=3)
           0,  0,  0,  0,  0   // N (i=4)
    };

    //------------------------------------------- 带宽估计
    inline int auto_band(int qlen, int tlen,
        double indel_rate = 0.1,
        int    margin = 200)           // 多一点保险
    {
        return margin + static_cast<int>(indel_rate * (qlen + tlen / 2));
    }

    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query);

    cigar::Cigar_t extendAlignKSW2(const std::string& ref, const std::string& query, int zdrop = 200);

    cigar::Cigar_t globalAlignWFA2(const std::string& ref, const std::string& query);

    // cigar::Cigar_t extendAlignWFA2(const std::string& ref, const std::string& query, int zdrop = 200);

    class RefAligner
    {
        public:
        // ------------------------------------------------------------------
        // 构造函数1：直接传入参数初始化
        // 参数：
        //   - consensus_string: 可选的共识序列字符串，如果提供则直接使用并保存到成员变量
        //   - keep_first_length: 是否保持第一条序列的长度不变
        //   - keep_all_length: 是否保持所有序列的长度不变
        // ------------------------------------------------------------------
        RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size = 21, int window_size = 10,
                    int sketch_size = 2000, bool noncanonical = true, std::string consensus_string = "",
                    bool keep_first_length = false, bool keep_all_length = false);

        // ------------------------------------------------------------------
        // 构造函数2：基于 Options 结构体初始化
        // 参数：
        //   - opt: Options 结构体（包含所有命令行参数和配置）
        //   - ref_fasta_path: 参考序列文件路径
        //   - consensus_string: 可选的共识序列字符串，如果提供则直接使用并保存到成员变量
        // 说明：
        //   - keep_first_length 和 keep_all_length 会从 opt 中自动提取
        // ------------------------------------------------------------------
        RefAligner(const Options& opt, const FilePath& ref_fasta_path, std::string consensus_string = "");

        // 说明：
        // - threads <= 0 表示使用 OpenMP 运行时默认线程数（例如由 OMP_NUM_THREADS 控制）
        // - batch_size 用于控制“流式读取”的批次大小，越大吞吐越高但占用内存更多
        void alignQueryToRef(const FilePath& qry_fasta_path, int threads = 0, std::size_t batch_size = 5120);

        void mergeAlignedResults(const FilePath& aligned_consensus_path);


        private:
        // ------------------------------------------------------------------
        // 函数：alignOneQueryToRef
        // 功能：对单个 query 执行比对并写入 SAM 文件
        // 参数：
        //   - q: 待比对的查询序列
        //   - out: 普通输出文件的 writer（无插入或二次比对后无插入的序列）
        //   - out_insertion: 插入序列输出文件的 writer（二次比对后仍有插入的序列）
        // 说明：
        //   - 该函数不修改共享 reference 数据结构（线程安全）
        //   - out 和 out_insertion 由调用者管理生命周期（每线程独立）
        // ------------------------------------------------------------------
        void alignOneQueryToRef(const seq_io::SeqRecord& q,
                               seq_io::SeqWriter& out,
                               seq_io::SeqWriter& out_insertion) const;

        // 辅助函数：根据 keep_first_length 标志选择参考序列名称
        std::string_view getRefNameForRecheck() const;

        // 辅助函数：执行二次比对（用于插入判断）
        cigar::Cigar_t performRecheckAlignment(const seq_io::SeqRecord& q) const;

        // 辅助函数：写入SAM记录（选择正确的参考名称和输出文件）
        void writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                           std::string_view ref_name, seq_io::SeqWriter& out) const;

        FilePath work_dir;
        seq_io::SeqRecords ref_sequences;
        mash::Sketches ref_sketch;
        std::vector<SeedHits> ref_minimizers;

        std::string consensus_seq;  // 共识序列（作为成员变量存储）

        int kmer_size = 21;
        int window_size = 10;
        int sketch_size = 2000;
        int random_seed = 42;

        bool keep_first_length = false;
        bool keep_all_length = false;
        bool noncanonical = true;

        std::vector<FilePath> outs_path;
        std::vector<FilePath> outs_with_insertion_path;



    };

} // namespace align

#endif //HALIGN4_ALIGN_H

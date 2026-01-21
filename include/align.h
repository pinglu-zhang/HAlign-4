// ==================================================================
// align.h - HAlign-4 序列比对模块核心头文件
// ==================================================================
// 功能概述：
// 本文件定义了 HAlign-4 项目的序列比对核心接口，包括：
// 1. CIGAR 操作的表示、转换与处理（cigar 命名空间）
// 2. 序列比对算法接口（KSW2、WFA2）
// 3. 参考序列比对器（RefAligner 类）：多序列比对（MSA）的核心组件

// ==================================================================

#ifndef HALIGN4_ALIGN_H
#define HALIGN4_ALIGN_H
#include "utils.h"
#include "mash.h"
#include "seed.h"

#include <unordered_map>
#include <filesystem>
#include <string>
#include <vector>
#include <functional>
#include "config.hpp"  // 包含 Options 结构体的完整定义
#include "consensus.h"
#include "preprocess.h"

// ==================================================================
// cigar 命名空间：CIGAR 操作的表示、解析与序列对齐
// ==================================================================
// 功能概述：
// 提供 CIGAR（Compact Idiosyncratic Gapped Alignment Report）操作的完整支持
// 包括：编码/解码、字符串转换、序列对齐（插入 gap）
//
// 核心设计：
// 1. **压缩编码**：使用 uint32_t 存储单个 CIGAR 操作（长度+操作符）
//    - 高 28 位：操作长度（0 到 2^28-1，约 2.68 亿）
//    - 低 4 位：操作类型编码（0-8，对应 M/I/D/N/S/H/P/=/X）
// 2. **零拷贝解析**：stringToCigar 直接从字符串解析为压缩格式
// 3. **原地对齐**：padQueryToRefByCigar 原地修改序列，避免内存分配
//
// 支持的 CIGAR 操作（SAM 标准）：
// - M (match/mismatch):    query 和 ref 都消耗，可能匹配或错配
// - I (insertion):         query 相对 ref 的插入，只消耗 query
// - D (deletion):          query 相对 ref 的缺失，只消耗 ref
// - N (skipped region):    ref 上跳过的区域（如剪接位点），只消耗 ref
// - S (soft clip):         query 中存在但未比对的部分，只消耗 query
// - H (hard clip):         query 中已被移除的部分，不消耗任何序列
// - P (padding):           silent deletion，不消耗序列（用于 MSA）
// - = (exact match):       精确匹配（扩展 CIGAR）
// - X (mismatch):          错配（扩展 CIGAR）
// ==================================================================
namespace cigar
{

    // ------------------------------------------------------------------
    // CIGAR 表示与转换
    // ------------------------------------------------------------------

    // 单个 CIGAR 操作的压缩编码（uint32_t）
    // 编码方式：
    //   - 高 28 位：操作长度（len），范围 0 到 2^28-1（约 268,435,455）
    //   - 低 4 位：操作类型（op），编码如下：
    //       0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
    //   - 示例：100M -> (100 << 4) | 0 = 0x640
    //           5I   -> (5 << 4) | 1 = 0x51
    // 内存占用：每个操作 4 字节（vs SAM 字符串每操作 3-10 字节）
    using CigarUnit = uint32_t;

    // 整个 CIGAR 操作序列（压缩形式）
    // 示例："100M5I95M" -> [0x640, 0x51, 0x5F0]
    // 性能：vector 存储，支持快速索引和迭代
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

    // ------------------------------------------------------------------
    // 函数：stringToCigar
    // 功能：将 SAM 格式的 CIGAR 字符串解析为压缩格式的 Cigar_t
    // 参数：cigar_str - SAM 格式的 CIGAR 字符串，例如 "100M5I95M"
    // 返回：Cigar_t（压缩的整数向量），每个元素编码一个操作（长度+操作符）
    // 异常：
    //   - 如果字符串格式无效（例如缺少数字、未知操作符、长度为0等），抛出 std::runtime_error
    //   - 如果长度超过 2^28-1，抛出 std::runtime_error
    // 性能：
    //   1. 复杂度 O(N)，N 为字符串长度
    //   2. 预分配结果向量空间，减少内存重新分配
    //   3. 使用原地解析，避免不必要的字符串拷贝
    // 示例：
    //   输入："100M5I95M"
    //   输出：[cigarToInt('M', 100), cigarToInt('I', 5), cigarToInt('M', 95)]
    //   输入："*"（未知 CIGAR）
    //   输出：空向量
    // 正确性说明：
    //   - 本函数是 cigarToString 的逆操作，满足：stringToCigar(cigarToString(c)) == c
    //   - 支持 SAM 标准的所有 CIGAR 操作符：M, I, D, N, S, H, P, =, X
    //   - 自动跳过空白字符（空格、制表符、换行符），容错性强
    // ------------------------------------------------------------------
    Cigar_t stringToCigar(const std::string& cigar_str);

    // ------------------------------------------------------------------
    // 函数：padQueryToRefByCigar
    // 功能：根据 CIGAR 操作对 query 序列插入 gap 字符（'-'），使其与参考序列对齐
    //
    // ------------------------------------------------------------------
    // 函数：padQueryToRefByCigar
    // 功能：根据 CIGAR 操作对 query 序列插入 gap 字符（'-'），使其与参考序列对齐
    //
    // **关键特性：保留原有 gap 字符**
    // - 输入 query 中已存在的 gap 字符（'-'）会被当作普通碱基处理（不删除、不特殊对待）
    // - 只根据 CIGAR 操作在适当位置插入新的 gap
    // - 这对于多次比对场景（MSA）或处理已经部分对齐的序列非常有用
    //
    // 参数：
    //   - query: 待插空的查询序列（会被原地修改）
    //            注意：如果输入序列已包含 gap 字符（'-'），会被保留并当作普通字符处理
    //   - cigar: CIGAR 操作序列（压缩形式）
    //
    // CIGAR 操作语义：
    //   - M/=/X (match/mismatch): query 和 ref 都消耗，拷贝原字符（**包括原有 gap**）
    //   - I (insertion): query 相对 ref 的插入，只消耗 query，拷贝原字符（**包括原有 gap**）
    //   - D (deletion): query 相对 ref 的缺失，只消耗 ref，插入新 gap
    //   - S (soft clip): query 中存在但未比对的部分，拷贝原字符（**包括原有 gap**）
    //   - H (hard clip): 已从 query 中移除，不处理
    //
    // 算法说明：
    //   1. 预计算对齐后的总长度（原长度 + 所有 D 操作的长度）
    //   2. 从后往前构建对齐序列，避免频繁插入导致的 O(N²) 复杂度
    //   3. 原有的 gap 字符（'-'）会被当作普通字符处理，与 A/C/G/T 等碱基无区别
    //
    // 性能：
    //   - 时间复杂度：O(M + N)，M 为 CIGAR 操作数，N 为 query 长度（**包括原有 gap**）
    //   - 空间复杂度：O(N)，需要临时存储原始 query
    //   - 内存分配次数：2 次（移动原序列 + 分配新空间）
    //   - 性能优化：移动语义、从后往前填充、单次内存分配、cache 友好
    //
    // 示例1（不含原有 gap，标准场景）：
    //   query = "ACGT", cigar = "2M1D2M"
    //   结果: query = "AC-GT"
    //
    // 示例2（包含原有 gap，MSA 场景）：
    //   query = "A-CG-T", cigar = "2M1D2M2M"
    //   说明：2M 消耗 "A-"（包括原有 gap），1D 插入新 gap，2M 消耗 "CG"，2M 消耗 "-T"
    //   结果: query = "A--CG-T"（长度 7 = 6 + 1个新gap）
    //
    // 异常：
    //   - 如果 CIGAR 操作超出 query 长度，会触发断言（Debug 模式）
    //   - 不支持的 CIGAR 操作符会抛出 std::runtime_error
    //
    // 使用场景：
    //   - 标准场景：输入序列不含 gap，根据 CIGAR 插入 gap
    //   - MSA 场景：输入序列已有 gap（来自之前的比对），需要根据新 CIGAR 调整对齐
    //   - 迭代比对：多次比对同一序列，每次都保留之前的 gap 并插入新 gap
    // ------------------------------------------------------------------
    void padQueryToRefByCigar(std::string& query, const Cigar_t& cigar);
}

// ==================================================================
// align 命名空间：序列比对算法与参考序列比对器
// ==================================================================
// 功能概述：
// 1. 提供多种序列比对算法接口（KSW2、WFA2）
// 2. RefAligner 类：高性能多序列比对（MSA）引擎
//    - 支持多线程并行处理
//    - 使用 MinHash + Minimizer 加速相似序列查找
//    - 支持插入序列的二次比对和 MSA 整合
// 3. 评分矩阵和配置结构体
//
// 设计理念：
// - **模块化比对器**：KSW2 和 WFA2 提供统一的 CIGAR 输出接口
// - **流式处理**：RefAligner 支持批量读取和处理，内存占用可控
// - **灵活配置**：通过 Options 结构体统一管理所有参数
// ==================================================================
namespace align {
    // ------------------------------------------------------------------
    // 类型别名：种子（Seed）与种子命中（Seed Hit）
    // ------------------------------------------------------------------
    // 说明：
    // - 种子（Seed）是序列中的短 k-mer 或 minimizer，用于快速定位潜在的同源区域
    // - 种子命中（Seed Hit）记录种子在参考序列和查询序列中的位置
    // - 当前实现使用 minimizer 作为种子策略（高效、内存友好）
    // ------------------------------------------------------------------
    using SeedHit = minimizer::MinimizerHit;   // 单个种子命中：(ref_pos, query_pos, hash)
    using SeedHits = std::vector<SeedHit>;     // 种子命中列表（用于锚点定位）
    static constexpr seed::SeedKind kSeedKind = seed::SeedKind::minimizer;  // 种子策略：minimizer

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

    // ==================================================================
    // KSW2AlignConfig：KSW2 比对算法的配置结构体
    // ==================================================================
    // 说明：
    // KSW2 是一个高性能的序列比对算法（SSE/AVX 加速），支持全局、局部和延伸比对
    // 本结构体封装了 KSW2 所需的所有参数
    //
    // 参数详解：
    // ------------------------------------------------------------------
    // 1. mat（替换矩阵）：
    //    - 类型：const int8_t* 指向一维数组（大小 alphabet_size^2）
    //    - 格式：一维展开的二维矩阵，mat[i*alphabet_size+j] 对应 mat[i][j]
    //    - 用途：评分 query[i] 与 ref[j] 的匹配/错配代价
    //    - 示例：dna5_simd_mat（5x5 矩阵，A/C/G/T/N）
    //
    // 2. alphabet_size（字母表大小）：
    //    - DNA 通常为 5（A/C/G/T/N）
    //    - 蛋白质为 20+1（20 种氨基酸 + 未知）
    //
    // 3. gap_open（gap 开启惩罚）：
    //    - 正整数，值越大越不倾向于引入新的 gap
    //    - 建议范围：4-10（DNA）、8-12（蛋白质）
    //    - 示例：gap_open=6 表示开启一个新 gap 的代价为 6 分
    //
    // 4. gap_extend（gap 延伸惩罚）：
    //    - 正整数，值越大越不倾向于延长现有 gap
    //    - 通常 gap_extend < gap_open（倾向于少量长 gap 而非多个短 gap）
    //    - 建议范围：1-3（DNA）、1-2（蛋白质）
    //    - 示例：gap_extend=2 表示延长 gap 每多 1bp 代价增加 2 分
    //
    // 5. end_bonus（末端奖励分）：
    //    - 用于局部比对或延伸比对，奖励到达序列末端
    //    - 全局比对通常设为 0
    //    - 局部比对可设为正值（如 5），鼓励比对延伸到序列末端
    //
    // 6. zdrop（Z-drop 剪枝参数）：
    //    - 默认值：100
    //    - 用于延伸比对（extension alignment）
    //    - 当比对得分下降超过 zdrop 时停止延伸（剪枝策略）
    //    - 值越大，比对越完整但计算越慢；值越小，速度越快但可能过早终止
    //
    // 7. band_width（带宽）：
    //    - 默认值：-1（表示使用全矩阵，无带宽限制）
    //    - 正值：限制 DP 矩阵的对角带宽度（加速长序列比对）
    //    - 建议：对于高相似度序列（如同源基因），可设为序列长度的 10%-20%
    //    - 示例：band_width=100 表示只计算对角线附近 ±100 的区域
    //
    // 8. flag（比对模式标志位）：
    //    - 默认值：0（使用完整替换矩阵）
    //    - 常用标志（按位或组合）：
    //        * 0：全局比对，使用完整替换矩阵（默认）
    //        * KSW_EZ_SCORE_ONLY：只计算得分，不生成 CIGAR
    //        * KSW_EZ_RIGHT：右对齐模式（用于延伸比对）
    //        * KSW_EZ_GENERIC_SC：必须设置此标志才能使用自定义替换矩阵
    //        * KSW_EZ_APPROX_MAX：近似最大得分（加速）
    //    - 注意：使用 dna5_simd_mat 时必须包含 KSW_EZ_GENERIC_SC 标志
    //
    // 使用示例：
    // ------------------------------------------------------------------
    // // 全局比对配置（DNA，5x5 矩阵）
    // KSW2AlignConfig cfg_global = {
    //     .mat = dna5_simd_mat,
    //     .alphabet_size = 5,
    //     .gap_open = 6,
    //     .gap_extend = 2,
    //     .end_bonus = 0,
    //     .zdrop = 100,
    //     .band_width = -1,  // 全矩阵
    //     .flag = KSW_EZ_GENERIC_SC
    // };
    //
    // // 延伸比对配置（带宽限制，高 zdrop）
    // KSW2AlignConfig cfg_extend = {
    //     .mat = dna5_simd_mat,
    //     .alphabet_size = 5,
    //     .gap_open = 6,
    //     .gap_extend = 2,
    //     .end_bonus = 5,       // 奖励到达末端
    //     .zdrop = 200,          // 更宽容的剪枝
    //     .band_width = 100,     // 限制带宽加速
    //     .flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT
    // };
    // ==================================================================
    struct KSW2AlignConfig {
        const int8_t* mat;                  // 一维替换矩阵 (flattened 5x5)
        int alphabet_size;                 // 通常为 5
        int gap_open;                      // gap open penalty (positive)
        int gap_extend;                    // gap extend penalty (positive)
        int end_bonus;                     // 末端奖励分
        int zdrop = 100;                   // Z-drop 剪枝参数
        int band_width = -1;               // -1 表示全矩阵
        int flag = 0;      // 默认使用全替换矩阵
    };

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

    // ==================================================================
    // 函数：auto_band - 自动估算 KSW2 带宽参数
    // ==================================================================
    // 功能：
    // 根据序列长度和预期的 indel（插入/删除）比率，自动估算合适的带宽
    // 用于加速长序列的 KSW2 比对（限制 DP 矩阵的计算范围）
    //
    // 参数：
    // @param qlen - query 序列长度
    // @param tlen - target（参考）序列长度
    // @param indel_rate - 预期的 indel 比率（默认 0.1 = 10%）
    //                     - 0.05：高相似度序列（如同种基因组内）
    //                     - 0.10：中等相似度（如近缘物种）
    //                     - 0.20：低相似度（如远缘物种）
    // @param margin - 安全边距（默认 200）
    //                 - 额外的缓冲区，防止真实 indel 超出带宽导致错误
    //
    // 返回：
    // 建议的带宽值（整数）
    //
    // 算法：
    // band_width = margin + indel_rate * (qlen + tlen / 2)
    // - 基本思想：带宽应覆盖"平均序列长度 * indel 比率"的区域
    // - tlen / 2：考虑 query 相对 target 的偏移（取中位数）
    // - margin：额外的安全边距，防止极端情况（如局部高 indel 区域）
    //
    // 使用建议：
    // - 如果序列长度差异很大（如 qlen << tlen），应使用 min(qlen, tlen)
    // - 如果已知序列高度相似，可减小 indel_rate（如 0.05）或 margin（如 100）
    // - 对于全基因组比对，建议设置上限（如 max_band = 5000）防止内存溢出
    //
    // 性能影响：
    // - 带宽越小，计算越快，但可能遗漏真实比对
    // - 带宽越大，结果越准确，但计算时间和内存消耗呈线性增长
    // - 经验值：带宽 = 序列长度的 10%-20% 时，速度与准确性平衡较好
    //
    // 示例：
    // auto_band(1000, 1000, 0.1, 200) = 200 + 0.1 * (1000 + 500) = 350
    // auto_band(10000, 10000, 0.05, 100) = 100 + 0.05 * (10000 + 5000) = 850
    // ==================================================================
    //------------------------------------------- 带宽估计
    inline int auto_band(int qlen, int tlen,
        double indel_rate = 0.1,
        int    margin = 200)           // 多一点保险
    {
        return margin + static_cast<int>(indel_rate * (qlen + tlen / 2));
    }

    // ==================================================================
    // 序列比对算法接口
    // ==================================================================
    // 说明：
    // 以下函数提供统一的序列比对接口，均返回 CIGAR 格式的比对结果
    // 输入序列可包含 A/C/G/T/N（大小写不敏感），其他字符视为 N
    // ==================================================================

    // ------------------------------------------------------------------
    // 函数：globalAlignKSW2 - KSW2 全局比对
    // ------------------------------------------------------------------
    // 功能：
    // 使用 KSW2 算法执行全局序列比对（Needleman-Wunsch 模式）
    // 适用场景：两条序列完整比对，考虑两端的所有碱基
    //
    // 参数：
    // @param ref - 参考序列（字符串，A/C/G/T/N）
    // @param query - 查询序列（字符串，A/C/G/T/N）
    //
    // 返回：
    // CIGAR 操作序列（cigar::Cigar_t），描述 query 如何比对到 ref
    //
    // 配置：
    // - 使用 dna5_simd_mat（5x5 替换矩阵）
    // - gap_open=6, gap_extend=2
    // - 全矩阵（无带宽限制）
    //
    // 性能：
    // - 时间复杂度：O(M * N)，M 和 N 分别为两序列长度
    // - 空间复杂度：O(M * N)（DP 矩阵）
    // - 加速：SSE/AVX 指令集加速（比普通 DP 快 4-8 倍）
    //
    // 示例：
    // globalAlignKSW2("ACGT", "AGT") -> "1M1D2M"（query 缺失了 ref 的 C）
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query);

    // ------------------------------------------------------------------
    // 函数：extendAlignKSW2 - KSW2 延伸比对
    // ------------------------------------------------------------------
    // 功能：
    // 使用 KSW2 算法执行延伸比对（extension alignment）
    // 适用场景：从种子（seed）位置向两侧延伸，直到得分下降过多或到达序列末端
    //
    // 参数：
    // @param ref - 参考序列（字符串，A/C/G/T/N）
    // @param query - 查询序列（字符串，A/C/G/T/N）
    // @param zdrop - Z-drop 剪枝阈值（默认 200）
    //                - 当比对得分下降超过 zdrop 时停止延伸
    //                - 值越大，延伸越完整但计算越慢
    //
    // 返回：
    // CIGAR 操作序列（cigar::Cigar_t），描述延伸区域的比对
    //
    // 配置：
    // - 使用 dna5_simd_mat（5x5 替换矩阵）
    // - gap_open=6, gap_extend=2
    // - end_bonus=5（奖励到达序列末端）
    // - 自动估算带宽（auto_band）
    //
    // 性能：
    // - 通常比全局比对快（因为剪枝和带宽限制）
    // - 适合长序列的局部比对（如基因组拼接中的重叠区域检测）
    //
    // 使用建议：
    // - 如果序列高度相似，可降低 zdrop（如 100）提高速度
    // - 如果需要完整比对，应使用 globalAlignKSW2 而非 extendAlignKSW2
    // ------------------------------------------------------------------
    cigar::Cigar_t extendAlignKSW2(const std::string& ref, const std::string& query, int zdrop = 200);

    // ------------------------------------------------------------------
    // 函数：globalAlignWFA2 - WFA2 全局比对
    // ------------------------------------------------------------------
    // 功能：
    // 使用 WFA2（Wavefront Alignment）算法执行全局序列比对
    // WFA2 特点：对于高相似度序列（indel 少），速度显著快于传统 DP 算法
    //
    // 参数：
    // @param ref - 参考序列（字符串，A/C/G/T/N）
    // @param query - 查询序列（字符串，A/C/G/T/N）
    //
    // 返回：
    // CIGAR 操作序列（cigar::Cigar_t），描述 query 如何比对到 ref
    //
    // 性能：
    // - 时间复杂度：O(s * N)，s 为编辑距离，N 为序列长度
    //   - 对于高相似度序列（s << N），WFA2 远快于 KSW2
    //   - 对于低相似度序列（s ≈ N），WFA2 可能比 KSW2 慢
    // - 空间复杂度：O(s^2)（波前矩阵）
    //
    // 使用建议：
    // - 推荐用于高相似度序列（如同种基因组、近缘物种）
    // - 对于低相似度或未知相似度的序列，建议先用 MinHash 估算相似度
    // - 如果估算相似度 > 90%，使用 WFA2；否则使用 KSW2
    //
    // 示例：
    // globalAlignWFA2("ACGT", "AGT") -> "1M1D2M"
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignWFA2(const std::string& ref, const std::string& query);

    // cigar::Cigar_t extendAlignWFA2(const std::string& ref, const std::string& query, int zdrop = 200);

    // ==================================================================
    // RefAligner 类：高性能参考序列比对器（多序列比对 MSA 引擎）
    // ==================================================================
    // 功能概述：
    // RefAligner 是 HAlign-4 的核心组件，负责将大量 query 序列比对到参考序列集合
    // 并生成最终的多序列比对（MSA）结果
    //
    // 核心特性：
    // 1. **多线程并行**：使用 OpenMP 并行处理 query 序列
    //    - 每个线程独立处理一批 query，写入独立的 SAM 文件
    //    - 无锁设计：共享数据（参考序列、sketch、minimizer）只读，线程安全
    //
    // 2. **快速相似度估计**：使用 MinHash sketch 快速筛选候选参考序列
    //    - 避免对每个 query 执行 O(N) 次全局比对（N = 参考序列数量）
    //    - 只对 Jaccard 相似度高的参考序列执行精确比对
    //
    // 3. **种子加速比对**：使用 Minimizer 定位同源区域
    //    - 快速找到潜在的比对锚点（anchor），减少全局比对的范围
    //    - 对于长序列（如全基因组），显著提高比对速度
    //
    // 4. **插入序列处理**：对无法比对到任何参考序列的 query 进行特殊处理
    //    - 第一次比对失败的序列记录到单独的 SAM 文件
    //    - 后续调用外部 MSA 工具（如 MAFFT）对插入序列进行多序列比对
    //    - 将插入序列的 MSA 结果整合到最终输出
    //
    // 5. **坐标系统一**：通过多级 CIGAR 投影，将所有序列投影到统一坐标系
    //    - query → ref：第一级投影（SAM CIGAR）
    //    - ref → consensus：第二级投影（参考序列 MSA 的 CIGAR）
    //    - consensus → insertion MSA：第三级投影（整合插入序列）
    //
    // 工作流程：
    // 1. 初始化：读取参考序列，生成 sketch 和 minimizer 索引
    // 2. alignQueryToRef：并行比对 query 序列，输出 SAM 文件
    // 3. mergeAlignedResults：合并所有 SAM 文件，生成最终 MSA FASTA 文件
    //
    // 性能优化：
    // - 批量读取：避免频繁 I/O
    // - 流式处理：内存占用与数据量无关（仅与批次大小相关）
    // - 缓冲写入：SeqWriter 使用大缓冲区，减少系统调用
    // - 编译期常量：评分矩阵、字符映射表等使用 constexpr
    //
    // 使用示例：
    // ------------------------------------------------------------------
    // // 1. 创建 RefAligner
    // Options opt = parseCommandLine(argc, argv);
    // RefAligner aligner(opt, "refs.fasta");
    //
    // // 2. 比对 query 序列
    // aligner.alignQueryToRef("queries.fasta", 5120);
    //
    // // 3. 合并结果生成 MSA
    // aligner.mergeAlignedResults("consensus_aligned.fasta", "mafft --auto");
    // ------------------------------------------------------------------
    // ==================================================================
    class RefAligner
    {
        public:
        // ------------------------------------------------------------------
        // 构造函数1：直接传入参数初始化
        // 参数：
        //   - threads: 线程数（用于共识序列生成）
        //   - msa_cmd: MSA 命令模板字符串（用于共识序列比对）
        //   - keep_first_length: 是否保持第一条序列的长度不变
        //   - keep_all_length: 是否保持所有序列的长度不变
        // 说明：
        //   1. 如果 keep_first_length == true，使用 ref_sequences[0] 作为参考
        //   2. 否则，调用 MSA 生成共识序列作为参考
        //   3. threads 和 msa_cmd 参数仅在生成共识序列时使用
        // ------------------------------------------------------------------
        RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path,
                   int kmer_size = 21, int window_size = 10,
                   int sketch_size = 2000, bool noncanonical = true,
                   int threads = 1, std::string msa_cmd = "",
                   bool keep_first_length = false, bool keep_all_length = false);

        // ------------------------------------------------------------------
        // 构造函数2：基于 Options 结构体初始化（推荐方式）
        // ------------------------------------------------------------------
        // 参数：
        //   @param opt - Options 结构体（包含所有命令行参数和配置）
        //                - 自动提取：threads, kmer_size, window_size, sketch_size,
        //                  noncanonical, msa_cmd, keep_first_length, keep_all_length
        //   @param ref_fasta_path - 参考序列 FASTA 文件路径
        //
        // 工作流程：
        //   1. 读取参考序列文件到 ref_sequences
        //   2. 为每个参考序列生成 MinHash sketch（用于快速相似度估计）
        //   3. 为每个参考序列生成 Minimizer 索引（用于快速定位同源区域）
        //   4. 如果 opt.consensus_string 非空，直接使用；否则生成共识序列
        //
        // 说明：
        //   - 本构造函数从 opt 中提取所有必要参数，无需手动指定
        //   - 推荐在实际应用中使用本构造函数（避免参数传递错误）
        //   - keep_first_length 和 keep_all_length 会从 opt 中自动提取
        //
        // 性能：
        //   - 时间复杂度：O(N * L)，N = 参考序列数量，L = 平均序列长度
        //   - 内存占用：O(N * (L + S + M))，S = sketch 大小，M = minimizer 数量
        //   - 建议：对于大规模参考序列（如全基因组），考虑减小 sketch_size
        //
        // 异常：
        //   - 如果参考序列文件不存在或格式错误，抛出 std::runtime_error
        //   - 如果内存不足，抛出 std::bad_alloc
        // ------------------------------------------------------------------
        RefAligner(const Options& opt, const FilePath& ref_fasta_path);

        // ==================================================================
        // 公共方法：比对与合并
        // ==================================================================

        // ------------------------------------------------------------------
        // 方法：alignQueryToRef - 并行比对 query 序列到参考序列
        // ------------------------------------------------------------------
        // 功能：
        // 读取 query FASTA 文件，将每条序列比对到最相似的参考序列
        // 输出多个 SAM 文件（每个线程一个），记录比对结果
        //
        // 参数：
        //   @param qry_fasta_path - query 序列 FASTA 文件路径
        //   @param batch_size - 批量读取的序列数量（默认 5120）
        //                       - 值越大，吞吐越高，但内存占用越多
        //                       - 建议：普通机器 1024-5120，高性能服务器 10000+
        //
        // 工作流程：
        //   1. 分批读取 query 序列（每批 batch_size 条）
        //   2. OpenMP 并行处理每批序列：
        //      a. 使用 MinHash 筛选候选参考序列（Jaccard 相似度）
        //      b. 使用 Minimizer 定位同源区域
        //      c. 执行精确比对（KSW2 或 WFA2）
        //      d. 检测 CIGAR 中是否有插入（'I'）
        //      e. 如果有插入，写入 out_insertion；否则写入 out
        //   3. 每个线程写入独立的 SAM 文件（避免锁竞争）
        // ------------------------------------------------------------------
        // 说明：
        // - threads <= 0 表示使用 OpenMP 运行时默认线程数（例如由 OMP_NUM_THREADS 控制）
        // - batch_size 用于控制"流式读取"的批次大小，越大吞吐越高但占用内存更多
        void alignQueryToRef(const FilePath& qry_fasta_path, std::size_t batch_size = 5120);

        // ------------------------------------------------------------------
        // 方法：mergeAlignedResults - 合并所有比对结果生成最终 MSA
        // ------------------------------------------------------------------
        // 功能：
        // 将多个线程产生的 SAM 文件合并为一个统一的多序列比对（MSA）FASTA 文件
        // 包括处理插入序列的二次比对和坐标系统一
        //
        // 参数：
        //   @param aligned_consensus_path - 已对齐的共识序列文件路径
        //                                   - 包含 consensus + 参考序列的 MSA 结果
        //   @param msa_cmd - 外部 MSA 工具命令模板（用于对齐插入序列）
        //                    - 示例："mafft --auto {input} > {output}"
        //                    - 支持的工具：MAFFT, Muscle, Clustal Omega
        //
        // 工作流程：
        //   1. 收集所有线程的插入序列 SAM 文件
        //   2. 合并插入序列 + 共识序列为一个 FASTA 文件
        //   3. 调用外部 MSA 工具对齐插入序列
        //   4. 解析 MSA 结果，生成 CIGAR 映射表（ref_aligned_map, insertion_aligned_map）
        //   5. 依次处理三类序列：
        //      a. 共识序列及其参考序列（来自 consensus_aligned_file）
        //      b. 插入序列（来自 aligned_insertion_fasta）
        //      c. 普通比对序列（来自各线程的 SAM 文件）
        //   6. 对每条序列应用多级 CIGAR 投影：
        //      - 第一级：query → ref（SAM CIGAR）
        //      - 第二级：ref → consensus（ref_aligned_map）
        //      - 第三级：consensus → insertion MSA（tmp_insertion_cigar）
        //   7. 移除冗余的 gap 列（可选，取决于 keep_first_length / keep_all_length）
        //   8. 长度一致性检测：确保所有序列等长
        //   9. 写入最终 MSA 文件：{work_dir}/results/final_aligned.fasta
        // ------------------------------------------------------------------
        void mergeAlignedResults(const FilePath& aligned_consensus_path, const std::string& msa_cmd);

        // ------------------------------------------------------------------
        // 静态方法：globalAlign - 全局序列比对（统一接口）
        // ------------------------------------------------------------------
        // 功能：
        // 根据两个序列的相似度自动选择合适的比对算法执行全局比对
        // 当前实现：直接使用 WFA2 算法（未来可扩展为相似度自适应策略）
        //
        // 参数：
        //   @param ref - 参考序列（字符串，A/C/G/T/N）
        //   @param query - 查询序列（字符串，A/C/G/T/N）
        //   @param similarity - 两个序列的相似度（0.0 到 1.0）
        //                       - 用于未来选择最优比对算法
        //                       - 当前版本未使用该参数，但保留接口兼容性
        //   @param ref_minimizer - 参考序列的 minimizer 索引（可选）
        //                          - 用于种子定位和锚点比对（未来扩展）
        //                          - 当前版本未使用，但保留接口用于优化
        //   @param query_minimizer - 查询序列的 minimizer 索引（可选）
        //                            - 用于种子定位和锚点比对（未来扩展）
        //                            - 当前版本未使用，但保留接口用于优化
        //
        // 返回：
        // CIGAR 操作序列（cigar::Cigar_t），描述 query 如何比对到 ref
        //
        // 算法选择策略（未来扩展）：
        // - 高相似度（similarity > 0.90）：使用 WFA2（速度优势明显）
        // - 中等相似度（0.70 < similarity <= 0.90）：使用 KSW2 延伸比对
        // - 低相似度（similarity <= 0.70）：使用 KSW2 全局比对（更稳健）
        //
        // 未来优化方向：
        // - 利用 minimizer 信息进行种子定位，减少比对范围
        // - 对于长序列，可先通过 minimizer 找到锚点，再分段比对
        // - 根据 minimizer 密度调整比对策略
        //
        // 当前实现：
        // - 直接调用 globalAlignWFA2（简化实现，保持一致性）
        // - 忽略 similarity、ref_minimizer、query_minimizer 参数（未来可用）
        //
        // 使用示例：
        // auto cigar = RefAligner::globalAlign(ref_seq, query_seq, 0.95, ref_mz, query_mz);
        //
        // 设计说明：
        // - 作为静态方法，不依赖 RefAligner 实例，可独立调用
        // - 提供统一的比对接口，方便未来扩展多算法支持
        // - 预留 minimizer 参数，避免后续修改所有调用点
        // ------------------------------------------------------------------
        static cigar::Cigar_t globalAlign(const std::string& ref,
                                          const std::string& query,
                                          double similarity,
                                          const SeedHits* ref_minimizer = nullptr,
                                          const SeedHits* query_minimizer = nullptr);

        // ------------------------------------------------------------------
        // 辅助函数：removeRefGapColumns
        // 功能：根据 ref_gap_pos 删除"参考为 gap 的那些列"（原地修改）
        //
        // 使用场景：
        // - 输入序列 seq 通常是"已经过 MSA 对齐"的序列（含 gap）
        // - ref_gap_pos 记录了参考（第一条序列）每一列是否为 gap
        // - 本函数会删除所有 ref_gap_pos[i]==true 的列，保留其余列
        //
        // 参数：
        //   - seq: 输入/输出序列（已对齐，包含 gap '-'）【原地修改】
        //   - ref_gap_pos: 参考（第一条序列）每一列是否为 gap；true 表示该列应被删除
        // ------------------------------------------------------------------
        static void removeRefGapColumns(
            std::string& seq,
            const std::vector<bool>& ref_gap_pos);


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


        // 辅助函数：写入SAM记录（选择正确的参考名称和输出文件）
        void writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                           std::string_view ref_name, seq_io::SeqWriter& out) const;

        // ------------------------------------------------------------------
        // 辅助函数：mergeConsensusAndSamToFasta
        // 功能：将共识序列和多个 SAM 文件合并为一个 FASTA 文件
        // 说明：
        // 1. 先写入共识序列（consensus_seq）到 FASTA 文件
        // 2. 然后逐个读取 SAM 文件，提取序列并追加写入
        // 3. 使用同一个 SeqWriter，避免追加模式的复杂性
        // 4. 流式处理，内存占用与文件数量和大小无关
        // 参数：
        //   - sam_paths: 输入 SAM 文件路径列表
        //   - fasta_path: 输出 FASTA 文件路径
        //   - line_width: FASTA 每行宽度（默认 80）
        // 返回：合并的总序列数（包括共识序列）
        // ------------------------------------------------------------------
        std::size_t mergeConsensusAndSamToFasta(
            const std::vector<FilePath>& sam_paths,
            const FilePath& fasta_path,
            std::size_t line_width = 80) const;

        // ------------------------------------------------------------------
        // 辅助函数：parseAlignedReferencesToCigar
        // 功能：读取 MSA 对齐后的参考序列文件，生成每个序列的 CIGAR（只包含 M 和 D）
        //
        // 重要变更（接口约定）：
        // 1) 不再通过返回值返回 map，而是通过参数输出（避免大对象返回/移动，调用端更明确）
        // 2) 新增 ref_gap_pos：标记“参考序列对齐后的每一列是否为 gap（'-'）”
        //    - 这里的“参考序列”指该对齐文件中的第一条序列（通常是 consensus 或中心序列）
        //    - ref_gap_pos[i] == true  表示第 i 列参考为 gap
        //    - ref_gap_pos[i] == false 表示第 i 列参考为碱基
        //
        // 说明：
        // - 对于每个后续序列（参考/插入序列），只根据其自身字符是否为 '-' 来编码：碱基 -> M，gap -> D
        // - 不做碱基一致性校验（例如与第一条序列的列一致性），保持当前逻辑等价
        // ------------------------------------------------------------------
        void parseAlignedReferencesToCigar(
            const FilePath& aligned_fasta_path,
            std::unordered_map<std::string, cigar::Cigar_t>& out_ref_aligned_map,
            std::vector<bool>& out_ref_gap_pos) const;


        // ==================================================================
        // 私有成员变量
        // ==================================================================

        // ------------------------------------------------------------------
        // 工作目录与数据路径
        // ------------------------------------------------------------------
        FilePath work_dir;  // 工作目录（存放中间文件和最终结果）

        // ------------------------------------------------------------------
        // 参考序列与索引
        // ------------------------------------------------------------------
        seq_io::SeqRecords ref_sequences;   // 参考序列集合（从 ref_fasta_path 读取）
        mash::Sketches ref_sketch;          // 每个参考序列的 MinHash sketch（用于快速相似度估计）
        std::vector<SeedHits> ref_minimizers;  // 每个参考序列的 Minimizer 索引（用于快速定位同源区域）

        // ------------------------------------------------------------------
        // 共识序列
        // ------------------------------------------------------------------
        // 说明：
        // - 共识序列是所有参考序列的"代表序列"
        // - 如果 keep_first_length=false，通过 MSA 工具生成（如 MAFFT）
        // - 如果 keep_first_length=true，直接使用 ref_sequences[0]
        // - 在 mergeAlignedResults 中作为坐标系的锚点
        // ------------------------------------------------------------------
        seq_io::SeqRecord consensus_seq;  // 共识序列（作为成员变量存储）
        mash::Sketch consensus_sketch;    // 共识序列的 MinHash sketch（用于二次比对时的相似度计算）
                                          // 说明：在构造函数中初始化一次，避免每次 alignOneQueryToRef 都重复计算
        SeedHits consensus_minimizer;     // 共识序列的 Minimizer 索引（用于二次比对时的种子定位）
                                          // 说明：在构造函数中初始化一次，避免每次 alignOneQueryToRef 都重复计算

        // ------------------------------------------------------------------
        // MinHash 和 Minimizer 参数
        // ------------------------------------------------------------------
        int kmer_size = 21;      // k-mer 大小（MinHash 和 Minimizer 共用）
        int window_size = 10;    // Minimizer 窗口大小
        int sketch_size = 2000;  // MinHash sketch 大小（哈希数量）
        int random_seed = 42;    // 随机种子（用于 MinHash 哈希函数）

        // ------------------------------------------------------------------
        // 线程与外部工具配置
        // ------------------------------------------------------------------
        int threads = 1;            // 线程数（用于并行比对）
        std::string msa_cmd;        // MSA 工具命令模板（用于共识序列生成和插入序列比对）

        // ------------------------------------------------------------------
        // MSA 输出选项
        // ------------------------------------------------------------------
        // 说明：
        // - keep_first_length：是否保持第一条序列（共识序列）的原始长度
        //   * true：移除共识序列为 gap 的所有列
        //   * false：保留所有列（包括共识序列的 gap）
        // - keep_all_length：是否保持所有序列的原始长度（优先级低于 keep_first_length）
        //   * true：移除插入 MSA 中共识序列为 gap 的列
        //   * false：保留所有列
        // ------------------------------------------------------------------
        bool keep_first_length = false;
        bool keep_all_length = false;

        // ------------------------------------------------------------------
        // MinHash 计算选项
        // ------------------------------------------------------------------
        // 说明：
        // - noncanonical：是否使用非规范化的 k-mer（不区分正负链）
        //   * true：考虑 k-mer 的反向互补（用于 DNA 序列）
        //   * false：只考虑原始 k-mer（用于蛋白质或单链 RNA）
        // ------------------------------------------------------------------
        bool noncanonical = true;

        // ------------------------------------------------------------------
        // 输出文件路径
        // ------------------------------------------------------------------
        // 说明：
        // - outs_path：每个线程的普通 SAM 输出文件路径
        //   * 长度 = threads
        //   * outs_path[i] = {work_dir}/data/aligned_thread_{i}.sam
        // - outs_with_insertion_path：每个线程的插入序列 SAM 输出文件路径
        //   * 长度 = threads
        //   * outs_with_insertion_path[i] = {work_dir}/data/aligned_insertion_thread_{i}.sam
        // ------------------------------------------------------------------
        std::vector<FilePath> outs_path;
        std::vector<FilePath> outs_with_insertion_path;
    };

} // namespace align

#endif //HALIGN4_ALIGN_H

#include "align.h"

#include <stdexcept>
#include <cctype>
#include <limits>
#include <string>

// ==================================================================
// cigar.cpp - CIGAR 编码/解码/转换与序列按 CIGAR 插空（投影）
// ==================================================================
// 本文件实现 namespace cigar 下的核心工具函数：
// 1) CIGAR 压缩编码（char+len <-> uint32_t）
// 2) SAM CIGAR 字符串 <-> 压缩 Cigar_t 的互转
// 3) 检测插入（hasInsertion）
// 4) padQueryToRefByCigar：按 CIGAR 将 query 投影到 ref 坐标系（插入 gap '-')
//
// 重要约定（正确性）：
// - 本项目的 CIGAR 压缩形式与 BAM/SAM 一致：unit=(len<<4)|op。
// - op 编码采用标准 BAM 编码：
//     0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
// - padQueryToRefByCigar 的语义：
//   "根据 CIGAR 在 query 中插入缺失列（D/N）对应的 gap '-'，从而让 query 与 ref 对齐"。
//   注意：M/I/S/=/X 都会消耗 query 的字符；D/N 不消耗 query，只在输出中放 '-'
//
// 性能设计点：
// - stringToCigar：单次线性扫描解析，预估 reserve 避免 vector 扩容
// - cigarToString：预估字符串容量，减少 realloc
// - padQueryToRefByCigar：采用"从后往前"填充，避免 std::string insert 的 O(N^2)
//
// 边界条件：
// - CIGAR 为 "*"（SAM 未知）-> 返回空 Cigar_t
// - 空字符串 -> 返回空 Cigar_t
// - 长度溢出（len >= 2^28）-> 抛异常
// - 非法操作符 -> 抛异常
// ==================================================================

namespace cigar
{
    // ------------------------------------------------------------------
    // 内部：操作符 <-> op 编码 映射
    // ------------------------------------------------------------------
    // 说明：
    // - BAM/SAM 的 op 编码是固定的（见上方约定）
    // - 这里用 switch 显式映射，避免查表引入额外静态数据
    // - 由于操作符集合很小，switch 在 O2 下通常会被编译器优化为跳转表
    // ------------------------------------------------------------------

    static inline uint32_t opCharToCode(const char op)
    {
        switch (op) {
        case 'M': return 0;
        case 'I': return 1;
        case 'D': return 2;
        case 'N': return 3;
        case 'S': return 4;
        case 'H': return 5;
        case 'P': return 6;
        case '=': return 7;
        case 'X': return 8;
        default:
            throw std::runtime_error(std::string("Unknown CIGAR op char: ") + op);
        }
    }

    static inline char opCodeToChar(const uint32_t code)
    {
        switch (code) {
        case 0: return 'M';
        case 1: return 'I';
        case 2: return 'D';
        case 3: return 'N';
        case 4: return 'S';
        case 5: return 'H';
        case 6: return 'P';
        case 7: return '=';
        case 8: return 'X';
        default:
            throw std::runtime_error("Unknown CIGAR op code: " + std::to_string(code));
        }
    }

    // ------------------------------------------------------------------
    // cigarToInt：编码 (operation,len) -> CigarUnit
    // ------------------------------------------------------------------
    // 编码布局：
    //   [len:28bits][op:4bits]
    // 其中 len 必须 < 2^28。
    // ------------------------------------------------------------------
    CigarUnit cigarToInt(char operation, uint32_t len)
    {
        // 防御式检查：len 超过 28bit 会覆盖 op 位并导致不可逆
        constexpr uint32_t kMaxLen = (1u << 28) - 1u;
        if (len == 0 || len > kMaxLen) {
            throw std::runtime_error("cigarToInt: invalid length=" + std::to_string(len));
        }
        const uint32_t op = opCharToCode(operation);
        return (len << 4) | (op & 0x0Fu);
    }

    // ------------------------------------------------------------------
    // intToCigar：解码 CigarUnit -> (operation,len)
    // ------------------------------------------------------------------
    // 说明：
    // - 低 4 位是 op，高位右移 4 得到 len
    // - 此函数不检查 len==0（因为理论上不应产生 0），保持与历史行为一致
    // ------------------------------------------------------------------
    void intToCigar(CigarUnit cigar, char& operation, uint32_t& len)
    {
        const uint32_t op = (cigar & 0x0Fu);
        len = (cigar >> 4);
        operation = opCodeToChar(op);
    }

    // ------------------------------------------------------------------
    // hasInsertion：检测是否存在 'I'
    // ------------------------------------------------------------------
    // 性能：O(#ops)，遇到第一个 I 立即返回（短路）
    // ------------------------------------------------------------------
    bool hasInsertion(const Cigar_t& cigar)
    {
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len;
            // 统一使用 intToCigar 解码以进行操作符检查
            intToCigar(cu, op_char, len);
            if (op_char == 'I') return true;
        }
        return false;
    }

    // ------------------------------------------------------------------
    // cigarToString：压缩 Cigar_t -> SAM CIGAR string
    // ------------------------------------------------------------------
    // 说明：
    // - 对每个 CigarUnit，解码为 len+op，再拼接到字符串
    // - 预分配字符串容量：粗略估计每个 op 平均需要 1-10 个数字 + 1 个 op 字符
    // ------------------------------------------------------------------
    std::string cigarToString(const Cigar_t& cigar)
    {
        if (cigar.empty()) return "";

        std::string out;
        // 粗略预估：每个操作平均 5 字符（如 10M/100M），避免频繁扩容
        out.reserve(cigar.size() * 5);

        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            // 使用已有函数解码，确保编码逻辑的一致性
            intToCigar(cu, op_char, len);
            out.append(std::to_string(len));
            out.push_back(op_char);
        }
        return out;
    }

    // ------------------------------------------------------------------
    // stringToCigar：SAM CIGAR string -> 压缩 Cigar_t
    // ------------------------------------------------------------------
    // 解析规则：
    // - 若输入为 "*"：表示未知，返回空
    // - 跳过空白字符（空格/\t/\n），增强容错性
    // - 每个操作必须是：<正整数><op_char>
    // - len==0、缺少数字、未知 op 都视为错误
    // ------------------------------------------------------------------
    Cigar_t stringToCigar(const std::string& cigar_str)
    {
        Cigar_t result;
        if (cigar_str.empty() || cigar_str == "*") return result;

        // 预估 reserve：操作数通常远小于字符串长度，取一个保守估计
        result.reserve(cigar_str.size() / 2 + 1);

        uint64_t len_acc = 0;
        bool has_number = false;

        for (std::size_t i = 0; i < cigar_str.size(); ++i) {
            const unsigned char c = static_cast<unsigned char>(cigar_str[i]);

            if (std::isspace(c)) {
                continue; // 容错：跳过空白
            }

            if (std::isdigit(c)) {
                has_number = true;
                len_acc = len_acc * 10 + (c - '0');

                // 溢出保护：len 不能超过 2^28-1
                if (len_acc > ((1ull << 28) - 1ull)) {
                    throw std::runtime_error("stringToCigar: op length overflow in '" + cigar_str + "'");
                }
                continue;
            }

            // 走到这里：c 不是数字也不是空白，认为是 op 字符
            if (!has_number || len_acc == 0) {
                throw std::runtime_error("stringToCigar: missing/invalid length before op in '" + cigar_str + "'");
            }

            const uint32_t len = static_cast<uint32_t>(len_acc);
            result.push_back(cigarToInt(static_cast<char>(c), len));

            // 重置解析状态，准备下一个 op
            len_acc = 0;
            has_number = false;
        }

        // 如果字符串以数字结尾（缺少 op），属于格式错误
        if (has_number) {
            throw std::runtime_error("stringToCigar: trailing number without op in '" + cigar_str + "'");
        }

        return result;
    }

    // ------------------------------------------------------------------
    // alignQueryToRef：按 CIGAR 将 query 投影到 ref 坐标系
    // ------------------------------------------------------------------
    // 目标：返回时 query 的长度等于 "对齐列数"。
    //
    // 关键策略：从后往前填充（reverse fill）
    // - 直接在 std::string 中间 insert 会导致 O(N^2)
    // - 我们先移动出原始 query，再根据 CIGAR 计算目标长度，分配新字符串
    // - 从末尾向前写入：
    //   * 消耗 query 的操作（M/I/S/=/X）从 old_query 取字符
    //   * 不消耗 query 的操作（D/N）写入 '-'
    //
    // 正确性不变量：
    // - 所有消耗 query 的操作总长度必须等于 old_query.size()
    // - 若不等，说明 CIGAR 与 query 长度不匹配（Debug 断言/运行时错误）
    // ------------------------------------------------------------------
    void padQueryToRefByCigar(std::string& query, const Cigar_t& cigar)
    {
        if (cigar.empty()) {
            return; // 空 CIGAR：不做任何调整（保持历史行为）
        }

        // 1) 统计结果长度：
        // - D/N 会在输出中产生 '-'（增加列数）
        // - 其它操作都会输出一个字符（来自 query）
        //   注意：这里"输出字符"的数量 = 对齐列数
        std::size_t out_len = 0;
        std::size_t consume_query = 0;

        for (const CigarUnit cu : cigar) {
            // 使用 intToCigar 统一解码 CigarUnit，增强可读性和一致性
            // 说明：intToCigar 将压缩整数还原为 (op_char, len)
            // 性能说明：intToCigar 内部本质上也是位运算和查表，开销极小；
            // 本次替换主要为了代码可维护性与避免重复低级位操作。
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            // D/N：ref 消耗，query 不消耗，但输出 '-' len 次
            if (op_char == 'D' || op_char == 'N') {
                out_len += len;
            } else if (op_char == 'H' || op_char == 'P') {
                // H/P：不消耗 ref，也不消耗 query，也不产生输出列（在 MSA 场景中可能出现）
                // 保持与历史行为：忽略
                continue;
            } else {
                // M/I/S/=/X：输出 len 个字符，且消耗 query len 个字符
                out_len += len;
                consume_query += len;
            }
        }

        // 2) 防御式检查：CIGAR 消耗的 query 长度必须与实际一致
        // 说明：release 构建中这里不强制 throw，以保持与历史行为一致；
        // 但在 Debug 下 assert 能更快发现上游 bug。
        // （如果你希望 release 也 fail-fast，可以改成 throw；本次按“逻辑不变”不改。）
        assert(consume_query == query.size());

        // 3) 备份原 query，并分配输出空间
        // 性能：move 避免拷贝原字符串；resize 一次性分配
        std::string old = std::move(query);
        query.assign(out_len, '-');

        // 4) 从后往前填充
        // w：写指针（指向 query 的当前位置）
        // r：读指针（指向 old 的当前位置）
        std::size_t w = out_len;
        std::size_t r = old.size();

        for (auto it = cigar.rbegin(); it != cigar.rend(); ++it) {
            // 统一通过 intToCigar 解码当前操作及长度
            char op_char;
            uint32_t len = 0;
            intToCigar(*it, op_char, len);

            if (op_char == 'D' || op_char == 'N') {
                // D/N：输出 '-' len 次，不消耗 old
                for (uint32_t i = 0; i < len; ++i) {
                    query[--w] = '-';
                }
                continue;
            }

            if (op_char == 'H' || op_char == 'P') {
                // H/P：忽略
                continue;
            }

            // M/I/S/=/X：从 old 拷贝 len 个字符
            for (uint32_t i = 0; i < len; ++i) {
                // Debug 保护：防止下溢（说明 CIGAR 与 query 不匹配）
                assert(r > 0);
                query[--w] = old[--r];
            }
        }

        // 5) Debug 保护：确保完全写满且完全消费 old
        assert(w == 0);
        assert(r == 0);
    }

    // ------------------------------------------------------------------
    // appendCigar：智能追加 CIGAR 并合并相邻同类型操作
    // ------------------------------------------------------------------
    // 设计目标：
    // - 合并相邻同类型操作（如 5M + 3M -> 8M），减少 CIGAR 长度
    // - 避免产生冗余操作（如 0M）
    //
    // 性能：O(cigar_to_add.size())，每个操作最多检查一次
    // ------------------------------------------------------------------
    void appendCigar(Cigar_t& result, const Cigar_t& cigar_to_add)
    {
        for (const CigarUnit cu : cigar_to_add) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            // 跳过长度为 0 的操作（防御性编程）
            if (len == 0) continue;

            if (result.empty()) {
                // 结果为空，直接追加
                result.push_back(cu);
            } else {
                // 检查是否可以与最后一个操作合并
                char last_op_char;
                uint32_t last_len = 0;
                intToCigar(result.back(), last_op_char, last_len);

                if (last_op_char == op_char) {
                    // 相同操作类型，合并长度
                    // 防溢出检查：确保合并后不超过 28 位
                    constexpr uint32_t kMaxLen = (1u << 28) - 1u;
                    if (static_cast<uint64_t>(last_len) + len > kMaxLen) {
                        throw std::runtime_error("appendCigar: merged length overflow");
                    }
                    result.back() = cigarToInt(op_char, last_len + len);
                } else {
                    // 不同操作类型，直接追加
                    result.push_back(cu);
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // getRefLength：计算 CIGAR 消耗的参考序列长度
    // ------------------------------------------------------------------
    // 说明：
    // - 消耗 ref 的操作：M, D, N, =, X
    // - 不消耗 ref 的操作：I, S, H, P
    //
    // 性能：O(cigar.size())，单次遍历
    // ------------------------------------------------------------------
    std::size_t getRefLength(const Cigar_t& cigar)
    {
        std::size_t total = 0;
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            // 消耗 ref 的操作
            if (op_char == 'M' || op_char == 'D' || op_char == 'N' ||
                op_char == '=' || op_char == 'X') {
                total += len;
            }
        }
        return total;
    }

    // ------------------------------------------------------------------
    // getQueryLength：计算 CIGAR 消耗的查询序列长度
    // ------------------------------------------------------------------
    // 说明：
    // - 消耗 query 的操作：M, I, S, =, X
    // - 不消耗 query 的操作：D, N, H, P
    //
    // 性能：O(cigar.size())，单次遍历
    // ------------------------------------------------------------------
    std::size_t getQueryLength(const Cigar_t& cigar)
    {
        std::size_t total = 0;
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            // 消耗 query 的操作
            if (op_char == 'M' || op_char == 'I' || op_char == 'S' ||
                op_char == '=' || op_char == 'X') {
                total += len;
            }
        }
        return total;
    }

} // namespace cigar


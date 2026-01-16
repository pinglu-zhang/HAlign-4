#include "align.h"

#include <cassert>
#include <cstdint>
#include <stdexcept>

namespace cigar {

// ------------------------------------------------------------------
// 约定：SAM CIGAR 在 BAM 中的压缩编码方式为：
//   cigar_unit = (len << 4) | op
// 其中 op 为 4-bit 操作码（0..15），len 为 28-bit（0..2^28-1）。
//
// 本文件提供：
// - cigarToInt : (char op, len) -> CigarUnit
// - intToCigar : CigarUnit -> (char op, len)
//
// 注意：
// - 这里的编码/解码只做“表示层”的转换，不对 len 做业务意义校验（例如 len==0 是否允许）。
// - 对于未知操作符，我们选择抛出异常（比静默映射成 '?' 更安全）。
// ------------------------------------------------------------------

namespace {

// SAM/BAM 的 CIGAR 操作码（4-bit）。
// 这套编号是行业约定（htslib/BAM spec），必须固定，否则与外部工具不兼容。
constexpr std::uint32_t kOpM  = 0;  // M: alignment match (match or mismatch)
constexpr std::uint32_t kOpI  = 1;  // I: insertion to the reference
constexpr std::uint32_t kOpD  = 2;  // D: deletion from the reference
constexpr std::uint32_t kOpN  = 3;  // N: skipped region from the reference
constexpr std::uint32_t kOpS  = 4;  // S: soft clipping (clipped sequences present in SEQ)
constexpr std::uint32_t kOpH  = 5;  // H: hard clipping (clipped sequences NOT present in SEQ)
constexpr std::uint32_t kOpP  = 6;  // P: padding (silent deletion from padded reference)
constexpr std::uint32_t kOpEq = 7;  // =: sequence match
constexpr std::uint32_t kOpX  = 8;  // X: sequence mismatch

constexpr std::uint32_t kOpMask = 0xFu;
constexpr std::uint32_t kLenBits = 28u;
constexpr std::uint32_t kMaxLen = (1u << kLenBits) - 1u;

// 把字符 op 映射为 4-bit opCode。
// 说明：用 switch 是最直接且零开销（编译器可生成跳表/比较链）。
constexpr std::uint32_t opCharToCode(const char operation)
{
    switch (operation) {
    case 'M': return kOpM;
    case 'I': return kOpI;
    case 'D': return kOpD;
    case 'N': return kOpN;
    case 'S': return kOpS;
    case 'H': return kOpH;
    case 'P': return kOpP;
    case '=': return kOpEq;
    case 'X': return kOpX;
    default:  return 0xFFFFFFFFu;
    }
}

// 把 4-bit opCode 映射回字符 op。
constexpr char opCodeToChar(const std::uint32_t op_code)
{
    switch (op_code) {
    case kOpM:  return 'M';
    case kOpI:  return 'I';
    case kOpD:  return 'D';
    case kOpN:  return 'N';
    case kOpS:  return 'S';
    case kOpH:  return 'H';
    case kOpP:  return 'P';
    case kOpEq: return '=';
    case kOpX:  return 'X';
    default:    return '?';
    }
}

} // namespace

CigarUnit cigarToInt(const char operation, const std::uint32_t len)
{
    // 关键正确性：len 只能占用 28bit。超出会破坏 opCode 位域。
    // 这里选择抛异常，避免 silent overflow。
    if (len > kMaxLen) {
        throw std::invalid_argument("cigarToInt: len exceeds 28-bit limit");
    }

    const std::uint32_t op_code = opCharToCode(operation);
    if (op_code == 0xFFFFFFFFu) {
        throw std::invalid_argument("cigarToInt: unknown CIGAR operation");
    }

    // (len << 4) | op
    return static_cast<CigarUnit>((len << 4) | (op_code & kOpMask));
}

void intToCigar(const CigarUnit cigar_unit, char& operation, std::uint32_t& len)
{
    // 低 4 位是 opCode，高 28 位是长度。
    const std::uint32_t op_code = static_cast<std::uint32_t>(cigar_unit) & kOpMask;
    len = static_cast<std::uint32_t>(cigar_unit) >> 4;

    operation = opCodeToChar(op_code);

    // 对于未知 opcode：这里用断言 + 保守返回 '?'。
    // - release 模式下不中断（避免解析外部脏数据时直接崩溃）
    // - debug 模式下尽早暴露问题
    assert(operation != '?');
}

// ------------------------------------------------------------------
// 函数：hasInsertion
// 功能：检测 CIGAR 序列中是否存在插入操作（'I'）
//
// 实现说明：
// 1. 性能优化：直接对 CigarUnit 的低 4 位进行位运算检查，避免完整解码
// 2. 短路优化：找到第一个 'I' 操作即返回 true，无需遍历整个序列
// 3. 复杂度：O(N)，最坏情况遍历所有操作；最好情况 O(1)（第一个就是 'I'）
// 4. 正确性：依赖 kOpI = 1 的编码约定（SAM/BAM 规范固定值）
//
// 参数：cigar - CIGAR 操作序列（压缩形式，每个元素为 CigarUnit）
// 返回：true 表示存在至少一个插入操作，false 表示不存在
// ------------------------------------------------------------------
bool hasInsertion(const Cigar_t& cigar)
{
    // 遍历所有 CIGAR 操作，检查低 4 位是否为 kOpI（插入操作码 = 1）
    for (const auto cigar_unit : cigar) {
        // 提取低 4 位的操作码
        const std::uint32_t op_code = static_cast<std::uint32_t>(cigar_unit) & kOpMask;

        // 如果是插入操作，立即返回 true（短路优化）
        if (op_code == kOpI) {
            return true;
        }
    }

    // 遍历完所有操作都没有找到插入，返回 false
    return false;
}

// ------------------------------------------------------------------
// 函数：cigarToString
// 功能：将 CIGAR 从压缩格式转换为 SAM 标准字符串格式
//
// 实现说明：
// 1. 性能优化：预分配字符串空间（cigar.size() * 5），避免多次内存重新分配
//    - 估算：平均每个 CIGAR 操作占用约 3-5 字符（如 "100M" = 4 字符）
//    - 对于大多数情况，reserve(cigar.size() * 5) 足够，避免了 string 的自动扩容
// 2. 使用 std::to_string 转换整数，编译器会优化为高效实现
// 3. 复杂度：O(N)，N 为 CIGAR 操作数量
// 4. 字符串拼接：使用 += 操作符，在预分配空间足够时为 O(1) 追加
//
// 参数：cigar - CIGAR 操作序列（压缩形式）
// 返回：SAM 格式的 CIGAR 字符串，例如 "100M5I95M"
//
// 示例：
//   输入：[cigarToInt('M', 100), cigarToInt('I', 5), cigarToInt('M', 95)]
//   输出："100M5I95M"
// ------------------------------------------------------------------
std::string cigarToString(const Cigar_t& cigar)
{
    // 预分配字符串空间，减少内存重新分配次数
    // 估算：每个 CIGAR 操作平均占用约 5 个字符（例如 "1000M" = 5 字符）
    std::string result;
    result.reserve(cigar.size() * 5);

    // 遍历所有 CIGAR 操作，逐个解码并拼接到结果字符串
    for (const auto cigar_unit : cigar) {
        char op;
        std::uint32_t len;
        intToCigar(cigar_unit, op, len);

        // 格式：<长度><操作符>，例如 "100M"
        result += std::to_string(len);
        result += op;
    }

    return result;
}

} // namespace cigar


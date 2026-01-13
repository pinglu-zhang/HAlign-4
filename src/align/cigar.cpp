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

} // namespace cigar


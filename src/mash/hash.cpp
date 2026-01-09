#include "hash.h"


// 将 code2bit 固化成“大端字节序”输入（跨平台一致）
// 在 little-endian 机器上做一次 bswap，使得内存字节序等价于原值的 big-endian bytes。
static inline std::uint64_t to_be64(std::uint64_t x) noexcept {
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
    return __builtin_bswap64(x);
#else
    return x;
#endif
}

hash_t getHash(const char* seq, int length, std::uint32_t seed)
{
    // XXH3_64bits_withSeed 的 seed 类型是 XXH64_hash_t（64-bit），直接扩展即可
    return static_cast<hash_t>(
        XXH3_64bits_withSeed(seq, static_cast<size_t>(length), static_cast<XXH64_hash_t>(seed))
    );
}

hash_t getHash2bit(std::uint64_t code2bit, std::uint32_t seed)
{
    const std::uint64_t be = to_be64(code2bit);
    return static_cast<hash_t>(
        XXH3_64bits_withSeed(&be, sizeof(be), static_cast<XXH64_hash_t>(seed))
    );
}

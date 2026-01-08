#include "hash.h"
#include "MurmurHash3.h"
#include <cstring>
#include <cstdint>

hash_t getHash(const char * seq, int length, std::uint32_t seed)
{
#ifdef ARCH_32
    // 32-bit: mimic Mash approach by concatenating two 32-bit hashes
    std::uint8_t out8[8];
    MurmurHash3_x86_32(seq, length > 16 ? 16 : length, seed, out8);
    MurmurHash3_x86_32(seq + (length > 16 ? 16 : length), length > 16 ? (length - 16) : 0, seed, out8 + 4);
    hash_t h = 0;
    std::memcpy(&h, out8, 8);
    return h;
#else
    std::uint8_t out16[16];
    MurmurHash3_x64_128(seq, length, seed, out16);
    hash_t h = 0;
    std::memcpy(&h, out16, 8);
    return h;
#endif
}

hash_t getHash2bit(std::uint64_t code2bit, std::uint32_t seed)
{
    // 固化输入字节序（big-endian），保证跨平台一致
    std::uint8_t key8[8];
    for (int i = 0; i < 8; ++i)
        key8[i] = static_cast<std::uint8_t>((code2bit >> (56 - 8 * i)) & 0xFF);

#ifdef ARCH_32
    std::uint8_t out8[8];
    // 低 32 + 高 32 拼接
    MurmurHash3_x86_32(key8, 8, seed, out8);
    MurmurHash3_x86_32(key8, 8, seed, out8 + 4);
    hash_t h = 0;
    std::memcpy(&h, out8, 8);
    return h;
#else
    std::uint8_t out16[16];
    MurmurHash3_x64_128(key8, 8, seed, out16);
    hash_t h = 0;
    std::memcpy(&h, out16, 8);
    return h;
#endif
}

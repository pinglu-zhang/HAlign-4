#ifndef HALIGN4_HASH_H
#define HALIGN4_HASH_H

#include <cstddef>
#include <cstdint>


using hash_t = std::uint64_t;

// Compute hash for an arbitrary byte sequence.
hash_t getHash(const char * seq, int length, std::uint32_t seed = 0);

// Compute hash for a 2-bit encoded k-mer code (canonicalized by caller if needed).
hash_t getHash2bit(std::uint64_t code2bit, std::uint32_t seed = 0);


#endif //HALIGN4_HASH_H
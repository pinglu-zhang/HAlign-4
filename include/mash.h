#ifndef HALIGN4_MASH_H
#define HALIGN4_MASH_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "hash.h"
#include "bloom_filter.hpp"

namespace mash
{
    inline constexpr std::uint8_t nt4_table[256] = {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,
        4,4,4,4,4,4,4,4,4,4,0,4,1,4,4,4,
        2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,
        3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
    };

    // ------------------------------------------------------------
    // Sketch (sorted unique bottom-k hashes)
    // ------------------------------------------------------------
    struct Sketch
    {
        std::size_t k = 0;
        bool noncanonical = true;
        std::vector<hash_t> hashes;

        std::size_t size() const noexcept { return hashes.size(); }
        bool empty() const noexcept { return hashes.empty(); }
    };

    // ------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------
    // 生成 MinHash sketch：对序列做 k-mer hash，然后取 bottom-k（sketch_size）并 sort+unique。
    // 注意：mash 模块与 seed/minimizer 无关；这里的 w 仅为兼容旧接口保留，当前实现不使用。
    Sketch sketchFromSequence(const std::string& seq,
                             std::size_t k,
                             std::size_t sketch_size,
                             bool noncanonical = true,
                             int seed = 0);

    bloom_filter filterFromSketch(const Sketch& sk, double false_positive_rate = 0.0001, int seed = 0xA5A5A5A5);


    // ------------------------------------------------------------
    // Similarity & distance
    // ------------------------------------------------------------

    // - both empty => 1.0
    // - one empty  => 0.0
    double jaccard(const Sketch& a, const Sketch& b);

    double jaccard(const bloom_filter& a, const Sketch& b);

    // d = - 1/k * ln( 2j/(1+j) )
    // - if j<=0 => +inf
    // - if j>=1 => 0
    double mashDistanceFromJaccard(double j, std::size_t k);

    // ANI ~= (2j/(1+j))^(1/k)
    // - clamp to [0,1]
    double aniFromJaccard(double j, std::size_t k);

    // ANI ~= exp(-d)
    // - clamp to [0,1]
    double aniFromMashDistance(double d);


    // ------------------------------------------------------------
    // Low-level helpers (sorted unique)
    // ------------------------------------------------------------

    std::size_t intersectionSizeSortedUnique(const std::vector<hash_t>& a,
                                             const std::vector<hash_t>& b) noexcept;

} // namespace mash

#endif //HALIGN4_MASH_H


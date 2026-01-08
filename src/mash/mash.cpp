// mash.cpp
// -----------------
// 该文件实现了简化版的 Mash 风格 sketch 与相关的度量函数（Jaccard、Mash 距离、ANI 推断等）。
// 我们在注释中说明设计考虑、参数意义、边界情况与实现细节，便于性能测试与维护。

#include "mash.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <queue>
#include <cstring> // for memcmp
#include <unordered_set>
#include "hash.h"
#include "robin_hood.h"

namespace mash
{
    // clamp01: 将 double 值限制在 [0,1] 范围内的辅助函数
    // 用于确保返回概率/比例型值的数值稳定性
    static inline double clamp01(double x) noexcept
    {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    }

    Sketch sketchFromSequence(const std::string& seq,
                              std::size_t k,
                              std::size_t sketch_size,
                              bool noncanonical,
                              int seed)
    {
        Sketch sk;
        sk.k = k;
        sk.noncanonical = noncanonical;

        if (k == 0 || sketch_size == 0 || seq.size() < k) return sk;
        if (k > 32) return sk;

        const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
        const std::uint64_t shift = 2ULL * (k - 1);

        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        std::size_t valid = 0;

        // max-heap (top is largest), keep bottom-k
        std::priority_queue<hash_t> maxHeap;
        robin_hood::unordered_set<hash_t> seen;
        seen.reserve(sketch_size * 2 + 1);

        for (std::size_t i = 0; i < seq.size(); ++i)
        {
            const uint8_t c = nt4_table[static_cast<unsigned char>(seq[i])];
            if (c >= 4)
            {
                fwd = rev = 0;
                valid = 0;
                continue;
            }

            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | (std::uint64_t(3U ^ c) << shift);

            if (valid < k) ++valid;
            if (valid < k) continue;

            const std::uint64_t code = noncanonical ? fwd : std::min(fwd, rev);
            const hash_t h = getHash2bit(code, static_cast<std::uint32_t>(seed));

            if (!seen.insert(h).second) continue;

            if (maxHeap.size() < sketch_size)
            {
                maxHeap.push(h);
            }
            else if (h < maxHeap.top())
            {
                // evict current largest
                seen.erase(maxHeap.top());
                maxHeap.pop();
                maxHeap.push(h);
            }
        }

        sk.hashes.reserve(maxHeap.size());
        while (!maxHeap.empty())
        {
            sk.hashes.push_back(maxHeap.top());
            maxHeap.pop();
        }

        std::sort(sk.hashes.begin(), sk.hashes.end());

        return sk;
    }

    // jaccard: 基于两个已排序且唯一化的 sketch 向量计算 Jaccard 相似度
    // 返回范围在 [0,1]，特殊 case：两个空集合返回 1.0（定义上的选择，表示完全相同的空集）
    double jaccard(const Sketch& a, const Sketch& b)
    {
        if (a.k != b.k) throw std::invalid_argument("mash::jaccard: mismatched k");

        if (a.hashes.empty() && b.hashes.empty()) return 1.0;
        if (a.hashes.empty() || b.hashes.empty()) return 0.0;

        const std::size_t inter = intersectionSizeSortedUnique(a.hashes, b.hashes);
        const std::size_t uni = a.hashes.size() + b.hashes.size() - inter;
        if (uni == 0) return 1.0;
        return static_cast<double>(inter) / static_cast<double>(uni);
    }

    // mashDistanceFromJaccard: 将 Jaccard 相似性转换为 Mash 距离（基于 Mash 的数学推导）
    // 公式：x = (2*j)/(1+j) ; distance = -ln(x) / k
    // 注意边界情况的处理：当 j<=0 或 x<=0 时返回 +inf；当 j>=1 时返回 0
    double mashDistanceFromJaccard(double j, std::size_t k)
    {
        if (k == 0) throw std::invalid_argument("mash::mashDistanceFromJaccard: k must be > 0");
        if (!(j > 0.0)) return std::numeric_limits<double>::infinity();
        if (j >= 1.0) return 0.0;

        const double x = (2.0 * j) / (1.0 + j);
        if (!(x > 0.0)) return std::numeric_limits<double>::infinity();
        return -std::log(x) / static_cast<double>(k);
    }

    // aniFromJaccard: 基于 Jaccard 估计平均核苷酸相似度（ANI），使用 Mash 的近似关系
    // 公式：x = (2*j)/(1+j) ; ANI ~ x^(1/k)
    // 返回值被 clamp 到 [0,1]
    double aniFromJaccard(double j, std::size_t k)
    {
        if (k == 0) throw std::invalid_argument("mash::aniFromJaccard: k must be > 0");
        if (!(j > 0.0)) return 0.0;
        if (j >= 1.0) return 1.0;

        const double x = (2.0 * j) / (1.0 + j);
        if (!(x > 0.0)) return 0.0;

        return clamp01(std::pow(x, 1.0 / static_cast<double>(k)));
    }

    // aniFromMashDistance: 从 Mash 距离反推 ANI：ANI ~ exp(-d)
    // 对无穷大或非有限输入做防护处理，输出限定在 [0,1]
    double aniFromMashDistance(double d)
    {
        if (!std::isfinite(d)) return 0.0;
        if (d <= 0.0) return 1.0;
        return clamp01(std::exp(-d));
    }

    // intersectionSizeSortedUnique
    // 计算两个已排序且唯一化（unique）的哈希向量的交集大小。算法为经典的双指针线性扫描。
    // 要求输入 a,b 已经满足升序且无重复；该函数不会修改输入，只返回交集元素数量（不返回集合本身）。
    // 时间复杂度：O(|a| + |b|)。适用于 sketch 的哈希集合比较。
    // 注意：由于 hash_u 是联合体，调用者需要确保两个向量使用相同的 use64 设置。
    // 这个函数主要用于向后兼容，建议直接使用 jaccard 函数。
    std::size_t intersectionSizeSortedUnique(const std::vector<hash_t>& a,
                                             const std::vector<hash_t>& b) noexcept
    {
        std::size_t i = 0, j = 0, inter = 0;
        while (i < a.size() && j < b.size()) {
            const auto av = a[i];
            const auto bv = b[j];
            if (av == bv) {
                ++inter;
                ++i;
                ++j;
            } else if (av < bv) {
                ++i;
            } else {
                ++j;
            }
        }
        return inter;
    }

} // namespace mash

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
#include "hash.h"
#include "robin_hood.h"
#include <unordered_map>
#include <unordered_set>
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

        // --- 参数/边界检查 ---
        // k=0 或 sketch_size=0：没有意义，直接返回空 sketch。
        // seq.size() < k：不足以形成任何 k-mer。
        if (k == 0 || sketch_size == 0 || seq.size() < k) return sk;
        // 这里实现使用 2-bit 编码将 k-mer 压缩到 64bit 中，因此限制 k<=32。
        if (k > 32) return sk;

        // mask：用于保留 rolling 编码的低 2*k 位。
        const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
        // shift：用于维护反向互补 rolling 编码时，把新碱基放到最高两位的位置。
        const std::uint64_t shift = 2ULL * (k - 1);

        // fwd：正向 rolling 2-bit 编码
        // rev：反向互补 rolling 2-bit 编码
        // valid：当前窗口累计的有效碱基数（遇到 N/非法字符会清零）
        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        std::size_t valid = 0;

        // max-heap（堆顶是当前集合里“最大”的 hash），用于维护 bottom-k（最小的 k 个 hash）
        // 复杂度：每个 k-mer 仅 O(log(sketch_size))。
        std::priority_queue<hash_t> maxHeap;

        // seen：用于去重，避免同一个 hash 重复进入 heap 造成 bias。
        // 这里用 robin_hood::unordered_set（通常比 std::unordered_set 更快、内存更紧凑）。
        // reserve *2：降低 rehash 次数。
        std::unordered_set<hash_t> seen;
        seen.reserve(sketch_size * 2 + 1);

        // 主循环：对序列做 rolling，生成每个位置的 k-mer（正向/反向互补），然后计算 hash。
        for (std::size_t i = 0; i < seq.size(); ++i)
        {
            // nt4_table：把碱基映射到 0/1/2/3；其他字符（例如 N）映射到 4。
            const uint8_t c = nt4_table[static_cast<unsigned char>(seq[i])];
            if (c >= 4)
            {
                // 遇到 N/非法字符：当前 rolling 窗口失效，需要从下一个字符重新累计。
                fwd = rev = 0;
                valid = 0;
                continue;
            }

            // forward rolling：左移 2 位并加入新碱基（低位），再用 mask 截断到 2*k 位。
            fwd = ((fwd << 2) | c) & mask;
            // reverse-complement rolling：右移 2 位，并把“互补碱基”放在最高位。
            // 互补关系：A(0)<->T(3), C(1)<->G(2)，对应 3^c。
            rev = (rev >> 2) | (std::uint64_t(3U ^ c) << shift);

            // valid 计数不足 k 时，还不能形成完整 k-mer。
            if (valid < k) ++valid;
            if (valid < k)
            {
                continue;
            }

            // canonical 选择：
            // - noncanonical=true ：只取正向
            // - noncanonical=false：取 min(fwd, rev) 作为 canonical k-mer
            const std::uint64_t code = noncanonical ? fwd : std::min(fwd, rev);

            // 把 2-bit 编码的 k-mer 直接 hash 成 64-bit。
            // seed 用于扰动 hash（不同 seed 产生不同 sketch）。
            const hash_t h = getHash2bit(code, static_cast<std::uint32_t>(seed));

            // 去重：如果之前见过这个 hash，则跳过。
            if (!seen.insert(h).second)
            {
                continue;
            }

            if (maxHeap.size() < sketch_size)
            {
                // heap 未满：直接插入
                maxHeap.push(h);
            }
            else if (h < maxHeap.top())
            {
                // heap 已满：只有当新 hash 更小，才有资格进入 bottom-k
                // 先把当前最大值弹出，同时从 seen 中移除；再插入新值。
                seen.erase(maxHeap.top());
                maxHeap.pop();
                maxHeap.push(h);
            }
            // else: 新 hash 不够小，直接丢弃（但注意：此时 seen 已经插入过，需要撤销）
            else
            {
                seen.erase(h);
            }
        }

        // 导出 heap 到向量（此时顺序不保证）
        sk.hashes.reserve(maxHeap.size());
        while (!maxHeap.empty())
        {
            sk.hashes.push_back(maxHeap.top());
            maxHeap.pop();
        }

        // 最后排序是必须的：
        // - intersection/jaccard 需要 sorted unique
        // - 其他地方也依赖 sketch 有序
        std::sort(sk.hashes.begin(), sk.hashes.end());

        return std::move(sk);
    }

    // jaccard: 基于两个已排序且唯一化的 sketch 向量计算 Jaccard 相似度
    // 返回范围在 [0,1]，特殊 case：两个空集合返回 1.0（定义上的选择，表示完全相同的空集）
    double jaccard(const Sketch& a, const Sketch& b)
    {
        if (a.k != b.k) throw std::invalid_argument("mash::jaccard: mismatched k");

        if (a.hashes.empty() && b.hashes.empty()) return 1.0;
        if (a.hashes.empty() || b.hashes.empty()) return 0.0;

        const std::size_t inter = intersectionSizeSortedUnique(a.hashes, b.hashes);
        const std::size_t uni = std::min(a.hashes.size(), b.hashes.size());
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
    // 注意：本项目中 hash_t = uint64_t，因此无需关心 use64 之类的分支。
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

    bloom_filter filterFromSketch(const Sketch& sk, double false_positive_rate, int seed)
    {
        // ------------------------------------------------------------
        // 基于 sketch 构建 Bloom Filter
        // ------------------------------------------------------------
        // 设计目标：
        // 1) 让后续的 contains(hash) 查询尽可能快；
        // 2) false_positive_rate 可控（近似）；
        // 3) 对空 sketch 友好：返回一个默认构造的 bloom_filter（operator!() 为 true）。
        //
        // bloom_filter.hpp 的参数模型：
        // - projected_element_count: 预计插入的元素数量
        // - false_positive_probability: 目标误判率
        // - random_seed: 影响盐值生成，从而影响 hash 函数族
        //
        // 注意：这里插入的是 hash_t（uint64_t）本身，也就是把“已经 hash 好的值”再当做 key。
        // 这和直接插入原始 k-mer 序列不同，但对于集合成员查询来说是合理的。
        // ------------------------------------------------------------

        bloom_filter bf;

        if (sk.hashes.empty())
        {
            return bf;
        }

        bloom_parameters p;
        p.projected_element_count = static_cast<unsigned long long>(sk.hashes.size());
        p.false_positive_probability = false_positive_rate;
        // 避免 random_seed 落到 bloom_filter.hpp 中的非法范围（0 或全 1）
        {
            const std::uint64_t s = static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed));
            p.random_seed = (s == 0ULL) ? 0xA5A5A5A55A5A5A5AULL : (0xA5A5A5A55A5A5A5AULL ^ (s * 0x9E3779B97F4A7C15ULL));
            if (p.random_seed == 0ULL) p.random_seed = 0x1ULL;
            if (p.random_seed == 0xFFFFFFFFFFFFFFFFULL) p.random_seed = 0xFFFFFFFFFFFFFFFEULL;
        }

        // 计算 bloom filter 的最优参数（bit 数与 hash 函数个数）
        if (!p.compute_optimal_parameters())
        {
            // 参数非法/计算失败时，退化为默认构造（空）
            return bloom_filter{};
        }

        bf = bloom_filter(p);

        // 插入所有 hash 值
        for (const hash_t hv : sk.hashes)
        {
            bf.insert(hv);
        }

        return bf;
    }

    double jaccard(const bloom_filter& a, const Sketch& b)
    {
        // ------------------------------------------------------------
        // BloomFilter vs Sketch 的 Jaccard（近似）
        // ------------------------------------------------------------
        // 我们在 Sketch-Sketch 的实现里使用： inter / min(|A|,|B|)
        // 这是 Mash/MinHash 场景常用的近似定义（当两边都是 bottom-k 集时）。
        //
        // 这里 a 是 BloomFilter（由某个 sketch 构建），它本身不包含“集合大小信息”，
        // 因此：
        // - 交集：通过对 b.hashes 中每个 hv 做 contains 查询来估计
        // - 并集大小：仍使用 min(|A|,|B|) 的约定，但 |A| 取 a.element_count()
        //   （构建时插入的元素数，注意它可能大于真实 unique 数，但本项目的 Sketch 本身是 unique）
        //
        // 边界：
        // - 两者都空 => 1.0
        // - 任一为空 => 0.0
        // ------------------------------------------------------------

        const std::size_t asz = static_cast<std::size_t>(a.element_count());
        const std::size_t bsz = b.hashes.size();

        if (asz == 0 && bsz == 0) return 1.0;
        if (asz == 0 || bsz == 0) return 0.0;

        std::size_t inter = 0;
        for (const hash_t hv : b.hashes)
        {
            if (a.contains(hv)) ++inter;
        }

        const std::size_t uni = std::min(asz, bsz);
        if (uni == 0) return 1.0;
        return static_cast<double>(inter) / static_cast<double>(uni);
    }

} // namespace mash

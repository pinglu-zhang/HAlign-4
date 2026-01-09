#include <doctest/doctest.h>

#include <mash.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// ------------------------------------------------------------
// test_jaccard.cpp
// ------------------------------------------------------------
// 该测试文件专门验证两种 Jaccard 计算方式的正确性与性能：
// 1) mash::jaccard(const Sketch&, const Sketch&)
//    - 这是精确的集合交集计数（要求 hashes 已经排序且 unique）。
// 2) mash::jaccard(const bloom_filter&, const Sketch&)
//    - 这是用 Bloom Filter 近似集合成员查询得到的“交集估计”，会受 false positive 影响。
//
// 你要求的性能评估：
// - 固定一个 ref（参考）sketch
// - 构造 10000 个 query sketch（每个大小 2k）
// - 评估：10000 次 "query vs ref" 的 jaccard 计算耗时
//
// 注意：性能用例默认跳过，避免 CI/普通测试耗时。
// 需要时设置环境变量：HALIGN4_RUN_PERF=1
// ------------------------------------------------------------

namespace {

// 生成随机 DNA 序列（A/C/G/T）。
static std::string random_dna(std::mt19937_64& rng, std::size_t len)
{
    static constexpr char bases[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(len);
    for (std::size_t i = 0; i < len; ++i)
        s.push_back(bases[rng() & 3ULL]);
    return s;
}

} // namespace

TEST_SUITE("jaccard") {
// 注意：某些环境下 doctest 的 TEST_SUITE 宏对换行/大括号位置较敏感，
// 这里采用与本项目其他测试文件一致的写法。

    TEST_CASE("Sketch-Sketch jaccard: basic correctness")
    {
        mash::Sketch a;
        mash::Sketch b;
        a.k = b.k = 21;

        SUBCASE("both empty")
        {
            CHECK(mash::jaccard(a, b) == doctest::Approx(1.0));
        }

        SUBCASE("one empty")
        {
            a.hashes = {1, 2, 3};
            CHECK(mash::jaccard(a, b) == doctest::Approx(0.0));
            CHECK(mash::jaccard(b, a) == doctest::Approx(0.0));
        }

        SUBCASE("identical")
        {
            a.hashes = {10, 20, 30};
            b.hashes = {10, 20, 30};
            CHECK(mash::jaccard(a, b) == doctest::Approx(1.0));
        }

        SUBCASE("partial overlap - uses min(|A|,|B|) denominator")
        {
            // 交集=2, min size=3 => 2/3
            a.hashes = {1, 2, 3};
            b.hashes = {2, 3, 4, 5};
            CHECK(mash::jaccard(a, b) == doctest::Approx(2.0 / 3.0));
        }
    }

    TEST_CASE("Sketch-Sketch jaccard: sketchFromSequence identical sequences should be 1")
    {
        // 如果你遇到“总是 0”的情况，最常见原因是：
        // 1) sketch.hashes 没有排序/去重，导致 intersectionSizeSortedUnique 算法失效；
        // 2) 两个 sketch 的 k 不一致（会抛异常）；
        // 3) seed/noncanonical 参数不一致。
        //
        // 这个用例直接用相同序列、相同参数生成两个 sketch，jaccard 必须接近 1。
        const std::size_t k = 21;
        const std::size_t sketch_size = 2000;
        const int seed = 42;

        const std::string s = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        auto a = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/true, seed);
        auto b = mash::sketchFromSequence(s, k, sketch_size, /*noncanonical*/true, seed);

        CHECK(a.k == k);
        CHECK(b.k == k);
        CHECK(std::is_sorted(a.hashes.begin(), a.hashes.end()));
        CHECK(std::is_sorted(b.hashes.begin(), b.hashes.end()));
        CHECK(std::adjacent_find(a.hashes.begin(), a.hashes.end()) == a.hashes.end());
        CHECK(std::adjacent_find(b.hashes.begin(), b.hashes.end()) == b.hashes.end());

        const auto inter = mash::intersectionSizeSortedUnique(a.hashes, b.hashes);
        const double j = mash::jaccard(a, b);
        MESSAGE("a.size=", a.size(), " b.size=", b.size(), " inter=", inter, " jaccard=", j);


        CHECK(j == doctest::Approx(1.0));
    }

    TEST_CASE("Sketch-Sketch jaccard: similar sequences should usually be > 0")
    {
        // 这个测试用来捕捉你描述的“经常只返回 0”。
        // 说明：对于随机序列，k=21 且 sketch_size=2000 时，两个不同序列的交集通常非常小，
        // 返回 0 其实是“正常现象”（底层 MinHash 交集为 0，表示估计的 Jaccard 很低）。
        // 因此这里用“相似但不完全相同”的序列，确保应该存在共享 k-mer，从而 jaccard 不应为 0。

        const std::size_t k = 21;
        const std::size_t sketch_size = 2000;
        const int seed = 11;

        std::string s1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        std::string s2 = s1;
        // 修改中间的一小段，保持大部分 k-mer 仍然相同
        if (s2.size() >= 10) {
            s2[5] = 'T';
            s2[6] = 'T';
            s2[7] = 'T';
        }

        auto a = mash::sketchFromSequence(s1, k, sketch_size, /*noncanonical*/true, seed);
        auto b = mash::sketchFromSequence(s2, k, sketch_size, /*noncanonical*/true, seed);

        const double j = mash::jaccard(a, b);
        MESSAGE("similar sequences jaccard=", j, " a.size=", a.size(), " b.size=", b.size());

        CHECK(j > 0.0);
    }

    TEST_CASE("BloomFilter-Sketch jaccard: exactness when FPP very small")
    {
        // 这里我们设置极低的 false positive rate，
        // 在测试规模较小的情况下，BloomFilter 查询几乎等价于精确集合包含。
        const std::size_t k = 21;

        mash::Sketch ref;
        ref.k = k;
        ref.hashes = {1, 2, 3, 100, 200, 500};

        mash::Sketch q1;
        q1.k = k;
        q1.hashes = {2, 3, 4, 5};

        // 构建 bloom filter（极低误判率）
        const double fpp = 1e-9;
        const int seed = 12345;
        auto bf = mash::filterFromSketch(ref, fpp, seed);

        // 由于 BloomFilter 会近似计算交集，所以我们只要求它“接近” Sketch-Sketch 的结果。
        const double j_exact = mash::jaccard(ref, q1);
        const double j_bf = mash::jaccard(bf, q1);

        // 在极小 fpp 下，结果应非常接近（允许极少量误差）。
        CHECK(j_bf == doctest::Approx(j_exact).epsilon(1e-6));
    }

    TEST_CASE("BloomFilter-Sketch jaccard: monotonic sanity (superset should not reduce) ")
    {
        // BloomFilter 的误判只会让 contains 偏大（不会把 true 变 false），
        // 所以对于同一个 BloomFilter：
        // - 如果 query 增加更多元素（原来已有的仍包含），估计交集不应减少。
        mash::Sketch ref;
        ref.k = 21;
        ref.hashes = {10, 20, 30, 40, 50};

        auto bf = mash::filterFromSketch(ref, /*fpp*/1e-6, /*seed*/7);

        mash::Sketch q_small;
        q_small.k = 21;
        q_small.hashes = {10, 999};

        mash::Sketch q_large;
        q_large.k = 21;
        q_large.hashes = {10, 20, 999, 888};

        const double j1 = mash::jaccard(bf, q_small);
        const double j2 = mash::jaccard(bf, q_large);

        CHECK(j2 + 1e-12 >= j1);
    }

    TEST_CASE("jaccard_perf: fixed ref vs 10000 queries (sketch size=2k)")
    {
        const char* env = std::getenv("HALIGN4_RUN_PERF");
        if (!env || std::string(env) != "1")
        {
            DOCTEST_INFO("jaccard_perf skipped; set HALIGN4_RUN_PERF=1 to enable");
            return;
        }

        // 固定参数：每个 sketch 2k
        const std::size_t k = 21;
        const std::size_t sketch_size = 2000;
        const std::size_t seq_len = 300000;       // 默认 3w
        const std::size_t num_queries = 10000;   // 1w
        const int seed = 42;

        std::mt19937_64 rng(123456);

        // 生成 ref
        const std::string ref_seq = random_dna(rng, seq_len);
        mash::Sketch ref = mash::sketchFromSequence(ref_seq, k, sketch_size, /*noncanonical*/true, seed);
        REQUIRE(ref.size() == ref.hashes.size());

        // ref 的 bloom filter（用于加速 query vs ref 的近似 jaccard）
        bloom_filter ref_bf = mash::filterFromSketch(ref, /*fpp*/1e-8, /*seed*/seed);

        // 生成 query sketches（先生成序列再 sketch，保持测试更贴近真实 workload）
        std::vector<mash::Sketch> queries;
        queries.reserve(num_queries);
        for (std::size_t i = 0; i < num_queries; ++i)
        {
            std::string qseq = random_dna(rng, seq_len);
            queries.emplace_back(mash::sketchFromSequence(qseq, k, sketch_size, /*noncanonical*/true, seed));
        }

        // --- 性能测试 1：Sketch-Sketch ---
        {
            volatile double sink = 0.0;
            auto t0 = std::chrono::steady_clock::now();
            for (const auto& q : queries)
            {
                sink += mash::jaccard(ref, q);
            }
            auto t1 = std::chrono::steady_clock::now();
            const double sec = std::chrono::duration<double>(t1 - t0).count();
            MESSAGE("jaccard_perf Sketch-Sketch: queries=" << num_queries << " sketch_size=" << sketch_size << " took " << sec << " s" << " sink=" << sink);
        }

        // --- 性能测试 2：BloomFilter-Sketch ---
        {
            volatile double sink = 0.0;
            auto t0 = std::chrono::steady_clock::now();
            for (const auto& q : queries)
            {
                sink += mash::jaccard(ref_bf, q);
            }
            auto t1 = std::chrono::steady_clock::now();
            const double sec = std::chrono::duration<double>(t1 - t0).count();
            MESSAGE("jaccard_perf BloomFilter-Sketch: queries=" << num_queries << " sketch_size=" << sketch_size << " took " << sec << " s" << " sink=" << sink);
        }

        CHECK(true);
    }
}

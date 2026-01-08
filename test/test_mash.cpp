#include <doctest/doctest.h>

#include <mash.h>

#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <iostream>

TEST_SUITE("mash")
{
    TEST_CASE("intersectionSizeSortedUnique")
    {
        std::vector<hash_t> a{1, 2, 3, 10};
        std::vector<hash_t> b{2, 3, 4, 5, 10, 11};
        CHECK(mash::intersectionSizeSortedUnique(a, b) == 3);
    }

    TEST_CASE("sketchFromSequence - empty & sketch_size=0")
    {
        const std::size_t k = 15;
        const std::size_t sketch_size = 0;

        auto sk = mash::sketchFromSequence(std::string("ACGTACGT"), k, sketch_size);
        CHECK(sk.empty());

        auto sk2 = mash::sketchFromSequence(std::string(""), k, 100);
        CHECK(sk2.empty());
    }

    TEST_CASE("jaccard - empty sets")
    {
        mash::Sketch a;
        a.k = 15;
        mash::Sketch b;
        b.k = 15;
        CHECK(mash::jaccard(a, b) == doctest::Approx(1.0));

        mash::Sketch c;
        c.k = 15;
        c.hashes = {1, 2, 3};
        CHECK(mash::jaccard(a, c) == doctest::Approx(0.0));
    }

    TEST_CASE("jaccard - identical sketches")
    {
        mash::Sketch a;
        a.k = 21;
        a.hashes = {1, 2, 3, 4};
        mash::Sketch b;
        b.k = 21;
        b.hashes = {1, 2, 3, 4};
        CHECK(mash::jaccard(a, b) == doctest::Approx(1.0));
        CHECK(mash::mashDistanceFromJaccard(1.0, 21) == doctest::Approx(0.0));
        CHECK(mash::aniFromJaccard(1.0, 21) == doctest::Approx(1.0));
    }

    TEST_CASE("jaccard - disjoint sketches")
    {
        mash::Sketch a;
        a.k = 21;
        a.hashes = {1, 2, 3};
        mash::Sketch b;
        b.k = 21;
        b.hashes = {4, 5, 6};
        CHECK(mash::jaccard(a, b) == doctest::Approx(0.0));
        CHECK(!std::isfinite(mash::mashDistanceFromJaccard(0.0, 21)));
        CHECK(mash::aniFromJaccard(0.0, 21) == doctest::Approx(0.0));
    }

    TEST_CASE("sketchFromSequence + jaccard smoke")
    {
        const std::size_t k = 15;
        const std::size_t sketch_size = 200;

        const std::string s1 = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        const std::string s2 = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        const std::string s3 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

        auto sk1 = mash::sketchFromSequence(s1, k, sketch_size);
        auto sk2 = mash::sketchFromSequence(s2, k, sketch_size);
        auto sk3 = mash::sketchFromSequence(s3, k, sketch_size);

        CHECK(sk1.k == k);
        CHECK(sk1.size() <= sketch_size);

        const double j12 = mash::jaccard(sk1, sk2);
        CHECK(j12 == doctest::Approx(1.0));

        const double j13 = mash::jaccard(sk1, sk3);
        CHECK(j13 >= 0.0);
        CHECK(j13 <= 1.0);
    }

    // Performance test: only run when HALIGN4_RUN_PERF=1 (set by test/run_tests.sh --perf)
    static std::string random_dna(std::mt19937_64 &rng, std::size_t len) {
        static const char bases[4] = {'A','C','G','T'};
        std::string s;
        s.reserve(len);
        for (std::size_t i = 0; i < len; ++i) s.push_back(bases[rng() & 3]);
        return s;
    }

    TEST_CASE("mash_perf") {
        const char* env = std::getenv("HALIGN4_RUN_PERF");
        if (!env || std::string(env) != "1") {
            DOCTEST_INFO("mash_perf skipped; run run_tests.sh --perf to enable");
            return;
        }

        const std::size_t N = []{
            const char* e = std::getenv("MASH_PERF_N");
            return e ? static_cast<std::size_t>(std::stoull(e)) : 30000ULL;
        }();
        const std::size_t L = []{
            const char* e = std::getenv("MASH_PERF_L");
            return e ? static_cast<std::size_t>(std::stoull(e)) : 30000ULL;
        }();
        const std::size_t k = 21;
        const std::size_t sketch_size = 2000;

        std::mt19937_64 rng(123456);
        std::vector<std::string> seqs;
        seqs.reserve(N);
        for (std::size_t i = 0; i < N; ++i) seqs.push_back(random_dna(rng, L));

        auto t0 = std::chrono::steady_clock::now();
        for (const auto &s : seqs) {
            volatile auto sk = mash::sketchFromSequence(s, k, sketch_size);
            (void)sk;
        }
        auto t1 = std::chrono::steady_clock::now();

        const double seconds = std::chrono::duration<double>(t1 - t0).count();
        MESSAGE("mash_perf: N=" << N << " L=" << L << " took " << seconds << "s");

        CHECK(true);
    }
}

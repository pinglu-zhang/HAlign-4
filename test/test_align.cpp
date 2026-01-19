#include <doctest/doctest.h>
#include <chrono>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "align.h"

// ------------------------------------------------------------------
// 辅助函数：生成随机 DNA 序列
// ------------------------------------------------------------------
static std::string generateRandomDNA(size_t length, unsigned seed = 42) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, 3);

    std::string seq;
    seq.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        seq += bases[dist(rng)];
    }
    return seq;
}

// ------------------------------------------------------------------
// 辅助函数：在序列中引入随机突变（SNP + Indel）
// ------------------------------------------------------------------
static std::string mutateSequence(const std::string& ref,
                                  double snp_rate = 0.01,     // 1% SNP
                                  double indel_rate = 0.005,  // 0.5% indel
                                  unsigned seed = 43) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_int_distribution<int> base_dist(0, 3);
    std::uniform_int_distribution<int> indel_len(1, 5);  // indel 长度 1-5bp

    std::string mutated;
    mutated.reserve(ref.size() * 1.1);  // 预留空间

    for (size_t i = 0; i < ref.size(); ++i) {
        double p = prob(rng);

        if (p < indel_rate) {
            // 随机插入或删除
            if (prob(rng) < 0.5) {
                // 插入
                int len = indel_len(rng);
                for (int j = 0; j < len; ++j) {
                    mutated += bases[base_dist(rng)];
                }
            } else {
                // 删除：跳过当前碱基
                continue;
            }
        } else if (p < indel_rate + snp_rate) {
            // SNP：替换为不同的碱基
            char original = ref[i];
            char replacement;
            do {
                replacement = bases[base_dist(rng)];
            } while (replacement == original);
            mutated += replacement;
        } else {
            // 保持不变
            mutated += ref[i];
        }
    }

    return mutated;
}

// ------------------------------------------------------------------
// 辅助函数：CIGAR 转字符串（便于调试）
// ------------------------------------------------------------------
static std::string cigarToString(const cigar::Cigar_t& cigar) {
    std::string result;
    for (auto unit : cigar) {
        char op;
        uint32_t len;
        cigar::intToCigar(unit, op, len);
        result += std::to_string(len) + op;
    }
    return result;
}

// ------------------------------------------------------------------
// 性能测试辅助：计时器
// ------------------------------------------------------------------
struct Timer {
    std::chrono::high_resolution_clock::time_point start;

    Timer() : start(std::chrono::high_resolution_clock::now()) {}

    double elapsedMs() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start).count();
    }
};

// ------------------------------------------------------------------
// 测试套件：正确性测试
// ------------------------------------------------------------------
TEST_SUITE("align") {

    TEST_CASE("globalAlignKSW2 - 精确匹配") {
        std::string seq = "ACGTACGTACGT";
        auto cigar = align::globalAlignKSW2(seq, seq);

        // 应该是完全匹配，CIGAR 应该是 12M
        REQUIRE(cigar.size() >= 1);
        char op;
        uint32_t len;
        cigar::intToCigar(cigar[0], op, len);
        CHECK(op == 'M');
        CHECK(len == 12);
    }

    TEST_CASE("globalAlignKSW2 - 单个错配") {
        std::string ref   = "ACGTACGTACGT";
        std::string query = "ACGTACCGTACGT";  // 第6位 T->C 错配

        auto cigar = align::globalAlignKSW2(ref, query);
        REQUIRE(cigar.size() > 0);

        // 验证 CIGAR 不为空
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("globalAlignKSW2 - 单个插入") {
        std::string ref   = "ACGTACGTACGT";
        std::string query = "ACGTAACGTACGT";  // 在第5位后插入 A

        auto cigar = align::globalAlignKSW2(ref, query);
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR: ", cigar_str);

        // 应该包含插入操作（I）
        bool has_insertion = cigar_str.find('I') != std::string::npos;
        CHECK(has_insertion);
    }

    TEST_CASE("extendAlignKSW2 - 基本延伸") {
        std::string ref = generateRandomDNA(1000, 100);
        std::string query = ref.substr(100, 500);  // 提取中间片段

        auto cigar = align::extendAlignKSW2(ref, query, 200);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("Extend CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("globalAlignWFA2 - 精确匹配") {
        std::string seq = "ACGTACGTACGT";
        auto cigar = align::globalAlignWFA2(seq, seq);

        REQUIRE(cigar.size() >= 1);
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("WFA2 CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("空序列边界测试") {
        std::string empty = "";
        std::string seq = "ACGT";

        // 测试空序列不崩溃（具体行为依赖于实现）
        CHECK_NOTHROW(align::globalAlignKSW2(empty, seq));
        CHECK_NOTHROW(align::globalAlignKSW2(seq, empty));
    }

    TEST_CASE("高相似度序列 - 99% 相似") {
        // 生成一个 1000bp 的序列，引入 1% 的错配
        std::string ref = generateRandomDNA(1000, 12345);
        std::string query = mutateSequence(ref, 0.01, 0.0, 12346);  // 只有 SNP，无 indel

        // 三种方法都应该能正确比对
        auto cigar_ksw2 = align::globalAlignKSW2(ref, query);
        auto cigar_extend = align::extendAlignKSW2(ref, query, 200);
        auto cigar_wfa2 = align::globalAlignWFA2(ref, query);

        // 验证 CIGAR 都不为空
        CHECK(cigar_ksw2.size() > 0);
        CHECK(cigar_extend.size() > 0);
        CHECK(cigar_wfa2.size() > 0);

        MESSAGE("KSW2 CIGAR: ", cigarToString(cigar_ksw2));
        MESSAGE("Extend CIGAR: ", cigarToString(cigar_extend));
        MESSAGE("WFA2 CIGAR: ", cigarToString(cigar_wfa2));
    }

    TEST_CASE("高相似度序列 - 98% 相似（含 indel）") {
        // 生成序列，引入 1% SNP + 1% indel
        std::string ref = generateRandomDNA(500, 54321);
        std::string query = mutateSequence(ref, 0.01, 0.01, 54322);

        auto cigar_ksw2 = align::globalAlignKSW2(ref, query);
        auto cigar_wfa2 = align::globalAlignWFA2(ref, query);

        // 验证 CIGAR 包含不同类型的操作
        std::string cigar_str_ksw2 = cigarToString(cigar_ksw2);
        std::string cigar_str_wfa2 = cigarToString(cigar_wfa2);

        MESSAGE("KSW2 CIGAR (with indels): ", cigar_str_ksw2);
        MESSAGE("WFA2 CIGAR (with indels): ", cigar_str_wfa2);

        CHECK(!cigar_str_ksw2.empty());
        CHECK(!cigar_str_wfa2.empty());

        // 应该包含匹配操作
        bool has_match_ksw2 = cigar_str_ksw2.find('M') != std::string::npos;
        bool has_match_wfa2 = cigar_str_wfa2.find('M') != std::string::npos;
        CHECK(has_match_ksw2);
        CHECK(has_match_wfa2);
    }

    TEST_CASE("极高相似度序列 - 99.9% 相似") {
        // 模拟测序错误：仅 0.1% 的错误率
        std::string ref = generateRandomDNA(10000, 99999);
        std::string query = mutateSequence(ref, 0.0005, 0.0005, 100000);

        // 在如此高的相似度下，所有方法都应该快速完成
        Timer timer;
        auto cigar = align::globalAlignKSW2(ref, query);
        double elapsed = timer.elapsedMs();

        CHECK(cigar.size() > 0);
        MESSAGE("极高相似度 10k 序列比对耗时: ", elapsed, " ms");

        // 对于 99.9% 相似的 10k 序列，应该在合理时间内完成（<100ms）
        CHECK(elapsed < 100.0);
    }
}

// ------------------------------------------------------------------
// 性能测试套件
// ------------------------------------------------------------------
TEST_SUITE("align_perf") {

    // ------------------------------------------------------------------
    // 性能测试：短序列（~100bp）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Short sequences (~100bp)") {
        constexpr int NUM_RUNS = 1000;
        constexpr size_t SEQ_LEN = 100;

        std::cout << "\n========== 短序列性能测试 (100bp, " << NUM_RUNS << " 次) ==========\n";

        // 生成测试数据
        std::vector<std::pair<std::string, std::string>> test_pairs;
        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);
        }

        // 测试 globalAlignKSW2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;  // 防止优化掉
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        // 测试 extendAlignKSW2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        // 测试 globalAlignWFA2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：中等长度序列（~1000bp）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Medium sequences (~1000bp)") {
        constexpr int NUM_RUNS = 100;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== 中等序列性能测试 (1000bp, " << NUM_RUNS << " 次) ==========\n";

        std::vector<std::pair<std::string, std::string>> test_pairs;
        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：长序列（~10000bp）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Long sequences (~10000bp)") {
        constexpr int NUM_RUNS = 10;
        constexpr size_t SEQ_LEN = 10000;

        std::cout << "\n========== 长序列性能测试 (10000bp, " << NUM_RUNS << " 次) ==========\n";

        std::vector<std::pair<std::string, std::string>> test_pairs;
        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2: " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：不同突变率下的表现
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Different mutation rates") {
        constexpr int NUM_RUNS = 100;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== 不同突变率性能测试 (1000bp, " << NUM_RUNS << " 次) ==========\n";

        std::vector<double> snp_rates = {0.0, 0.01, 0.05, 0.10};

        for (double snp_rate : snp_rates) {
            std::cout << "\n  SNP率: " << (snp_rate * 100) << "%\n";

            std::vector<std::pair<std::string, std::string>> test_pairs;
            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
                std::string query = mutateSequence(ref, snp_rate, 0.005, i * 2 + 1);
                test_pairs.emplace_back(ref, query);
            }

            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                std::cout << "    KSW2:  " << std::fixed << std::setprecision(2)
                          << (elapsed / NUM_RUNS) << " ms/次\n";
            }

            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                std::cout << "    WFA2:  " << std::fixed << std::setprecision(2)
                          << (elapsed / NUM_RUNS) << " ms/次\n";
            }
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：高相似度序列（95%-99.9%，真实测序场景）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - High similarity sequences (real-world)") {
        constexpr int NUM_RUNS = 200;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== 高相似度序列性能测试 (1000bp, " << NUM_RUNS << " 次) ==========\n";
        std::cout << "说明：模拟真实基因组测序场景，序列相似度 >95%\n\n";

        // 测试不同相似度等级（通过控制 SNP + indel 率）
        struct SimilarityLevel {
            double snp_rate;
            double indel_rate;
            const char* desc;
            double similarity;  // 预期相似度
        };

        std::vector<SimilarityLevel> levels = {
            {0.0001, 0.0001, "极高相似度 (>99.9%)", 99.98},
            {0.001,  0.0005, "很高相似度 (~99%)",   99.0},
            {0.005,  0.002,  "高相似度 (~98%)",     98.0},
            {0.01,   0.005,  "中等相似度 (~97%)",   97.0},
            {0.03,   0.01,   "较低相似度 (~95%)",   95.0}
        };

        for (const auto& level : levels) {
            std::cout << "---------- " << level.desc << " ----------\n";

            // 生成测试数据
            std::vector<std::pair<std::string, std::string>> test_pairs;
            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
                std::string query = mutateSequence(ref, level.snp_rate, level.indel_rate, i * 2 + 1);
                test_pairs.emplace_back(ref, query);
            }

            // KSW2 测试
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  KSW2 全局比对:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/次  (吞吐: " << (1000.0 / avg_ms)
                          << " 次/秒)\n";
            }

            // KSW2 延伸模式测试
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  KSW2 延伸模式:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/次  (吞吐: " << (1000.0 / avg_ms)
                          << " 次/秒)\n";
            }

            // WFA2 测试
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  WFA2 全局比对:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/次  (吞吐: " << (1000.0 / avg_ms)
                          << " 次/秒)\n";
            }

            std::cout << "\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：长度敏感性（高相似度，不同长度）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Length scaling (high similarity)") {
        constexpr int NUM_RUNS = 50;
        constexpr double SNP_RATE = 0.01;    // 1% SNP（~98% 相似度）
        constexpr double INDEL_RATE = 0.005; // 0.5% indel

        std::cout << "\n========== 长度扩展性测试 (相似度 ~98%, " << NUM_RUNS << " 次) ==========\n";

        std::vector<size_t> lengths = {100, 500, 1000, 5000, 10000};

        std::cout << std::setw(10) << "长度(bp)"
                  << std::setw(15) << "KSW2(ms)"
                  << std::setw(15) << "Extend(ms)"
                  << std::setw(15) << "WFA2(ms)" << "\n";
        std::cout << std::string(55, '-') << "\n";

        for (size_t len : lengths) {
            // 生成测试数据
            std::vector<std::pair<std::string, std::string>> test_pairs;
            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(len, i * 2);
                std::string query = mutateSequence(ref, SNP_RATE, INDEL_RATE, i * 2 + 1);
                test_pairs.emplace_back(ref, query);
            }

            double ksw2_time = 0.0;
            double extend_time = 0.0;
            double wfa2_time = 0.0;

            // KSW2
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                ksw2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // Extend
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                extend_time = timer.elapsedMs() / NUM_RUNS;
            }

            // WFA2
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                wfa2_time = timer.elapsedMs() / NUM_RUNS;
            }

            std::cout << std::setw(10) << len
                      << std::setw(15) << std::fixed << std::setprecision(3) << ksw2_time
                      << std::setw(15) << extend_time
                      << std::setw(15) << wfa2_time << "\n";
        }

        std::cout << "========================================================\n\n";
    }

    TEST_CASE("removeRefGapColumns - drop ref gap columns from aligned seq") {
        // 说明：测试"按 ref_gap_pos 删列"的纯过滤功能（原地修改）。
        // 输入序列已经是对齐后的（含 gap），本函数只负责删除参考为 gap 的列。

        align::RefAligner aligner(".", "example/data/covid-ref.fasta", 21, 10, 2000, true, 1, "", true, false);

        // 输入：已对齐序列 "AC-GT"（长度 5）
        std::string seq = "AC-GT";

        // ref_gap_pos: true 表示参考该列是 gap，需要删除该列。
        // 我们删除第2列(0-based==2)这一列，输出应为 "ACGT"
        const std::vector<bool> ref_gap_pos = {false, false, true, false, false};

        aligner.removeRefGapColumns(seq, ref_gap_pos);
        CHECK(seq == "ACGT");
    }

    TEST_CASE("removeRefGapColumns - keep existing '-' as base when not in ref gap pos") {
        // 说明：输入序列本身含有 '-'，但只要 ref_gap_pos 该列为 false，就保留。
        align::RefAligner aligner(".", "example/data/covid-ref.fasta", 21, 10, 2000, true, 1, "", true, false);

        std::string seq = "A-CG"; // 包含 1 个原生 '-'
        const std::vector<bool> ref_gap_pos = {false, false, false, false}; // 所有列都保留

        aligner.removeRefGapColumns(seq, ref_gap_pos);
        CHECK(seq == "A-CG");
    }

    TEST_CASE("removeRefGapColumns - empty ref_gap_pos means no change") {
        align::RefAligner aligner(".", "example/data/covid-ref.fasta", 21, 10, 2000, true, 1, "", true, false);

        std::string seq = "A--A";

        // ref_gap_pos 为空时，seq 保持不变
        aligner.removeRefGapColumns(seq, {});
        CHECK(seq == "A--A");
    }
}

// ------------------------------------------------------------------
// 一键运行性能测试的辅助说明
// ------------------------------------------------------------------
// 说明：
// 1. 正确性测试默认运行，用于验证基本功能
// 2. 性能测试默认跳过（doctest::skip(true)），需要手动启用
// 3. 运行性能测试：
//    ./halign4_tests -tc="*Performance*" --no-skip
// 4. 只运行特定性能测试：
//    ./halign4_tests -tc="*Short*" --no-skip
// ------------------------------------------------------------------

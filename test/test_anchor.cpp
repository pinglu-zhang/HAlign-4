// ==================================================================
// test_anchor.cpp - anchor 模块单元测试
// ==================================================================
//
// 测试覆盖：
// 1. collect_anchors：锚点收集功能（基于 minimizer hits）
// 2. chainAnchors：链化算法（DP 动态规划）
// 3. 过滤参数：q_occ_frac, f_top_frac, sample_every_bp
// ==================================================================

#include <doctest/doctest.h>
#include "seed.h"
#include "anchor.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

using hash_t = std::uint64_t;

// ==================================================================
// 辅助函数
// ==================================================================

// 创建 MinimizerHit（用于测试 collect_anchors）
static minimizer::MinimizerHit makeHit(hash_t hash56, std::uint32_t pos,
                                       std::uint32_t rid = 0, bool strand = false, std::uint32_t span = 15)
{
    // 使用便捷构造函数：MinimizerHit(hash56, pos, rid, strand, span)
    return minimizer::MinimizerHit(hash56, pos, rid, strand, static_cast<std::uint8_t>(span));
}

// 创建 Anchor（用于测试 chainAnchors）
static anchor::Anchor makeAnchor(hash_t hash, std::uint32_t pos_ref, std::uint32_t pos_qry,
                                  std::uint32_t rid_ref = 0, std::uint32_t rid_qry = 0,
                                  bool is_rev = false, std::uint32_t span = 15)
{
    anchor::Anchor a;
    a.hash = hash;
    a.rid_ref = rid_ref;
    a.pos_ref = pos_ref;
    a.rid_qry = rid_qry;
    a.pos_qry = pos_qry;
    a.span = span;
    a.is_rev = is_rev;
    return a;
}

// ==================================================================
// TEST SUITE: collect_anchors - 锚点收集测试
// ==================================================================

TEST_SUITE("anchor")
{
    // ------------------------------------------------------------------
    // 基础功能测试
    // ------------------------------------------------------------------

    TEST_CASE("collect_anchors - 空输入返回空")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - ref 为空时无锚点")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;
        qry_hits.push_back(makeHit(0x111111, 100));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - qry 为空时无锚点")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;
        ref_hits.push_back(makeHit(0x111111, 100));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - 单一完美匹配")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0x123456, 100, 0, false, 20));
        qry_hits.push_back(makeHit(0x123456, 50, 0, false, 20));

        // 禁用所有过滤
        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].hash == 0x123456);
        CHECK(anchors[0].pos_ref == 100);
        CHECK(anchors[0].pos_qry == 50);
        CHECK(anchors[0].span == 20);
        CHECK(anchors[0].is_rev == false);
    }

    TEST_CASE("collect_anchors - 反向匹配检测")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0xABCDEF, 200, 0, false, 15));
        qry_hits.push_back(makeHit(0xABCDEF, 80, 0, true, 15));  // 反向

        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].is_rev == true);  // ref XOR qry = false XOR true = true
    }

    TEST_CASE("collect_anchors - 无共同 hash 时返回空")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0x111111, 100));
        ref_hits.push_back(makeHit(0x222222, 200));
        qry_hits.push_back(makeHit(0x333333, 50));
        qry_hits.push_back(makeHit(0x444444, 80));

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits);
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - 一对多展开（occurrence expansion）")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        // ref 端同一个 hash 出现 3 次
        ref_hits.push_back(makeHit(0x555555, 100));
        ref_hits.push_back(makeHit(0x555555, 200));
        ref_hits.push_back(makeHit(0x555555, 300));

        // qry 端该 hash 出现 1 次
        qry_hits.push_back(makeHit(0x555555, 50));

        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        // 应该生成 3 个锚点（1 qry × 3 ref）
        REQUIRE(anchors.size() == 3);

        // 验证所有锚点的 qry 位置相同，ref 位置不同
        CHECK(anchors[0].pos_qry == 50);
        CHECK(anchors[1].pos_qry == 50);
        CHECK(anchors[2].pos_qry == 50);

        // ref 位置应该是 100, 200, 300（顺序可能不同，因为排序）
        std::vector<std::uint32_t> ref_positions;
        for (const auto& a : anchors) {
            ref_positions.push_back(a.pos_ref);
        }
        std::sort(ref_positions.begin(), ref_positions.end());
        CHECK(ref_positions[0] == 100);
        CHECK(ref_positions[1] == 200);
        CHECK(ref_positions[2] == 300);
    }

    // ------------------------------------------------------------------
    // 过滤参数测试
    // ------------------------------------------------------------------

    TEST_CASE("collect_anchors - q_occ_frac 过滤高频 hash")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0x888888, 100));

        // qry 端有 100 个 hit，其中同一个 hash 出现 50 次
        for (std::uint32_t i = 0; i < 50; ++i) {
            qry_hits.push_back(makeHit(0x888888, i));  // 高频 hash
        }
        for (std::uint32_t i = 0; i < 50; ++i) {
            qry_hits.push_back(makeHit(0x999900 + i, i));  // 不同的 hash
        }

        // 设置 q_occ_frac = 10%（即 10 次），超过则丢弃
        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.10;  // 10% of 100 = 10
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        // 高频 hash (0x888888) 应该被过滤掉，因为 50 > 10
        CHECK(anchors.empty());
    }

    TEST_CASE("collect_anchors - span 取 min(ref.span, qry.span)")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0xAAAAAA, 100, 0, false, 30));  // span=30
        qry_hits.push_back(makeHit(0xAAAAAA, 50, 0, false, 20));   // span=20

        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].span == 20);  // min(30, 20) = 20
    }

    TEST_CASE("collect_anchors - 多序列 rid 正确传递")
    {
        minimizer::MinimizerHits ref_hits;
        minimizer::MinimizerHits qry_hits;

        ref_hits.push_back(makeHit(0xBBBBBB, 100, 2, false, 15));  // rid_ref=2
        qry_hits.push_back(makeHit(0xBBBBBB, 50, 3, false, 15));   // rid_qry=3

        anchor::SeedFilterParams params;
        params.q_occ_frac = 0.0;
        params.f_top_frac = 0.0;

        auto anchors = minimizer::collect_anchors(ref_hits, qry_hits, params);

        REQUIRE(anchors.size() == 1);
        CHECK(anchors[0].rid_ref == 2);
        CHECK(anchors[0].rid_qry == 3);
    }
}

// ==================================================================
// TEST SUITE: chainAnchors - 链化算法测试
// ==================================================================

TEST_SUITE("anchor")
{
    // ------------------------------------------------------------------
    // 基础功能测试
    // ------------------------------------------------------------------

    TEST_CASE("chainAnchors - 空输入返回空链")
    {
        anchor::Anchors anchors;
        auto best_chain = anchor::chainAnchors(anchors);
        CHECK(best_chain.empty());
    }

    TEST_CASE("chainAnchors - 单个锚点形成单链")
    {
        anchor::Anchors anchors;
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));

        anchor::ChainParams params;
        params.min_cnt = 1;  // 允许单锚点链
        params.min_score = 10;

        auto best_chain = anchor::chainAnchors(anchors, params);

        REQUIRE(!best_chain.empty());
        CHECK(best_chain.size() == 1);
        CHECK(best_chain[0].span == 20);
    }

    TEST_CASE("chainAnchors - 两个可链接锚点形成单链")
    {
        anchor::Anchors anchors;

        // 两个锚点，位置递增，gap 适中
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 20));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 30;

        auto best_chain = anchor::chainAnchors(anchors, params);

        REQUIRE(!best_chain.empty());
        CHECK(best_chain.size() == 2);
        // 验证锚点按位置排序
        CHECK(best_chain[0].pos_ref < best_chain[1].pos_ref);
    }

    TEST_CASE("chainAnchors - 不同参考序列的锚点只返回最佳链")
    {
        anchor::Anchors anchors;

        // 链 A：rid_ref=0，得分较高
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x555555, 200, 150, 0, 0, false, 20));

        // 链 B：rid_ref=1，得分较低
        anchors.push_back(makeAnchor(0x333333, 200, 150, 1, 0, false, 15));
        anchors.push_back(makeAnchor(0x444444, 250, 200, 1, 0, false, 15));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 20;

        auto best_chain = anchor::chainAnchors(anchors, params);

        // 应该返回得分更高的链 A（rid_ref=0）
        REQUIRE(!best_chain.empty());
        CHECK(best_chain[0].rid_ref == 0);
        CHECK(best_chain.size() >= 2);
    }

    TEST_CASE("chainAnchors - 正向和反向锚点只返回最佳链")
    {
        anchor::Anchors anchors;

        // 正向链（得分更高）
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x666666, 200, 150, 0, 0, false, 20));

        // 反向链（得分较低）
        anchors.push_back(makeAnchor(0x333333, 300, 250, 0, 0, true, 15));
        anchors.push_back(makeAnchor(0x444444, 350, 300, 0, 0, true, 15));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 20;

        auto best_chain = anchor::chainAnchors(anchors, params);

        // 应该返回得分更高的正向链
        REQUIRE(!best_chain.empty());
        CHECK(best_chain[0].is_rev == false);
        CHECK(best_chain.size() >= 2);
    }

    TEST_CASE("chainAnchors - 距离过远的锚点不链接")
    {
        anchor::Anchors anchors;

        // 两个锚点，ref 距离超过 max_dist_x (默认 5000)
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 6000, 5500, 0, 0, false, 20));

        anchor::ChainParams params;
        params.min_cnt = 1;  // 允许单锚点链
        params.min_score = 10;
        params.max_dist_x = 5000;

        auto best_chain = anchor::chainAnchors(anchors, params);

        // 应该只返回一个锚点（距离过远无法链接）
        REQUIRE(!best_chain.empty());
        CHECK(best_chain.size() == 1);
    }

    TEST_CASE("chainAnchors - 对角线偏移超过带宽时不链接")
    {
        anchor::Anchors anchors;

        // 两个锚点，对角线偏移 = |dr - dq| = |150 - 100| = 50
        // 如果 bw < 50，应该不链接
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 250, 150, 0, 0, false, 20));
        // dr = 250 - 100 = 150, dq = 150 - 50 = 100, dd = |150 - 100| = 50

        anchor::ChainParams params;
        params.min_cnt = 1;
        params.min_score = 10;
        params.bw = 30;  // 带宽 < 50，应该不链接

        auto best_chain = anchor::chainAnchors(anchors, params);

        // 应该只返回一个锚点（对角线偏移超过带宽）
        REQUIRE(!best_chain.empty());
        CHECK(best_chain.size() == 1);
    }

    TEST_CASE("chainAnchors - 多条链按得分降序排列")
    {
        anchor::Anchors anchors;

        // 链 A：3 个锚点（得分更高）
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x333333, 200, 150, 0, 0, false, 20));

        // 链 B：2 个锚点（得分较低，且距离远离��� A）
        anchors.push_back(makeAnchor(0x444444, 10000, 5000, 0, 0, false, 15));
        anchors.push_back(makeAnchor(0x555555, 10100, 5100, 0, 0, false, 15));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 20;

        auto best_chain = anchor::chainAnchors(anchors, params);

        REQUIRE(!best_chain.empty());

        // 应该返回得分更高的链 A（3 个锚点）
        CHECK(best_chain.size() == 3);
    }

    // ------------------------------------------------------------------
    // 参数过滤测试
    // ------------------------------------------------------------------

    TEST_CASE("chainAnchors - min_cnt 过滤短链")
    {
        anchor::Anchors anchors;

        // 只有 2 个锚点
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 15));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 15));

        // 要求至少 3 个锚点
        anchor::ChainParams params;
        params.min_cnt = 3;
        params.min_score = 10;

        auto best_chain = anchor::chainAnchors(anchors, params);

        // 应该没有链（不满足 min_cnt）
        CHECK(best_chain.empty());
    }

    TEST_CASE("chainAnchors - min_score 过滤低分链")
    {
        anchor::Anchors anchors;

        // 2 个锚点，span 很小，得分低
        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 5));
        anchors.push_back(makeAnchor(0x222222, 150, 100, 0, 0, false, 5));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 50;  // 要求得分 >= 50

        auto best_chain = anchor::chainAnchors(anchors, params);

        // 得分太低，应该被过滤
        CHECK(best_chain.empty());
    }

    TEST_CASE("chainAnchors - 返回锚点按位置排序")
    {
        anchor::Anchors anchors;

        anchors.push_back(makeAnchor(0x111111, 100, 50, 0, 0, false, 20));
        anchors.push_back(makeAnchor(0x222222, 200, 150, 0, 0, false, 25));

        anchor::ChainParams params;
        params.min_cnt = 2;
        params.min_score = 30;

        auto best_chain = anchor::chainAnchors(anchors, params);

        REQUIRE(!best_chain.empty());
        REQUIRE(best_chain.size() == 2);

        // 验证锚点按位置顺序排列
        CHECK(best_chain[0].pos_ref == 100);
        CHECK(best_chain[1].pos_ref == 200);
        CHECK(best_chain[0].pos_qry == 50);
        CHECK(best_chain[1].pos_qry == 150);
    }
}



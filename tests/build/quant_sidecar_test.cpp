/*
 * Round-trip tests for the `.qtx` sidecar format (quant_sidecar::Writer /
 * quant_sidecar::Reader).
 *
 * These tests treat the Writer + Reader as a black box: write a known
 * set of records, reopen the file, and verify header fields + record
 * contents via both `lookup` (binary search) and `for_each` (linear
 * scan). No build pipeline is involved — the goal is to catch format
 * regressions in isolation before they propagate into end-to-end
 * workflows once sidecars are wired into the build path.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <vector>

#include "quant_sidecar.hpp"

namespace fs = std::filesystem;

class QuantSidecarTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Per-test unique directory so a failed test leaves an
        // inspectable artifact behind without colliding with siblings.
        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_qtx_test_";
        dir_name += info ? info->name() : "unknown";
        tmp_dir = fs::temp_directory_path() / dir_name;

        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
        fs::create_directories(tmp_dir, ec);
        ASSERT_FALSE(ec) << "Failed to create tmp_dir: " << ec.message();
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
    }

    fs::path tmp_dir;
};

// ── Round trip: header + records ────────────────────────────────────

TEST_F(QuantSidecarTest, EmptyFileRoundTrip) {
    fs::path path = tmp_dir / "empty.qtx";

    {
        quant_sidecar::Writer w(path, /*expr_type=*/4, /*grove_id=*/0xabcdef12);
        EXPECT_TRUE(w.empty());
        w.finalize();
    }
    ASSERT_TRUE(fs::exists(path));

    quant_sidecar::Reader r(path);
    EXPECT_EQ(r.size(), 0u);
    EXPECT_TRUE(r.empty());
    EXPECT_EQ(r.header().count, 0u);
    EXPECT_EQ(r.header().expr_type, 4u);
    EXPECT_EQ(r.header().grove_id, 0xabcdef12u);
    EXPECT_EQ(r.header().version, quant_sidecar::QTX_VERSION);
    EXPECT_EQ(std::memcmp(r.header().magic,
                          quant_sidecar::MAGIC,
                          sizeof(quant_sidecar::MAGIC)), 0);
    EXPECT_TRUE(r.grove_id_matches(0xabcdef12));
    EXPECT_FALSE(r.grove_id_matches(0));
}

TEST_F(QuantSidecarTest, SmallFileRoundTrip) {
    fs::path path = tmp_dir / "small.qtx";

    // Append in unsorted order to verify finalize() sorts.
    {
        quant_sidecar::Writer w(path, /*expr_type=*/1, /*grove_id=*/0xdeadbeef);
        w.append(7,   0.5f);
        w.append(1,   42.0f);
        w.append(42, 123.75f);
        w.append(3,   0.0f);   // zero is a valid value
        EXPECT_EQ(w.size(), 4u);
        EXPECT_FALSE(w.empty());
        w.finalize();
    }

    quant_sidecar::Reader r(path);
    ASSERT_EQ(r.size(), 4u);
    EXPECT_EQ(r.header().count, 4u);
    EXPECT_EQ(r.header().expr_type, 1u);
    EXPECT_EQ(r.header().grove_id, 0xdeadbeefu);

    // Records come back sorted by segment_index regardless of
    // insertion order.
    auto& recs = r.records();
    ASSERT_EQ(recs.size(), 4u);
    EXPECT_EQ(recs[0].segment_index, 1u);
    EXPECT_FLOAT_EQ(recs[0].value, 42.0f);
    EXPECT_EQ(recs[1].segment_index, 3u);
    EXPECT_FLOAT_EQ(recs[1].value, 0.0f);
    EXPECT_EQ(recs[2].segment_index, 7u);
    EXPECT_FLOAT_EQ(recs[2].value, 0.5f);
    EXPECT_EQ(recs[3].segment_index, 42u);
    EXPECT_FLOAT_EQ(recs[3].value, 123.75f);
}

// ── lookup() hits, misses, bounds ───────────────────────────────────

TEST_F(QuantSidecarTest, LookupHitAndMiss) {
    fs::path path = tmp_dir / "lookup.qtx";

    {
        quant_sidecar::Writer w(path, /*expr_type=*/2, /*grove_id=*/0);
        w.append(100, 1.0f);
        w.append(200, 2.0f);
        w.append(300, 3.0f);
        w.append(400, 4.0f);
        w.append(500, 5.0f);
        w.finalize();
    }

    quant_sidecar::Reader r(path);

    // Hits
    auto hit_200 = r.lookup(200);
    ASSERT_TRUE(hit_200.has_value());
    EXPECT_FLOAT_EQ(*hit_200, 2.0f);

    auto hit_500 = r.lookup(500);
    ASSERT_TRUE(hit_500.has_value());
    EXPECT_FLOAT_EQ(*hit_500, 5.0f);

    // Boundary hit (first record)
    auto hit_first = r.lookup(100);
    ASSERT_TRUE(hit_first.has_value());
    EXPECT_FLOAT_EQ(*hit_first, 1.0f);

    // Misses
    EXPECT_FALSE(r.lookup(0).has_value());         // below range
    EXPECT_FALSE(r.lookup(150).has_value());       // between records
    EXPECT_FALSE(r.lookup(600).has_value());       // above range
    EXPECT_FALSE(r.lookup(0xffffffffull).has_value());  // way above
}

// ── for_each() ordering ─────────────────────────────────────────────

TEST_F(QuantSidecarTest, ForEachIteratesSorted) {
    fs::path path = tmp_dir / "for_each.qtx";

    {
        quant_sidecar::Writer w(path, /*expr_type=*/3, /*grove_id=*/99);
        // Shuffled insertion order
        w.append(50,   5.0f);
        w.append(10,   1.0f);
        w.append(40,   4.0f);
        w.append(20,   2.0f);
        w.append(30,   3.0f);
        w.finalize();
    }

    quant_sidecar::Reader r(path);

    std::vector<std::pair<uint64_t, float>> seen;
    r.for_each([&](uint64_t sid, float val) {
        seen.emplace_back(sid, val);
    });

    ASSERT_EQ(seen.size(), 5u);
    // Strictly increasing by segment_index
    for (size_t i = 1; i < seen.size(); ++i) {
        EXPECT_LT(seen[i - 1].first, seen[i].first)
            << "for_each yielded records out of order at index " << i;
    }
    EXPECT_EQ(seen[0].first, 10u);
    EXPECT_FLOAT_EQ(seen[0].second, 1.0f);
    EXPECT_EQ(seen[4].first, 50u);
    EXPECT_FLOAT_EQ(seen[4].second, 5.0f);
}

// ── Atomicity: no lingering .tmp file after successful finalize ────

TEST_F(QuantSidecarTest, NoLeftoverTmpAfterFinalize) {
    fs::path path = tmp_dir / "atomic.qtx";

    {
        quant_sidecar::Writer w(path, 0, 0);
        w.append(1, 1.0f);
        w.finalize();
    }

    EXPECT_TRUE(fs::exists(path));
    fs::path tmp_path = path;
    tmp_path += ".tmp";
    EXPECT_FALSE(fs::exists(tmp_path))
        << "Writer left behind .tmp after atomic rename";
}

// ── Destructor flushes when finalize not called explicitly ─────────

TEST_F(QuantSidecarTest, DestructorFinalizesIfForgotten) {
    fs::path path = tmp_dir / "dtor.qtx";

    {
        quant_sidecar::Writer w(path, 0, 0);
        w.append(1, 10.0f);
        w.append(2, 20.0f);
        // No explicit finalize() — destructor should handle it.
    }

    ASSERT_TRUE(fs::exists(path));
    quant_sidecar::Reader r(path);
    EXPECT_EQ(r.size(), 2u);
    auto v1 = r.lookup(1);
    auto v2 = r.lookup(2);
    ASSERT_TRUE(v1.has_value());
    ASSERT_TRUE(v2.has_value());
    EXPECT_FLOAT_EQ(*v1, 10.0f);
    EXPECT_FLOAT_EQ(*v2, 20.0f);
}

// ── Idempotent finalize() ───────────────────────────────────────────

TEST_F(QuantSidecarTest, RepeatFinalizeIsNoOp) {
    fs::path path = tmp_dir / "idempotent.qtx";

    quant_sidecar::Writer w(path, 0, 0);
    w.append(5, 5.5f);
    w.finalize();
    w.finalize();  // must not throw or corrupt
    w.finalize();

    quant_sidecar::Reader r(path);
    EXPECT_EQ(r.size(), 1u);
}

// ── Reader error handling ───────────────────────────────────────────

// Note on the brace-init pattern below (`Reader{path}` instead of
// `Reader(path)`): when `path` is a local identifier in scope, C++'s
// most-vexing-parse rule interprets `quant_sidecar::Reader(path)` as
// declaring a default-constructed Reader named `path`, which is a
// compile error because Reader has no default constructor. Brace
// initialization is unambiguous — it creates a temporary.
TEST_F(QuantSidecarTest, ReaderRejectsMissingFile) {
    EXPECT_THROW(
        quant_sidecar::Reader{tmp_dir / "does_not_exist.qtx"},
        std::runtime_error);
}

TEST_F(QuantSidecarTest, ReaderRejectsBadMagic) {
    fs::path path = tmp_dir / "bad_magic.qtx";

    // Write a 32-byte header with wrong magic.
    {
        std::ofstream out(path, std::ios::binary);
        quant_sidecar::Header hdr{};
        const char junk[4] = {'N', 'O', 'P', 'E'};
        std::memcpy(hdr.magic, junk, 4);
        hdr.version = quant_sidecar::QTX_VERSION;
        hdr.count = 0;
        out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    }

    EXPECT_THROW(quant_sidecar::Reader{path}, std::runtime_error);
}

TEST_F(QuantSidecarTest, ReaderRejectsFutureVersion) {
    fs::path path = tmp_dir / "bad_version.qtx";

    {
        std::ofstream out(path, std::ios::binary);
        quant_sidecar::Header hdr{};
        std::memcpy(hdr.magic, quant_sidecar::MAGIC, sizeof(quant_sidecar::MAGIC));
        hdr.version = quant_sidecar::QTX_VERSION + 99;  // far future
        hdr.count = 0;
        out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    }

    EXPECT_THROW(quant_sidecar::Reader{path}, std::runtime_error);
}

TEST_F(QuantSidecarTest, ReaderRejectsTruncatedRecords) {
    fs::path path = tmp_dir / "truncated.qtx";

    // Header claims 5 records, but only write 2 records' worth of bytes.
    {
        std::ofstream out(path, std::ios::binary);
        quant_sidecar::Header hdr{};
        std::memcpy(hdr.magic, quant_sidecar::MAGIC, sizeof(quant_sidecar::MAGIC));
        hdr.version = quant_sidecar::QTX_VERSION;
        hdr.count = 5;
        out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

        quant_sidecar::Record recs[2] = {
            {1, 1.0f},
            {2, 2.0f},
        };
        out.write(reinterpret_cast<const char*>(recs), sizeof(recs));
    }

    EXPECT_THROW(quant_sidecar::Reader{path}, std::runtime_error);
}

// ── Stable sort: insertion order preserved on duplicate keys ────────
//
// Not an invariant the build path relies on (duplicates shouldn't
// reach the Writer), but since we explicitly use stable_sort, we want
// a regression test to catch any accidental switch back to std::sort.

TEST_F(QuantSidecarTest, StableSortPreservesInsertionOrderOnDupKeys) {
    fs::path path = tmp_dir / "stable.qtx";

    {
        quant_sidecar::Writer w(path, 0, 0);
        w.append(5, 1.0f);   // first-inserted for key 5
        w.append(3, 10.0f);
        w.append(5, 2.0f);   // second-inserted for key 5
        w.append(1, 100.0f);
        w.append(5, 3.0f);   // third-inserted for key 5
        w.finalize();
    }

    quant_sidecar::Reader r(path);
    ASSERT_EQ(r.size(), 5u);
    auto& recs = r.records();

    // Sorted by key, ties resolved by insertion order.
    EXPECT_EQ(recs[0].segment_index, 1u);
    EXPECT_EQ(recs[1].segment_index, 3u);
    EXPECT_EQ(recs[2].segment_index, 5u);
    EXPECT_FLOAT_EQ(recs[2].value, 1.0f);
    EXPECT_EQ(recs[3].segment_index, 5u);
    EXPECT_FLOAT_EQ(recs[3].value, 2.0f);
    EXPECT_EQ(recs[4].segment_index, 5u);
    EXPECT_FLOAT_EQ(recs[4].value, 3.0f);
}

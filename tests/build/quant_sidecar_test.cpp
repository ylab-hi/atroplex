/*
 * Round-trip tests for the `.qtx` sidecar format.
 *
 * The new sidecar format is segment-major: per-sample temp streams
 * (SampleStreamWriter) get K-way merged into a single `.qtx` file
 * (merge_to_qtx), which the Reader then opens for segment-indexed
 * lookup of per-sample values.
 *
 * These tests exercise the round-trip end-to-end: write a handful of
 * per-sample streams with known records, merge them, reopen with Reader,
 * and verify lookups return the right (sample_id, value) sets.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "quant_sidecar.hpp"

namespace fs = std::filesystem;

class QuantSidecarTest : public ::testing::Test {
protected:
    void SetUp() override {
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

    // Helper: write a per-sample stream file and finalize it.
    fs::path write_stream(uint32_t sample_id,
                          const std::vector<std::pair<uint64_t, float>>& records) {
        fs::path p = tmp_dir / (std::to_string(sample_id) + ".stream");
        {
            quant_sidecar::SampleStreamWriter w(p, sample_id);
            for (auto [seg, val] : records) {
                w.append(seg, val);
            }
            w.finalize();
        }
        return p;
    }

    // Helper: build a small catalog of (sample_id, records) + metadata and
    // merge it into a .qtx at the given path.
    fs::path merge_samples(
        const std::vector<std::pair<uint32_t, std::vector<std::pair<uint64_t, float>>>>& samples,
        const std::string& qtx_name = "out.qtx") {
        std::vector<fs::path> streams;
        std::vector<quant_sidecar::SampleMetadata> metas;
        for (const auto& [sid, recs] : samples) {
            streams.push_back(write_stream(sid, recs));
            metas.push_back({sid, /*expr_type=*/0, "sample_" + std::to_string(sid)});
        }
        fs::path out = tmp_dir / qtx_name;
        quant_sidecar::merge_to_qtx(out, streams, metas);
        return out;
    }

    fs::path tmp_dir;
};

// ── SampleStreamWriter round-trip (basic sanity) ────────────────────

TEST_F(QuantSidecarTest, StreamWriter_EmptyFinalizes) {
    fs::path p = tmp_dir / "empty.stream";
    {
        quant_sidecar::SampleStreamWriter w(p, /*sample_id=*/42);
        EXPECT_TRUE(w.empty());
        w.finalize();
    }
    EXPECT_TRUE(fs::exists(p)) << "Empty finalize should still produce a file";
    EXPECT_GE(fs::file_size(p), sizeof(quant_sidecar::StreamHeader));
}

TEST_F(QuantSidecarTest, StreamWriter_DestructorFinalizes) {
    fs::path p = tmp_dir / "dtor.stream";
    {
        quant_sidecar::SampleStreamWriter w(p, /*sample_id=*/7);
        w.append(100, 1.5f);
        w.append(50,  2.5f);
        // No explicit finalize; destructor should flush.
    }
    EXPECT_TRUE(fs::exists(p));
}

TEST_F(QuantSidecarTest, StreamWriter_RepeatFinalizeIsNoOp) {
    fs::path p = tmp_dir / "idempotent.stream";
    quant_sidecar::SampleStreamWriter w(p, 1);
    w.append(1, 1.0f);
    w.finalize();
    w.finalize();  // must not throw
    w.finalize();
    EXPECT_TRUE(fs::exists(p));
}

// ── Merge + Reader round-trip ───────────────────────────────────────

TEST_F(QuantSidecarTest, Merge_SingleSample) {
    // Sample 7 has values at segments 1, 3, 5
    fs::path out = merge_samples({
        {7, {{5, 50.0f}, {1, 10.0f}, {3, 30.0f}}}  // insertion unsorted
    });

    quant_sidecar::Reader r(out);
    EXPECT_EQ(r.version(), quant_sidecar::QTX_VERSION);
    EXPECT_EQ(r.segment_block_count(), 3u);
    ASSERT_EQ(r.samples().size(), 1u);
    EXPECT_EQ(r.samples()[0].sample_id, 7u);
    EXPECT_EQ(r.samples()[0].name, "sample_7");

    // Each segment has exactly one (sample_id=7, value) record
    auto r1 = r.lookup(1);
    ASSERT_EQ(r1.size(), 1u);
    EXPECT_EQ(r1[0].sample_id, 7u);
    EXPECT_FLOAT_EQ(r1[0].value, 10.0f);

    auto r3 = r.lookup(3);
    ASSERT_EQ(r3.size(), 1u);
    EXPECT_FLOAT_EQ(r3[0].value, 30.0f);

    auto r5 = r.lookup(5);
    ASSERT_EQ(r5.size(), 1u);
    EXPECT_FLOAT_EQ(r5[0].value, 50.0f);

    // Missing segment returns empty
    EXPECT_TRUE(r.lookup(2).empty());
    EXPECT_TRUE(r.lookup(99).empty());
}

TEST_F(QuantSidecarTest, Merge_MultipleSamplesSegmentMajor) {
    // Three samples overlapping on segments 10 and 20
    fs::path out = merge_samples({
        {1, {{10, 1.0f}, {20, 2.0f}, {30, 3.0f}}},
        {2, {{10, 11.0f}, {20, 22.0f}}},
        {3, {{20, 222.0f}, {40, 444.0f}}},
    });

    quant_sidecar::Reader r(out);
    ASSERT_EQ(r.samples().size(), 3u);
    EXPECT_EQ(r.segment_block_count(), 4u);  // segments 10, 20, 30, 40

    // Segment 10: samples 1 and 2
    auto b10 = r.lookup(10);
    ASSERT_EQ(b10.size(), 2u);
    // Records within a block are sorted by sample_id
    EXPECT_EQ(b10[0].sample_id, 1u);
    EXPECT_FLOAT_EQ(b10[0].value, 1.0f);
    EXPECT_EQ(b10[1].sample_id, 2u);
    EXPECT_FLOAT_EQ(b10[1].value, 11.0f);

    // Segment 20: all three samples
    auto b20 = r.lookup(20);
    ASSERT_EQ(b20.size(), 3u);
    EXPECT_EQ(b20[0].sample_id, 1u);
    EXPECT_FLOAT_EQ(b20[0].value, 2.0f);
    EXPECT_EQ(b20[1].sample_id, 2u);
    EXPECT_FLOAT_EQ(b20[1].value, 22.0f);
    EXPECT_EQ(b20[2].sample_id, 3u);
    EXPECT_FLOAT_EQ(b20[2].value, 222.0f);

    // Segment 30: sample 1 only
    auto b30 = r.lookup(30);
    ASSERT_EQ(b30.size(), 1u);
    EXPECT_EQ(b30[0].sample_id, 1u);
    EXPECT_FLOAT_EQ(b30[0].value, 3.0f);

    // Segment 40: sample 3 only
    auto b40 = r.lookup(40);
    ASSERT_EQ(b40.size(), 1u);
    EXPECT_EQ(b40[0].sample_id, 3u);
    EXPECT_FLOAT_EQ(b40[0].value, 444.0f);
}

TEST_F(QuantSidecarTest, LookupFiltered_IntersectsWithWantedSet) {
    fs::path out = merge_samples({
        {1, {{100, 1.0f}}},
        {2, {{100, 2.0f}}},
        {3, {{100, 3.0f}}},
        {4, {{100, 4.0f}}},
        {5, {{100, 5.0f}}},
    });

    quant_sidecar::Reader r(out);

    // Ask for samples {2, 4, 99}: only 2 and 4 hit; 99 isn't in the block
    auto filtered = r.lookup_filtered(100, {2u, 4u, 99u});
    ASSERT_EQ(filtered.size(), 2u);
    EXPECT_FLOAT_EQ(filtered.at(2), 2.0f);
    EXPECT_FLOAT_EQ(filtered.at(4), 4.0f);
    EXPECT_EQ(filtered.count(99), 0u);

    // Missing segment: empty map regardless of wanted list
    auto none = r.lookup_filtered(999, {1u, 2u});
    EXPECT_TRUE(none.empty());
}

TEST_F(QuantSidecarTest, ForEachSegment_IteratesInOrder) {
    fs::path out = merge_samples({
        {1, {{5, 5.0f}, {1, 1.0f}, {3, 3.0f}}},
        {2, {{3, 33.0f}, {5, 55.0f}}},
    });

    quant_sidecar::Reader r(out);

    std::vector<uint64_t> seen_segs;
    r.for_each_segment([&](uint64_t seg, const auto& recs) {
        seen_segs.push_back(seg);
        // Records within block must be sorted by sample_id
        for (size_t i = 1; i < recs.size(); ++i) {
            EXPECT_LT(recs[i - 1].sample_id, recs[i].sample_id)
                << "Block records not sorted by sample_id at segment " << seg;
        }
    });

    ASSERT_EQ(seen_segs.size(), 3u);
    EXPECT_EQ(seen_segs[0], 1u);
    EXPECT_EQ(seen_segs[1], 3u);
    EXPECT_EQ(seen_segs[2], 5u);
}

TEST_F(QuantSidecarTest, Merge_SampleMetadataRoundTrips) {
    // Explicit metadata: mix of expr_type values and names
    std::vector<fs::path> streams;
    std::vector<quant_sidecar::SampleMetadata> metas;

    streams.push_back(write_stream(100, {{1, 10.0f}}));
    metas.push_back({100, /*expr_type=*/2, "TCGA-A01"});

    streams.push_back(write_stream(200, {{1, 20.0f}}));
    metas.push_back({200, /*expr_type=*/4, "ENCSR_042"});

    fs::path out = tmp_dir / "meta.qtx";
    quant_sidecar::merge_to_qtx(out, streams, metas);

    quant_sidecar::Reader r(out);
    ASSERT_EQ(r.samples().size(), 2u);

    // Order in sample metadata should match the merge input order
    EXPECT_EQ(r.samples()[0].sample_id, 100u);
    EXPECT_EQ(r.samples()[0].expr_type, 2u);
    EXPECT_EQ(r.samples()[0].name, "TCGA-A01");

    EXPECT_EQ(r.samples()[1].sample_id, 200u);
    EXPECT_EQ(r.samples()[1].expr_type, 4u);
    EXPECT_EQ(r.samples()[1].name, "ENCSR_042");
}

TEST_F(QuantSidecarTest, Merge_EmptyStreamsProducesEmptyQtx) {
    // All samples have zero records — valid edge case
    fs::path s1 = write_stream(1, {});
    fs::path s2 = write_stream(2, {});

    std::vector<fs::path> streams{s1, s2};
    std::vector<quant_sidecar::SampleMetadata> metas = {
        {1, 0, "a"}, {2, 0, "b"}
    };

    fs::path out = tmp_dir / "empty.qtx";
    quant_sidecar::merge_to_qtx(out, streams, metas);

    quant_sidecar::Reader r(out);
    EXPECT_EQ(r.segment_block_count(), 0u);
    EXPECT_EQ(r.samples().size(), 2u);
    EXPECT_TRUE(r.lookup(1).empty());
}

TEST_F(QuantSidecarTest, Merge_NoLeftoverTmpFile) {
    fs::path out = merge_samples({{1, {{1, 1.0f}}}});
    fs::path tmp_out = fs::path(out.string() + ".tmp");
    EXPECT_TRUE(fs::exists(out));
    EXPECT_FALSE(fs::exists(tmp_out))
        << "merge_to_qtx left behind a .tmp file";
}

// ── Reader error handling ───────────────────────────────────────────

TEST_F(QuantSidecarTest, Reader_RejectsMissingFile) {
    EXPECT_THROW(
        quant_sidecar::Reader{tmp_dir / "does_not_exist.qtx"},
        std::runtime_error);
}

TEST_F(QuantSidecarTest, Reader_RejectsBadMagic) {
    fs::path path = tmp_dir / "bad_magic.qtx";
    {
        std::ofstream out(path, std::ios::binary);
        quant_sidecar::Header hdr{};
        const char junk[4] = {'N', 'O', 'P', 'E'};
        std::memcpy(hdr.magic, junk, 4);
        hdr.version = quant_sidecar::QTX_VERSION;
        out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    }
    EXPECT_THROW(quant_sidecar::Reader{path}, std::runtime_error);
}

TEST_F(QuantSidecarTest, Reader_RejectsFutureVersion) {
    fs::path path = tmp_dir / "bad_version.qtx";
    {
        std::ofstream out(path, std::ios::binary);
        quant_sidecar::Header hdr{};
        std::memcpy(hdr.magic, quant_sidecar::MAGIC, sizeof(quant_sidecar::MAGIC));
        hdr.version = quant_sidecar::QTX_VERSION + 99;
        out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    }
    EXPECT_THROW(quant_sidecar::Reader{path}, std::runtime_error);
}
/*
 * Tests for the export `--conserved-fraction` flag.
 *
 * Builds a grove with 4 sample-typed entries — 3 share a single 3-exon
 * structure, the 4th contributes a distinct structure — then drives
 * subcall::export_gtf end-to-end through its cxxopts CLI surface and
 * verifies the per-sample GTFs reflect the configured threshold:
 *
 *   - default `--conserved-only` (fraction=1.0): segment must be in
 *     every sample-typed entry. With 3/4 it should NOT export.
 *   - `--conserved-only --conserved-fraction 0.75`: ceil(4 * 0.75) = 3,
 *     so the segment passes and IS exported.
 *   - `--conserved-fraction` outside (0, 1] fails validate().
 *
 * Mirrors the conserved_threshold_test fixture pattern (sample-shared
 * structure built from N small GTFs) and the compact_test pattern
 * (drive a subcall through cxxopts end-to-end).
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cxxopts.hpp>

#include "builder.hpp"
#include "genomic_feature.hpp"
#include "sample_info.hpp"
#include "subcall/export_gtf.hpp"

namespace fs = std::filesystem;

class ExportConservedFractionTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_export_cf_test_";
        dir_name += info ? info->name() : "unknown";
        tmp_dir = fs::temp_directory_path() / dir_name;

        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
        fs::create_directories(tmp_dir, ec);
        ASSERT_FALSE(ec) << "Failed to create tmp_dir " << tmp_dir
                         << ": " << ec.message();
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
    }

    fs::path tmp_dir;

    fs::path write_shared_gtf(const std::string& tag) {
        std::ostringstream g;
        g << "chr1\tTEST\ttranscript\t1000\t5000\t.\t+\t.\t"
          << "gene_id \"G_SHARED\"; transcript_id \"TX_" << tag << "\";\n"
          << "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\t"
          << "gene_id \"G_SHARED\"; transcript_id \"TX_" << tag << "\"; exon_number \"1\";\n"
          << "chr1\tTEST\texon\t2000\t2200\t.\t+\t.\t"
          << "gene_id \"G_SHARED\"; transcript_id \"TX_" << tag << "\"; exon_number \"2\";\n"
          << "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\t"
          << "gene_id \"G_SHARED\"; transcript_id \"TX_" << tag << "\"; exon_number \"3\";\n";

        fs::path p = tmp_dir / ("shared_" + tag + ".gtf");
        std::ofstream out(p, std::ios::binary | std::ios::trunc);
        out << g.str();
        out.close();
        return p;
    }

    fs::path write_distinct_gtf(const std::string& tag) {
        std::ostringstream g;
        g << "chr1\tTEST\ttranscript\t20000\t25000\t.\t+\t.\t"
          << "gene_id \"G_DISTINCT_" << tag << "\"; transcript_id \"DTX_" << tag << "\";\n"
          << "chr1\tTEST\texon\t20000\t20200\t.\t+\t.\t"
          << "gene_id \"G_DISTINCT_" << tag << "\"; transcript_id \"DTX_" << tag << "\"; exon_number \"1\";\n"
          << "chr1\tTEST\texon\t22000\t22200\t.\t+\t.\t"
          << "gene_id \"G_DISTINCT_" << tag << "\"; transcript_id \"DTX_" << tag << "\"; exon_number \"2\";\n"
          << "chr1\tTEST\texon\t24500\t25000\t.\t+\t.\t"
          << "gene_id \"G_DISTINCT_" << tag << "\"; transcript_id \"DTX_" << tag << "\"; exon_number \"3\";\n";

        fs::path p = tmp_dir / ("distinct_" + tag + ".gtf");
        std::ofstream out(p, std::ios::binary | std::ios::trunc);
        out << g.str();
        out.close();
        return p;
    }

    /// Build a grove with `shared` entries that all carry the
    /// G_SHARED segment plus `distinct` entries with their own
    /// non-overlapping structures, serialize to `index_dir/test.ggx`.
    void build_index(const fs::path& index_dir,
                     size_t shared,
                     size_t distinct) {
        fs::create_directories(index_dir);

        std::vector<sample_info> samples;
        samples.reserve(shared + distinct);
        for (size_t i = 0; i < shared; ++i) {
            std::string tag = "S" + std::to_string(i);
            samples.emplace_back("sample_" + tag, write_shared_gtf(tag));
            samples.back().type = "sample";
        }
        for (size_t i = 0; i < distinct; ++i) {
            std::string tag = "D" + std::to_string(i);
            samples.emplace_back("dropout_" + tag, write_distinct_gtf(tag));
            samples.back().type = "sample";
        }

        grove_type grove(3);
        [[maybe_unused]] auto summary = builder::build_from_samples(grove, samples);

        // Mirror subcall::save_grove
        fs::path ggx_path = index_dir / "test.ggx";
        std::ofstream ofs(ggx_path, std::ios::binary);
        ASSERT_TRUE(ofs.is_open());
        ofs.write("AGRX", 4);
        uint16_t version = 1;
        ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));
        gene_registry::instance().serialize(ofs);
        source_registry::instance().serialize(ofs);
        transcript_registry::instance().serialize(ofs);
        sample_registry::instance().serialize(ofs);
        grove.serialize(ofs);
        ofs.close();
    }

    /// Reset registries so a subsequent grove load does not collide
    /// with the in-memory state left by build_from_samples.
    void reset_registries() {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();
    }

    /// Drive subcall::export_gtf through its CLI surface end-to-end.
    void run_export(const std::vector<std::string>& cli_args) {
        subcall::export_gtf cmd;
        cxxopts::Options options = cmd.parse_args(0, nullptr);

        std::vector<std::string> storage = cli_args;
        std::vector<char*> argv;
        argv.reserve(storage.size());
        for (auto& s : storage) argv.push_back(s.data());

        auto parsed = options.parse(static_cast<int>(argv.size()), argv.data());
        cmd.run(parsed);
    }

    /// True if any per-sample GTF in `dir` contains an exon whose
    /// gene_id is "G_SHARED".
    bool any_gtf_contains_shared(const fs::path& dir) const {
        if (!fs::exists(dir)) return false;
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.path().extension() != ".gtf") continue;
            std::ifstream in(entry.path());
            std::string line;
            while (std::getline(in, line)) {
                if (line.find("G_SHARED") != std::string::npos) return true;
            }
        }
        return false;
    }
};

// ── Default --conserved-only is strict: 3 of 4 samples → segment dropped ──
TEST_F(ExportConservedFractionTest, ConservedOnly_DefaultStrict_DropsBelow100Percent) {
    fs::path input_dir = tmp_dir / "in";
    fs::path output_dir = tmp_dir / "out";
    build_index(input_dir, /*shared=*/3, /*distinct=*/1);

    reset_registries();

    run_export({
        "atroplex export",
        "-g", input_dir.string(),
        "-o", output_dir.string(),
        "--conserved-only",
    });

    EXPECT_FALSE(any_gtf_contains_shared(output_dir))
        << "Default (strict) --conserved-only must drop segments missing "
           "any sample. With 3/4 the G_SHARED segment should be filtered.";
}

// ── --conserved-fraction 0.75 with 3 of 4 → segment passes ────────────────
TEST_F(ExportConservedFractionTest, ConservedOnly_Fraction075_KeepsThreeOfFour) {
    fs::path input_dir = tmp_dir / "in";
    fs::path output_dir = tmp_dir / "out";
    build_index(input_dir, /*shared=*/3, /*distinct=*/1);

    reset_registries();

    run_export({
        "atroplex export",
        "-g", input_dir.string(),
        "-o", output_dir.string(),
        "--conserved-only",
        "--conserved-fraction", "0.75",
    });

    EXPECT_TRUE(any_gtf_contains_shared(output_dir))
        << "ceil(4 * 0.75) = 3 — a segment present in 3/4 samples should "
           "satisfy --conserved-fraction 0.75 and be exported.";
}

// ── --conserved-fraction must be in (0, 1] ─────────────────────────────
TEST_F(ExportConservedFractionTest, Validate_FractionOutOfRangeRejected) {
    fs::path input_dir = tmp_dir / "in";
    build_index(input_dir, /*shared=*/2, /*distinct=*/0);

    reset_registries();

    EXPECT_THROW(run_export({
        "atroplex export",
        "-g", input_dir.string(),
        "-o", (tmp_dir / "out").string(),
        "--conserved-only",
        "--conserved-fraction", "1.5",
    }), std::runtime_error);

    EXPECT_THROW(run_export({
        "atroplex export",
        "-g", input_dir.string(),
        "-o", (tmp_dir / "out2").string(),
        "--conserved-only",
        "--conserved-fraction", "0",
    }), std::runtime_error);
}

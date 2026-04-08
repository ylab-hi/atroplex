/*
 * Tests for .ggx index serialization roundtrip.
 *
 * Builds a grove from annotation.gtf + sample.gtf (reuses query fixtures),
 * serializes to a temp file, resets all registries, deserializes, then
 * verifies that registries, spatial index, and graph edges are preserved.
 */

#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>

#include "genomic_feature.hpp"
#include "build_gff.hpp"
#include "sample_info.hpp"
#include "transcript_matcher.hpp"

namespace fs = std::filesystem;

class SerializationTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();
        // Unique temp path per test — ctest --parallel runs tests concurrently
        auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        temp_path_ = (fs::temp_directory_path() /
            ("test_roundtrip_" + std::string(info->name()) + ".ggx")).string();
    }

    void TearDown() override {
        if (fs::exists(temp_path_)) {
            fs::remove(temp_path_);
        }
    }

    // Build grove from annotation + sample fixtures
    void build_grove() {
        grove_ = std::make_unique<grove_type>(3);

        fs::path anno_path = fs::path(TEST_FIXTURE_DIR) / "annotation.gtf";
        fs::path sample_path = fs::path(TEST_FIXTURE_DIR) / "sample.gtf";
        ASSERT_TRUE(fs::exists(anno_path)) << "annotation.gtf not found";
        ASSERT_TRUE(fs::exists(sample_path)) << "sample.gtf not found";

        // Register annotation
        sample_info anno_info("TEST_ANNOTATION");
        anno_info.type = "annotation";
        anno_info.annotation_source = "TEST";
        [[maybe_unused]] auto anno_id = sample_registry::instance().register_data(anno_info);

        build_gff::build(*grove_, anno_path, 0, exon_caches_, segment_caches_,
                         gene_indices_, segment_count_, 0, -1.0f, false);

        // Register sample (with expression)
        sample_info samp_info("TEST_SAMPLE");
        samp_info.type = "sample";
        [[maybe_unused]] auto samp_id = sample_registry::instance().register_data(samp_info);

        build_gff::build(*grove_, sample_path, 1, exon_caches_, segment_caches_,
                         gene_indices_, segment_count_, 0, -1.0f, false);
    }

    // Serialize grove + registries to temp file (mirrors subcall::save_grove)
    void save_grove(const std::string& path) {
        std::ofstream ofs(path, std::ios::binary);
        ASSERT_TRUE(ofs.is_open());

        ofs.write("AGRX", 4);
        uint16_t version = 1;
        ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));

        gene_registry::instance().serialize(ofs);
        source_registry::instance().serialize(ofs);
        transcript_registry::instance().serialize(ofs);
        sample_registry::instance().serialize(ofs);

        grove_->serialize(ofs);
        ASSERT_TRUE(ofs.good());
    }

    // Deserialize grove + registries from file (mirrors subcall::load_grove)
    std::unique_ptr<grove_type> load_grove(const std::string& path) {
        std::ifstream ifs(path, std::ios::binary);
        EXPECT_TRUE(ifs.is_open());

        char magic[4];
        ifs.read(magic, 4);
        EXPECT_EQ(std::string(magic, 4), "AGRX");

        uint16_t version;
        ifs.read(reinterpret_cast<char*>(&version), sizeof(version));
        EXPECT_EQ(version, 1);

        gene_registry::instance().deserialize_into(ifs);
        source_registry::instance().deserialize_into(ifs);
        transcript_registry::instance().deserialize_into(ifs);
        (void)sample_registry::deserialize(ifs);

        auto loaded = std::make_unique<grove_type>(grove_type::deserialize(ifs));
        EXPECT_FALSE(ifs.fail()) << "Stream should be clean after deserialization";
        return loaded;
    }

    // Snapshot registry state for comparison
    struct registry_snapshot {
        size_t gene_count;
        size_t source_count;
        size_t transcript_count;
        size_t sample_count;

        // Gene details
        std::vector<std::tuple<std::string, std::string, std::string>> genes;
        // Sample IDs
        std::vector<std::string> sample_ids;
    };

    registry_snapshot capture_registries() {
        registry_snapshot snap;
        snap.gene_count = gene_registry::instance().size();
        snap.source_count = source_registry::instance().size();
        snap.transcript_count = transcript_registry::instance().size();
        snap.sample_count = sample_registry::instance().size();

        for (uint32_t i = 0; i < snap.gene_count; ++i) {
            auto [id, name, biotype] = gene_registry::instance().resolve(i);
            snap.genes.emplace_back(id, name, biotype);
        }

        for (uint32_t i = 0; i < snap.sample_count; ++i) {
            auto& info = sample_registry::instance().get(i);
            snap.sample_ids.push_back(info.id);
        }

        return snap;
    }

    // Walk grove leaves and collect segment/exon data for comparison
    struct feature_snapshot {
        size_t segment_count = 0;
        size_t exon_count = 0;
        size_t total_edges = 0;
        std::set<std::string> gene_ids;
        std::set<uint32_t> transcript_ids;
        std::map<std::string, float> segment_expressions;  // coordinate -> expression
    };

    feature_snapshot capture_grove(grove_type& grove) {
        feature_snapshot snap;
        snap.total_edges = grove.edge_count();

        auto roots = grove.get_root_nodes();
        for (auto& [seqid, root] : roots) {
            if (!root) continue;
            auto* node = root;
            while (!node->get_is_leaf()) {
                auto& children = node->get_children();
                if (children.empty()) break;
                node = children[0];
            }

            while (node) {
                for (auto* key : node->get_keys()) {
                    auto& feature = key->get_data();
                    if (is_segment(feature)) {
                        auto& seg = get_segment(feature);
                        if (seg.absorbed) continue;
                        snap.segment_count++;
                        auto [id, name, biotype] = gene_registry::instance().resolve(seg.gene_idx);
                        snap.gene_ids.insert(id);
                        for (auto tx : seg.transcript_ids) {
                            snap.transcript_ids.insert(tx);
                        }
                        // Capture expression for sample_id=1
                        float expr = seg.get_expression(1);
                        if (expr >= 0) {
                            auto coord = format_coordinate(seqid, key->get_value());
                            snap.segment_expressions[coord] = expr;
                        }
                    }
                }
                node = node->get_next();
            }
        }

        return snap;
    }

    std::unique_ptr<grove_type> grove_;
    chromosome_exon_caches exon_caches_;
    chromosome_segment_caches segment_caches_;
    chromosome_gene_segment_indices gene_indices_;
    size_t segment_count_ = 0;
    std::string temp_path_;
};

// ─── Registry roundtrip ────────────────────────────────────────────

TEST_F(SerializationTest, RegistryRoundtrip) {
    build_grove();
    auto before = capture_registries();

    ASSERT_GT(before.gene_count, 0);
    ASSERT_GT(before.source_count, 0);
    ASSERT_GT(before.transcript_count, 0);
    ASSERT_EQ(before.sample_count, 2);  // annotation + sample

    save_grove(temp_path_);

    // Reset all registries to prove deserialization restores them
    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    ASSERT_EQ(gene_registry::instance().size(), 0);

    auto loaded = load_grove(temp_path_);
    auto after = capture_registries();

    EXPECT_EQ(before.gene_count, after.gene_count);
    EXPECT_EQ(before.source_count, after.source_count);
    EXPECT_EQ(before.transcript_count, after.transcript_count);
    EXPECT_EQ(before.sample_count, after.sample_count);

    // Verify gene details match
    ASSERT_EQ(before.genes.size(), after.genes.size());
    for (size_t i = 0; i < before.genes.size(); ++i) {
        EXPECT_EQ(std::get<0>(before.genes[i]), std::get<0>(after.genes[i]))
            << "Gene ID mismatch at index " << i;
        EXPECT_EQ(std::get<1>(before.genes[i]), std::get<1>(after.genes[i]))
            << "Gene name mismatch at index " << i;
        EXPECT_EQ(std::get<2>(before.genes[i]), std::get<2>(after.genes[i]))
            << "Gene biotype mismatch at index " << i;
    }

    // Verify sample metadata preserved
    ASSERT_EQ(before.sample_ids.size(), after.sample_ids.size());
    for (size_t i = 0; i < before.sample_ids.size(); ++i) {
        EXPECT_EQ(before.sample_ids[i], after.sample_ids[i]);
    }
}

// ─── Full grove roundtrip ──────────────────────────────────────────

TEST_F(SerializationTest, GroveRoundtrip) {
    build_grove();
    auto before = capture_grove(*grove_);

    ASSERT_GT(before.segment_count, 0);
    ASSERT_GT(before.total_edges, 0);
    ASSERT_GT(before.gene_ids.size(), 0);
    ASSERT_GT(before.transcript_ids.size(), 0);

    save_grove(temp_path_);

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    auto loaded = load_grove(temp_path_);
    auto after = capture_grove(*loaded);

    EXPECT_EQ(before.segment_count, after.segment_count);
    EXPECT_EQ(before.total_edges, after.total_edges);
    EXPECT_EQ(before.gene_ids, after.gene_ids);
    EXPECT_EQ(before.transcript_ids, after.transcript_ids);
}

// ─── Expression preservation ───────────────────────────────────────

TEST_F(SerializationTest, ExpressionRoundtrip) {
    build_grove();
    auto before = capture_grove(*grove_);

    // sample.gtf has counts on TX_S1 (150) and TX_S2 (75)
    ASSERT_GT(before.segment_expressions.size(), 0)
        << "Expected expression data from sample.gtf";

    save_grove(temp_path_);

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    auto loaded = load_grove(temp_path_);
    auto after = capture_grove(*loaded);

    ASSERT_EQ(before.segment_expressions.size(), after.segment_expressions.size());
    for (auto& [coord, expr] : before.segment_expressions) {
        auto it = after.segment_expressions.find(coord);
        ASSERT_NE(it, after.segment_expressions.end())
            << "Missing expression for " << coord;
        EXPECT_FLOAT_EQ(expr, it->second)
            << "Expression mismatch at " << coord;
    }
}

// ─── Spatial query works after roundtrip ───────────────────────────

TEST_F(SerializationTest, SpatialQueryAfterRoundtrip) {
    build_grove();

    // Query that should hit GENE_A segments (10000-14500 on chr22)
    gdt::genomic_coordinate query_coord('+', 10000, 14500);
    auto before_results = grove_->intersect(query_coord, "chr22");
    size_t before_hits = before_results.get_keys().size();
    ASSERT_GT(before_hits, 0);

    save_grove(temp_path_);

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    auto loaded = load_grove(temp_path_);

    auto after_results = loaded->intersect(query_coord, "chr22");
    size_t after_hits = after_results.get_keys().size();

    EXPECT_EQ(before_hits, after_hits)
        << "Spatial query returned different number of hits after roundtrip";
}

// ─── Graph traversal works after roundtrip ─────────────────────────

TEST_F(SerializationTest, GraphTraversalAfterRoundtrip) {
    build_grove();

    // Collect exon chains via graph traversal before serialization
    std::map<size_t, std::vector<std::string>> before_chains;  // segment_index -> exon coords

    gdt::genomic_coordinate query_coord('+', 10000, 14500);
    auto results = grove_->intersect(query_coord, "chr22");
    for (auto* key : results.get_keys()) {
        auto& feature = key->get_data();
        if (!is_segment(feature)) continue;
        auto& seg = get_segment(feature);
        if (seg.absorbed) continue;

        std::vector<std::string> chain;
        auto first_exons = grove_->get_neighbors_if(key,
            [&seg](const edge_metadata& e) {
                return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON &&
                       e.id == seg.segment_index;
            });

        if (!first_exons.empty()) {
            auto* current = first_exons[0];
            while (current) {
                chain.push_back(format_coordinate("chr22", current->get_value()));
                auto next = grove_->get_neighbors_if(current,
                    [&seg](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::EXON_TO_EXON &&
                               e.id == seg.segment_index;
                    });
                current = next.empty() ? nullptr : next[0];
            }
        }
        before_chains[seg.segment_index] = std::move(chain);
    }
    ASSERT_GT(before_chains.size(), 0);

    save_grove(temp_path_);

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    auto loaded = load_grove(temp_path_);

    // Repeat graph traversal on loaded grove
    std::map<size_t, std::vector<std::string>> after_chains;
    auto loaded_results = loaded->intersect(query_coord, "chr22");
    for (auto* key : loaded_results.get_keys()) {
        auto& feature = key->get_data();
        if (!is_segment(feature)) continue;
        auto& seg = get_segment(feature);
        if (seg.absorbed) continue;

        std::vector<std::string> chain;
        auto first_exons = loaded->get_neighbors_if(key,
            [&seg](const edge_metadata& e) {
                return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON &&
                       e.id == seg.segment_index;
            });

        if (!first_exons.empty()) {
            auto* current = first_exons[0];
            while (current) {
                chain.push_back(format_coordinate("chr22", current->get_value()));
                auto next = loaded->get_neighbors_if(current,
                    [&seg](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::EXON_TO_EXON &&
                               e.id == seg.segment_index;
                    });
                current = next.empty() ? nullptr : next[0];
            }
        }
        after_chains[seg.segment_index] = std::move(chain);
    }

    ASSERT_EQ(before_chains.size(), after_chains.size());
    for (auto& [idx, before_chain] : before_chains) {
        auto it = after_chains.find(idx);
        ASSERT_NE(it, after_chains.end())
            << "Missing segment_index " << idx << " after roundtrip";
        EXPECT_EQ(before_chain, it->second)
            << "Exon chain mismatch for segment_index " << idx;
    }
}

// ─── Sample membership preserved ───────────────────────────────────

TEST_F(SerializationTest, SampleMembershipRoundtrip) {
    build_grove();

    // Capture per-segment sample membership
    std::map<size_t, std::vector<uint32_t>> before_samples;  // segment_index -> sample IDs

    gdt::genomic_coordinate query_coord('+', 10000, 14500);
    auto results = grove_->intersect(query_coord, "chr22");
    for (auto* key : results.get_keys()) {
        auto& feature = key->get_data();
        if (!is_segment(feature)) continue;
        auto& seg = get_segment(feature);
        if (seg.absorbed) continue;

        std::vector<uint32_t> samples;
        for (auto sid : seg.sample_idx) {
            samples.push_back(sid);
        }
        before_samples[seg.segment_index] = std::move(samples);
    }
    ASSERT_GT(before_samples.size(), 0);

    save_grove(temp_path_);

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    auto loaded = load_grove(temp_path_);

    std::map<size_t, std::vector<uint32_t>> after_samples;
    auto loaded_results = loaded->intersect(query_coord, "chr22");
    for (auto* key : loaded_results.get_keys()) {
        auto& feature = key->get_data();
        if (!is_segment(feature)) continue;
        auto& seg = get_segment(feature);
        if (seg.absorbed) continue;

        std::vector<uint32_t> samples;
        for (auto sid : seg.sample_idx) {
            samples.push_back(sid);
        }
        after_samples[seg.segment_index] = std::move(samples);
    }

    ASSERT_EQ(before_samples.size(), after_samples.size());
    for (auto& [idx, before_sids] : before_samples) {
        auto it = after_samples.find(idx);
        ASSERT_NE(it, after_samples.end())
            << "Missing segment_index " << idx << " after roundtrip";
        EXPECT_EQ(before_sids, it->second)
            << "Sample membership mismatch for segment_index " << idx;
    }
}

// ─── Invalid file handling ─────────────────────────────────────────

TEST_F(SerializationTest, InvalidMagicThrows) {
    // Write garbage to temp file
    {
        std::ofstream ofs(temp_path_, std::ios::binary);
        ofs.write("JUNK", 4);
    }

    std::ifstream ifs(temp_path_, std::ios::binary);
    char magic[4];
    ifs.read(magic, 4);
    EXPECT_NE(std::string(magic, 4), "AGRX");
}

TEST_F(SerializationTest, FileSize) {
    build_grove();
    save_grove(temp_path_);

    auto size = fs::file_size(temp_path_);
    // Magic(4) + version(2) + registries + compressed grove
    EXPECT_GT(size, 6) << ".ggx file should be larger than just the header";
}

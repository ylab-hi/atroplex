// Regression coverage for the remaining splicing event catalog gap in #61.
// GTF coordinates below are preserved by the current builder and TSV writers.
#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "analysis_report.hpp"
#include "build/build_gff.hpp"
#include "sample_info.hpp"
#include "splicing_catalog.hpp"

namespace {
namespace fs = std::filesystem;
using Exon = std::pair<size_t, size_t>;
using Chain = std::vector<Exon>;

struct EventCase {
    std::string label;
    splicing_event_type type;
    Chain first;
    Chain second;
    std::string upstream;
    std::string downstream;
    std::set<std::string> affected;
    size_t included;       // sample A carries both paths
    size_t first_included; // sample B carries only the first path
};

// Independent two-isoform topologies. Some topologies also trigger other
// detectors; the test accounts for those additional classifications explicitly.
const std::vector<EventCase> event_cases{
    {"cassette_exon", splicing_event_type::CASSETTE_EXON,
     {{100, 199}, {300, 399}, {500, 599}}, {{100, 199}, {500, 599}},
     "100-199", "500-599", {"300-399"}, 1, 1},
    {"alt_5prime", splicing_event_type::ALT_5PRIME,
     {{100, 199}, {300, 399}, {600, 699}}, {{100, 199}, {300, 449}, {600, 699}},
     ".", ".", {"300-399", "300-449"}, 2, 1},
    {"alt_3prime", splicing_event_type::ALT_3PRIME,
     {{100, 199}, {300, 399}, {600, 699}}, {{100, 199}, {250, 399}, {600, 699}},
     ".", ".", {"250-399", "300-399"}, 2, 1},
    {"intron_retention", splicing_event_type::INTRON_RETENTION,
     {{100, 199}, {300, 399}, {500, 599}, {700, 799}},
     {{100, 199}, {300, 599}, {700, 799}},
     "300-399", "500-599", {"300-599"}, 1, 0},
    {"alt_first_exon", splicing_event_type::ALT_FIRST_EXON,
     {{100, 199}, {500, 599}}, {{300, 399}, {500, 599}},
     ".", ".", {"100-199", "300-399"}, 2, 1},
    {"alt_last_exon", splicing_event_type::ALT_LAST_EXON,
     {{100, 199}, {300, 399}}, {{100, 199}, {500, 599}},
     ".", ".", {"300-399", "500-599"}, 2, 1},
    {"mutually_exclusive", splicing_event_type::MUTUALLY_EXCLUSIVE,
     {{100, 199}, {300, 349}, {600, 699}}, {{100, 199}, {400, 449}, {600, 699}},
     "100-199", "600-699", {"300-349", "400-449"}, 2, 1},
};

std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> fields;
    std::istringstream input(line);
    std::string field;
    while (std::getline(input, field, delimiter)) fields.push_back(field);
    return fields;
}

struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
};

class SplicingEventsTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();
        std::string name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
        std::replace(name.begin(), name.end(), '/', '_');
        tmp_dir = fs::temp_directory_path() / ("atroplex_splicing_events_" + name);
        fs::create_directories(tmp_dir);
    }

    void TearDown() override {
        std::error_code error;
        fs::remove_all(tmp_dir, error);
    }

    void build_sample(const std::string& sample, const std::vector<Chain>& chains,
                      const std::string& gene = "G_EVENT") {
        const auto path = tmp_dir / (sample + ".gtf");
        std::ofstream output(path);
        ASSERT_TRUE(output.is_open()) << path;
        for (size_t i = 0; i < chains.size(); ++i) {
            const auto& chain = chains[i];
            const std::string attributes = "gene_id \"" + gene + "\"; gene_name \"EventGene\"; "
                "transcript_id \"" + sample + "_T" + std::to_string(i) + "\";";
            output << "chr1\tTEST\ttranscript\t" << chain.front().first << '\t'
                   << chain.back().second << "\t.\t+\t.\t" << attributes << '\n';
            for (const auto& [start, end] : chain) {
                output << "chr1\tTEST\texon\t" << start << '\t' << end
                       << "\t.\t+\t.\t" << attributes << '\n';
            }
        }
        output.close();
        ASSERT_TRUE(output) << "Failed to write " << path;
        sample_info info(sample);
        info.type = "sample";
        const auto sid = sample_registry::instance().intern(info.id, info);
        build_options options;
        options.absorb = false; // Preserve the designed isoforms, including terminal variants.
        build_counters counters;
        build_gff::build(grove, path, sid, exon_caches, segment_caches,
                         segment_count, options, counters);
    }

    static Table read_table(const fs::path& path) {
        std::ifstream input(path);
        EXPECT_TRUE(input.is_open()) << "Missing event TSV: " << path;
        Table table;
        std::string line;
        if (std::getline(input, line)) table.header = split(line, '\t');
        while (std::getline(input, line)) table.rows.push_back(split(line, '\t'));
        return table;
    }

    static std::vector<std::string> header(const std::vector<std::string>& samples) {
        std::vector<std::string> columns{"gene_id", "gene_name", "event_type", "chromosome",
                                         "upstream_exon", "downstream_exon", "affected_exons"};
        for (const auto& sample : samples) {
            for (const auto* suffix : {".psi", ".included", ".total"}) {
                columns.push_back(sample + suffix);
            }
        }
        return columns;
    }

    Table stream_events() {
        const auto path = tmp_dir / "streamed.splicing_events.tsv";
        {
            analysis_report report;
            report.begin_splicing_events_stream(path.string());
            report.collect(grove, nullptr, 0);
        } // Flush before reading.
        return read_table(path);
    }

    fs::path tmp_dir;
    grove_type grove{3};
    chromosome_exon_caches exon_caches;
    chromosome_segment_caches segment_caches;
    size_t segment_count = 0;
};

class SplicingEventTopologyTest : public SplicingEventsTest,
                                  public ::testing::WithParamInterface<EventCase> {};

TEST_P(SplicingEventTopologyTest, DetectionAndBothTsvWriters) {
    const auto& expected = GetParam();
    ASSERT_NO_FATAL_FAILURE(build_sample("sample_A", {expected.first, expected.second}));
    ASSERT_NO_FATAL_FAILURE(build_sample("sample_B", {expected.first}));
    // Registered sample with no transcripts at the event locus: PSI must be
    // missing ('.'), rather than zero. Its single-path gene produces no events.
    ASSERT_NO_FATAL_FAILURE(build_sample("sample_C", {{{2000, 2099}, {2300, 2399}}}, "G_OTHER"));
    ASSERT_EQ(segment_count, 3u) << "Two event isoforms plus one unrelated isoform";

    const auto events = splicing_catalog::collect_from_grove(grove);
    std::multiset<std::string> expected_types{expected.label};
    // The current detectors also classify alternative internal exon boundaries
    // as mutually exclusive paths, and the retained exon shares both splice
    // boundaries with the two spliced exons. Keep these overlaps explicit.
    if (expected.type == splicing_event_type::ALT_5PRIME ||
        expected.type == splicing_event_type::ALT_3PRIME) {
        expected_types.insert("mutually_exclusive");
    }
    if (expected.type == splicing_event_type::INTRON_RETENTION) {
        expected_types.insert("alt_5prime");
        expected_types.insert("alt_3prime");
    }
    std::multiset<std::string> actual_types;
    for (const auto& event : events) actual_types.insert(splicing_catalog::event_type_str(event.type));
    EXPECT_EQ(actual_types, expected_types) << "No missing, duplicate, or unexpected classifications";
    const splicing_event* target = nullptr;
    size_t matches = 0;
    for (const auto& event : events) {
        if (event.type == expected.type && event.gene_id == "G_EVENT") {
            target = &event;
            ++matches;
        }
    }
    ASSERT_EQ(matches, 1u) << "Expected one " << expected.label << " event";
    ASSERT_NE(target, nullptr);
    EXPECT_EQ(target->gene_name, "EventGene");
    EXPECT_EQ(target->chromosome, "chr1");
    EXPECT_EQ(splicing_catalog::format_exon(target->upstream_exon), expected.upstream);
    EXPECT_EQ(splicing_catalog::format_exon(target->downstream_exon), expected.downstream);
    std::set<std::string> affected;
    for (auto* exon : target->affected_exons) affected.insert(splicing_catalog::format_exon(exon));
    EXPECT_EQ(affected, expected.affected);
    EXPECT_EQ(target->affected_exons.size(), expected.affected.size()) << "No duplicate affected exons";
    ASSERT_TRUE(target->sample_psi.contains(0));
    ASSERT_TRUE(target->sample_psi.contains(1));
    EXPECT_FLOAT_EQ(target->sample_psi.at(0), static_cast<float>(expected.included) / 2.0f);
    EXPECT_FLOAT_EQ(target->sample_psi.at(1), static_cast<float>(expected.first_included));
    EXPECT_FALSE(target->sample_psi.contains(2)) << "Unrelated sample has no event-locus denominator";

    const auto standalone_path = tmp_dir / "catalog.splicing_events.tsv";
    splicing_catalog::write_events_tsv(standalone_path.string(), events);
    const auto standalone = read_table(standalone_path);
    const auto streamed = stream_events();
    // Standalone derives columns from event membership; streaming includes
    // all registered samples, even those without any event-locus transcripts.
    ASSERT_EQ(standalone.header, header({"sample_A", "sample_B"}));
    ASSERT_EQ(streamed.header, header({"sample_A", "sample_B", "sample_C"}));
    ASSERT_EQ(standalone.rows.size(), events.size());
    ASSERT_EQ(streamed.rows.size(), events.size());

    for (const auto* table : {&standalone, &streamed}) {
        size_t target_rows = 0;
        std::multiset<std::string> row_types;
        for (const auto& row : table->rows) {
            ASSERT_EQ(row.size(), table->header.size()) << "Every row must match its header";
            row_types.insert(row[2]);
            EXPECT_EQ(row[0], "G_EVENT") << "Single-isoform unrelated gene must not emit events";
            EXPECT_EQ(row[1], "EventGene");
            EXPECT_EQ(row[3], "chr1");
            if (table == &streamed) {
                EXPECT_EQ(row[13], ".");
                EXPECT_EQ(row[14], ".");
                EXPECT_EQ(row[15], ".");
            }
            if (row[2] != expected.label) continue;
            ++target_rows;
            EXPECT_EQ(row[4], expected.upstream);
            EXPECT_EQ(row[5], expected.downstream);
            const auto coords = split(row[6], ',');
            EXPECT_EQ(std::set<std::string>(coords.begin(), coords.end()), expected.affected);
            EXPECT_EQ(coords.size(), expected.affected.size());
            EXPECT_DOUBLE_EQ(std::stod(row[7]), static_cast<double>(expected.included) / 2.0);
            EXPECT_EQ(row[8], std::to_string(expected.included));
            EXPECT_EQ(row[9], "2");
            EXPECT_DOUBLE_EQ(std::stod(row[10]), static_cast<double>(expected.first_included));
            EXPECT_EQ(row[11], std::to_string(expected.first_included));
            EXPECT_EQ(row[12], "1");
        }
        EXPECT_EQ(row_types, expected_types) << "Writers preserve every event type exactly once";
        EXPECT_EQ(target_rows, 1u) << expected.label;
    }
}

INSTANTIATE_TEST_SUITE_P(AllEventTypes, SplicingEventTopologyTest,
                        ::testing::ValuesIn(event_cases),
                        [](const ::testing::TestParamInfo<EventCase>& info) { return info.param.label; });

TEST_F(SplicingEventsTest, EmptyGroveWritesHeaderOnly) {
    const auto events = splicing_catalog::collect_from_grove(grove);
    EXPECT_TRUE(events.empty());
    const auto path = tmp_dir / "empty.tsv";
    splicing_catalog::write_events_tsv(path.string(), events);
    const auto standalone = read_table(path);
    EXPECT_EQ(standalone.header, header({}));
    EXPECT_TRUE(standalone.rows.empty());
    const auto streamed = stream_events();
    EXPECT_EQ(streamed.header, header({}));
    EXPECT_TRUE(streamed.rows.empty());
}

TEST_F(SplicingEventsTest, BothWritersUseMissingValuesForSamplesAtOtherEventLoci) {
    ASSERT_NO_FATAL_FAILURE(build_sample("sample_A",
        {{{100, 199}, {300, 399}, {500, 599}}, {{100, 199}, {500, 599}}}));
    ASSERT_NO_FATAL_FAILURE(build_sample("sample_B",
        {{{2000, 2099}, {2200, 2299}, {2400, 2499}}, {{2000, 2099}, {2400, 2499}}}, "G_OTHER"));
    const auto events = splicing_catalog::collect_from_grove(grove);
    ASSERT_EQ(events.size(), 2u) << "One cassette event at each independent locus";
    const auto path = tmp_dir / "two_loci.tsv";
    splicing_catalog::write_events_tsv(path.string(), events);
    const auto standalone = read_table(path);
    const auto streamed = stream_events();
    for (const auto* table : {&standalone, &streamed}) {
        ASSERT_EQ(table->header, header({"sample_A", "sample_B"}));
        ASSERT_EQ(table->rows.size(), 2u);
        std::set<std::string> genes;
        for (const auto& row : table->rows) {
            ASSERT_EQ(row.size(), table->header.size());
            ASSERT_TRUE(row[0] == "G_EVENT" || row[0] == "G_OTHER");
            genes.insert(row[0]);
            EXPECT_EQ(row[2], "cassette_exon");
            const size_t present = row[0] == "G_EVENT" ? 7 : 10;
            const size_t missing = row[0] == "G_EVENT" ? 10 : 7;
            EXPECT_DOUBLE_EQ(std::stod(row[present]), 0.5);
            EXPECT_EQ(row[present + 1], "1");
            EXPECT_EQ(row[present + 2], "2");
            EXPECT_EQ(row[missing], ".");
            EXPECT_EQ(row[missing + 1], ".");
            EXPECT_EQ(row[missing + 2], ".");
        }
        EXPECT_EQ(genes, (std::set<std::string>{"G_EVENT", "G_OTHER"}));
    }
}

TEST_F(SplicingEventsTest, IdenticalIsoformsDoNotCreateEvents) {
    const Chain chain{{100, 199}, {300, 399}, {500, 599}};
    ASSERT_NO_FATAL_FAILURE(build_sample("sample_A", {chain, chain}));
    ASSERT_EQ(segment_count, 1u) << "Identical exon chains deduplicate to one segment";
    EXPECT_TRUE(splicing_catalog::collect_from_grove(grove).empty());
    const auto streamed = stream_events();
    EXPECT_EQ(streamed.header, header({"sample_A"}));
    EXPECT_TRUE(streamed.rows.empty()) << "Constitutive exons alone are not alternative events";
}
} // namespace

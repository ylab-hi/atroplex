/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/inspect.hpp"

#include <filesystem>

#include "analysis_report.hpp"
#include "utility.hpp"

namespace subcall {

cxxopts::Options inspect::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex inspect",
        "Pan-transcriptome inspection: overview, sharing, splicing hubs and events");

    add_common_options(options);

    options.add_options("Filtering")
        ("min-samples", "Only include segments present in >= N samples. "
            "Applied non-destructively during traversal — the grove is not modified.",
            cxxopts::value<size_t>()->default_value("0"))
        ("conserved-fraction", "Fraction of sample-typed entries (0,1] a segment "
            "or exon must appear in to be classified as conserved. Default 1.0 = "
            "strict (must be in every sample). Relaxing (e.g. 0.95) yields a "
            "dropout-tolerant 'conserved core' for both per-sample counts and "
            "the conserved_segments / conserved_exons TSVs.",
            cxxopts::value<double>()->default_value("1.0"))
        ;

    options.add_options("Output")
        ("events", "Write per-gene splicing event catalog (cassette, alt-5'/3', IR, "
            "alt-terminal, mutually exclusive). Off by default — the file can be very "
            "large on cohort-scale builds. Hub analysis always runs regardless.")
        ;

    return options;
}

void inspect::validate(const cxxopts::ParseResult& args) {
    if (!args.count("build-from") && !args.count("manifest") && !args.count("genogrove")) {
        throw std::runtime_error(
            "No input specified. Use -m/--manifest, -b/--build-from, or -g/--genogrove");
    }

    if (args.count("manifest")) {
        std::string manifest_path = args["manifest"].as<std::string>();
        if (!std::filesystem::exists(manifest_path)) {
            throw std::runtime_error("Manifest file not found: " + manifest_path);
        }
    }

    if (args.count("build-from")) {
        auto files = args["build-from"].as<std::vector<std::string>>();
        for (const auto& f : files) {
            if (!std::filesystem::exists(f)) {
                throw std::runtime_error("Annotation file not found: " + f);
            }
        }
    }

    if (args.count("conserved-fraction")) {
        double f = args["conserved-fraction"].as<double>();
        if (!(f > 0.0 && f <= 1.0)) {
            throw std::runtime_error(
                "--conserved-fraction must be in (0, 1]; got " + std::to_string(f));
        }
    }
}

void inspect::execute(const cxxopts::ParseResult& args) {
    if (!grove) {
        throw std::runtime_error("Grove not available for inspection");
    }

    // Determine output path
    std::string fallback;
    if (args.count("manifest")) {
        fallback = args["manifest"].as<std::string>();
    } else if (args.count("build-from")) {
        fallback = args["build-from"].as<std::vector<std::string>>().front();
    }

    // Outputs land directly under -o, grouped by category subfolders
    // (overview / sharing / splicing_hubs / splicing_events). No extra
    // wrapper folder — matches the convention of the other subcommands.
    auto out_dir = resolve_output_dir(args, fallback);
    std::filesystem::create_directories(out_dir);

    std::string basename = resolve_prefix(args);

    // ── Streaming analysis report (Phase 8) ──────────────────────────
    // Single-pass overview stats — scalable to 20K+ samples
    auto overview_dir = out_dir / "overview";
    std::filesystem::create_directories(overview_dir);

    analysis_report report;
    report.set_conserved_fraction(args["conserved-fraction"].as<double>());

    // Whether to include per-sample expression columns in hub / branch
    // / conserved-exon / sample_stats outputs. Enabled when a .qtx
    // sidecar was loaded alongside the grove in setup_grove.
    const bool emit_expression_cols = qtx_reader_ptr() != nullptr;

    // Phase 8.5: open the conserved-exons stream BEFORE collect() so
    // rows are written inline during traversal (no accumulator).
    auto sharing_dir = out_dir / "sharing";
    std::filesystem::create_directories(sharing_dir);
    report.begin_conserved_exon_stream(
        (sharing_dir / (basename + ".conserved_exons.tsv")).string(),
        emit_expression_cols);
    report.begin_conserved_segment_stream(
        (sharing_dir / (basename + ".conserved_segments.tsv")).string(),
        emit_expression_cols);

    // Phase 8.3: open splicing hub streams before collect() — rows are
    // emitted at each gene's finalization (not at the end of the grove)
    // so peak memory stays at O(num_samples), not O(num_hubs × num_samples).
    auto hubs_dir = out_dir / "splicing_hubs";
    std::filesystem::create_directories(hubs_dir);
    report.begin_splicing_hub_streams(
        (hubs_dir / (basename + ".splicing_hubs.tsv")).string(),
        (hubs_dir / (basename + ".branch_details.tsv")).string(),
        emit_expression_cols);

    // Phase 8.6: splicing events are opt-in (--events) because the
    // per-event TSV can be very large on cohort-scale builds. When
    // disabled, the exon chain capture and detect_gene_events() call
    // are both skipped (gated by splicing_events_stream == nullptr).
    if (args.count("events")) {
        auto events_dir = out_dir / "splicing_events";
        std::filesystem::create_directories(events_dir);
        report.begin_splicing_events_stream(
            (events_dir / (basename + ".splicing_events.tsv")).string());
    }

    size_t min_samples = args["min-samples"].as<size_t>();
    report.collect(*grove, qtx_reader_ptr(), min_samples);

    report.write_overview((overview_dir / (basename + ".overview.tsv")).string());
    report.write_per_sample((overview_dir / (basename + ".per_sample.tsv")).string());
    // Phase 8.7a: per-source aggregation (HAVANA / ENSEMBL / TALON / ...).
    // Bounded by source_registry capacity (16 sources max), no extra
    // traversal — counters are populated inline during collect().
    report.write_per_source((overview_dir / (basename + ".per_source.tsv")).string());
    // Phase 8.7a: biotype breakdown (long-form gene + transcript rows
    // with global + per-sample columns). Bounded by ~30 biotype strings;
    // counters are populated inline during collect().
    report.write_biotype((overview_dir / (basename + ".biotype.tsv")).string());

    // Phase 8.2: sharing TSVs — metrics × samples transposition of
    // per_sample counters, no extra memory
    report.write_exon_sharing((sharing_dir / (basename + ".exon_sharing.tsv")).string());
    report.write_segment_sharing((sharing_dir / (basename + ".segment_sharing.tsv")).string());

    logging::info("Inspection written to: " + out_dir.string());
}

} // namespace subcall
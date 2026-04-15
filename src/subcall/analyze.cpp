/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/analyze.hpp"

#include <filesystem>

#include "analysis_report.hpp"
#include "utility.hpp"

namespace subcall {

cxxopts::Options analyze::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex analyze",
        "Full pan-transcriptome analysis with per-sample statistics");

    add_common_options(options);

    return options;
}

void analyze::validate(const cxxopts::ParseResult& args) {
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
}

void analyze::execute(const cxxopts::ParseResult& args) {
    if (!grove) {
        throw std::runtime_error("Grove not available for analysis");
    }

    // Determine output path
    std::string fallback;
    if (args.count("manifest")) {
        fallback = args["manifest"].as<std::string>();
    } else if (args.count("build-from")) {
        fallback = args["build-from"].as<std::vector<std::string>>().front();
    }

    auto out_dir = resolve_output_dir(args, fallback);
    auto analysis_dir = out_dir / "analysis";
    std::filesystem::create_directories(analysis_dir);

    std::string basename = resolve_prefix(args);

    // ── Streaming analysis report (Phase 8) ──────────────────────────
    // Single-pass overview stats — scalable to 20K+ samples
    {
        auto overview_dir = analysis_dir / "overview";
        std::filesystem::create_directories(overview_dir);

        analysis_report report;

        // Phase 8.5: open the conserved-exons stream BEFORE collect() so
        // rows are written inline during traversal (no accumulator).
        auto sharing_dir = analysis_dir / "sharing";
        std::filesystem::create_directories(sharing_dir);
        report.begin_conserved_exon_stream(
            (sharing_dir / (basename + ".conserved_exons.tsv")).string());

        // Phase 8.3: open splicing hub streams before collect() — rows are
        // emitted at each gene's finalization (not at the end of the grove)
        // so peak memory stays at O(num_samples), not O(num_hubs × num_samples).
        auto hubs_dir = analysis_dir / "splicing_hubs";
        std::filesystem::create_directories(hubs_dir);
        report.begin_splicing_hub_streams(
            (hubs_dir / (basename + ".splicing_hubs.tsv")).string(),
            (hubs_dir / (basename + ".branch_details.tsv")).string());

        // Phase 8.6: open splicing events stream before collect() — events
        // are detected per-gene at gene finalization and rows are streamed
        // inline. Per-gene segment chains are captured during the existing
        // exon walk and dropped at chromosome boundary.
        auto events_dir = analysis_dir / "splicing_events";
        std::filesystem::create_directories(events_dir);
        report.begin_splicing_events_stream(
            (events_dir / (basename + ".splicing_events.tsv")).string());

        report.collect(*grove);

        report.write_overview((overview_dir / (basename + ".overview.tsv")).string());
        report.write_per_sample((overview_dir / (basename + ".per_sample.tsv")).string());
        // Phase 8.7a: per-source aggregation (HAVANA / ENSEMBL / TALON / ...).
        // Bounded by source_registry capacity (16 sources max), no extra
        // traversal — counters are populated inline during collect().
        report.write_per_source((overview_dir / (basename + ".per_source.tsv")).string());

        // Phase 8.2: sharing TSVs — metrics × samples transposition of
        // per_sample counters, no extra memory
        report.write_exon_sharing((sharing_dir / (basename + ".exon_sharing.tsv")).string());
        report.write_segment_sharing((sharing_dir / (basename + ".segment_sharing.tsv")).string());
    }

    // TODO(Phase 8): Legacy analysis disabled — OOMs at 1000+ samples.
    // Will be replaced incrementally by analysis_report phases 8.2-8.7.
    // Disabled code covers: sharing stats, splicing hubs, splicing events,
    // diversity metrics, per-sample CSV, full text report.
    //
    // To re-enable for small groves, uncomment the block below.

    /*
    index_stats::collect_options opts;
    opts.detailed = true;
    opts.output_dir = analysis_dir.string();
    opts.basename = basename;
    auto stats = index_stats::collect(*grove, opts);

    stats.write((analysis_dir / (basename + ".analysis.txt")).string());
    stats.write_overview_tsv((analysis_dir / (basename + ".overview.tsv")).string());
    if (!stats.per_chromosome.empty())
        stats.write_per_chromosome_tsv((analysis_dir / (basename + ".per_chromosome.tsv")).string());
    if (!stats.per_source.empty())
        stats.write_per_source_tsv((analysis_dir / (basename + ".per_source.tsv")).string());
    if (!stats.per_sample.empty()) {
        stats.write_per_sample_tsv((analysis_dir / (basename + ".per_sample.tsv")).string());
        auto sharing_dir = analysis_dir / "sharing";
        std::filesystem::create_directories(sharing_dir);
        stats.write_exon_sharing_tsv((sharing_dir / (basename + ".exon_sharing.tsv")).string());
        stats.write_segment_sharing_tsv((sharing_dir / (basename + ".segment_sharing.tsv")).string());
    }
    if (!stats.genes_by_biotype.empty() || !stats.transcripts_by_biotype.empty())
        stats.write_biotype_tsv((analysis_dir / (basename + ".biotype.tsv")).string());
    if (!stats.splicing_hubs.empty()) {
        auto hubs_dir = analysis_dir / "splicing_hubs";
        std::filesystem::create_directories(hubs_dir);
        stats.write_splicing_hubs_tsv((hubs_dir / (basename + ".splicing_hubs.tsv")).string());
    }
    {
        auto events_dir = analysis_dir / "splicing_events";
        std::filesystem::create_directories(events_dir);
        auto events = gene_indices_.empty()
            ? splicing_catalog::collect_from_grove(*grove)
            : splicing_catalog::collect(gene_indices_, *grove);
        splicing_catalog::write_events_tsv((events_dir / (basename + ".splicing_events.tsv")).string(), events);
        splicing_catalog::write_summary((events_dir / (basename + ".splicing_events.summary.txt")).string(), events);
    }
    */

    logging::info("Analysis written to: " + analysis_dir.string());
}

} // namespace subcall
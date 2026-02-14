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

#include "index_stats.hpp"
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
    std::string first_input;
    if (args.count("manifest")) {
        first_input = args["manifest"].as<std::string>();
    } else if (args.count("build-from")) {
        first_input = args["build-from"].as<std::vector<std::string>>().front();
    }

    auto out_dir = resolve_output_dir(args, first_input);
    auto analysis_dir = out_dir / "analysis";
    std::filesystem::create_directories(analysis_dir);

    std::string basename = first_input.empty()
        ? "atroplex"
        : std::filesystem::path(first_input).stem().string();

    // Full collection including Jaccard diversity (Phase 4b)
    // Stream large outputs (conserved exons, branch details) directly to files
    logging::info("Running pan-transcriptome analysis...");
    index_stats::collect_options opts;
    opts.detailed = true;
    opts.output_dir = analysis_dir.string();
    opts.basename = basename;
    auto stats = index_stats::collect(*grove, opts);

    // Write full text report
    std::string stats_path = (analysis_dir / (basename + ".analysis.txt")).string();
    stats.write(stats_path);

    // Write overview TSV (global stats)
    std::string overview_path = (analysis_dir / (basename + ".overview.tsv")).string();
    stats.write_overview_tsv(overview_path);

    // Write per-chromosome TSV
    if (!stats.per_chromosome.empty()) {
        std::string chr_path = (analysis_dir / (basename + ".per_chromosome.tsv")).string();
        stats.write_per_chromosome_tsv(chr_path);
    }

    // Write per-source TSV
    if (!stats.per_source.empty()) {
        std::string source_path = (analysis_dir / (basename + ".per_source.tsv")).string();
        stats.write_per_source_tsv(source_path);
    }

    // Write per-sample TSV
    if (!stats.per_sample.empty()) {
        std::string sample_path = (analysis_dir / (basename + ".per_sample.tsv")).string();
        stats.write_per_sample_tsv(sample_path);
    }



    // Write exon and segment sharing into subfolder
    // (conserved_exons.tsv is already streamed by collect())
    if (!stats.per_sample.empty()) {
        auto sharing_dir = analysis_dir / "sharing";
        std::filesystem::create_directories(sharing_dir);

        std::string exon_path = (sharing_dir / (basename + ".exon_sharing.tsv")).string();
        stats.write_exon_sharing_tsv(exon_path);

        std::string seg_path = (sharing_dir / (basename + ".segment_sharing.tsv")).string();
        stats.write_segment_sharing_tsv(seg_path);
    }

    // Write biotype breakdown
    if (!stats.genes_by_biotype.empty() || !stats.transcripts_by_biotype.empty()) {
        std::string biotype_path = (analysis_dir / (basename + ".biotype.tsv")).string();
        stats.write_biotype_tsv(biotype_path);
    }

    // Write splicing hubs summary (branch_details.tsv is already streamed by collect())
    if (!stats.splicing_hubs.empty()) {
        auto hubs_dir = analysis_dir / "splicing_hubs";
        std::filesystem::create_directories(hubs_dir);

        std::string hubs_path = (hubs_dir / (basename + ".splicing_hubs.tsv")).string();
        stats.write_splicing_hubs_tsv(hubs_path);
    }

    logging::info("Analysis written to: " + analysis_dir.string());
}

} // namespace subcall
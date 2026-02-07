/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/discover.hpp"

#include <filesystem>

#include "utility.hpp"

namespace gio = genogrove::io;

namespace subcall {

cxxopts::Options discover::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex discover",
        "Detect transcripts from long-read sequencing data");

    options.add_options("Input/Output")
        ("i,input", "Input BAM/SAM file with aligned reads",
            cxxopts::value<std::string>())
        ;

    options.add_options("Clustering")
        ("junction-bin", "Bin size for fuzzy junction clustering (bp)",
            cxxopts::value<int>()->default_value("10"))
        ("junction-tolerance", "Max junction position difference within cluster (bp)",
            cxxopts::value<int>()->default_value("5"))
        ("min-mapq", "Minimum mapping quality",
            cxxopts::value<int>()->default_value("20"))
        ;

    options.add_options("Matching")
        ("min-junction-score", "Minimum junction match score (0-1)",
            cxxopts::value<double>()->default_value("0.8"))
        ("min-overlap", "Minimum overlap with reference (bp)",
            cxxopts::value<int>()->default_value("50"))
        ;

    options.add_options("Novel transcript discovery")
        ("add-novel", "Add novel transcripts to the grove structure")
        ("min-novel-support", "Minimum read support to create a novel transcript",
            cxxopts::value<uint32_t>()->default_value("3"))
        ("min-novel-junctions", "Minimum novel junctions to classify as novel transcript",
            cxxopts::value<uint32_t>()->default_value("1"))
        ;

    add_common_options(options);

    return options;
}

void discover::validate(const cxxopts::ParseResult& args) {
    // Must have a grove source
    if (!args.count("genogrove") && !args.count("build-from") && !args.count("manifest")) {
        throw std::runtime_error("Must provide -g/--genogrove, -m/--manifest, or -b/--build-from");
    }

    // If input provided, check it exists
    if (args.count("input")) {
        std::string input = args["input"].as<std::string>();
        if (!std::filesystem::exists(input)) {
            throw std::runtime_error("Input file not found: " + input);
        }
    }

    // Check manifest file exists
    if (args.count("manifest")) {
        std::string manifest_path = args["manifest"].as<std::string>();
        if (!std::filesystem::exists(manifest_path)) {
            throw std::runtime_error("Manifest file not found: " + manifest_path);
        }
    }

    // Check build-from files exist
    if (args.count("build-from")) {
        auto files = args["build-from"].as<std::vector<std::string>>();
        for (const auto& f : files) {
            if (!std::filesystem::exists(f)) {
                throw std::runtime_error("Annotation file not found: " + f);
            }
        }
    }
}

void discover::execute(const cxxopts::ParseResult& args) {
    logging::info("Starting transcript discovery...");

    // Process reads if input provided
    if (args.count("input")) {
        process_reads(args);
    } else {
        logging::warning("No input file specified - grove built but no analysis performed");
    }

    logging::info("Discovery pipeline complete");
}

void discover::process_reads(const cxxopts::ParseResult& args) {
    std::string input_path = args["input"].as<std::string>();

    // Detect input file type
    gio::filetype_detector detector;
    auto [detected_type, detected_compression] = detector.detect_filetype(input_path);

    if (detected_type != gio::filetype::BAM && detected_type != gio::filetype::SAM) {
        logging::error("Input must be BAM/SAM file: " + input_path);
        return;
    }

    logging::info("Processing reads from: " + input_path);

    // Get parameters
    int junction_bin_size = args["junction-bin"].as<int>();
    int junction_tolerance = args["junction-tolerance"].as<int>();
    int min_mapq = args["min-mapq"].as<int>();
    double min_junction_score = args["min-junction-score"].as<double>();
    int min_overlap_bp = args["min-overlap"].as<int>();

    // Configure BAM reader
    gio::bam_reader_options opts = gio::bam_reader_options::primary_only();
    opts.min_mapq = min_mapq;
    gio::bam_reader reader(input_path, opts);

    // Step 1: Cluster reads by splice junction signature
    logging::info("Clustering reads by splice junction signature...");
    read_clusterer::config cluster_cfg;
    cluster_cfg.junction_bin_size = junction_bin_size;
    cluster_cfg.junction_tolerance = junction_tolerance;
    cluster_cfg.min_mapq = min_mapq;

    read_clusterer clusterer(cluster_cfg);
    auto clusters = clusterer.cluster_reads(reader);

    const auto& cluster_stats = clusterer.get_stats();
    logging::info("Clustering complete:");
    logging::info("  Total reads: " + std::to_string(cluster_stats.total_reads));
    logging::info("  Filtered (low mapq): " + std::to_string(cluster_stats.filtered_reads));
    logging::info("  Single-exon reads: " + std::to_string(cluster_stats.single_exon_reads));
    logging::info("  Multi-exon reads: " + std::to_string(cluster_stats.multi_exon_reads));
    logging::info("  Total clusters: " + std::to_string(cluster_stats.total_clusters));
    logging::info("  Mean cluster size: " + std::to_string(cluster_stats.mean_cluster_size));

    // Step 2: Match clusters to reference transcripts
    logging::info("Matching clusters to reference transcripts...");
    transcript_matcher::config match_cfg;
    match_cfg.junction_tolerance = junction_tolerance;
    match_cfg.min_junction_score = min_junction_score;
    match_cfg.min_overlap_bp = min_overlap_bp;

    transcript_matcher matcher(*grove, match_cfg);
    auto results = matcher.match_batch(clusters);

    // Step 3: Update grove with read support
    bool add_novel = args.count("add-novel") > 0;
    uint32_t min_novel_support = args["min-novel-support"].as<uint32_t>();

    if (add_novel) {
        logging::info("Updating grove with read support (adding novel transcripts with >= " +
                      std::to_string(min_novel_support) + " reads)...");
    } else {
        logging::info("Updating grove with read support (novel transcripts not added)...");
    }

    for (size_t i = 0; i < clusters.size(); ++i) {
        // Only add novel transcripts if flag is set and cluster has sufficient support
        bool should_add_novel = add_novel && clusters[i].members.size() >= min_novel_support;
        matcher.update_grove(clusters[i], results[i], should_add_novel);
    }

    const auto& match_stats = matcher.get_stats();
    logging::info("Matching complete:");
    logging::info("  Total matches: " + std::to_string(match_stats.total_matches));
    logging::info("  FSM (Full Splice Match): " + std::to_string(match_stats.fsm_matches));
    logging::info("  ISM (Incomplete Splice Match): " + std::to_string(match_stats.ism_matches));
    logging::info("  NIC (Novel In Catalog): " + std::to_string(match_stats.nic_matches));
    logging::info("  NNC (Novel Not in Catalog): " + std::to_string(match_stats.nnc_matches));
    logging::info("  Genic (intron/genomic): " + std::to_string(match_stats.genic_intron_matches + match_stats.genic_genomic_matches));
    logging::info("  Intergenic: " + std::to_string(match_stats.intergenic_matches));
    logging::info("  Ambiguous: " + std::to_string(match_stats.ambiguous_matches));
    logging::info("  Segments updated: " + std::to_string(match_stats.segments_updated));
    logging::info("  Segments created: " + std::to_string(match_stats.segments_created));

    // Step 4: Write output files
    auto out_dir = resolve_output_dir(args, input_path);
    std::string basename = std::filesystem::path(input_path).stem().string();
    std::string results_path = (out_dir / (basename + ".atroplex.tsv")).string();
    std::string summary_path = (out_dir / (basename + ".atroplex.summary.txt")).string();

    transcript_matcher::write_results(results_path, clusters, results);
    matcher.write_summary(summary_path);

    logging::info("Results written to: " + results_path);
}

} // namespace subcall
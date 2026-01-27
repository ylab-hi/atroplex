#include "base.hpp"

// standard
#include <fstream>
#include <iostream>

// class
#include "utility.hpp"
#include "builder.hpp"

// genogrove
#include <genogrove/io/gff_reader.hpp>
#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/filetype_detector.hpp>

namespace gdt = genogrove::data_type;
namespace ggs = genogrove::structure;
namespace gio = genogrove::io;

base::base(const cxxopts::ParseResult& params)
    : params{params}, ftype{gio::filetype::UNKNOWN}, compression{false} {

    // Check if a genogrove file is provided - if so load genogrove structure
    if(params.count("genogrove")) {
        std::string gg_path = params["genogrove"].as<std::string>();
        if(std::filesystem::exists(gg_path)) {
            logging::info("Loading existing genogrove from: " + gg_path);
            load_genogrove(gg_path);
        } else {
            logging::warning("Genogrove file not found: " + gg_path);
        }
    }
    if(params.count("build-from")) {
        auto build_files = params["build-from"].as<std::vector<std::string>>();
        build(build_files);
    }
}

base::~base() = default;

void base::validate(const cxxopts::ParseResult& args){
    // check if the input file exists
    if(!std::filesystem::exists(args["input"].as<std::string>())) {
        logging::error("Input file does not exist: " + args["input"].as<std::string>());
        throw std::runtime_error("Input file not found");
    }
}

void base::process() {
    logging::info("Starting atroplex pipeline...");
    start();

    // Process reads if input file is provided
    if (params.count("input")) {
        process_reads();
    }
}

void base::start() {
    logging::info("Initializing genogrove structure...");

    // If no grove loaded and build-from files are provided, create new grove
    if(!grove && params.count("build-from")) {
        auto build_files = params["build-from"].as<std::vector<std::string>>();
        build(build_files);
    }

    if(!grove) {
        logging::warning("No grove loaded or created");
    } else {
        logging::info("Grove ready with spatial index and graph structure");
    }
}

void base::build(const std::vector<std::string>& build_files) {
    int order = params["order"].as<int>();
    uint32_t threads = params["threads"].as<uint32_t>();

    logging::info("Creating grove with order: " + std::to_string(order));

    // Base creates and owns the grove
    grove = std::make_unique<grove_type>(order);

    // Builder populates the grove from annotation files
    builder::build_from_files(*grove, build_files, threads);
}

void base::load_genogrove(const std::string& gg_path) {
    // TODO: Implement genogrove deserialization
    logging::error("Genogrove loading not yet implemented");
}

void base::align_reads() {
    // TODO: Implement alignment for FASTQ input
    logging::warning("Alignment not yet implemented");
}

void base::process_reads() {
    if (!grove) {
        logging::warning("No grove available for read processing");
        return;
    }

    std::string input_path = params["input"].as<std::string>();

    // Detect input file type
    gio::filetype_detector detector;
    auto [detected_type, detected_compression] = detector.detect_filetype(input_path);

    if (detected_type != gio::filetype::BAM && detected_type != gio::filetype::SAM) {
        logging::info("Input is not BAM/SAM - skipping read processing");
        return;
    }

    logging::info("Processing reads from: " + input_path);

    // Configure bam_reader for primary alignments with quality filter
    gio::bam_reader_options opts = gio::bam_reader_options::primary_only();
    opts.min_mapq = 20;

    gio::bam_reader reader(input_path, opts);

    // Step 1: Cluster reads by splice junction signature
    logging::info("Clustering reads by splice junction signature...");
    read_clusterer::config cluster_cfg;
    cluster_cfg.junction_bin_size = 10;
    cluster_cfg.junction_tolerance = 5;
    cluster_cfg.min_mapq = opts.min_mapq;

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
    match_cfg.junction_tolerance = 5;
    match_cfg.min_junction_score = 0.8;
    match_cfg.min_overlap_bp = 50;

    transcript_matcher matcher(*grove, match_cfg);
    auto results = matcher.match_batch(clusters);

    // Step 3: Update grove with read support
    logging::info("Updating grove with read support...");
    for (size_t i = 0; i < clusters.size(); ++i) {
        matcher.update_grove(clusters[i], results[i]);
    }

    const auto& match_stats = matcher.get_stats();
    logging::info("Matching complete:");
    logging::info("  Total matches: " + std::to_string(match_stats.total_matches));
    logging::info("  Exact: " + std::to_string(match_stats.exact_matches));
    logging::info("  Compatible: " + std::to_string(match_stats.compatible_matches));
    logging::info("  Novel junction: " + std::to_string(match_stats.novel_junction_matches));
    logging::info("  Novel exon: " + std::to_string(match_stats.novel_exon_matches));
    logging::info("  Intergenic: " + std::to_string(match_stats.intergenic_matches));
    logging::info("  Ambiguous: " + std::to_string(match_stats.ambiguous_matches));
    logging::info("  Segments updated: " + std::to_string(match_stats.segments_updated));
    logging::info("  Segments created: " + std::to_string(match_stats.segments_created));

    // Step 4: Write output files
    std::filesystem::path input_file(input_path);
    std::string base_name = input_file.stem().string();
    std::string output_dir = input_file.parent_path().string();
    if (output_dir.empty()) output_dir = ".";

    std::string results_path = output_dir + "/" + base_name + ".atroplex.tsv";
    std::string summary_path = output_dir + "/" + base_name + ".atroplex.summary.txt";

    transcript_matcher::write_results(results_path, clusters, results);
    matcher.write_summary(summary_path);

    logging::info("Read processing complete");
}
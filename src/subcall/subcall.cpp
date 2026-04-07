/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/subcall.hpp"

#include <filesystem>
#include <fstream>

#include "utility.hpp"
#include "builder.hpp"
#include "build_gff.hpp"
#include "sample_manifest.hpp"

namespace subcall {

void subcall::add_common_options(cxxopts::Options& options) {
    options.add_options("Common")
        ("o,output-dir", "Output directory for results (default: input file directory)",
            cxxopts::value<std::string>())
        ("p,prefix", "Output file prefix (default: derived from manifest or first input file)",
            cxxopts::value<std::string>())
        ("t,threads", "Number of threads (0 = auto-detect)",
            cxxopts::value<uint32_t>()->default_value("1"))
        ("progress", "Show progress output")
        ("h,help", "Show help message")
        ;

    options.add_options("Grove")
        ("g,genogrove", "Pre-built genogrove index (.ggx)",
            cxxopts::value<std::string>())
        ("m,manifest", "Sample manifest file (TSV with metadata)",
            cxxopts::value<std::string>())
        ("b,build-from", "Build grove from file(s) without metadata (GFF/GTF)",
            cxxopts::value<std::vector<std::string>>())
        ("k,order", "Genogrove tree order",
            cxxopts::value<int>()->default_value("3"))
        ("min-expression", "Minimum expression value to include a transcript (filters low-count novels)",
            cxxopts::value<float>()->default_value("-1"))
        ("no-absorb", "Disable ISM (Incomplete Splice Match) segment absorption into longer parents")
        ("fuzzy-tolerance", "Max bp difference for fuzzy exon boundary matching during absorption (0 = exact only)",
            cxxopts::value<size_t>()->default_value("5"))
        ("min-replicates", "Merge biological replicates within groups; require features in >= N replicates (0 = no merge)",
            cxxopts::value<int>()->default_value("0"))
        ;
}

void subcall::apply_common_options(const cxxopts::ParseResult& args) {
    if (args.count("progress")) {
        logging::set_progress_enabled(true);
    }
}

std::filesystem::path subcall::resolve_output_dir(const cxxopts::ParseResult& args,
                                                   const std::string& fallback_input_path) const {
    std::filesystem::path dir;

    if (args.count("output-dir")) {
        dir = args["output-dir"].as<std::string>();
    } else if (!fallback_input_path.empty()) {
        dir = std::filesystem::path(fallback_input_path).parent_path();
    }

    if (dir.empty()) {
        dir = std::filesystem::current_path();
    }

    std::filesystem::create_directories(dir);
    return dir;
}

void subcall::run(const cxxopts::ParseResult& args) {
    validate(args);
    apply_common_options(args);
    setup_grove(args);
    write_build_summary(args);
    execute(args);
}

void subcall::setup_grove(const cxxopts::ParseResult& args) {
    if (args.count("genogrove")) {
        std::string gg_path = args["genogrove"].as<std::string>();
        logging::info("Loading grove from: " + gg_path);
        load_grove(gg_path);
        return;  // Grove loaded from file, skip building
    }

    int order = args["order"].as<int>();
    uint32_t threads = args["threads"].as<uint32_t>();

    // Collect samples from both sources
    std::vector<sample_info> all_samples;

    // From manifest (metadata in TSV)
    if (args.count("manifest")) {
        std::string manifest_path = args["manifest"].as<std::string>();
        logging::info("Loading manifest: " + manifest_path);
        sample_manifest manifest(manifest_path);
        logging::info("Found " + std::to_string(manifest.size()) + " sample(s) in manifest");

        for (const auto& info : manifest) {
            all_samples.push_back(info);
        }
    }

    // From build-from (metadata parsed from GFF headers)
    if (args.count("build-from")) {
        auto build_files = args["build-from"].as<std::vector<std::string>>();
        logging::info("Parsing headers from " + std::to_string(build_files.size()) + " file(s)");

        for (const auto& filepath : build_files) {
            sample_info info = build_gff::parse_header(filepath);
            all_samples.push_back(std::move(info));
        }
    }

    // Build grove if we have samples
    if (!all_samples.empty()) {
        float min_expr = args["min-expression"].as<float>();
        bool absorb = !args.count("no-absorb");
        size_t fuzzy_tol = args["fuzzy-tolerance"].as<size_t>();
        int min_reps = args["min-replicates"].as<int>();
        logging::info("Creating grove with order: " + std::to_string(order));
        if (min_expr >= 0) {
            logging::info("Filtering transcripts with expression < " + std::to_string(min_expr));
        }
        if (!absorb) {
            logging::info("ISM segment absorption disabled");
        } else if (fuzzy_tol > 0) {
            logging::info("Fuzzy absorption tolerance: " + std::to_string(fuzzy_tol) + "bp");
        }
        if (min_reps > 0) {
            logging::info("Replicate merging enabled: min_replicates = " + std::to_string(min_reps));
        }
        grove = std::make_unique<grove_type>(order);
        auto build_start = std::chrono::steady_clock::now();
        build_stats = builder::build_from_samples(*grove, all_samples, threads, min_expr, absorb, min_reps, fuzzy_tol, &exon_caches_);
        auto build_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - build_start).count();
        build_stats->build_time_seconds = build_elapsed;
        logging::info("Grove ready with spatial index and graph structure");
    }
}

std::string subcall::resolve_prefix(const cxxopts::ParseResult& args) const {
    if (args.count("prefix")) {
        return args["prefix"].as<std::string>();
    }

    // Derive from manifest filename
    if (args.count("manifest")) {
        return std::filesystem::path(args["manifest"].as<std::string>()).stem().string();
    }

    // Derive from first build-from file
    if (args.count("build-from")) {
        return std::filesystem::path(
            args["build-from"].as<std::vector<std::string>>().front()).stem().string();
    }

    // Derive from genogrove index (strip .ggx extension)
    if (args.count("genogrove")) {
        return std::filesystem::path(args["genogrove"].as<std::string>()).stem().string();
    }

    return "atroplex";
}

void subcall::load_grove(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open()) {
        throw std::runtime_error("Cannot open index file: " + path);
    }

    // Read magic bytes + version
    char magic[4];
    ifs.read(magic, 4);
    if (std::string(magic, 4) != "AGRX") {
        throw std::runtime_error("Invalid .ggx file (bad magic): " + path);
    }
    uint16_t version;
    ifs.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != 1) {
        throw std::runtime_error("Unsupported .ggx version: " + std::to_string(version));
    }

    // Deserialize registries (order must match save_grove)
    gene_registry::instance().deserialize_into(ifs);
    source_registry::instance().deserialize_into(ifs);
    transcript_registry::instance().deserialize_into(ifs);
    (void)sample_registry::deserialize(ifs);

    // Deserialize grove (handles its own zlib decompression)
    grove = std::make_unique<grove_type>(grove_type::deserialize(ifs));

    logging::info("Loaded index from: " + path);
    logging::info("  Genes: " + std::to_string(gene_registry::instance().size()) +
                  ", Sources: " + std::to_string(source_registry::instance().size()) +
                  ", Transcripts: " + std::to_string(transcript_registry::instance().size()) +
                  ", Samples: " + std::to_string(sample_registry::instance().size()));
}

void subcall::build_grove(const std::vector<std::string>& files, int order, uint32_t threads) {
    logging::info("Creating grove with order: " + std::to_string(order));
    grove = std::make_unique<grove_type>(order);
    auto build_start = std::chrono::steady_clock::now();
    build_stats = builder::build_from_files(*grove, files, threads);
    auto build_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - build_start).count();
    build_stats->build_time_seconds = build_elapsed;
    logging::info("Grove ready with spatial index and graph structure");
}

void subcall::save_grove(const std::string& path) {
    if (!grove) {
        logging::warning("No grove to save");
        return;
    }

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.is_open()) {
        throw std::runtime_error("Cannot create index file: " + path);
    }

    // Write magic bytes + version
    ofs.write("AGRX", 4);  // Atroplex GRove indeX
    uint16_t version = 1;
    ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));

    // Serialize registries (order must match load_grove)
    gene_registry::instance().serialize(ofs);
    source_registry::instance().serialize(ofs);
    transcript_registry::instance().serialize(ofs);
    sample_registry::instance().serialize(ofs);

    // Serialize grove (handles its own zlib compression)
    grove->serialize(ofs);

    logging::info("Index saved to: " + path + " (" +
                  std::to_string(std::filesystem::file_size(path) / 1024) + " KB)");
}

void subcall::write_build_summary(const cxxopts::ParseResult& args) {
    if (!build_stats) return;

    std::string prefix = resolve_prefix(args);
    std::string fallback;
    if (args.count("manifest")) {
        fallback = args["manifest"].as<std::string>();
    } else if (args.count("build-from")) {
        fallback = args["build-from"].as<std::vector<std::string>>().front();
    }

    auto out_dir = resolve_output_dir(args, fallback);
    std::string summary_path = (out_dir / (prefix + ".ggx.summary")).string();
    build_stats->write_summary(summary_path);
}

} // namespace subcall
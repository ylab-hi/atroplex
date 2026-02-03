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

#include "utility.hpp"
#include "builder.hpp"
#include "index_stats.hpp"

namespace subcall {

void subcall::add_common_options(cxxopts::Options& options) {
    options.add_options("Common")
        ("o,output-dir", "Output directory for results (default: input file directory)",
            cxxopts::value<std::string>())
        ("t,threads", "Number of threads (0 = auto-detect)",
            cxxopts::value<uint32_t>()->default_value("1"))
        ("s,stats", "Write index statistics to output directory")
        ("progress", "Show progress output")
        ("h,help", "Show help message")
        ;

    options.add_options("Grove")
        ("g,genogrove", "Pre-built genogrove index (.gg)",
            cxxopts::value<std::string>())
        ("b,build-from", "Build grove from annotation file(s) (GFF/GTF)",
            cxxopts::value<std::vector<std::string>>())
        ("k,order", "Genogrove tree order",
            cxxopts::value<int>()->default_value("3"))
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
    execute(args);
    write_index_stats(args);
}

void subcall::setup_grove(const cxxopts::ParseResult& args) {
    if (args.count("genogrove")) {
        std::string gg_path = args["genogrove"].as<std::string>();
        logging::info("Loading grove from: " + gg_path);
        load_grove(gg_path);
    }

    if (args.count("build-from")) {
        auto build_files = args["build-from"].as<std::vector<std::string>>();
        int order = args["order"].as<int>();
        uint32_t threads = args["threads"].as<uint32_t>();
        build_grove(build_files, order, threads);
    }
}

void subcall::load_grove(const std::string& path) {
    // TODO: Implement genogrove deserialization
    logging::error("Genogrove loading not yet implemented: " + path);
    throw std::runtime_error("Genogrove deserialization not implemented");
}

void subcall::build_grove(const std::vector<std::string>& files, int order, uint32_t threads) {
    logging::info("Creating grove with order: " + std::to_string(order));
    grove = std::make_unique<grove_type>(order);
    builder::build_from_files(*grove, files, threads);
    logging::info("Grove ready with spatial index and graph structure");
}

void subcall::save_grove(const std::string& path) {
    // TODO: Implement genogrove serialization
    logging::error("Genogrove saving not yet implemented: " + path);
    throw std::runtime_error("Genogrove serialization not implemented");
}

void subcall::write_index_stats(const cxxopts::ParseResult& args) {
    if (!args.count("stats") || !grove) return;

    logging::info("Collecting index statistics...");
    auto stats = index_stats::collect(*grove);

    // Determine output path: use a fallback from available args
    std::string fallback;
    if (args.count("input")) {
        fallback = args["input"].as<std::string>();
    } else if (args.count("build-from")) {
        fallback = args["build-from"].as<std::vector<std::string>>().front();
    } else if (args.count("manifest")) {
        fallback = args["manifest"].as<std::string>();
    }

    auto out_dir = resolve_output_dir(args, fallback);
    std::string basename = fallback.empty()
        ? "atroplex"
        : std::filesystem::path(fallback).stem().string();

    std::string stats_path = (out_dir / (basename + ".index_stats.txt")).string();
    stats.write(stats_path);
}

} // namespace subcall
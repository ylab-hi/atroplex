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

#include "utility.hpp"
#include "builder.hpp"

namespace subcall {

void subcall::add_common_options(cxxopts::Options& options) {
    options.add_options("Common")
        ("t,threads", "Number of threads (0 = auto-detect)",
            cxxopts::value<uint32_t>()->default_value("1"))
        ("progress", "Show progress output")
        ("h,help", "Show help message")
        ;
}

void subcall::apply_common_options(const cxxopts::ParseResult& args) {
    if (args.count("progress")) {
        logging::set_progress_enabled(true);
    }
    // Thread count is accessed directly by subcommands when needed
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

} // namespace subcall
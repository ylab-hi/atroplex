/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/build.hpp"

#include <filesystem>

#include "utility.hpp"
#include "builder.hpp"

namespace subcall {

cxxopts::Options build::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex build",
        "Build genogrove index from annotation files");

    options.add_options("Input/Output")
        ("i,input", "Input annotation file(s) (GFF/GTF) - can specify multiple",
            cxxopts::value<std::vector<std::string>>())
        ("o,output", "Output genogrove index file (.gg)",
            cxxopts::value<std::string>())
        ;

    options.add_options("Grove")
        ("k,order", "Genogrove tree order (higher = more memory, faster queries)",
            cxxopts::value<int>()->default_value("3"))
        ;

    add_common_options(options);

    return options;
}

void build::validate(const cxxopts::ParseResult& args) {
    if (!args.count("input")) {
        throw std::runtime_error("No input files specified. Use -i/--input");
    }

    if (!args.count("output")) {
        throw std::runtime_error("No output path specified. Use -o/--output");
    }

    // Check input files exist
    auto files = args["input"].as<std::vector<std::string>>();
    for (const auto& f : files) {
        if (!std::filesystem::exists(f)) {
            throw std::runtime_error("Input file not found: " + f);
        }
    }
}

void build::execute(const cxxopts::ParseResult& args) {
    apply_common_options(args);

    auto input_files = args["input"].as<std::vector<std::string>>();
    std::string output_path = args["output"].as<std::string>();
    int order = args["order"].as<int>();
    uint32_t threads = args["threads"].as<uint32_t>();

    logging::info("Building genogrove index from " +
                  std::to_string(input_files.size()) + " file(s)");

    // Build grove
    build_grove(input_files, order, threads);

    // Save to output
    logging::info("Saving grove to: " + output_path);
    save_grove(output_path);

    logging::info("Build complete");
}

} // namespace subcall
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

namespace subcall {

cxxopts::Options build::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex build",
        "Build genogrove index from annotation files");

    add_common_options(options);

    return options;
}

void build::validate(const cxxopts::ParseResult& args) {
    if (!args.count("build-from") && !args.count("manifest")) {
        throw std::runtime_error("No input specified. Use -m/--manifest or -b/--build-from");
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

void build::execute(const cxxopts::ParseResult& args) {
    // Determine output path from first input source
    std::string first_input;
    if (args.count("manifest")) {
        first_input = args["manifest"].as<std::string>();
    } else if (args.count("build-from")) {
        first_input = args["build-from"].as<std::vector<std::string>>().front();
    }

    auto out_dir = resolve_output_dir(args, first_input);
    std::string basename = std::filesystem::path(first_input).stem().string();
    std::string output_path = (out_dir / (basename + ".gg")).string();

    // Grove already built by setup_grove()
    logging::info("Saving grove to: " + output_path);
    save_grove(output_path);

    logging::info("Build complete");
}

} // namespace subcall
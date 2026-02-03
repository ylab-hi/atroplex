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
    if (!args.count("build-from")) {
        throw std::runtime_error("No annotation files specified. Use -b/--build-from");
    }

    // Check build-from files exist
    auto files = args["build-from"].as<std::vector<std::string>>();
    for (const auto& f : files) {
        if (!std::filesystem::exists(f)) {
            throw std::runtime_error("Annotation file not found: " + f);
        }
    }
}

void build::execute(const cxxopts::ParseResult& args) {
    auto build_files = args["build-from"].as<std::vector<std::string>>();
    auto out_dir = resolve_output_dir(args, build_files.front());
    std::string basename = std::filesystem::path(build_files.front()).stem().string();
    std::string output_path = (out_dir / (basename + ".gg")).string();

    // Grove already built by setup_grove()
    logging::info("Saving grove to: " + output_path);
    save_grove(output_path);

    logging::info("Build complete");
}

} // namespace subcall
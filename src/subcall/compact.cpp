/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/compact.hpp"

#include <filesystem>
#include <stdexcept>

#include "builder.hpp"
#include "utility.hpp"

namespace subcall {

namespace {

/// Locate the single `.ggx` file inside `grove_dir`. Mirrors the scan
/// performed by subcall::setup_grove so validate() can produce the same
/// error messages before the heavy load step runs.
std::filesystem::path find_ggx(const std::filesystem::path& grove_dir) {
    if (!std::filesystem::is_directory(grove_dir)) {
        throw std::runtime_error(
            "-g/--genogrove expects a directory, not a file: " + grove_dir.string());
    }
    std::filesystem::path gg_path;
    for (const auto& entry : std::filesystem::directory_iterator(grove_dir)) {
        if (entry.path().extension() == ".ggx") {
            if (!gg_path.empty()) {
                throw std::runtime_error(
                    "Multiple .ggx files in directory: " + grove_dir.string());
            }
            gg_path = entry.path();
        }
    }
    if (gg_path.empty()) {
        throw std::runtime_error("No .ggx file found in directory: " + grove_dir.string());
    }
    return gg_path;
}

} // namespace

cxxopts::Options compact::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex compact",
        "Compact a genogrove index by physically removing absorbed segments");

    add_common_options(options);

    options.add_options("Compact")
        ("no-qtx", "Allow compact to proceed without a .qtx sidecar alongside "
            "the input .ggx. Default: off — compact fails fast if no .qtx is "
            "found, since a structure-only output would silently desynchronise "
            "from any expression sidecar that does exist.")
        ;

    return options;
}

void compact::validate(const cxxopts::ParseResult& args) {
    if (!args.count("genogrove")) {
        throw std::runtime_error(
            "compact requires -g/--genogrove pointing at an existing index directory");
    }
    if (!args.count("output-dir")) {
        throw std::runtime_error(
            "compact requires -o/--output-dir for the compacted index. "
            "Pick a directory other than the input to keep the operation reversible.");
    }

    std::filesystem::path grove_dir = args["genogrove"].as<std::string>();
    std::filesystem::path gg_path = find_ggx(grove_dir);

    // Compact never builds — reject the build-side options so users don't
    // assume `-m` or `-b` mean "rebuild then compact".
    if (args.count("manifest") || args.count("build-from")) {
        throw std::runtime_error(
            "compact operates on an existing index only; remove -m/--manifest "
            "and -b/--build-from");
    }

    bool no_qtx = args.count("no-qtx") > 0;
    std::filesystem::path qtx_path = gg_path;
    qtx_path.replace_extension(".qtx");
    bool qtx_exists = std::filesystem::exists(qtx_path);

    if (!qtx_exists && !no_qtx) {
        throw std::runtime_error(
            "No .qtx sidecar found alongside " + gg_path.string()
            + " (expected " + qtx_path.string() + "). "
            "Pass --no-qtx if the index is intentionally structure-only.");
    }
    if (qtx_exists && no_qtx) {
        logging::warning("--no-qtx set but a .qtx sidecar exists at "
            + qtx_path.string()
            + " — it will NOT be copied to the output directory. "
            "The compacted index will be structure-only.");
    }

    std::filesystem::path out_dir = args["output-dir"].as<std::string>();
    if (std::filesystem::weakly_canonical(out_dir) ==
        std::filesystem::weakly_canonical(grove_dir)) {
        throw std::runtime_error(
            "compact refuses to write to the input directory. "
            "Choose a different -o/--output-dir.");
    }
}

void compact::execute(const cxxopts::ParseResult& args) {
    std::filesystem::path grove_dir = args["genogrove"].as<std::string>();
    std::filesystem::path src_ggx = find_ggx(grove_dir);
    std::string stem = src_ggx.stem().string();

    bool no_qtx = args.count("no-qtx") > 0;

    auto out_dir = resolve_output_dir(args, src_ggx.string());
    std::filesystem::path dst_ggx = out_dir / (stem + ".ggx");
    std::filesystem::path dst_qtx = out_dir / (stem + ".qtx");
    std::filesystem::path dst_summary = out_dir / (stem + ".ggx.summary");

    // setup_grove already loaded the grove (and tried the qtx reader) for
    // us via run(). Walk the live grove and unlink every absorbed segment.
    chromosome_segment_caches empty_caches;
    size_t removed = builder::remove_tombstones(*grove, empty_caches, /*physical=*/true);

    if (removed == 0) {
        logging::info("Index already compact — no absorbed segments to remove");
    } else {
        logging::info("Compacted index: " + std::to_string(removed)
            + " absorbed segment(s) physically removed");
    }

    save_grove(dst_ggx.string());

    // The .qtx written by build is already remapped against live segments
    // (merge_to_qtx consumes the tombstone remap at build time), so a
    // straight copy keeps it in sync with the compacted .ggx. Same for
    // the build summary, which is timestamped with the original build's
    // counters (input_transcripts, scaffold_filtered, etc.) that are not
    // recoverable from the loaded grove.
    if (!no_qtx) {
        std::filesystem::path src_qtx = src_ggx;
        src_qtx.replace_extension(".qtx");
        if (std::filesystem::exists(src_qtx)) {
            std::filesystem::copy_file(src_qtx, dst_qtx,
                std::filesystem::copy_options::overwrite_existing);
            logging::info("Copied .qtx sidecar to: " + dst_qtx.string());
        }
    }

    std::filesystem::path src_summary = src_ggx;
    src_summary.replace_extension(".ggx.summary");
    if (std::filesystem::exists(src_summary)) {
        std::filesystem::copy_file(src_summary, dst_summary,
            std::filesystem::copy_options::overwrite_existing);
        logging::info("Copied .ggx.summary to: " + dst_summary.string());
    }

    logging::info("Compact complete: " + dst_ggx.string());
}

} // namespace subcall
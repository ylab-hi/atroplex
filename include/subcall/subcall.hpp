/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_HPP
#define ATROPLEX_SUBCALL_HPP

#include <filesystem>
#include <memory>
#include <string>

#include <cxxopts.hpp>

#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>

#include "genomic_feature.hpp"

namespace subcall {

/**
 * Abstract base class for all atroplex subcommands.
 */
class subcall {
public:
    virtual ~subcall() = default;

    /**
     * Parse command-line arguments and return the options object.
     * The subclass defines its own options here.
     * Should call add_common_options() to include shared options.
     */
    virtual cxxopts::Options parse_args(int argc, char** argv) = 0;

    /**
     * Validate parsed arguments. Throws on invalid input.
     */
    virtual void validate(const cxxopts::ParseResult& args) = 0;

    /**
     * Execute the subcommand (called after grove setup).
     */
    virtual void execute(const cxxopts::ParseResult& args) = 0;

    /**
     * Template method: validate → apply_common_options → setup_grove → execute.
     */
    void run(const cxxopts::ParseResult& args);

    /**
     * Add common options shared across all subcommands.
     * Call this in parse_args() implementations.
     */
    static void add_common_options(cxxopts::Options& options);

    /**
     * Get the subcommand name (for help text).
     */
    virtual std::string name() const = 0;

    /**
     * Get brief description (for help text).
     */
    virtual std::string description() const = 0;

protected:
    std::unique_ptr<grove_type> grove;
    std::filesystem::path output_dir;

    /**
     * Resolve the output directory from --output-dir or the parent of a fallback path.
     * Creates the directory if it doesn't exist.
     */
    std::filesystem::path resolve_output_dir(const cxxopts::ParseResult& args,
                                             const std::string& fallback_input_path) const;

    /**
     * Load or build grove from common options (--genogrove, --build-from, --order).
     */
    void setup_grove(const cxxopts::ParseResult& args);

    /**
     * Load grove from .gg file
     */
    void load_grove(const std::string& path);

    /**
     * Build grove from annotation files
     */
    void build_grove(const std::vector<std::string>& files, int order, uint32_t threads);

    /**
     * Save grove to .gg file
     */
    void save_grove(const std::string& path);

private:
    /**
     * Apply common options (threads, progress, etc.)
     */
    static void apply_common_options(const cxxopts::ParseResult& args);

    /**
     * Collect and write index statistics if --stats is set and grove exists.
     */
    void write_index_stats(const cxxopts::ParseResult& args);
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_HPP
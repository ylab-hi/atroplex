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
     * Execute the subcommand.
     */
    virtual void execute(const cxxopts::ParseResult& args) = 0;

    /**
     * Add common options shared across all subcommands.
     * Call this in parse_args() implementations.
     */
    static void add_common_options(cxxopts::Options& options);

    /**
     * Apply common options (threads, progress, etc.)
     * Call this at the start of execute() implementations.
     */
    static void apply_common_options(const cxxopts::ParseResult& args);

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
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_HPP
/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_BASE_HPP
#define ATROPLEX_BASE_HPP

// standard
#include <filesystem>
#include <memory>

// cxxopts
#include <cxxopts.hpp>

// class
#include "filetype_detector.hpp"

/**
 * base handles the initial file type detection and routing
 * to appropriate processing stages (alignment for FASTQ, structure creation, detection)
 */

class base {
public:
    explicit base(const cxxopts::ParseResult& params);
    ~base();

    // main processing pipeline
    void process();

    // Getter

private:
    cxxopts::ParseResult params;
    std::filesystem::path input_path;
    filetype detected_filetype;
};

#endif //ATROPLEX_BASE_HPP

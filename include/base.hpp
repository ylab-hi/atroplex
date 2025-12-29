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
#include <string>

// cxxopts
#include <cxxopts.hpp>

// genogrove
#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/gff_reader.hpp>
#include <genogrove/io/filetype_detector.hpp>

#include "genomic_feature.hpp"

// class
#include "genogrove_builder.hpp"

namespace gdt = genogrove::data_type;
namespace gst = genogrove::structure;
namespace gio = genogrove::io;

/**
 * base handles the initial file type detection and routing
 * to appropriate processing stages (alignment for FASTQ, structure creation, detection)
 */

class base {
public:
    explicit base(const cxxopts::ParseResult& params);
    ~base();

    void validate(const cxxopts::ParseResult& args);

    // main processing pipeline
    void process();

    // Getters
    genomic_structures* get_structures() const { return structures; }

private:
    cxxopts::ParseResult params;
    gst::grove<gdt::interval, genomic_feature> grove;

    gio::filetype ftype;
    gio::compression_type compression;
    bool gzipped;

    // Genomic data structures (grove + transcript graph)
    genomic_structures* structures = nullptr;

    // steps for the pipeline
    void start();
    void detect_input_filetype();
    void align_reads();
    void create_genogrove(const std::vector<std::string>& build_files);
    void load_genogrove(const std::string& gg_path);

};

#endif //ATROPLEX_BASE_HPP

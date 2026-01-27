/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
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
#include <genogrove/io/bam_reader.hpp>

#include "genomic_feature.hpp"
#include "read_cluster.hpp"
#include "transcript_matcher.hpp"

// class
#include "builder.hpp"

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
    grove_type* get_grove() const { return grove.get(); }

private:
    cxxopts::ParseResult params;
    std::unique_ptr<grove_type> grove;

    gio::filetype ftype;
    gio::compression_type compression;
    bool gzipped;

    // steps for the pipeline
    void start();
    void detect_input_filetype();
    void align_reads();
    void build(const std::vector<std::string>& build_files);
    void load_genogrove(const std::string& gg_path);
    void process_reads();

};

#endif //ATROPLEX_BASE_HPP

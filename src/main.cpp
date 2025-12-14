/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

// standard
#include <iostream>
#include <string>

// cxxopts
#include <cxxopts.hpp>

// class
#include "config.hpp"
#include "utility.hpp"
#include "base.hpp"

void showVersion(std::ostream& _str) {
    _str << "morrigan v" << atroplex_VERSION_MAJOR;
    _str << "." << atroplex_VERSION_MINOR << ".";
    _str << atroplex_VERSION_PATCH << " - ";
    _str << "Detect full-length transcripts from long-read ";
    _str << "sequencing data";
    _str << std::endl;
}

int main(int argc, char** argv) {
    try {
        cxxopts::Options options("morrigan",
            "Detect full-length transcripts from long-read sequencing data");

        options.add_options("General")
            ("i,input", "Input file (can be BAM or fastq.(gz))",
                cxxopts::value<std::string>())
            ("g,genogrove", "Genogrove file in .gg format (optional)",
                cxxopts::value<std::string>())
            ("a,annotation", "Annotation file in GFF/GTF format (optional)",
                cxxopts::value<std::string>())
            ("k,order", "Order of genogrove (only used when genogrove is newly created)",
                cxxopts::value<int>()->default_value("3"))
            ;

        options.add_options("Other")
            ("h,help", "Print help message")
            ("v,version", "Print version number")
            ;

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        if (result.count("version")) {
            showVersion(std::cout);
            return 0;
        }

        // Check if input is provided
        if (!result.count("input")) {
            logging::error("Please specify an input file with -i/--input");
            std::cout << options.help() << std::endl;
            return 1;
        }

        // Create base processor which handles file type detection and routing
        base processor(result);
        processor.process();

    } catch(const cxxopts::exceptions::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}


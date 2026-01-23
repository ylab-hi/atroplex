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
    _str << "atroplex v" << atroplex_VERSION_MAJOR;
    _str << "." << atroplex_VERSION_MINOR << ".";
    _str << atroplex_VERSION_PATCH << " - ";
    _str << "Detect full-length transcripts from long-read ";
    _str << "sequencing data";
    _str << std::endl;
}

int main(int argc, char** argv) {
    try {
        cxxopts::Options options("atroplex",
            "Detect full-length transcripts from long-read sequencing data");

        options.add_options("General")
            ("i,input", "Input file (can be BAM or fastq.(gz))",
                cxxopts::value<std::string>())
            ("g,genogrove", "Genogrove file in .gg format (optional)",
                cxxopts::value<std::string>())
            ("b,build-from", "Build the genogrove structure from file(s) (can be specified multiple times)",
                cxxopts::value<std::vector<std::string>>())
            ("o,output", "Output file path for built genogrove index (.gg)",
                cxxopts::value<std::string>())
            ("k,order", "Order of genogrove (only used when genogrove is newly created)",
                cxxopts::value<int>()->default_value("3"))
            ("t,threads", "Number of threads for parallel processing (0 = auto-detect)",
                cxxopts::value<uint32_t>()->default_value("1"))
            ("build-only", "Only build the genogrove index, skip analysis")
            ;

        options.add_options("Pan-transcriptome")
            ("sample-manifest", "TSV file with sample metadata (sample_id, file, tissue, condition, ...)",
                cxxopts::value<std::string>())
            ("sample-id", "Sample identifier for single-sample build",
                cxxopts::value<std::string>())
            ("tissue", "Tissue type for single-sample build",
                cxxopts::value<std::string>())
            ("condition", "Condition for single-sample build",
                cxxopts::value<std::string>())
            ;

        options.add_options("Other")
            ("h,help", "Print help message")
            ("v,version", "Print version number")
            ("progress", "Show progress with carriage return updates (for interactive use)")
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

        // Enable progress output if requested
        if (result.count("progress")) {
            logging::set_progress_enabled(true);
        }

        // Validate arguments based on mode
        bool build_only = result.count("build-only") > 0;
        bool has_manifest = result.count("sample-manifest") > 0;
        bool has_build_files = result.count("build-from") > 0;

        // Build-only mode requires either manifest or build-from files
        if (build_only) {
            if (!has_manifest && !has_build_files) {
                logging::error("Build-only mode requires --sample-manifest or --build-from files");
                std::cout << options.help() << std::endl;
                return 1;
            }
        } else {
            // Analysis mode requires input file (for debugging disables)
            // if (!result.count("input")) {
            //     logging::error("Please specify an input file with -i/--input");
            //     std::cout << options.help() << std::endl;
            //     return 1;
            // }
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


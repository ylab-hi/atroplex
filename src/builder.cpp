/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

// standard
#include <filesystem>

// class
#include "utility.hpp"
#include "builder.hpp"
#include "build_gff.hpp"

void builder::build_from_files(grove_type& grove, const std::vector<std::string>& files) {
    if (files.empty()) {
        logging::warning("No files provided to build genogrove");
        return;
    }

    logging::info("Populating grove from " + std::to_string(files.size()) + " file(s)");

    // Process each file
    for (const auto& filepath : files) {
        if (!std::filesystem::exists(filepath)) {
            logging::error("File not found: " + filepath);
            continue;
        }

        // Detect file type
        gio::filetype_detector detector;
        auto [ftype, is_gzipped] = detector.detect_filetype(filepath);

        // Dispatch to appropriate builder based on file type
        if (ftype == gio::filetype::GFF || ftype == gio::filetype::GTF) {
            logging::info("Processing GFF/GTF file: " + filepath);
            build_gff::build(grove, filepath);
        } else {
            logging::warning("Unsupported file type for: " + filepath);
        }

        // Future format support:
        // else if (ftype == gio::filetype::BED) {
        //     build_bed::build(grove, filepath);
        // }
        // else if (ftype == gio::filetype::BAM) {
        //     build_bam::build(grove, filepath);
        // }
    }

    logging::info("Grove construction complete");
}
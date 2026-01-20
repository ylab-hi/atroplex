/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_UTILITY_HPP
#define ATROPLEX_UTILITY_HPP

// standard
#include <chrono>
#include <string>

namespace logging {
    void info(const std::string& message);
    void warning(const std::string& message);
    void error(const std::string& message);

    /**
     * Display progress with carriage return (overwrites current line)
     * Shows line count with comma formatting and lines/sec rate
     * @param lines Number of lines read
     * @param prefix Prefix message (e.g., "Processing gencode.gtf")
     */
    void progress(size_t lines, const std::string& prefix = "Processing");

    /**
     * Clear the progress line and print final summary
     * @param segments Total segments created
     * @param prefix Prefix message
     */
    void progress_done(size_t segments, const std::string& prefix = "Processed");

    /**
     * Reset the progress timer (call before starting a new progress sequence)
     */
    void progress_start();
}

#endif //ATROPLEX_UTILITY_HPP

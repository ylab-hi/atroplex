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
}

#endif //ATROPLEX_UTILITY_HPP

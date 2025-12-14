/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_ALIGNMENT_HPP
#define ATROPLEX_ALIGNMENT_HPP

// cxxopts
#include <cxxopts.hpp>

class alignment {
    public:
        alignment(const cxxopts::ParseResult& params);

    private:
        cxxopts::ParseResult params;
};

#endif //ATROPLEX_ALIGNMENT_HPP

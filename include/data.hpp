/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_DATA_HPP
#define ATROPLEX_DATA_HPP

#include <string>

struct exon {
    std::string exon_id;
    std::string gene_id;

    std::string gene_name;
    std::string gene_type;
};

struct segment {
    std::string segment_id;
};

struct transcript {
    std::string transcript_id;
};

#endif //ATROPLEX_DATA_HPP

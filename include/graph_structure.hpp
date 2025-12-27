/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_GRAPH_STRUCTURE_HPP
#define ATROPLEX_GRAPH_STRUCTURE_HPP

// standard

// class
#include "segment_identifier.hpp"

// genogrove
#include <genogrove/data_type/genomic_coordinate.hpp>

namespace gdt = genogrove::data_type;


struct segment {
    segment_identifier id;
    gdt::genomic_coordinate genomic_coordinate;

};

struct exon {
    int num;
};



#endif //ATROPLEX_GRAPH_STRUCTURE_HPP

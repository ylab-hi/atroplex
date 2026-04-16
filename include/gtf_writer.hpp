/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_GTF_WRITER_HPP
#define ATROPLEX_GTF_WRITER_HPP

#include <ostream>
#include <string>
#include <vector>

#include "sample_info.hpp"

namespace gtf_writer {

struct gtf_gene {
    std::string seqid;
    std::string source;
    char strand;
    size_t start;
    size_t end;
    std::string gene_id;
    std::string gene_name;
    std::string gene_biotype;
};

struct gtf_transcript {
    std::string seqid;
    std::string source;
    char strand;
    size_t start;
    size_t end;
    std::string gene_id;
    std::string gene_name;
    std::string transcript_id;
    std::string transcript_biotype;
    int exon_count;
};

struct gtf_exon {
    std::string seqid;
    std::string source;
    char strand;
    size_t start;
    size_t end;
    std::string gene_id;
    std::string gene_name;
    std::string transcript_id;
    int exon_number;
};

void write_header(std::ostream& os, const sample_info& info);
void write_gene(std::ostream& os, const gtf_gene& g);
void write_transcript(std::ostream& os, const gtf_transcript& t);
void write_exon(std::ostream& os, const gtf_exon& e);

} // namespace gtf_writer

#endif // ATROPLEX_GTF_WRITER_HPP
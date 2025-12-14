/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_DATA_ENTRIES_HPP
#define ATROPLEX_DATA_ENTRIES_HPP

#include <string>
#include <optional>
#include <unordered_map>

// represents a single GFF/GTF annotation entry
struct gff_entry {
    std::string seqid;
    std::string source;
    std::string type;
    int start;
    int end;
    float score;
    char strand;
    int8_t phase;
    std::unordered_map<std::string, std::string> attributes;

    gff_entry() = default;
    gff_entry(std::string seqid, std::string source, std::string type, int start, int end,
        float score, char strand, int8_t phase,
        std::unordered_map<std::string, std::string> attributes)
        : seqid{seqid}, source{source}, type{type}, start{start}, end{end}, score{score},
        strand{strand}, phase{phase}, attributes{attributes} {}
};

// represents a sequencing read from FASTQ file
struct fastq_entry {
    std::string seqid;
    std::string seq;
    std::string qual;

    fastq_entry() = default;
    fastq_entry(std::string seqid, std::string seq, std::string qual)
        : seqid{seqid}, seq{seq}, qual{qual} {}
};

struct alignment_entry {
    std::string qname;
    std::string rname;
    int pos;
    int mapq;
    std::string cigar;
    std::string seq;
    std::string qual;
    int flag;

    // optional fields
    std::optional<int> tlen;

    alignment_entry() = default;
    alignment_entry(
        std::string qname, std::string rname, int pos, int mapq, std::string cigar,
        std::string seq, std::string qual, int flag)
        : qname{qname}, rname{rname}, pos{pos}, mapq{mapq},
            cigar{cigar}, seq{seq}, qual{qual}, flag{flag} {}
};

#endif //ATROPLEX_DATA_ENTRIES_HPP

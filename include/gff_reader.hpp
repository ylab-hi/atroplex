/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */
#ifndef ATROPLEX_GFF_READER_HPP
#define ATROPLEX_GFF_READER_HPP

// standard
#include <string>
#include <filesystem>
#include <fstream>

// class
#include "file_reader.hpp"
#include "file_entries.hpp"

class gff_reader : public file_reader<gff_entry> {
public:
    explicit gff_reader(const std::filesystem::path& filepath);
    ~gff_reader();

    // Read next entry
    bool read_next(gff_entry& entry);

    // Check if more entries available
    bool has_next() const { return !eof_reached; }

    // Get current line number (for error reporting)
    size_t get_current_line() const { return line_num; }

private:
    std::ifstream file;
    size_t line_num;
    bool eof_reached;
    bool is_gtf;  // true for GTF, false for GFF3

    // Parse a single line into gff_entry
    bool parse_line(const std::string& line, gff_entry& entry);

    // Parse attributes field (format differs between GFF3 and GTF)
    void parse_attributes(const std::string& attr_str, gff_entry& entry);
    void parse_gff3_attributes(const std::string& attr_str, gff_entry& entry);
    void parse_gtf_attributes(const std::string& attr_str, gff_entry& entry);
};

#endif //ATROPLEX_GFF_READER_HPP

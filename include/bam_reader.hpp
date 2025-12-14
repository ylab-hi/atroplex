/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_BAM_READER_HPP
#define ATROPLEX_BAM_READER_HPP

// standard
#include <string>
#include <filesystem>

// class
#include "file_reader.hpp"

// htslib for reading BAM/SAM files

class bam_reader : public file_reader<alignment_entry> {
    public:
        explicit bam_reader(const std::filesystem::path& path);
        bool read_next(alignment_entry& entry) override;
        bool has_next() override;
        std::string get_error_message() override;
        size_t get_current_line() override;
        ~bam_reader();

    private:
        samFile* file;
        bam_hdr_t* header;
        bam1_t* record;
        size_t line_num;
        std::string error_message;
        bool eof_reached;
};

#endif //ATROPLEX_BAM_READER_HPP

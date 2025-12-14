/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_FILETYPE_DETECTOR_HPP
#define ATROPLEX_FILETYPE_DETECTOR_HPP

// standard
#include <tuple>
#include <filesystem>

enum class filetype {
    FASTQ, BAM, GFF, GTF, GG, UNKNOWN
};

class filetype_detector {
public:
    std::tuple<filetype, bool> detect_filetype(const std::filesystem::path& filepath);

private:
    std::tuple<filetype, bool> detect_plain_filetype(const char* buffer, std::streamsize size);
    std::tuple<filetype, bool> detect_gzipped_filetype(const std::filesystem::path& filepath);
};

#endif //ATROPLEX_FILETYPE_DETECTOR_HPP

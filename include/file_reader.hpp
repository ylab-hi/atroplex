/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_FILE_READER_HPP
#define ATROPLEX_FILE_READER_HPP

#include <string>

// Base class for all file readers
class file_reader_base {
    public:
        virtual bool has_next() = 0;
        virtual std::string get_error_message() = 0;
        virtual size_t get_current_line() = 0;
        virtual ~file_reader_base() = default;
};

// Templated derived class for type-specific reading
template<typename EntryType>
class file_reader : public file_reader_base {
    public:
        virtual bool read_next(EntryType& entry) = 0;
};

#endif //ATROPLEX_FILE_READER_HPP

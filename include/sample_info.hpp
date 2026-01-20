/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SAMPLE_INFO_HPP
#define ATROPLEX_SAMPLE_INFO_HPP

#include <string>
#include <unordered_map>
#include <filesystem>
#include <ostream>
#include <istream>
#include <type_traits>

// genogrove
#include <genogrove/data_type/data_registry.hpp>

namespace gdt = genogrove::data_type;

/**
 * Sample metadata for pan-transcriptome tracking
 *
 * Stores information about each sample/annotation source.
 * Designed to be stored in genogrove's data_registry.
 */
struct sample_info {
    // Expression quantification types
    enum class expression_type {
        UNKNOWN,    // Not specified
        TPM,        // Transcripts per million
        FPKM,       // Fragments per kilobase per million
        RPKM,       // Reads per kilobase per million
        COUNTS,     // Raw read counts
        CPM         // Counts per million
    };

    // Core identifiers
    std::string id;                     // Unique sample identifier
    std::filesystem::path source_file;  // Original annotation file path

    // Biological metadata
    std::string tissue;                 // e.g., "brain", "liver", "heart"
    std::string condition;              // e.g., "tumor", "healthy", "treated"
    std::string species;                // e.g., "human", "mouse"

    // Technical metadata
    std::string annotation_source;      // e.g., "GENCODE", "Ensembl", "RefSeq"
    std::string annotation_version;     // e.g., "v44", "release 110"

    // Expression metadata
    std::filesystem::path expression_file;  // TSV with feature_id -> expression
    expression_type expr_type = expression_type::UNKNOWN;  // Quantification unit

    // Extensible key-value attributes
    std::unordered_map<std::string, std::string> attributes;

    // Constructors
    sample_info() = default;

    explicit sample_info(std::string sample_id)
        : id(std::move(sample_id)) {}

    sample_info(std::string sample_id, std::filesystem::path file)
        : id(std::move(sample_id)), source_file(std::move(file)) {}

    // Builder-style setters for fluent API
    sample_info& with_tissue(std::string t) {
        tissue = std::move(t);
        return *this;
    }

    sample_info& with_condition(std::string c) {
        condition = std::move(c);
        return *this;
    }

    sample_info& with_species(std::string s) {
        species = std::move(s);
        return *this;
    }

    sample_info& with_annotation_source(std::string src) {
        annotation_source = std::move(src);
        return *this;
    }

    sample_info& with_annotation_version(std::string ver) {
        annotation_version = std::move(ver);
        return *this;
    }

    sample_info& with_expression_file(std::filesystem::path path) {
        expression_file = std::move(path);
        return *this;
    }

    sample_info& with_expression_type(expression_type type) {
        expr_type = type;
        return *this;
    }

    sample_info& with_attribute(const std::string& key, std::string value) {
        attributes[key] = std::move(value);
        return *this;
    }

    // Accessors
    bool has_expression_data() const {
        return !expression_file.empty();
    }

    std::string get_attribute(const std::string& key,
                              const std::string& default_value = "") const {
        auto it = attributes.find(key);
        return it != attributes.end() ? it->second : default_value;
    }

    // --- Serialization for data_registry ---

    void serialize(std::ostream& os) const {
        // Helper to write strings
        auto write_string = [&os](const std::string& s) {
            size_t len = s.size();
            os.write(reinterpret_cast<const char*>(&len), sizeof(len));
            os.write(s.data(), len);
        };

        auto write_path = [&os, &write_string](const std::filesystem::path& p) {
            write_string(p.string());
        };

        write_string(id);
        write_path(source_file);
        write_string(tissue);
        write_string(condition);
        write_string(species);
        write_string(annotation_source);
        write_string(annotation_version);
        write_path(expression_file);

        // Expression type (as underlying int)
        auto expr_int = static_cast<std::underlying_type_t<expression_type>>(expr_type);
        os.write(reinterpret_cast<const char*>(&expr_int), sizeof(expr_int));

        // Attributes map
        size_t attr_count = attributes.size();
        os.write(reinterpret_cast<const char*>(&attr_count), sizeof(attr_count));
        for (const auto& [key, value] : attributes) {
            write_string(key);
            write_string(value);
        }
    }

    static sample_info deserialize(std::istream& is) {
        // Helper to read strings
        auto read_string = [&is]() -> std::string {
            size_t len;
            is.read(reinterpret_cast<char*>(&len), sizeof(len));
            std::string s(len, '\0');
            is.read(s.data(), len);
            return s;
        };

        auto read_path = [&read_string]() -> std::filesystem::path {
            return std::filesystem::path(read_string());
        };

        sample_info info;
        info.id = read_string();
        info.source_file = read_path();
        info.tissue = read_string();
        info.condition = read_string();
        info.species = read_string();
        info.annotation_source = read_string();
        info.annotation_version = read_string();
        info.expression_file = read_path();

        // Expression type
        std::underlying_type_t<expression_type> expr_int;
        is.read(reinterpret_cast<char*>(&expr_int), sizeof(expr_int));
        info.expr_type = static_cast<expression_type>(expr_int);

        // Attributes map
        size_t attr_count;
        is.read(reinterpret_cast<char*>(&attr_count), sizeof(attr_count));
        for (size_t i = 0; i < attr_count; ++i) {
            std::string key = read_string();
            std::string value = read_string();
            info.attributes[key] = std::move(value);
        }

        return info;
    }
};

// Type alias for the sample registry singleton
using sample_registry = gdt::data_registry<sample_info>;

#endif //ATROPLEX_SAMPLE_INFO_HPP
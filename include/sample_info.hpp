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
 * Stores information about each input file (annotation or biological sample).
 * In the pan-transcriptome context, both reference annotations and sample
 * assemblies are treated as "samples" of the transcriptome.
 *
 * Designed to be stored in genogrove's data_registry.
 *
 * Field descriptions aligned with ENCODE standards:
 * https://www.encodeproject.org/data-standards/
 *
 * Core fields:
 *   - id: Unique identifier (auto-generated from file metadata or user-provided)
 *   - description: Free-text description of the sample/annotation
 *   - source_file: Path to the input file
 *
 * Biological metadata (ENCODE-aligned):
 *   - assay: Experimental assay (e.g., "RNA-seq", "CAGE", "polyA RNA-seq")
 *   - biosample_type: Type of sample (e.g., "tissue", "cell line", "primary cell")
 *   - biosample: Specific biosample (e.g., "brain", "HeLa", "GM12878")
 *   - condition: Experimental condition (e.g., "tumor", "healthy", "treated")
 *   - treatment: Treatment applied (e.g., "dexamethasone 100nM 1h")
 *   - species: Species (e.g., "Homo sapiens", "Mus musculus")
 *   - replication_type: Replicate type (e.g., "biological", "technical", "isogenic")
 *
 * Technical/processing metadata:
 *   - platform: Sequencing platform (e.g., "Illumina NovaSeq", "PacBio Sequel II")
 *   - pipeline: Analysis pipeline (e.g., "GENCODE", "StringTie", "ENCODE long-read")
 *   - pipeline_version: Pipeline version (e.g., "v3.2.1")
 *
 * Annotation metadata (for reference annotations):
 *   - annotation_source: Source database (e.g., "GENCODE", "Ensembl", "RefSeq")
 *   - annotation_version: Version string (e.g., "v44", "release 110")
 *
 * Links:
 *   - source_url: URL to original data (e.g., ENCODE portal, GEO accession)
 *   - publication: Publication DOI or reference
 *
 * Expression (optional override for GFF attribute parsing):
 *   - expr_type: Which expression attribute to parse from GFF (TPM, FPKM, etc.)
 *                If UNKNOWN, auto-detects from available GFF attributes
 *
 * Flexible key-value storage:
 *   - attributes: Map for any additional metadata
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
    std::string id;                     // Unique identifier
    std::string type = "sample";        // Entry type: "sample" or "annotation"
    std::string description;            // Free-text description
    std::filesystem::path source_file;  // Original input file path

    // Biological metadata (ENCODE-aligned)
    std::string assay;                  // e.g., "RNA-seq", "CAGE", "polyA RNA-seq"
    std::string biosample_type;         // e.g., "tissue", "cell line", "primary cell"
    std::string biosample;              // e.g., "brain", "HeLa", "GM12878"
    std::string condition;              // e.g., "tumor", "healthy", "treated"
    std::string treatment;              // e.g., "dexamethasone 100nM 1h"
    std::string species;                // e.g., "Homo sapiens", "Mus musculus"
    std::string replication_type;       // e.g., "biological", "technical", "isogenic"

    // Technical/processing metadata
    std::string platform;               // e.g., "Illumina NovaSeq", "PacBio Sequel II"
    std::string pipeline;               // e.g., "GENCODE", "StringTie", "ENCODE long-read"
    std::string pipeline_version;       // e.g., "v3.2.1"

    // Annotation metadata (for reference annotations)
    std::string annotation_source;      // e.g., "GENCODE", "Ensembl", "RefSeq"
    std::string annotation_version;     // e.g., "v44", "release 110"

    // Links
    std::string source_url;             // URL to original data
    std::string publication;            // Publication DOI or reference

    // Expression type override (optional)
    // If UNKNOWN, auto-detect from GFF attributes (TPM, FPKM, cov, etc.)
    // If set, only parse that specific attribute from GFF
    expression_type expr_type = expression_type::UNKNOWN;

    // Extensible key-value attributes
    std::unordered_map<std::string, std::string> attributes;

    // Constructors
    sample_info() = default;

    explicit sample_info(std::string sample_id)
        : id(std::move(sample_id)) {}

    sample_info(std::string sample_id, std::filesystem::path file)
        : id(std::move(sample_id)), source_file(std::move(file)) {}

    // Builder-style setters for fluent API
    sample_info& with_description(std::string d) {
        description = std::move(d);
        return *this;
    }

    sample_info& with_assay(std::string a) {
        assay = std::move(a);
        return *this;
    }

    sample_info& with_biosample_type(std::string bt) {
        biosample_type = std::move(bt);
        return *this;
    }

    sample_info& with_biosample(std::string b) {
        biosample = std::move(b);
        return *this;
    }

    sample_info& with_condition(std::string c) {
        condition = std::move(c);
        return *this;
    }

    sample_info& with_treatment(std::string t) {
        treatment = std::move(t);
        return *this;
    }

    sample_info& with_species(std::string s) {
        species = std::move(s);
        return *this;
    }

    sample_info& with_replication_type(std::string rt) {
        replication_type = std::move(rt);
        return *this;
    }

    sample_info& with_platform(std::string p) {
        platform = std::move(p);
        return *this;
    }

    sample_info& with_pipeline(std::string p) {
        pipeline = std::move(p);
        return *this;
    }

    sample_info& with_pipeline_version(std::string pv) {
        pipeline_version = std::move(pv);
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

    sample_info& with_source_url(std::string url) {
        source_url = std::move(url);
        return *this;
    }

    sample_info& with_publication(std::string pub) {
        publication = std::move(pub);
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
    bool has_expression_type() const {
        return expr_type != expression_type::UNKNOWN;
    }

    bool is_annotation() const {
        return !annotation_source.empty() || !annotation_version.empty();
    }

    bool is_sample() const {
        return !biosample.empty() || !assay.empty();
    }

    std::string get_attribute(const std::string& key,
                              const std::string& default_value = "") const {
        auto it = attributes.find(key);
        return it != attributes.end() ? it->second : default_value;
    }

    // --- Serialization for data_registry ---

    void serialize(std::ostream& os) const {
        auto write_string = [&os](const std::string& s) {
            size_t len = s.size();
            os.write(reinterpret_cast<const char*>(&len), sizeof(len));
            os.write(s.data(), static_cast<std::streamsize>(len));
        };

        auto write_path = [&write_string](const std::filesystem::path& p) {
            write_string(p.string());
        };

        // Core
        write_string(id);
        write_string(type);
        write_string(description);
        write_path(source_file);

        // Biological
        write_string(assay);
        write_string(biosample_type);
        write_string(biosample);
        write_string(condition);
        write_string(treatment);
        write_string(species);
        write_string(replication_type);

        // Technical
        write_string(platform);
        write_string(pipeline);
        write_string(pipeline_version);

        // Annotation
        write_string(annotation_source);
        write_string(annotation_version);

        // Links
        write_string(source_url);
        write_string(publication);

        // Expression type
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
        auto read_string = [&is]() -> std::string {
            size_t len;
            is.read(reinterpret_cast<char*>(&len), sizeof(len));
            std::string s(len, '\0');
            is.read(s.data(), static_cast<std::streamsize>(len));
            return s;
        };

        auto read_path = [&read_string]() -> std::filesystem::path {
            return std::filesystem::path(read_string());
        };

        sample_info info;

        // Core
        info.id = read_string();
        info.type = read_string();
        info.description = read_string();
        info.source_file = read_path();

        // Biological
        info.assay = read_string();
        info.biosample_type = read_string();
        info.biosample = read_string();
        info.condition = read_string();
        info.treatment = read_string();
        info.species = read_string();
        info.replication_type = read_string();

        // Technical
        info.platform = read_string();
        info.pipeline = read_string();
        info.pipeline_version = read_string();

        // Annotation
        info.annotation_source = read_string();
        info.annotation_version = read_string();

        // Links
        info.source_url = read_string();
        info.publication = read_string();

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
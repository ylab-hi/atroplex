/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SEQUENCE_VARIANT_HPP
#define ATROPLEX_SEQUENCE_VARIANT_HPP

#include <cstdint>
#include <string>
#include <unordered_set>
#include <ostream>
#include <istream>
#include <type_traits>

// genogrove
#include <genogrove/data_type/data_registry.hpp>

namespace gdt = genogrove::data_type;

/**
 * Sequence variant for pan-transcriptome provenance tracking
 *
 * Stores SNPs, indels, and other sequence-level variants with
 * provenance information (which samples contain this variant).
 * Designed for data_registry storage with lightweight ID references.
 *
 * Usage:
 *   auto& registry = variant_registry::instance();
 *   uint32_t var_id = registry.register_data(
 *       sequence_variant("chr1", 12345, "A", "G")
 *           .with_id("rs12345")
 *           .with_sample(gencode_sample_id)
 *   );
 */
struct sequence_variant {
    // Variant effect on transcript
    enum class effect_type {
        UNKNOWN,
        SYNONYMOUS,      // No amino acid change
        MISSENSE,        // Amino acid substitution
        NONSENSE,        // Premature stop codon
        FRAMESHIFT,      // Indel not multiple of 3
        SPLICE_SITE,     // Affects splicing
        UTR,             // In untranslated region
        INTRONIC         // In intron (if retained)
    };

    // Core variant identifiers
    std::string id;              // e.g., rs12345 or chr1:12345:A>G
    std::string seqid;           // Chromosome/contig
    uint64_t position;           // 0-based genomic position
    std::string ref;             // Reference allele
    std::string alt;             // Alternate allele

    // Provenance: sample_registry IDs where this variant was observed
    std::unordered_set<uint32_t> sample_ids;

    // Optional annotation
    effect_type effect = effect_type::UNKNOWN;
    std::string gene_id;         // Affected gene (if known)
    std::string transcript_id;   // Affected transcript (if known)

    // Constructors
    sequence_variant() : position(0) {}

    sequence_variant(std::string chr, uint64_t pos, std::string reference, std::string alternate)
        : seqid(std::move(chr))
        , position(pos)
        , ref(std::move(reference))
        , alt(std::move(alternate)) {
        // Generate default ID if not provided
        id = seqid + ":" + std::to_string(position) + ":" + ref + ">" + alt;
    }

    // Builder-style setters
    sequence_variant& with_id(std::string var_id) {
        id = std::move(var_id);
        return *this;
    }

    sequence_variant& with_sample(uint32_t sample_id) {
        sample_ids.insert(sample_id);
        return *this;
    }

    sequence_variant& with_effect(effect_type eff) {
        effect = eff;
        return *this;
    }

    sequence_variant& with_gene(std::string gid) {
        gene_id = std::move(gid);
        return *this;
    }

    sequence_variant& with_transcript(std::string tid) {
        transcript_id = std::move(tid);
        return *this;
    }

    // --- Provenance methods ---

    // Add a sample that has this variant
    void add_sample(uint32_t sample_id) {
        sample_ids.insert(sample_id);
    }

    // Check if variant is present in a specific sample
    bool in_sample(uint32_t sample_id) const {
        return sample_ids.contains(sample_id);
    }

    // Check if variant is novel (not in any sample in pan-transcriptome)
    bool is_novel() const {
        return sample_ids.empty();
    }

    // Get number of samples with this variant
    size_t sample_count() const {
        return sample_ids.size();
    }

    // Calculate allele frequency in pan-transcriptome
    float frequency(size_t total_samples) const {
        if (total_samples == 0) return 0.0f;
        return static_cast<float>(sample_ids.size()) / static_cast<float>(total_samples);
    }

    // Check if variant is rare (below threshold frequency)
    bool is_rare(size_t total_samples, float threshold = 0.01f) const {
        return frequency(total_samples) < threshold;
    }

    // Check if variant is common (above threshold frequency)
    bool is_common(size_t total_samples, float threshold = 0.05f) const {
        return frequency(total_samples) >= threshold;
    }

    // --- Serialization for data_registry ---

    void serialize(std::ostream& os) const {
        auto write_string = [&os](const std::string& s) {
            size_t len = s.size();
            os.write(reinterpret_cast<const char*>(&len), sizeof(len));
            os.write(s.data(), static_cast<std::streamsize>(len));
        };

        write_string(id);
        write_string(seqid);
        os.write(reinterpret_cast<const char*>(&position), sizeof(position));
        write_string(ref);
        write_string(alt);

        // Sample IDs
        size_t sample_count = sample_ids.size();
        os.write(reinterpret_cast<const char*>(&sample_count), sizeof(sample_count));
        for (uint32_t sid : sample_ids) {
            os.write(reinterpret_cast<const char*>(&sid), sizeof(sid));
        }

        // Effect type
        auto effect_int = static_cast<std::underlying_type_t<effect_type>>(effect);
        os.write(reinterpret_cast<const char*>(&effect_int), sizeof(effect_int));

        write_string(gene_id);
        write_string(transcript_id);
    }

    static sequence_variant deserialize(std::istream& is) {
        auto read_string = [&is]() -> std::string {
            size_t len;
            is.read(reinterpret_cast<char*>(&len), sizeof(len));
            std::string s(len, '\0');
            is.read(s.data(), static_cast<std::streamsize>(len));
            return s;
        };

        sequence_variant var;
        var.id = read_string();
        var.seqid = read_string();
        is.read(reinterpret_cast<char*>(&var.position), sizeof(var.position));
        var.ref = read_string();
        var.alt = read_string();

        // Sample IDs
        size_t sample_count;
        is.read(reinterpret_cast<char*>(&sample_count), sizeof(sample_count));
        for (size_t i = 0; i < sample_count; ++i) {
            uint32_t sid;
            is.read(reinterpret_cast<char*>(&sid), sizeof(sid));
            var.sample_ids.insert(sid);
        }

        // Effect type
        std::underlying_type_t<effect_type> effect_int;
        is.read(reinterpret_cast<char*>(&effect_int), sizeof(effect_int));
        var.effect = static_cast<effect_type>(effect_int);

        var.gene_id = read_string();
        var.transcript_id = read_string();

        return var;
    }
};

// Type alias for the variant registry singleton
using variant_registry = gdt::data_registry<sequence_variant>;

#endif //ATROPLEX_SEQUENCE_VARIANT_HPP
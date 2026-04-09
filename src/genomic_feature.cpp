/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "genomic_feature.hpp"

#include <sstream>
#include <stdexcept>

/**
 * Helper to extract attribute from map, return empty string if not found
 */
static std::string get_attribute(
    const std::map<std::string, std::string, std::less<>>& attributes,
    const std::string& key
) {
    auto it = attributes.find(key);
    return (it != attributes.end()) ? it->second : "";
}

exon_feature exon_feature::from_gff_entry(
    const std::map<std::string, std::string, std::less<>>& attributes,
    const std::string& seqid,
    const gdt::interval& interval,
    char strand
) {
    exon_feature feature;

    // Extract common attributes
    feature.id = get_attribute(attributes, "ID");
    if (feature.id.empty()) {
        feature.id = get_attribute(attributes, "exon_id");
    }

    std::string gid = get_attribute(attributes, "gene_id");
    std::string gname = get_attribute(attributes, "gene_name");
    std::string gbiotype = get_attribute(attributes, "gene_type");
    if (gbiotype.empty()) {
        gbiotype = get_attribute(attributes, "gene_biotype");
    }
    feature.gene_idx = gene_registry::instance().intern(gid, gname, gbiotype);

    // Transcript information
    std::string transcript_id = get_attribute(attributes, "transcript_id");
    if (transcript_id.empty()) {
        transcript_id = get_attribute(attributes, "Parent");
    }
    if (!transcript_id.empty()) {
        feature.transcript_ids.insert(transcript_registry::instance().intern(transcript_id));
    }

    // Note: GFF source (column 2) is added separately via add_source() in build_gff

    // Exon number
    std::string exon_num_str = get_attribute(attributes, "exon_number");
    if (!exon_num_str.empty()) {
        try {
            feature.exon_number = std::stoi(exon_num_str);
        } catch (...) {
            feature.exon_number = -1;
        }
    }

    return feature;
}

// ============================================================================
// Binary serialization helpers
// ============================================================================

namespace {

template<typename T>
void write_pod(std::ostream& os, const T& value) {
    os.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T>
T read_pod(std::istream& is) {
    T value;
    is.read(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
}

void write_string(std::ostream& os, const std::string& s) {
    uint32_t len = static_cast<uint32_t>(s.size());
    write_pod(os, len);
    os.write(s.data(), len);
}

std::string read_string(std::istream& is) {
    uint32_t len = read_pod<uint32_t>(is);
    std::string s(len, '\0');
    is.read(s.data(), len);
    return s;
}

} // anonymous namespace

// ============================================================================
// sample_bitset serialization
// ============================================================================

void sample_bitset::serialize(std::ostream& os) const {
    uint32_t n = static_cast<uint32_t>(words_.size());
    write_pod(os, n);
    if (n > 0) {
        os.write(reinterpret_cast<const char*>(words_.data()),
                 static_cast<std::streamsize>(n * sizeof(uint64_t)));
    }
}

sample_bitset sample_bitset::deserialize(std::istream& is) {
    sample_bitset bs;
    uint32_t n = read_pod<uint32_t>(is);
    bs.words_.resize(n);
    if (n > 0) {
        is.read(reinterpret_cast<char*>(bs.words_.data()),
                static_cast<std::streamsize>(n * sizeof(uint64_t)));
    }
    return bs;
}

// ============================================================================
// expression_store serialization
// ============================================================================

void expression_store::serialize(std::ostream& os) const {
    if (!data_) {
        write_pod<uint32_t>(os, 0);
        return;
    }
    uint32_t n = static_cast<uint32_t>(data_->size());
    write_pod(os, n);
    if (n > 0) {
        os.write(reinterpret_cast<const char*>(data_->data()),
                 static_cast<std::streamsize>(n * sizeof(float)));
    }
}

expression_store expression_store::deserialize(std::istream& is) {
    expression_store store;
    uint32_t n = read_pod<uint32_t>(is);
    if (n > 0) {
        store.data_ = std::make_unique<std::vector<float>>(n);
        is.read(reinterpret_cast<char*>(store.data_->data()),
                static_cast<std::streamsize>(n * sizeof(float)));
    }
    return store;
}

// ============================================================================
// sorted_vec serialization
// ============================================================================

void sorted_vec::serialize(std::ostream& os) const {
    uint32_t n = static_cast<uint32_t>(data_.size());
    write_pod(os, n);
    if (n > 0) {
        os.write(reinterpret_cast<const char*>(data_.data()),
                 static_cast<std::streamsize>(n * sizeof(uint32_t)));
    }
}

sorted_vec sorted_vec::deserialize(std::istream& is) {
    sorted_vec sv;
    uint32_t n = read_pod<uint32_t>(is);
    sv.data_.resize(n);
    if (n > 0) {
        is.read(reinterpret_cast<char*>(sv.data_.data()),
                static_cast<std::streamsize>(n * sizeof(uint32_t)));
    }
    return sv;
}

// ============================================================================
// edge_metadata serialization
// ============================================================================

void edge_metadata::serialize(std::ostream& os) const {
    write_pod(os, id);
    write_pod(os, static_cast<uint8_t>(type));
}

edge_metadata edge_metadata::deserialize(std::istream& is) {
    edge_metadata em;
    em.id = read_pod<size_t>(is);
    em.type = static_cast<edge_type>(read_pod<uint8_t>(is));
    return em;
}

// ============================================================================
// exon_feature serialization
// ============================================================================

void exon_feature::serialize(std::ostream& os) const {
    write_string(os, id);
    write_pod(os, gene_idx);
    transcript_ids.serialize(os);
    write_pod(os, sources);
    sample_idx.serialize(os);
    expression.serialize(os);
    write_pod(os, exon_number);
}

exon_feature exon_feature::deserialize(std::istream& is) {
    exon_feature ef;
    ef.id = read_string(is);
    ef.gene_idx = read_pod<uint32_t>(is);
    ef.transcript_ids = sorted_vec::deserialize(is);
    ef.sources = read_pod<uint16_t>(is);
    ef.sample_idx = sample_bitset::deserialize(is);
    ef.expression = expression_store::deserialize(is);
    ef.exon_number = read_pod<int>(is);
    return ef;
}

// ============================================================================
// segment_feature serialization
// ============================================================================

void segment_feature::serialize(std::ostream& os) const {
    write_pod(os, gene_idx);
    transcript_ids.serialize(os);
    write_pod(os, sources);
    sample_idx.serialize(os);
    expression.serialize(os);

    // transcript_biotypes map
    uint32_t bt_count = static_cast<uint32_t>(transcript_biotypes.size());
    write_pod(os, bt_count);
    for (const auto& [tx_id, biotype] : transcript_biotypes) {
        write_pod(os, tx_id);
        write_string(os, biotype);
    }

    write_pod(os, segment_index);
    write_pod(os, exon_count);
    write_pod(os, read_coverage);
    write_pod(os, absorbed);
    write_pod(os, absorbed_count);
}

segment_feature segment_feature::deserialize(std::istream& is) {
    segment_feature sf;
    sf.gene_idx = read_pod<uint32_t>(is);
    sf.transcript_ids = sorted_vec::deserialize(is);
    sf.sources = read_pod<uint16_t>(is);
    sf.sample_idx = sample_bitset::deserialize(is);
    sf.expression = expression_store::deserialize(is);

    uint32_t bt_count = read_pod<uint32_t>(is);
    for (uint32_t i = 0; i < bt_count; ++i) {
        uint32_t tx_id = read_pod<uint32_t>(is);
        std::string biotype = read_string(is);
        sf.transcript_biotypes[tx_id] = std::move(biotype);
    }

    sf.segment_index = read_pod<size_t>(is);
    sf.exon_count = read_pod<int>(is);
    sf.read_coverage = read_pod<size_t>(is);
    sf.absorbed = read_pod<bool>(is);
    sf.absorbed_count = read_pod<size_t>(is);
    return sf;
}

// ============================================================================
// genomic_feature variant serialization
// ============================================================================

namespace genogrove::data_type {

void serialization_traits<genomic_feature>::serialize(std::ostream& os,
                                                      const genomic_feature& feature) {
    uint8_t tag = static_cast<uint8_t>(feature.index());  // 0 = exon, 1 = segment
    os.write(reinterpret_cast<const char*>(&tag), sizeof(tag));
    std::visit([&os](const auto& f) { f.serialize(os); }, feature);
}

genomic_feature serialization_traits<genomic_feature>::deserialize(std::istream& is) {
    uint8_t tag;
    is.read(reinterpret_cast<char*>(&tag), sizeof(tag));
    if (tag == 0) {
        return exon_feature::deserialize(is);
    } else if (tag == 1) {
        return segment_feature::deserialize(is);
    }
    throw std::runtime_error("Invalid genomic_feature variant tag: " + std::to_string(tag));
}

} // namespace genogrove::data_type

// ============================================================================
// Registry serialization
// ============================================================================

void gene_registry::serialize(std::ostream& os) const {
    uint32_t n = static_cast<uint32_t>(entries_.size());
    write_pod(os, n);
    for (const auto& entry : entries_) {
        write_string(os, entry.gene_id);
        write_string(os, entry.gene_name);
        write_string(os, entry.gene_biotype);
    }
}

void gene_registry::deserialize_into(std::istream& is) {
    uint32_t n = read_pod<uint32_t>(is);
    for (uint32_t i = 0; i < n; ++i) {
        std::string gene_id = read_string(is);
        std::string gene_name = read_string(is);
        std::string gene_biotype = read_string(is);
        intern(gene_id, gene_name, gene_biotype);
    }
}

void source_registry::serialize(std::ostream& os) const {
    uint32_t n = static_cast<uint32_t>(bit_to_str_.size());
    write_pod(os, n);
    for (const auto& s : bit_to_str_) {
        write_string(os, s);
    }
}

void source_registry::deserialize_into(std::istream& is) {
    uint32_t n = read_pod<uint32_t>(is);
    for (uint32_t i = 0; i < n; ++i) {
        std::string s = read_string(is);
        intern(s);
    }
}

void transcript_registry::serialize(std::ostream& os) const {
    uint32_t n = static_cast<uint32_t>(id_to_str_.size());
    write_pod(os, n);
    for (const auto& s : id_to_str_) {
        write_string(os, s);
    }
}

void transcript_registry::deserialize_into(std::istream& is) {
    uint32_t n = read_pod<uint32_t>(is);
    for (uint32_t i = 0; i < n; ++i) {
        std::string s = read_string(is);
        intern(s);
    }
}
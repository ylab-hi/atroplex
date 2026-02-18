/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_GENOMIC_FEATURE_HPP
#define ATROPLEX_GENOMIC_FEATURE_HPP

// standard
#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <variant>
#include <map>

// genogrove
#include <genogrove/data_type/interval.hpp>
#include <genogrove/data_type/genomic_coordinate.hpp>
#include <genogrove/structure/grove/grove.hpp>

namespace gdt = genogrove::data_type;
namespace gst = genogrove::structure;

/**
 * Dynamic bitset for tracking sample membership.
 * Stores a vector<uint64_t> where each bit corresponds to a sample registry ID.
 * Grows automatically when a sample ID beyond the current capacity is set.
 *
 * Typical sizes: 128 samples = 2 words (16 bytes data), 1024 = 16 words (128 bytes).
 * Much more compact than unordered_set<uint32_t> which costs ~56 bytes base + 12 per entry.
 */
class sample_bitset {
    std::vector<uint64_t> words_;

    static size_t word_index(uint32_t id) { return id / 64; }
    static uint64_t bit_mask(uint32_t id) { return uint64_t{1} << (id % 64); }

    void ensure_capacity(uint32_t id) {
        size_t needed = word_index(id) + 1;
        if (needed > words_.size()) {
            words_.resize(needed, 0);
        }
    }

public:
    sample_bitset() = default;

    void set(uint32_t id) {
        ensure_capacity(id);
        words_[word_index(id)] |= bit_mask(id);
    }

    bool test(uint32_t id) const {
        size_t wi = word_index(id);
        if (wi >= words_.size()) return false;
        return (words_[wi] & bit_mask(id)) != 0;
    }

    size_t count() const {
        size_t n = 0;
        for (uint64_t w : words_) n += static_cast<size_t>(__builtin_popcountll(w));
        return n;
    }

    bool empty() const {
        for (uint64_t w : words_) if (w) return false;
        return true;
    }

    /// Merge another bitset into this one (bitwise OR)
    void merge(const sample_bitset& other) {
        if (other.words_.size() > words_.size()) {
            words_.resize(other.words_.size(), 0);
        }
        for (size_t i = 0; i < other.words_.size(); ++i) {
            words_[i] |= other.words_[i];
        }
    }

    /// Iterate over set bits, calling fn(uint32_t sample_id) for each
    template<typename Fn>
    void for_each(Fn&& fn) const {
        for (size_t wi = 0; wi < words_.size(); ++wi) {
            uint64_t w = words_[wi];
            while (w) {
                uint32_t bit = static_cast<uint32_t>(__builtin_ctzll(w));
                fn(static_cast<uint32_t>(wi * 64 + bit));
                w &= w - 1;  // clear lowest set bit
            }
        }
    }

    /// STL-style iterator for range-based for loops
    class iterator {
        const sample_bitset* bs_;
        size_t wi_;
        uint64_t remaining_;
        uint32_t current_;

        void advance() {
            while (remaining_ == 0) {
                wi_++;
                if (wi_ >= bs_->words_.size()) return;
                remaining_ = bs_->words_[wi_];
            }
            uint32_t bit = static_cast<uint32_t>(__builtin_ctzll(remaining_));
            current_ = static_cast<uint32_t>(wi_ * 64 + bit);
            remaining_ &= remaining_ - 1;
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = uint32_t;
        using difference_type = std::ptrdiff_t;
        using pointer = const uint32_t*;
        using reference = uint32_t;

        // end sentinel
        iterator() : bs_(nullptr), wi_(0), remaining_(0), current_(0) {}

        // begin
        explicit iterator(const sample_bitset& bs)
            : bs_(&bs), wi_(0), remaining_(0), current_(0) {
            if (!bs_->words_.empty()) {
                remaining_ = bs_->words_[0];
                advance();
            } else {
                bs_ = nullptr; // empty → equal to end
            }
        }

        uint32_t operator*() const { return current_; }
        iterator& operator++() { advance(); return *this; }
        iterator operator++(int) { auto tmp = *this; advance(); return tmp; }

        bool operator==(const iterator& o) const {
            if (bs_ == nullptr && o.bs_ == nullptr) return true;
            if (bs_ == nullptr || o.bs_ == nullptr) return false;
            return wi_ == o.wi_ && remaining_ == o.remaining_ && current_ == o.current_;
        }
        bool operator!=(const iterator& o) const { return !(*this == o); }
    };

    iterator begin() const { return iterator(*this); }
    iterator end() const { return iterator(); }

    /// Number of 64-bit words allocated (for memory estimation)
    size_t word_count() const { return words_.size(); }
};

/**
 * Lazily allocated flat expression array.
 * Stores per-sample expression values as a vector<float> indexed by sample registry ID.
 * Uses -1.0f as sentinel for "no data". Only allocates when the first value is set.
 * Null pointer (no allocation) when no expression data exists — 8 bytes per feature.
 */
class expression_store {
    static constexpr float SENTINEL = -1.0f;
    std::unique_ptr<std::vector<float>> data_;

    void ensure(uint32_t sid) {
        if (!data_) {
            data_ = std::make_unique<std::vector<float>>(sid + 1, SENTINEL);
        } else if (sid >= data_->size()) {
            data_->resize(sid + 1, SENTINEL);
        }
    }

public:
    expression_store() = default;

    // Copy support (unique_ptr needs explicit copy)
    expression_store(const expression_store& other)
        : data_(other.data_ ? std::make_unique<std::vector<float>>(*other.data_) : nullptr) {}
    expression_store& operator=(const expression_store& other) {
        if (this != &other) {
            data_ = other.data_ ? std::make_unique<std::vector<float>>(*other.data_) : nullptr;
        }
        return *this;
    }
    expression_store(expression_store&&) noexcept = default;
    expression_store& operator=(expression_store&&) noexcept = default;

    void set(uint32_t sid, float value) {
        ensure(sid);
        (*data_)[sid] = value;
    }

    float get(uint32_t sid) const {
        if (!data_ || sid >= data_->size()) return 0.0f;
        float v = (*data_)[sid];
        return v > SENTINEL ? v : 0.0f;
    }

    bool has(uint32_t sid) const {
        if (!data_ || sid >= data_->size()) return false;
        return (*data_)[sid] > SENTINEL;
    }

    bool empty() const { return !data_; }

    /// Accumulate: add value to existing (or set if not present)
    void accumulate(uint32_t sid, float value) {
        float current = has(sid) ? get(sid) : 0.0f;
        set(sid, current + value);
    }

    /// Merge another store into this one (accumulating values)
    void merge(const expression_store& other) {
        if (!other.data_) return;
        for (size_t i = 0; i < other.data_->size(); ++i) {
            float v = (*other.data_)[i];
            if (v > SENTINEL) {
                accumulate(static_cast<uint32_t>(i), v);
            }
        }
    }

    /// Iterate over (sample_id, value) pairs where expression exists
    template<typename Fn>
    void for_each(Fn&& fn) const {
        if (!data_) return;
        for (size_t i = 0; i < data_->size(); ++i) {
            float v = (*data_)[i];
            if (v > SENTINEL) {
                fn(static_cast<uint32_t>(i), v);
            }
        }
    }

    /// Convert to unordered_map for output struct compatibility
    std::unordered_map<uint32_t, float> to_map() const {
        std::unordered_map<uint32_t, float> m;
        for_each([&](uint32_t sid, float v) { m[sid] = v; });
        return m;
    }

    /// Bytes of heap data allocated (for memory estimation)
    size_t data_bytes() const {
        return data_ ? data_->capacity() * sizeof(float) : 0;
    }
};

/**
 * Sorted flat vector for compact set storage.
 * Drop-in replacement for unordered_set<uint32_t> with much lower memory overhead:
 * unordered_set costs ~40-56 bytes base + ~12 bytes per entry (hash node + bucket pointer).
 * sorted_vec costs ~24 bytes base + 4 bytes per entry (just the uint32_t).
 *
 * Supports insert, count, size, empty, begin/end, and iteration.
 * Insertions maintain sorted order via binary search + insert.
 */
class sorted_vec {
    std::vector<uint32_t> data_;

public:
    sorted_vec() = default;

    void insert(uint32_t val) {
        auto it = std::lower_bound(data_.begin(), data_.end(), val);
        if (it == data_.end() || *it != val) {
            data_.insert(it, val);
        }
    }

    size_t count(uint32_t val) const {
        return std::binary_search(data_.begin(), data_.end(), val) ? 1 : 0;
    }

    size_t size() const { return data_.size(); }
    bool empty() const { return data_.empty(); }

    auto begin() const { return data_.begin(); }
    auto end() const { return data_.end(); }
    auto begin() { return data_.begin(); }
    auto end() { return data_.end(); }

    /// Reserve capacity (useful when approximate count is known)
    void reserve(size_t n) { data_.reserve(n); }

    /// Bytes of heap data (for memory estimation)
    size_t data_bytes() const { return data_.capacity() * sizeof(uint32_t); }
};

/**
 * Transcript ID string pool (intern strings → uint32_t)
 * Singleton; all transcript ID strings are stored once and referenced by index.
 */
class transcript_registry {
    std::unordered_map<std::string, uint32_t> str_to_id_;
    std::vector<std::string> id_to_str_;
public:
    static transcript_registry& instance() {
        static transcript_registry reg;
        return reg;
    }

    uint32_t intern(const std::string& name) {
        auto it = str_to_id_.find(name);
        if (it != str_to_id_.end()) return it->second;
        uint32_t id = static_cast<uint32_t>(id_to_str_.size());
        str_to_id_[name] = id;
        id_to_str_.push_back(name);
        return id;
    }

    const std::string& resolve(uint32_t id) const {
        if (id >= id_to_str_.size()) {
            throw std::out_of_range("transcript_registry: invalid id");
        }
        return id_to_str_[id];
    }

    size_t size() const { return id_to_str_.size(); }
};

/**
 * Gene metadata (stored once per unique gene_id in the registry)
 */
struct gene_info {
    std::string gene_id;
    std::string gene_name;
    std::string gene_biotype;
};

/**
 * Gene string pool (intern gene_id/gene_name/gene_biotype triplet → uint32_t)
 * Singleton; gene metadata is stored once per unique gene_id and referenced by index.
 * All exon/segment features carry a single uint32_t gene_idx instead of three strings.
 */
class gene_registry {
    std::unordered_map<std::string, uint32_t> id_to_idx_;
    std::vector<gene_info> entries_;
public:
    static gene_registry& instance() {
        static gene_registry reg;
        return reg;
    }

    uint32_t intern(const std::string& gene_id,
                    const std::string& gene_name,
                    const std::string& gene_biotype) {
        auto it = id_to_idx_.find(gene_id);
        if (it != id_to_idx_.end()) return it->second;
        uint32_t idx = static_cast<uint32_t>(entries_.size());
        id_to_idx_[gene_id] = idx;
        entries_.push_back({gene_id, gene_name, gene_biotype});
        return idx;
    }

    const gene_info& resolve(uint32_t idx) const {
        if (idx >= entries_.size()) {
            throw std::out_of_range("gene_registry: invalid idx");
        }
        return entries_[idx];
    }

    size_t size() const { return entries_.size(); }
};

/**
 * Source string pool (intern GFF column 2 values → bit positions in a uint16_t)
 * Singleton; supports up to 16 unique source strings (HAVANA, ENSEMBL, TALON, etc.).
 * Features store a uint16_t bitfield instead of unordered_set<string>.
 */
class source_registry {
    std::unordered_map<std::string, uint8_t> str_to_bit_;
    std::vector<std::string> bit_to_str_;
public:
    static source_registry& instance() {
        static source_registry reg;
        return reg;
    }

    /// Returns the bit position for this source string (allocates new if unseen)
    uint8_t intern(const std::string& source) {
        auto it = str_to_bit_.find(source);
        if (it != str_to_bit_.end()) return it->second;
        uint8_t bit = static_cast<uint8_t>(bit_to_str_.size());
        if (bit >= 16) {
            throw std::overflow_error(
                "source_registry: max 16 unique sources exceeded (trying to add '" + source +
                "'). If you need more, widen the sources bitfield from uint16_t to uint32_t.");
        }
        str_to_bit_[source] = bit;
        bit_to_str_.push_back(source);
        return bit;
    }

    /// Returns the bitmask (single bit set) for a source string
    uint16_t mask(const std::string& source) {
        return static_cast<uint16_t>(1u << intern(source));
    }

    const std::string& resolve(uint8_t bit) const {
        if (bit >= bit_to_str_.size()) {
            throw std::out_of_range("source_registry: invalid bit");
        }
        return bit_to_str_[bit];
    }

    size_t size() const { return bit_to_str_.size(); }

    /// Iterate over set bits in a bitfield, calling fn(const string&) for each
    template<typename Fn>
    void for_each(uint16_t bitfield, Fn&& fn) const {
        uint16_t bits = bitfield;
        while (bits) {
            uint8_t bit = static_cast<uint8_t>(__builtin_ctz(bits));
            fn(bit_to_str_[bit]);
            bits &= bits - 1;  // clear lowest set bit
        }
    }

    /// Count number of set bits
    static size_t count(uint16_t bitfield) {
        return static_cast<size_t>(__builtin_popcount(bitfield));
    }
};

/**
 * Edge metadata for grove graph
 * Tracks which transcript/segment an edge belongs to
 */
struct edge_metadata {
    // Numeric segment index identifying which segment this edge belongs to
    size_t id;

    // Type of edge
    enum class edge_type {
        EXON_TO_EXON,      // Exon chain within a segment
        SEGMENT_TO_EXON,   // Segment links to its first exon
        SEGMENT_TO_SEGMENT // Transcript path through segments
    } type;

    edge_metadata() : id(0), type(edge_type::EXON_TO_EXON) {}

    edge_metadata(size_t id, edge_type t)
        : id(id), type(t) {}
};

/// Format a coordinate string from seqid + genomic_coordinate: "chr1:+:100-200"
inline std::string format_coordinate(const std::string& seqid,
                                     const gdt::genomic_coordinate& coord) {
    std::string s;
    s.reserve(seqid.size() + 24);
    s += seqid;
    s += ':';
    s += coord.get_strand();
    s += ':';
    s += std::to_string(coord.get_start());
    s += '-';
    s += std::to_string(coord.get_end());
    return s;
}

/**
 * Exon feature: spatial interval with biological annotations
 * Used for read overlap queries and CDS/UTR disambiguation
 */
struct exon_feature {
    // Core identifiers
    std::string id;                  // Exon ID
    uint32_t gene_idx;               // Index into gene_registry (gene_id, gene_name, gene_biotype)

    // Gene accessors (resolve through gene_registry)
    const std::string& gene_id() const { return gene_registry::instance().resolve(gene_idx).gene_id; }
    const std::string& gene_name() const { return gene_registry::instance().resolve(gene_idx).gene_name; }
    const std::string& gene_biotype() const { return gene_registry::instance().resolve(gene_idx).gene_biotype; }

    // Transcript associations (interned via transcript_registry)
    sorted_vec transcript_ids;  // Reference transcripts using this exon

    // Source tracking (GFF column 2) — bitfield indexed by source_registry
    uint16_t sources = 0;

    // Sample tracking — dynamic bitset indexed by sample_registry IDs
    sample_bitset sample_idx;

    // Expression quantification (lazily allocated, -1.0f sentinel for missing)
    expression_store expression;

    // Additional annotations
    int exon_number;                 // Exon number within transcript (-1 for unknown)

    // Constructors
    exon_feature() : gene_idx(0), sources(0), exon_number(-1) {}

    // Create exon from GFF entry
    static exon_feature from_gff_entry(
        const std::map<std::string, std::string>& attributes,
        const std::string& seqid,
        const gdt::interval& interval,
        char strand
    );

    // --- Sample tracking methods ---

    // Add sample reference (registry ID from sample_registry)
    void add_sample(uint32_t sample_id) {
        sample_idx.set(sample_id);
    }

    // Add source (GFF column 2 value like HAVANA, ENSEMBL, StringTie)
    void add_source(const std::string& source) {
        sources |= source_registry::instance().mask(source);
    }

    // Merge sources from another feature's bitfield
    void merge_sources(uint16_t other) { sources |= other; }

    // Check if feature exists in a specific sample
    bool in_sample(uint32_t sample_id) const {
        return sample_idx.test(sample_id);
    }

    // Check if feature has a specific source
    bool has_source(const std::string& source) const {
        return (sources & source_registry::instance().mask(source)) != 0;
    }

    // Get number of sources
    size_t source_count() const { return source_registry::count(sources); }

    // Get number of samples containing this feature
    size_t sample_count() const {
        return sample_idx.count();
    }

    // --- Expression methods ---

    void set_expression(uint32_t sample_id, float value) { expression.set(sample_id, value); }
    float get_expression(uint32_t sample_id) const { return expression.get(sample_id); }
    bool has_expression(uint32_t sample_id) const { return expression.has(sample_id); }

    // --- Pan-transcriptome statistics ---

    float frequency(size_t total_samples) const {
        if (total_samples == 0) return 0.0f;
        return static_cast<float>(sample_idx.count()) / static_cast<float>(total_samples);
    }

    bool is_conserved(size_t total_samples) const {
        return sample_idx.count() == total_samples;
    }

    bool is_sample_specific() const {
        return sample_idx.count() == 1;
    }
};

/**
 * Segment feature: spans consecutive exons, represents transcript paths
 * Used for graph traversal and transcript assignment
 */
struct segment_feature {
    // Core identifiers
    uint32_t gene_idx;               // Index into gene_registry (gene_id, gene_name, gene_biotype)

    // Gene accessors (resolve through gene_registry)
    const std::string& gene_id() const { return gene_registry::instance().resolve(gene_idx).gene_id; }
    const std::string& gene_name() const { return gene_registry::instance().resolve(gene_idx).gene_name; }
    const std::string& gene_biotype() const { return gene_registry::instance().resolve(gene_idx).gene_biotype; }

    // Transcript associations (interned via transcript_registry)
    sorted_vec transcript_ids;  // Transcripts using this segment

    // Source tracking (GFF column 2) — bitfield indexed by source_registry
    uint16_t sources = 0;

    // Sample tracking — dynamic bitset indexed by sample_registry IDs
    sample_bitset sample_idx;

    // Expression quantification (lazily allocated, -1.0f sentinel for missing)
    expression_store expression;

    // Transcript biotype mapping (interned transcript_id -> biotype)
    // e.g., protein_coding, retained_intron, nonsense_mediated_decay
    std::unordered_map<uint32_t, std::string> transcript_biotypes;

    // Segment structure
    size_t segment_index;            // Unique numeric index (used as edge ID for traversal)
    int exon_count;                  // Number of exons in this segment

    // Read support metadata (populated during discovery phase)
    size_t read_coverage;            // Number of reads supporting this segment
    std::vector<std::string> supporting_reads;  // Read IDs that support this segment

    // Absorption tracking (ISM segments absorbed into longer parents)
    bool absorbed = false;           // Tombstone: this segment was absorbed into a longer parent
    size_t absorbed_count = 0;       // Number of ISM segments absorbed into this one

    // Constructors
    segment_feature()
        : gene_idx(0), sources(0), segment_index(0), exon_count(0), read_coverage(0) {}

    // --- Sample tracking methods ---

    // Add sample reference (registry ID from sample_registry)
    void add_sample(uint32_t sample_id) {
        sample_idx.set(sample_id);
    }

    // Add source (GFF column 2 value like HAVANA, ENSEMBL, StringTie)
    void add_source(const std::string& source) {
        sources |= source_registry::instance().mask(source);
    }

    // Merge sources from another feature's bitfield
    void merge_sources(uint16_t other) { sources |= other; }

    // Check if feature exists in a specific sample
    bool in_sample(uint32_t sample_id) const {
        return sample_idx.test(sample_id);
    }

    // Check if feature has a specific source
    bool has_source(const std::string& source) const {
        return (sources & source_registry::instance().mask(source)) != 0;
    }

    // Get number of sources
    size_t source_count() const { return source_registry::count(sources); }

    // Get number of samples containing this feature
    size_t sample_count() const {
        return sample_idx.count();
    }

    // --- Expression methods ---

    void set_expression(uint32_t sample_id, float value) { expression.set(sample_id, value); }
    float get_expression(uint32_t sample_id) const { return expression.get(sample_id); }
    bool has_expression(uint32_t sample_id) const { return expression.has(sample_id); }

    // --- Read support methods ---

    // Add read support
    void add_read_support(const std::string& read_id) {
        supporting_reads.push_back(read_id);
        read_coverage++;
    }

    // --- Pan-transcriptome statistics ---

    // Calculate frequency across samples
    float frequency(size_t total_samples) const {
        if (total_samples == 0) return 0.0f;
        return static_cast<float>(sample_idx.count()) / static_cast<float>(total_samples);
    }

    // Check if feature is conserved (present in all samples)
    bool is_conserved(size_t total_samples) const {
        return sample_idx.count() == total_samples;
    }

    // Check if feature is sample-specific (only in one sample)
    bool is_sample_specific() const {
        return sample_idx.count() == 1;
    }
};

/**
 * Genomic feature variant: can be either exon or segment
 * Use std::holds_alternative<T> to check type
 * Use std::get<T> to access specific type
 */
using genomic_feature = std::variant<exon_feature, segment_feature>;

// Helper functions to check variant type
inline bool is_exon(const genomic_feature& feature) {
    return std::holds_alternative<exon_feature>(feature);
}

inline bool is_segment(const genomic_feature& feature) {
    return std::holds_alternative<segment_feature>(feature);
}

// Helper functions to get variant value
inline exon_feature& get_exon(genomic_feature& feature) {
    return std::get<exon_feature>(feature);
}

inline const exon_feature& get_exon(const genomic_feature& feature) {
    return std::get<exon_feature>(feature);
}

inline segment_feature& get_segment(genomic_feature& feature) {
    return std::get<segment_feature>(feature);
}

inline const segment_feature& get_segment(const genomic_feature& feature) {
    return std::get<segment_feature>(feature);
}

// ========== Grove type aliases ==========
// Canonical definitions used throughout the codebase

using grove_type = gst::grove<gdt::genomic_coordinate, genomic_feature, edge_metadata>;
using key_ptr = gdt::key<gdt::genomic_coordinate, genomic_feature>*;

// Cache types for deduplication across files (pan-transcriptome)
using exon_cache_type = std::map<gdt::genomic_coordinate, key_ptr>;
using segment_cache_type = std::unordered_map<std::string, key_ptr>;
using chromosome_exon_caches = std::map<std::string, exon_cache_type>;
using chromosome_segment_caches = std::map<std::string, segment_cache_type>;

// ISM absorption: per-gene index of segments with their exon chains
struct segment_chain_entry {
    key_ptr segment;
    std::vector<key_ptr> exon_chain;
    std::string structure_key;
};

using gene_segment_index_type = std::unordered_map<std::string, std::vector<segment_chain_entry>>;
using chromosome_gene_segment_indices = std::map<std::string, gene_segment_index_type>;

#endif //ATROPLEX_GENOMIC_FEATURE_HPP
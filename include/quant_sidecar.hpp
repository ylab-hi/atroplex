/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_QUANT_SIDECAR_HPP
#define ATROPLEX_QUANT_SIDECAR_HPP

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

/**
 * Per-sample quantification sidecar (.qtx) format and I/O.
 *
 * A .qtx file stores (segment_index, value) records for one sample,
 * sorted by segment_index. Values come from GFF/GTF attributes during
 * build (counts / TPM / FPKM / RPKM / cov, after the expression filter
 * passes) and are written to disk instead of being stored on features
 * in the grove. This keeps the in-memory grove purely structural —
 * membership lives on `sample_bitset`, quantification lives in `.qtx`.
 *
 * Layout:
 *   Header (32 bytes, little-endian, packed):
 *     char[4]  magic        "AQTX"
 *     u32      version      QTX_VERSION
 *     u32      count        number of records
 *     u8       expr_type    sample_info::expression_type code
 *     u8[3]    pad          zeros
 *     u64      grove_id     build-run identifier shared with .ggx header
 *     u64      reserved     zero
 *
 *   Records (count × 12 bytes, sorted by segment_index ascending):
 *     u64      segment_index   segment_feature::segment_index
 *     f32      value           the filter-passing attribute value
 *
 * Write path (build time): records are buffered in memory by a single
 * Writer per sample, sorted and flushed atomically at finalize() via
 * write-to-tmp + rename. Typical per-sample size is a few K to a few
 * tens of K records, so in-memory accumulation is cheap.
 *
 * Read path (query time): Reader loads the full file into memory. For
 * point lookups `lookup(segment_index)` binary-searches the sorted
 * records; for linear scans `for_each` iterates in segment_index order.
 *
 * The `grove_id` field should match the corresponding `.ggx` header so
 * the reader can refuse to operate on mismatched pairs (sidecar vs.
 * grove built in different runs, where segment_index values would be
 * stale). Populating it is the caller's responsibility — typically a
 * build-run identifier generated at the start of the build and passed
 * to both the grove serializer and every sidecar writer.
 */
namespace quant_sidecar {

constexpr char     MAGIC[4]    = {'A', 'Q', 'T', 'X'};
constexpr uint32_t QTX_VERSION = 1;

#pragma pack(push, 1)
struct Record {
    uint64_t segment_index;
    float    value;
};
#pragma pack(pop)
static_assert(sizeof(Record) == 12, "quant_sidecar::Record must be packed 12 bytes");

#pragma pack(push, 1)
struct Header {
    char     magic[4];
    uint32_t version;
    uint32_t count;
    uint8_t  expr_type;
    uint8_t  pad[3];
    uint64_t grove_id;
    uint64_t reserved;
};
#pragma pack(pop)
static_assert(sizeof(Header) == 32, "quant_sidecar::Header must be 32 bytes");

/**
 * Writer for building a .qtx file during the build phase.
 *
 * Typical usage (one writer per sample, owned by the builder):
 *
 *     Writer w(path, static_cast<uint8_t>(info.expr_type), build_id);
 *     for (every passing transcript in this sample) {
 *         w.append(segment_index, value);
 *     }
 *     w.finalize();  // sorts + flushes; or let destructor call it
 *
 * `append` is O(1) amortized. `finalize` sorts by segment_index and
 * writes the file atomically via tmp + rename. Calling `finalize` more
 * than once is a no-op; the destructor invokes it if it has not been
 * called already, but destructor-driven flushes swallow exceptions,
 * so explicit finalization is strongly recommended.
 *
 * Not thread-safe. Callers that want parallel builds must create one
 * Writer per sample and ensure each is touched by only one thread.
 */
class Writer {
public:
    Writer(std::filesystem::path path, uint8_t expr_type, uint64_t grove_id);
    Writer(const Writer&) = delete;
    Writer& operator=(const Writer&) = delete;
    Writer(Writer&&) noexcept = default;
    Writer& operator=(Writer&&) noexcept = default;
    ~Writer();

    /// Append one (segment_index, value) record. Ordering is not
    /// enforced — `finalize` sorts before writing.
    void append(uint64_t segment_index, float value);

    /// Current in-memory record count (mostly for diagnostics/tests).
    size_t size() const { return records_.size(); }

    /// True if no records have been appended.
    bool empty() const { return records_.empty(); }

    /// Sort and flush to disk atomically. Safe to call multiple times;
    /// subsequent calls are no-ops. Throws on I/O failure.
    void finalize();

    /// Final path the file will be written to (after finalize()).
    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
    uint8_t               expr_type_ = 0;
    uint64_t              grove_id_  = 0;
    std::vector<Record>   records_;
    bool                  finalized_ = false;
};

/**
 * Reader for .qtx files.
 *
 * Loads the full file into memory at construction — sidecar files are
 * small (fixed 12-byte records, typical sizes well under 1 MB), so
 * slurping is cheaper than mmap bookkeeping. Supports both point
 * lookups (binary search) and linear scans.
 *
 * Throws on magic/version mismatch or truncated input.
 */
class Reader {
public:
    /// Open and load a .qtx file. Throws std::runtime_error on any
    /// I/O or format failure.
    explicit Reader(const std::filesystem::path& path);

    /// Optionally enforce a grove_id match — useful when the caller
    /// knows the grove's build-run identifier and wants to reject
    /// mismatched sidecars up front. Returns false if the header's
    /// grove_id differs from `expected`.
    bool grove_id_matches(uint64_t expected) const {
        return header_.grove_id == expected;
    }

    const Header& header() const { return header_; }

    size_t size() const { return records_.size(); }
    bool   empty() const { return records_.empty(); }

    /// O(log N) point lookup by segment_index. Returns the associated
    /// value if present, or std::nullopt otherwise.
    std::optional<float> lookup(uint64_t segment_index) const;

    /// Linear scan over all records in segment_index order, calling
    /// fn(segment_index, value) for each.
    template <typename Fn>
    void for_each(Fn&& fn) const {
        for (const auto& r : records_) {
            fn(r.segment_index, r.value);
        }
    }

    /// Direct zero-copy access to the underlying records (still sorted).
    const std::vector<Record>& records() const { return records_; }

private:
    Header               header_{};
    std::vector<Record>  records_;
};

} // namespace quant_sidecar

#endif // ATROPLEX_QUANT_SIDECAR_HPP

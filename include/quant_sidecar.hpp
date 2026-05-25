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

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * Quantification sidecar — segment-major monolithic .qtx format.
 *
 * At rest, one .qtx file stores expression values for every sample in a
 * grove, grouped by segment. This replaces the old per-sample sidecar
 * design: instead of N small files (one per sample, sorted by
 * segment_index), we produce a single file with a per-segment block for
 * every segment that has at least one recorded value. Within a block,
 * records are sorted by sample_id ascending so that downstream readers
 * can efficiently intersect a "samples of interest" list.
 *
 * At build time, we still accumulate per-sample streams — each sample
 * gets its own temp file — and then run a K-way merge at end of build to
 * produce the final segment-major file. Per-sample temp streams are tiny
 * and don't require holding an N-sample × M-segment matrix in RAM.
 *
 * On-disk layout (final .qtx, little-endian, all POD structs are packed):
 *
 *   Header (40 bytes):
 *     char[4]  magic                 "AQTX"
 *     u32      version               QTX_VERSION (2)
 *     u64      segment_block_count   number of segments with >=1 record
 *     u64      toc_offset            file offset where the TOC starts
 *     u64      sample_meta_offset    file offset where sample metadata starts
 *     u64      reserved              0
 *
 *   Sample metadata section (starts at sample_meta_offset):
 *     u32 num_samples
 *     For each sample:
 *       u32  sample_id
 *       u8   types_mask        (bitfield: bit N set ↔ sample has ≥1 record
 *                                of expression_type N)
 *       u8   name_len
 *       char[name_len] name
 *
 *   Segment blocks (one per segment with >=1 record, in segment_index order):
 *     u64 segment_index
 *     u32 num_records         (count of packed records, ie. distinct
 *                              sample_ids that have at least one value
 *                              at this segment)
 *     For each record (sorted by sample_id ascending):
 *       u32 sample_id
 *       u8  type_mask          (bit N set ↔ this (seg, sample) carries a
 *                                value of expression_type N; see
 *                                expression_type_label())
 *       f32 values[popcount(type_mask)]  (values in bit-position order:
 *                                lower bits first)
 *
 *   Segment TOC (at toc_offset, sorted by segment_index ascending):
 *     For each block:
 *       u64 segment_index
 *       u64 offset            (absolute file offset of the block)
 *       u32 length            (total bytes of the block, including the
 *                              u64+u32 segment_index/num_records header)
 *
 * Build-time per-sample temp stream format (internal, not user-facing):
 *
 *   Header (16 bytes):
 *     char[4] magic "AQTS"
 *     u32     version  2
 *     u32     sample_id       (representative; 0 for batched intermediates)
 *     u32     num_records
 *
 *   Records (sorted by (segment_index, sample_id, type)):
 *     u64 segment_index
 *     u32 sample_id           (authoritative — same format used for both
 *                              per-sample streams and batched intermediates)
 *     u8  type                (expression_type enum byte)
 *     u8  pad[3]
 *     f32 value
 *
 * Each record is 20 bytes. The sample_id in the header is redundant for
 * single-sample streams but simplifies the merge code since intermediate
 * batched files have mixed sample_ids and the header slot becomes 0.
 *
 * Multi-value storage: one sample × one segment may carry up to N records,
 * one per expression_type the source GFF row carried (e.g. TPM + FPKM +
 * counts). Filtering at build time is CLI-driven (--min-X / --filter-
 * precedence), but the index stores everything the row provided so callers
 * can re-query in any unit later.
 */
namespace quant_sidecar {

// ── Constants ──────────────────────────────────────────────────────────

constexpr char     MAGIC[4]    = {'A', 'Q', 'T', 'X'};
constexpr uint32_t QTX_VERSION = 3;

// Temp-stream magic and version (used by SampleStreamWriter and the
// intermediate files produced by batched merge passes).
constexpr char     STREAM_MAGIC[4]    = {'A', 'Q', 'T', 'S'};
constexpr uint32_t STREAM_VERSION     = 2;

// ── Packed POD structs ────────────────────────────────────────────────

#pragma pack(push, 1)

/// Final .qtx file header (40 bytes).
struct Header {
    char     magic[4];
    uint32_t version;
    uint64_t segment_block_count;
    uint64_t toc_offset;
    uint64_t sample_meta_offset;
    uint64_t reserved;
};

/// Table-of-contents entry in the final .qtx (20 bytes).
struct TOCEntry {
    uint64_t segment_index;
    uint64_t offset;
    uint32_t length;
};

/// Temp-stream header for per-sample temp files and batched intermediates.
struct StreamHeader {
    char     magic[4];
    uint32_t version;
    uint32_t sample_id;    // representative; 0 for batched intermediates
    uint32_t num_records;
};

/// On-disk record in the temp-stream format (20 bytes). The `type` byte
/// carries an `expression_type` enum code (see expression_type_label()).
/// `_pad` is reserved and must be written as zero.
struct StreamRecord {
    uint64_t segment_index;
    uint32_t sample_id;
    uint8_t  type;
    uint8_t  _pad[3];
    float    value;
};

#pragma pack(pop)

static_assert(sizeof(Header) == 40,        "quant_sidecar::Header must be 40 bytes");
static_assert(sizeof(TOCEntry) == 20,      "quant_sidecar::TOCEntry must be 20 bytes");
static_assert(sizeof(StreamHeader) == 16,  "quant_sidecar::StreamHeader must be 16 bytes");
static_assert(sizeof(StreamRecord) == 20,  "quant_sidecar::StreamRecord must be 20 bytes");

// ── Expression type labels ─────────────────────────────────────────────

/**
 * Return the canonical string label for an expression type byte. The
 * type byte values match the `sample_info::expression_type` enum:
 *   0=UNKNOWN, 1=TPM, 2=FPKM, 3=RPKM, 4=counts, 5=CPM, 6=cov.
 * Used to build TSV column suffixes (e.g. "sampleA.TPM") and to format
 * diagnostics. Returns "unknown" for any unrecognized byte.
 */
inline const char* expression_type_label(uint8_t type_byte) {
    switch (type_byte) {
        case 1: return "TPM";
        case 2: return "FPKM";
        case 3: return "RPKM";
        case 4: return "counts";
        case 5: return "CPM";
        case 6: return "cov";
        default: return "unknown";
    }
}

// ── Sample metadata ────────────────────────────────────────────────────

/**
 * Describes one sample in a .qtx file. Kept small and plain so it can
 * be moved/copied freely between build and merge steps.
 *
 * `types_mask` is a bitfield indicating which expression_type codes the
 * sample carries at least one record for: bit N set ↔ at least one
 * record with type=N exists for this sample. Populated by merge_to_qtx
 * as records flow through; zero before merge.
 */
struct SampleMetadata {
    uint32_t    sample_id;
    uint8_t     types_mask;  ///< bit N set ↔ sample has ≥1 record of type N
    std::string name;        ///< the sample's info.id string
};

// ── SampleStreamWriter ─────────────────────────────────────────────────

/**
 * Per-sample temp stream writer.
 *
 * Used during grove build. Accumulates (segment_index, type, value)
 * records in memory; at finalize() it sorts by (segment_index, type)
 * and writes them to disk atomically (tmp + rename). The resulting
 * file carries its sample_id in the header and in every record, so it
 * can be fed verbatim into the K-way merge.
 *
 * Typical usage (one writer per sample, owned by the builder):
 *
 *     SampleStreamWriter w(tmp_path, sample_id);
 *     for (every passing transcript in this sample) {
 *         for (auto [type, value] : transcript_expression_values) {
 *             w.append(segment_index, type, value);
 *         }
 *     }
 *     w.finalize();  // or let the destructor do it (exceptions swallowed)
 *
 * Not thread-safe. Callers that want parallel builds must create one
 * writer per sample and ensure each is touched by only one thread.
 */
class SampleStreamWriter {
public:
    SampleStreamWriter(std::filesystem::path path, uint32_t sample_id);
    SampleStreamWriter(const SampleStreamWriter&) = delete;
    SampleStreamWriter& operator=(const SampleStreamWriter&) = delete;
    SampleStreamWriter(SampleStreamWriter&&) noexcept = default;
    SampleStreamWriter& operator=(SampleStreamWriter&&) noexcept = default;
    ~SampleStreamWriter();

    /// Append one (segment_index, type, value) record. Ordering is not
    /// enforced — finalize() sorts before writing. `type` is an
    /// expression_type enum byte (see expression_type_label()).
    void append(uint64_t segment_index, uint8_t type, float value);

    size_t size() const { return records_.size(); }
    bool   empty() const { return records_.empty(); }
    uint32_t sample_id() const { return sample_id_; }
    const std::filesystem::path& path() const { return path_; }

    /// Sort and flush to disk atomically. Idempotent.
    void finalize();

private:
    std::filesystem::path     path_;
    uint32_t                  sample_id_ = 0;
    std::vector<StreamRecord> records_;
    bool                      finalized_ = false;
};

// ── K-way merge into final .qtx ────────────────────────────────────────

/**
 * K-way merge per-sample stream files into a segment-major monolithic .qtx.
 *
 * Each stream file must have been produced by SampleStreamWriter (or by
 * a prior batched-merge pass) and contain records sorted by
 * (segment_index, sample_id). In the output, records are grouped into
 * per-segment blocks, sorted first by segment_index then by sample_id
 * within each block.
 *
 * To stay under `ulimit -n`, the merge is done in batches: at most
 * `max_fds_per_pass` stream files are open simultaneously. If the
 * number of streams exceeds that, an intermediate pass produces
 * batched output files (using the same temp-stream format) that are
 * then merged in a second pass. Batched intermediates are written
 * next to the final output (same parent directory) with a
 * `.qtx.batch.<i>` suffix and are cleaned up on success.
 *
 * The TOC is streamed to disk progressively — entries are appended to
 * a temp TOC file as blocks are written. At end of merge, the TOC is
 * concatenated to the main file and the header is patched with the
 * correct offsets and block count.
 *
 * Atomic: the output is written to `output_path + ".tmp"` and renamed
 * at end. Caller is responsible for deleting the input stream files.
 *
 * @param output_path       final .qtx destination
 * @param streams           per-sample stream file paths (must exist, sorted)
 * @param samples           sample metadata matching the streams (same length,
 *                          same order)
 * @param max_fds_per_pass  max simultaneous open streams (default 256)
 * @param excluded_segments Optional set of segment_index values whose
 *                          records should be omitted from the output.
 *                          Used by callers that know about tombstoned
 *                          segments so they don't end up as dead
 *                          records in the final file. Empty set means
 *                          include all records.
 * @param segment_remap     Optional map {tombstoned_idx → live_parent_idx}.
 *                          When a record's segment_index is in this map,
 *                          its index is rewritten to the parent and the
 *                          record is re-inserted into the merge heap so
 *                          it sorts into the parent's block. Entries
 *                          override `excluded_segments` — a record in
 *                          the remap is never dropped, even if its
 *                          source index also appears in the drop set.
 *                          The caller must have resolved transitive
 *                          chains (path compression) so every value is
 *                          a live segment_index. Issue #34.
 */
void merge_to_qtx(
    const std::filesystem::path& output_path,
    const std::vector<std::filesystem::path>& streams,
    const std::vector<SampleMetadata>& samples,
    size_t max_fds_per_pass = 256,
    const std::unordered_set<uint64_t>& excluded_segments = {},
    const std::unordered_map<uint64_t, uint64_t>& segment_remap = {});

// ── Reader ────────────────────────────────────────────────────────────

/**
 * Reader for segment-major .qtx files.
 *
 * Opens the file once, loads the TOC and sample metadata into memory
 * (both are small — a few MB at cohort scale) and caches the file
 * handle for on-demand block reads. Point lookups are O(log N) in the
 * number of segments via binary search on the TOC, plus one block
 * read.
 *
 * Throws std::runtime_error on any I/O or format failure.
 */
class Reader {
public:
    /// In-memory record carrying every dimension of a stored value.
    /// Returned by lookup_all() and for_each_segment(). Records within a
    /// block are sorted by (sample_id, type) ascending.
    struct TypedValueRecord {
        uint32_t sample_id;
        uint8_t  type;
        float    value;
    };

    /// In-memory record from a single-type lookup. The type axis is
    /// implicit (caller passed it as an argument).
    struct ValueRecord {
        uint32_t sample_id;
        float    value;
    };

    explicit Reader(const std::filesystem::path& path);

    /// Look up every (sample, type, value) record for a segment. Returns
    /// an empty vector when the segment has no records. Records are
    /// sorted by (sample_id, type) ascending.
    [[nodiscard]] std::vector<TypedValueRecord> lookup_all(uint64_t segment_index) const;

    /// Look up every record matching a single expression type for a
    /// segment. Returns an empty vector when no record matches. Records
    /// are sorted by sample_id ascending.
    [[nodiscard]] std::vector<ValueRecord> lookup(uint64_t segment_index, uint8_t type) const;

    /// Look up a segment, filter to a single expression type, and
    /// intersect with a sample-id allowlist. Returns a map of sample_id
    /// → value for samples that appear in both the segment's block (with
    /// the requested type) and `wanted_samples`.
    [[nodiscard]] std::unordered_map<uint32_t, float> lookup_filtered(
        uint64_t segment_index,
        uint8_t type,
        const std::vector<uint32_t>& wanted_samples) const;

    /// Bitmask of expression types this sample carries at least one
    /// record for. Returns 0 if the sample is unknown to this reader.
    /// Driven by `SampleMetadata::types_mask` (populated at merge time).
    [[nodiscard]] uint8_t types_for_sample(uint32_t sample_id) const;

    /// Iterate over every segment block in segment_index order. `fn` is
    /// called with (segment_index, vector<TypedValueRecord>&) for each
    /// block. Intended for streaming analyses that want a single linear
    /// pass through the file.
    template <typename Fn>
    void for_each_segment(Fn&& fn) const {
        std::vector<TypedValueRecord> buf;
        for (const auto& entry : toc_) {
            read_block(entry, buf);
            fn(entry.segment_index, buf);
        }
    }

    const std::vector<SampleMetadata>& samples() const { return samples_; }
    size_t   segment_block_count() const { return toc_.size(); }
    uint32_t version() const { return header_.version; }
    const std::filesystem::path& path() const { return path_; }

private:
    /// Read the block described by `entry` into `out` (cleared first).
    /// Factored out so lookup() and for_each_segment() share the same
    /// on-disk parsing logic.
    void read_block(const TOCEntry& entry,
                    std::vector<TypedValueRecord>& out) const;

    std::filesystem::path       path_;
    mutable std::ifstream       file_;    // mutable: lookup() is const
    Header                      header_{};
    std::vector<SampleMetadata> samples_;
    std::vector<TOCEntry>       toc_;     // sorted by segment_index
};

} // namespace quant_sidecar

#endif // ATROPLEX_QUANT_SIDECAR_HPP
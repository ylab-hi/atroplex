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
 *       u8   expr_type
 *       u8   name_len
 *       char[name_len] name
 *
 *   Segment blocks (one per segment with >=1 record, in segment_index order):
 *     u64 segment_index
 *     u32 num_records
 *     For each record (sorted by sample_id ascending):
 *       u32 sample_id
 *       f32 value
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
 *     u32     version  1
 *     u32     sample_id       (representative; 0 for batched intermediates)
 *     u32     num_records
 *
 *   Records (sorted by (segment_index, sample_id)):
 *     u64 segment_index
 *     u32 sample_id           (authoritative — same format used for both
 *                              per-sample streams and batched intermediates)
 *     f32 value
 *
 * Each record is 16 bytes. The sample_id in the header is redundant for
 * single-sample streams but simplifies the merge code since intermediate
 * batched files have mixed sample_ids and the header slot becomes 0.
 */
namespace quant_sidecar {

// ── Constants ──────────────────────────────────────────────────────────

constexpr char     MAGIC[4]    = {'A', 'Q', 'T', 'X'};
constexpr uint32_t QTX_VERSION = 2;

// Temp-stream magic and version (used by SampleStreamWriter and the
// intermediate files produced by batched merge passes).
constexpr char     STREAM_MAGIC[4]    = {'A', 'Q', 'T', 'S'};
constexpr uint32_t STREAM_VERSION     = 1;

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

/// On-disk record in the temp-stream format (16 bytes).
struct StreamRecord {
    uint64_t segment_index;
    uint32_t sample_id;
    float    value;
};

#pragma pack(pop)

static_assert(sizeof(Header) == 40,        "quant_sidecar::Header must be 40 bytes");
static_assert(sizeof(TOCEntry) == 20,      "quant_sidecar::TOCEntry must be 20 bytes");
static_assert(sizeof(StreamHeader) == 16,  "quant_sidecar::StreamHeader must be 16 bytes");
static_assert(sizeof(StreamRecord) == 16,  "quant_sidecar::StreamRecord must be 16 bytes");

// ── Sample metadata ────────────────────────────────────────────────────

/**
 * Describes one sample in a .qtx file. Kept small and plain so it can
 * be moved/copied freely between build and merge steps.
 */
struct SampleMetadata {
    uint32_t    sample_id;
    uint8_t     expr_type;   ///< sample_info::expression_type code
    std::string name;        ///< the sample's info.id string
};

// ── SampleStreamWriter ─────────────────────────────────────────────────

/**
 * Per-sample temp stream writer.
 *
 * Used during grove build. Accumulates (segment_index, value) records in
 * memory; at finalize() it sorts by segment_index and writes them to
 * disk atomically (tmp + rename). The resulting file carries its
 * sample_id in the header and in every record, so it can be fed
 * verbatim into the K-way merge.
 *
 * Typical usage (one writer per sample, owned by the builder):
 *
 *     SampleStreamWriter w(tmp_path, sample_id);
 *     for (every passing transcript in this sample) {
 *         w.append(segment_index, value);
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

    /// Append one (segment_index, value) record. Ordering is not
    /// enforced — finalize() sorts before writing.
    void append(uint64_t segment_index, float value);

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
 */
void merge_to_qtx(
    const std::filesystem::path& output_path,
    const std::vector<std::filesystem::path>& streams,
    const std::vector<SampleMetadata>& samples,
    size_t max_fds_per_pass = 256);

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
    struct ValueRecord {
        uint32_t sample_id;
        float    value;
    };

    explicit Reader(const std::filesystem::path& path);

    /// Look up all (sample_id, value) records for a segment. Returns an
    /// empty vector if the segment has no records in the file.
    std::vector<ValueRecord> lookup(uint64_t segment_index) const;

    /// Look up a segment and filter to a set of samples of interest.
    /// Returns a map of sample_id → value for samples that appear in
    /// both the segment's block and `wanted_samples`. Efficient because
    /// records are already sorted by sample_id within a block.
    std::unordered_map<uint32_t, float> lookup_filtered(
        uint64_t segment_index,
        const std::vector<uint32_t>& wanted_samples) const;

    /// Iterate over every segment block in segment_index order.
    /// `fn` is called with (segment_index, vector<ValueRecord>&) for
    /// each block. Intended for streaming analyses that want a single
    /// linear pass through the file.
    template <typename Fn>
    void for_each_segment(Fn&& fn) const {
        std::vector<ValueRecord> buf;
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
                    std::vector<ValueRecord>& out) const;

    std::filesystem::path       path_;
    mutable std::ifstream       file_;    // mutable: lookup() is const
    Header                      header_{};
    std::vector<SampleMetadata> samples_;
    std::vector<TOCEntry>       toc_;     // sorted by segment_index
};

} // namespace quant_sidecar

#endif // ATROPLEX_QUANT_SIDECAR_HPP
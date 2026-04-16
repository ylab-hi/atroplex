/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "quant_sidecar.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <memory>
#include <queue>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

namespace quant_sidecar {

namespace {

// ── Small I/O helpers ──────────────────────────────────────────────────
//
// These wrap std::istream/std::ostream binary read/write with
// explicit size checks and descriptive error messages. We keep them
// local to this TU so they don't leak into the public header.

/// Throws if the read didn't return the expected number of bytes.
void read_exact(std::istream& in, void* dst, std::streamsize n,
                const std::filesystem::path& path, const char* what) {
    in.read(reinterpret_cast<char*>(dst), n);
    if (!in || in.gcount() != n) {
        throw std::runtime_error(
            std::string("quant_sidecar: truncated ") + what + " in " +
            path.string());
    }
}

/// Throws if the write fails.
void write_exact(std::ostream& out, const void* src, std::streamsize n,
                 const std::filesystem::path& path, const char* what) {
    out.write(reinterpret_cast<const char*>(src), n);
    if (!out) {
        throw std::runtime_error(
            std::string("quant_sidecar: failed to write ") + what + " to " +
            path.string());
    }
}

/// Append `.tmp` to a path so it can be renamed in place atomically.
std::filesystem::path tmp_suffixed(const std::filesystem::path& path) {
    std::filesystem::path out = path;
    out += ".tmp";
    return out;
}

// ── Temp-stream readers used by the merge algorithm ────────────────────

/**
 * Streaming reader for a single temp stream file. Reads records one at
 * a time on demand — we never materialize a full stream in memory
 * during merge, which is what lets the merge scale regardless of how
 * much data any individual stream holds.
 */
class StreamReader {
public:
    StreamReader(std::filesystem::path path)
        : path_(std::move(path)) {
        in_.open(path_, std::ios::binary);
        if (!in_) {
            throw std::runtime_error(
                "quant_sidecar: cannot open stream file " + path_.string());
        }

        StreamHeader hdr{};
        read_exact(in_, &hdr, sizeof(hdr), path_, "stream header");

        if (std::memcmp(hdr.magic, STREAM_MAGIC, sizeof(STREAM_MAGIC)) != 0) {
            throw std::runtime_error(
                "quant_sidecar: bad magic in stream file " + path_.string());
        }
        if (hdr.version != STREAM_VERSION) {
            throw std::runtime_error(
                "quant_sidecar: unsupported stream version " +
                std::to_string(hdr.version) + " in " + path_.string() +
                " (expected " + std::to_string(STREAM_VERSION) + ")");
        }

        remaining_ = hdr.num_records;
        // Prime the reader so current() / has_current() are meaningful
        // immediately after construction.
        advance();
    }

    bool has_current() const { return has_current_; }
    const StreamRecord& current() const { return current_; }
    const std::filesystem::path& path() const { return path_; }

    /// Step to the next record. After this returns, has_current() is
    /// true iff another record was available.
    void advance() {
        if (remaining_ == 0) {
            has_current_ = false;
            return;
        }
        read_exact(in_, &current_, sizeof(current_), path_, "stream record");
        --remaining_;
        has_current_ = true;
    }

private:
    std::filesystem::path path_;
    std::ifstream         in_;
    StreamRecord          current_{};
    bool                  has_current_ = false;
    uint32_t              remaining_   = 0;
};

/**
 * Streaming writer for a temp stream file. Used to produce
 * intermediate batched files during the merge. Records are written
 * as-is; the caller is responsible for pushing them in
 * (segment_index, sample_id) order. The header is patched at close()
 * with the final record count.
 */
class BatchWriter {
public:
    BatchWriter(std::filesystem::path path, uint32_t header_sample_id)
        : path_(std::move(path)), header_sample_id_(header_sample_id) {
        out_.open(path_, std::ios::binary | std::ios::trunc);
        if (!out_) {
            throw std::runtime_error(
                "quant_sidecar: cannot open intermediate stream " +
                path_.string() + " for write");
        }

        // Reserve header space with a placeholder; patched on close().
        StreamHeader hdr{};
        std::memcpy(hdr.magic, STREAM_MAGIC, sizeof(STREAM_MAGIC));
        hdr.version     = STREAM_VERSION;
        hdr.sample_id   = header_sample_id_;
        hdr.num_records = 0;
        write_exact(out_, &hdr, sizeof(hdr), path_, "placeholder header");
    }

    void append(const StreamRecord& rec) {
        write_exact(out_, &rec, sizeof(rec), path_, "stream record");
        ++count_;
    }

    /// Flushes, patches the header's record count, and closes.
    void close() {
        if (closed_) return;
        closed_ = true;

        out_.flush();
        if (!out_) {
            throw std::runtime_error(
                "quant_sidecar: failed to flush intermediate stream " +
                path_.string());
        }

        // Patch num_records in place. Seek to the num_records field
        // which lives at offset (sizeof(StreamHeader) - sizeof(u32)).
        out_.seekp(static_cast<std::streamoff>(offsetof(StreamHeader,
                                                         num_records)));
        if (!out_) {
            throw std::runtime_error(
                "quant_sidecar: seekp failed while patching header for " +
                path_.string());
        }
        write_exact(out_, &count_, sizeof(count_), path_, "num_records patch");

        out_.close();
        if (!out_) {
            throw std::runtime_error(
                "quant_sidecar: close failed for " + path_.string());
        }
    }

    const std::filesystem::path& path() const { return path_; }

    ~BatchWriter() {
        try {
            close();
        } catch (...) {
            // intentionally swallowed — destructor-driven cleanup
        }
    }

private:
    std::filesystem::path path_;
    std::ofstream         out_;
    uint32_t              header_sample_id_ = 0;
    uint32_t              count_    = 0;
    bool                  closed_   = false;
};

// ── K-way merge heap node ──────────────────────────────────────────────

/// Min-heap entry keyed by (segment_index, sample_id, stream_index).
/// The stream_index tiebreaker keeps the ordering total and
/// deterministic even if two streams somehow carry the same sample_id.
struct HeapNode {
    uint64_t segment_index;
    uint32_t sample_id;
    float    value;
    size_t   stream_index;

    bool operator>(const HeapNode& other) const {
        if (segment_index != other.segment_index)
            return segment_index > other.segment_index;
        if (sample_id != other.sample_id)
            return sample_id > other.sample_id;
        return stream_index > other.stream_index;
    }
};

/// Pull the next record from a given reader and push it onto the heap.
/// No-op if the reader has no more records.
void push_from_reader(std::priority_queue<HeapNode,
                                          std::vector<HeapNode>,
                                          std::greater<HeapNode>>& heap,
                      StreamReader& reader,
                      size_t stream_index) {
    if (!reader.has_current()) return;
    const auto& r = reader.current();
    heap.push(HeapNode{r.segment_index, r.sample_id, r.value, stream_index});
    reader.advance();
}

// ── Batched merge → intermediate temp-stream file ─────────────────────

/**
 * Merge a batch of stream files into a single intermediate temp-stream
 * file. Records are emitted in (segment_index, sample_id) order.
 * Used when the total stream count exceeds max_fds_per_pass.
 */
void merge_batch_to_stream(const std::vector<std::filesystem::path>& inputs,
                           const std::filesystem::path& output_path) {
    std::vector<std::unique_ptr<StreamReader>> readers;
    readers.reserve(inputs.size());
    for (const auto& p : inputs) {
        readers.emplace_back(std::make_unique<StreamReader>(p));
    }

    std::priority_queue<HeapNode, std::vector<HeapNode>,
                        std::greater<HeapNode>> heap;
    for (size_t i = 0; i < readers.size(); ++i) {
        push_from_reader(heap, *readers[i], i);
    }

    // header_sample_id = 0 for intermediates (mixed sample_ids).
    BatchWriter out(output_path, /*header_sample_id=*/0);
    while (!heap.empty()) {
        HeapNode node = heap.top();
        heap.pop();
        out.append(StreamRecord{node.segment_index, node.sample_id, node.value});
        push_from_reader(heap, *readers[node.stream_index], node.stream_index);
    }
    out.close();
}

// ── Final-pass merge → segment-major .qtx ─────────────────────────────

/**
 * Write the sample metadata section to `out` at the current position.
 * Returns the offset at which the section started (to be stored in the
 * header).
 */
uint64_t write_sample_metadata(std::ofstream& out,
                               const std::vector<SampleMetadata>& samples,
                               const std::filesystem::path& path) {
    const auto offset = static_cast<uint64_t>(out.tellp());

    const uint32_t num_samples = static_cast<uint32_t>(samples.size());
    write_exact(out, &num_samples, sizeof(num_samples), path, "num_samples");

    for (const auto& s : samples) {
        write_exact(out, &s.sample_id, sizeof(s.sample_id), path, "sample_id");
        write_exact(out, &s.expr_type, sizeof(s.expr_type), path, "expr_type");

        if (s.name.size() > 255) {
            throw std::runtime_error(
                "quant_sidecar: sample name exceeds 255 bytes: " + s.name);
        }
        const uint8_t name_len = static_cast<uint8_t>(s.name.size());
        write_exact(out, &name_len, sizeof(name_len), path, "name_len");
        if (name_len > 0) {
            write_exact(out, s.name.data(),
                        static_cast<std::streamsize>(name_len),
                        path, "sample name");
        }
    }

    return offset;
}

/**
 * Drain the K-way merge and write out one segment-block per distinct
 * segment_index as we go. Returns a vector of TOCEntry records
 * (ordered by segment_index, which is the natural output order).
 *
 * We buffer per-segment records in a small vector so we know how many
 * records there are before writing the block header. For a normal
 * workload this buffer stays tiny (O(num_samples with a value at that
 * segment)) so there's no memory concern.
 */
std::vector<TOCEntry> run_final_merge(
    std::ofstream& out,
    const std::filesystem::path& path,
    std::vector<std::unique_ptr<StreamReader>>& readers,
    const std::unordered_set<uint64_t>& excluded_segments) {
    std::priority_queue<HeapNode, std::vector<HeapNode>,
                        std::greater<HeapNode>> heap;
    for (size_t i = 0; i < readers.size(); ++i) {
        push_from_reader(heap, *readers[i], i);
    }

    std::vector<TOCEntry> toc;
    std::vector<std::pair<uint32_t, float>> block_records;  // (sample_id, value)

    // Helper: flush the current block_records buffer to disk and
    // append a matching TOCEntry. `current_seg` is the segment_index
    // we've been accumulating for.
    auto flush_block = [&](uint64_t current_seg) {
        if (block_records.empty()) return;
        const auto block_offset = static_cast<uint64_t>(out.tellp());

        const uint32_t num_records =
            static_cast<uint32_t>(block_records.size());

        write_exact(out, &current_seg, sizeof(current_seg), path,
                    "block segment_index");
        write_exact(out, &num_records, sizeof(num_records), path,
                    "block num_records");
        for (const auto& [sid, val] : block_records) {
            write_exact(out, &sid, sizeof(sid), path, "block sample_id");
            write_exact(out, &val, sizeof(val), path, "block value");
        }

        const auto block_end = static_cast<uint64_t>(out.tellp());
        TOCEntry entry{};
        entry.segment_index = current_seg;
        entry.offset        = block_offset;
        entry.length        = static_cast<uint32_t>(block_end - block_offset);
        toc.push_back(entry);

        block_records.clear();
    };

    bool     have_current = false;
    uint64_t current_seg  = 0;

    while (!heap.empty()) {
        HeapNode node = heap.top();
        heap.pop();

        // Advance the stream this node came from immediately so the
        // heap always has the next record queued — this lets the skip
        // branch below `continue` without stalling the merge.
        push_from_reader(heap, *readers[node.stream_index], node.stream_index);

        // Skip records for excluded segments (e.g., tombstoned by
        // reverse absorption). We still advanced the reader above so
        // the merge progresses normally.
        if (!excluded_segments.empty() &&
            excluded_segments.count(node.segment_index) > 0) {
            continue;
        }

        if (!have_current) {
            current_seg  = node.segment_index;
            have_current = true;
        } else if (node.segment_index != current_seg) {
            flush_block(current_seg);
            current_seg = node.segment_index;
        }
        block_records.emplace_back(node.sample_id, node.value);
    }
    if (have_current) {
        flush_block(current_seg);
    }

    return toc;
}

} // namespace

// ── SampleStreamWriter ─────────────────────────────────────────────────

SampleStreamWriter::SampleStreamWriter(std::filesystem::path path,
                                       uint32_t sample_id)
    : path_(std::move(path)), sample_id_(sample_id) {
    // Small initial reservation; vector doubling handles any tail.
    records_.reserve(8 * 1024);
}

SampleStreamWriter::~SampleStreamWriter() {
    if (!finalized_) {
        try {
            finalize();
        } catch (...) {
            // Destructor-driven flush must not throw.
        }
    }
}

void SampleStreamWriter::append(uint64_t segment_index, float value) {
    records_.push_back(StreamRecord{segment_index, sample_id_, value});
}

void SampleStreamWriter::finalize() {
    if (finalized_) return;
    finalized_ = true;

    // Sort by segment_index; sample_id is constant within a single-sample
    // stream, so it doesn't affect ordering. We use stable_sort so that
    // if duplicate segment_index values ever appear, insertion order wins
    // deterministically rather than being implementation-defined.
    std::stable_sort(records_.begin(), records_.end(),
                     [](const StreamRecord& a, const StreamRecord& b) {
                         return a.segment_index < b.segment_index;
                     });

    const auto tmp_path = tmp_suffixed(path_);

    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error(
            "quant_sidecar::SampleStreamWriter: cannot open " +
            tmp_path.string() + " for write");
    }

    StreamHeader hdr{};
    std::memcpy(hdr.magic, STREAM_MAGIC, sizeof(STREAM_MAGIC));
    hdr.version     = STREAM_VERSION;
    hdr.sample_id   = sample_id_;
    hdr.num_records = static_cast<uint32_t>(records_.size());
    write_exact(out, &hdr, sizeof(hdr), tmp_path, "stream header");

    if (!records_.empty()) {
        const auto bytes = static_cast<std::streamsize>(
            records_.size() * sizeof(StreamRecord));
        write_exact(out, records_.data(), bytes, tmp_path, "stream records");
    }

    out.flush();
    out.close();
    if (!out) {
        throw std::runtime_error(
            "quant_sidecar::SampleStreamWriter: flush/close failed for " +
            tmp_path.string());
    }

    std::error_code ec;
    std::filesystem::rename(tmp_path, path_, ec);
    if (ec) {
        throw std::runtime_error(
            "quant_sidecar::SampleStreamWriter: rename " +
            tmp_path.string() + " -> " + path_.string() + " failed: " +
            ec.message());
    }

    // Drop in-memory state.
    records_.clear();
    records_.shrink_to_fit();
}

// ── merge_to_qtx ───────────────────────────────────────────────────────

void merge_to_qtx(
    const std::filesystem::path& output_path,
    const std::vector<std::filesystem::path>& streams,
    const std::vector<SampleMetadata>& samples,
    size_t max_fds_per_pass,
    const std::unordered_set<uint64_t>& excluded_segments) {
    if (streams.size() != samples.size()) {
        throw std::runtime_error(
            "quant_sidecar::merge_to_qtx: streams.size() (" +
            std::to_string(streams.size()) + ") != samples.size() (" +
            std::to_string(samples.size()) + ")");
    }
    if (max_fds_per_pass < 2) {
        throw std::runtime_error(
            "quant_sidecar::merge_to_qtx: max_fds_per_pass must be >= 2");
    }

    // ── Step 1: if we have more streams than max_fds_per_pass, do
    //    a batched pre-merge pass into intermediate files. Each
    //    intermediate still has the temp-stream format, so the final
    //    pass doesn't care whether its inputs are per-sample streams
    //    or batched intermediates.
    //
    //    We cascade batches if the number of intermediates is still
    //    too large for one final pass. In practice one or two
    //    cascaded passes is plenty: 256^2 = 65K streams.

    std::vector<std::filesystem::path> current_inputs = streams;
    std::vector<std::filesystem::path> owned_intermediates;  // to clean up

    const auto parent_dir = output_path.parent_path();
    size_t batch_counter = 0;

    while (current_inputs.size() > max_fds_per_pass) {
        std::vector<std::filesystem::path> next_round;
        for (size_t i = 0; i < current_inputs.size(); i += max_fds_per_pass) {
            const size_t end = std::min(i + max_fds_per_pass,
                                        current_inputs.size());
            std::vector<std::filesystem::path> batch(
                current_inputs.begin() + static_cast<std::ptrdiff_t>(i),
                current_inputs.begin() + static_cast<std::ptrdiff_t>(end));

            std::filesystem::path batch_path =
                parent_dir / (output_path.filename().string() +
                              ".batch." + std::to_string(batch_counter++));

            merge_batch_to_stream(batch, batch_path);
            owned_intermediates.push_back(batch_path);
            next_round.push_back(batch_path);
        }
        current_inputs = std::move(next_round);
    }

    // ── Step 2: open all final-pass readers. These are at most
    //    max_fds_per_pass, so we're under ulimit by design.

    std::vector<std::unique_ptr<StreamReader>> readers;
    readers.reserve(current_inputs.size());
    for (const auto& p : current_inputs) {
        readers.emplace_back(std::make_unique<StreamReader>(p));
    }

    // ── Step 3: open the output file. Reserve space for the header
    //    by writing a zeroed placeholder; we'll seek back and patch
    //    it once we know the final offsets.

    const auto tmp_path = tmp_suffixed(output_path);
    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error(
            "quant_sidecar::merge_to_qtx: cannot open " + tmp_path.string() +
            " for write");
    }

    Header placeholder{};
    std::memcpy(placeholder.magic, MAGIC, sizeof(MAGIC));
    placeholder.version             = QTX_VERSION;
    placeholder.segment_block_count = 0;
    placeholder.toc_offset          = 0;
    placeholder.sample_meta_offset  = 0;
    placeholder.reserved            = 0;
    write_exact(out, &placeholder, sizeof(placeholder), tmp_path,
                "placeholder header");

    // ── Step 4: write sample metadata right after the header. Doing
    //    it up front keeps readers from needing to seek to the end
    //    for names when opening the file.

    const uint64_t sample_meta_offset =
        write_sample_metadata(out, samples, tmp_path);

    // ── Step 5: run the K-way merge, emitting per-segment blocks and
    //    accumulating a TOC as we go. The TOC is held in memory
    //    briefly — segment count at cohort scale is ~10M at the high
    //    end, × 20 bytes = 200 MB, which is tolerable. If this ever
    //    becomes a concern, switch to streaming the TOC to a temp
    //    file and concatenating at the end.

    const std::vector<TOCEntry> toc =
        run_final_merge(out, tmp_path, readers, excluded_segments);

    // ── Step 6: write the TOC at the current offset.

    const auto toc_offset = static_cast<uint64_t>(out.tellp());
    for (const auto& entry : toc) {
        write_exact(out, &entry, sizeof(entry), tmp_path, "TOC entry");
    }

    // ── Step 7: seek back and patch the real header with the final
    //    offsets and block count.

    out.seekp(0);
    if (!out) {
        throw std::runtime_error(
            "quant_sidecar::merge_to_qtx: seekp(0) failed for " +
            tmp_path.string());
    }

    Header hdr{};
    std::memcpy(hdr.magic, MAGIC, sizeof(MAGIC));
    hdr.version             = QTX_VERSION;
    hdr.segment_block_count = static_cast<uint64_t>(toc.size());
    hdr.toc_offset          = toc_offset;
    hdr.sample_meta_offset  = sample_meta_offset;
    hdr.reserved            = 0;
    write_exact(out, &hdr, sizeof(hdr), tmp_path, "final header");

    out.flush();
    out.close();
    if (!out) {
        throw std::runtime_error(
            "quant_sidecar::merge_to_qtx: flush/close failed for " +
            tmp_path.string());
    }

    // ── Step 8: atomic rename into place.

    std::error_code ec;
    std::filesystem::rename(tmp_path, output_path, ec);
    if (ec) {
        throw std::runtime_error(
            "quant_sidecar::merge_to_qtx: rename " + tmp_path.string() +
            " -> " + output_path.string() + " failed: " + ec.message());
    }

    // ── Step 9: clean up our own batched intermediates. (We only
    //    touch files we created here; per-sample input streams are
    //    the caller's responsibility.)

    for (const auto& p : owned_intermediates) {
        std::error_code rm_ec;
        std::filesystem::remove(p, rm_ec);
        // Swallow failure — stale intermediates are a housekeeping
        // annoyance, not a correctness problem, and we've already
        // committed the output file.
    }
}

// ── Reader ────────────────────────────────────────────────────────────

Reader::Reader(const std::filesystem::path& path)
    : path_(path) {
    file_.open(path_, std::ios::binary);
    if (!file_) {
        throw std::runtime_error(
            "quant_sidecar::Reader: cannot open " + path_.string());
    }

    // Load and validate the header.
    read_exact(file_, &header_, sizeof(header_), path_, "header");
    if (std::memcmp(header_.magic, MAGIC, sizeof(MAGIC)) != 0) {
        throw std::runtime_error(
            "quant_sidecar::Reader: bad magic in " + path_.string() +
            " (not a .qtx file)");
    }
    if (header_.version != QTX_VERSION) {
        throw std::runtime_error(
            "quant_sidecar::Reader: unsupported version " +
            std::to_string(header_.version) + " in " + path_.string() +
            " (expected " + std::to_string(QTX_VERSION) + ")");
    }

    // Load sample metadata (seek to the offset recorded in the header).
    file_.seekg(static_cast<std::streamoff>(header_.sample_meta_offset));
    if (!file_) {
        throw std::runtime_error(
            "quant_sidecar::Reader: seek to sample metadata failed in " +
            path_.string());
    }

    uint32_t num_samples = 0;
    read_exact(file_, &num_samples, sizeof(num_samples), path_, "num_samples");
    samples_.reserve(num_samples);
    for (uint32_t i = 0; i < num_samples; ++i) {
        SampleMetadata s;
        read_exact(file_, &s.sample_id, sizeof(s.sample_id), path_, "sample_id");
        read_exact(file_, &s.expr_type, sizeof(s.expr_type), path_, "expr_type");

        uint8_t name_len = 0;
        read_exact(file_, &name_len, sizeof(name_len), path_, "name_len");
        if (name_len > 0) {
            s.name.resize(name_len);
            read_exact(file_, s.name.data(),
                       static_cast<std::streamsize>(name_len),
                       path_, "sample name");
        }
        samples_.push_back(std::move(s));
    }

    // Load the TOC (seek to header_.toc_offset).
    file_.seekg(static_cast<std::streamoff>(header_.toc_offset));
    if (!file_) {
        throw std::runtime_error(
            "quant_sidecar::Reader: seek to TOC failed in " + path_.string());
    }
    toc_.resize(header_.segment_block_count);
    if (header_.segment_block_count > 0) {
        const auto bytes = static_cast<std::streamsize>(
            header_.segment_block_count * sizeof(TOCEntry));
        read_exact(file_, toc_.data(), bytes, path_, "TOC");
    }
}

void Reader::read_block(const TOCEntry& entry,
                        std::vector<ValueRecord>& out) const {
    out.clear();

    file_.seekg(static_cast<std::streamoff>(entry.offset));
    if (!file_) {
        throw std::runtime_error(
            "quant_sidecar::Reader: seek to block offset " +
            std::to_string(entry.offset) + " failed in " + path_.string());
    }

    uint64_t block_seg  = 0;
    uint32_t num_records = 0;
    read_exact(file_, &block_seg, sizeof(block_seg), path_,
               "block segment_index");
    read_exact(file_, &num_records, sizeof(num_records), path_,
               "block num_records");

    if (block_seg != entry.segment_index) {
        throw std::runtime_error(
            "quant_sidecar::Reader: TOC/block segment_index mismatch at offset " +
            std::to_string(entry.offset) + " in " + path_.string() +
            " (TOC says " + std::to_string(entry.segment_index) +
            ", block says " + std::to_string(block_seg) + ")");
    }

    out.resize(num_records);
    for (uint32_t i = 0; i < num_records; ++i) {
        read_exact(file_, &out[i].sample_id, sizeof(out[i].sample_id),
                   path_, "block sample_id");
        read_exact(file_, &out[i].value, sizeof(out[i].value),
                   path_, "block value");
    }
}

std::vector<Reader::ValueRecord> Reader::lookup(uint64_t segment_index) const {
    // Binary search on the TOC. We compare using a custom predicate
    // rather than constructing a sentinel TOCEntry so we don't depend
    // on operator< being defined for the struct.
    auto it = std::lower_bound(
        toc_.begin(), toc_.end(), segment_index,
        [](const TOCEntry& e, uint64_t v) { return e.segment_index < v; });
    if (it == toc_.end() || it->segment_index != segment_index) {
        return {};
    }

    std::vector<ValueRecord> out;
    read_block(*it, out);
    return out;
}

std::unordered_map<uint32_t, float> Reader::lookup_filtered(
    uint64_t segment_index,
    const std::vector<uint32_t>& wanted_samples) const {
    std::unordered_map<uint32_t, float> result;
    if (wanted_samples.empty()) return result;

    auto it = std::lower_bound(
        toc_.begin(), toc_.end(), segment_index,
        [](const TOCEntry& e, uint64_t v) { return e.segment_index < v; });
    if (it == toc_.end() || it->segment_index != segment_index) {
        return result;
    }

    std::vector<ValueRecord> block;
    read_block(*it, block);

    // Two-pointer intersection: block is sorted by sample_id by
    // construction; we sort a local copy of wanted_samples so this is
    // O(n + m) rather than O(n * m). Leaving wanted_samples itself
    // untouched keeps the API surface non-surprising.
    std::vector<uint32_t> wanted_sorted = wanted_samples;
    std::sort(wanted_sorted.begin(), wanted_sorted.end());

    size_t bi = 0, wi = 0;
    while (bi < block.size() && wi < wanted_sorted.size()) {
        if (block[bi].sample_id < wanted_sorted[wi]) {
            ++bi;
        } else if (block[bi].sample_id > wanted_sorted[wi]) {
            ++wi;
        } else {
            result.emplace(block[bi].sample_id, block[bi].value);
            ++bi;
            ++wi;
        }
    }
    return result;
}

} // namespace quant_sidecar
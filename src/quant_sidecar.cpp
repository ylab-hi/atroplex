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
#include <fstream>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace quant_sidecar {

// ── Writer ──────────────────────────────────────────────────────────

Writer::Writer(std::filesystem::path path, uint8_t expr_type, uint64_t grove_id)
    : path_(std::move(path)), expr_type_(expr_type), grove_id_(grove_id) {
    // Reserve a modest buffer. For TCGA StringTie samples this typically
    // ends up around ~10K records; small initial reservation keeps
    // per-sample overhead low while the vector's doubling handles the
    // tail growth.
    records_.reserve(8 * 1024);
}

Writer::~Writer() {
    if (!finalized_) {
        // Destructor must not throw — if the caller never called
        // finalize() explicitly, do our best and swallow any error.
        try {
            finalize();
        } catch (...) {
            // intentionally swallowed
        }
    }
}

void Writer::append(uint64_t segment_index, float value) {
    records_.push_back(Record{segment_index, value});
}

void Writer::finalize() {
    if (finalized_) return;
    finalized_ = true;

    // Sort by segment_index ascending. Duplicates are not expected —
    // the build path accumulates per-segment expression in segment_builder
    // before a value reaches us — but we don't dedup here, since the
    // build invariant is enforced upstream.
    std::sort(records_.begin(), records_.end(),
              [](const Record& a, const Record& b) {
                  return a.segment_index < b.segment_index;
              });

    // Atomic write: stage to path_ + ".tmp", then rename. This keeps
    // half-written sidecars from being mistaken for real ones if the
    // build crashes mid-flush.
    std::filesystem::path tmp_path = path_;
    tmp_path += ".tmp";

    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("quant_sidecar::Writer: cannot open "
                                 + tmp_path.string() + " for write");
    }

    Header hdr{};
    std::memcpy(hdr.magic, MAGIC, sizeof(MAGIC));
    hdr.version   = QTX_VERSION;
    hdr.count     = static_cast<uint32_t>(records_.size());
    hdr.expr_type = expr_type_;
    hdr.pad[0] = hdr.pad[1] = hdr.pad[2] = 0;
    hdr.grove_id  = grove_id_;
    hdr.reserved  = 0;

    out.write(reinterpret_cast<const char*>(&hdr), sizeof(Header));
    if (!out) {
        throw std::runtime_error("quant_sidecar::Writer: header write failed for "
                                 + tmp_path.string());
    }

    if (!records_.empty()) {
        const auto bytes = static_cast<std::streamsize>(
            records_.size() * sizeof(Record));
        out.write(reinterpret_cast<const char*>(records_.data()), bytes);
        if (!out) {
            throw std::runtime_error("quant_sidecar::Writer: record write failed for "
                                     + tmp_path.string());
        }
    }

    out.flush();
    out.close();
    if (!out) {
        throw std::runtime_error("quant_sidecar::Writer: flush/close failed for "
                                 + tmp_path.string());
    }

    std::error_code ec;
    std::filesystem::rename(tmp_path, path_, ec);
    if (ec) {
        throw std::runtime_error("quant_sidecar::Writer: rename "
                                 + tmp_path.string() + " -> " + path_.string()
                                 + " failed: " + ec.message());
    }

    // Drop in-memory state — the Writer is now at end of life.
    records_.clear();
    records_.shrink_to_fit();
}

// ── Reader ──────────────────────────────────────────────────────────

Reader::Reader(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("quant_sidecar::Reader: cannot open "
                                 + path.string());
    }

    in.read(reinterpret_cast<char*>(&header_), sizeof(Header));
    if (!in) {
        throw std::runtime_error("quant_sidecar::Reader: truncated header in "
                                 + path.string());
    }

    if (std::memcmp(header_.magic, MAGIC, sizeof(MAGIC)) != 0) {
        throw std::runtime_error("quant_sidecar::Reader: bad magic in "
                                 + path.string() + " (not a .qtx file)");
    }
    if (header_.version != QTX_VERSION) {
        throw std::runtime_error("quant_sidecar::Reader: unsupported version "
                                 + std::to_string(header_.version) + " in "
                                 + path.string() + " (expected "
                                 + std::to_string(QTX_VERSION) + ")");
    }

    records_.resize(header_.count);
    if (header_.count > 0) {
        const auto bytes = static_cast<std::streamsize>(
            static_cast<size_t>(header_.count) * sizeof(Record));
        in.read(reinterpret_cast<char*>(records_.data()), bytes);
        if (!in) {
            throw std::runtime_error("quant_sidecar::Reader: truncated records in "
                                     + path.string());
        }
    }
}

std::optional<float> Reader::lookup(uint64_t segment_index) const {
    auto it = std::lower_bound(
        records_.begin(), records_.end(), segment_index,
        [](const Record& r, uint64_t v) { return r.segment_index < v; });
    if (it != records_.end() && it->segment_index == segment_index) {
        return it->value;
    }
    return std::nullopt;
}

} // namespace quant_sidecar

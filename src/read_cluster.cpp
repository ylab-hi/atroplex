/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "read_cluster.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>

// ============================================================================
// processed_read implementation
// ============================================================================

processed_read processed_read::from_sam_entry(const gio::sam_entry& entry) {
    processed_read read;
    read.read_id = entry.qname;
    read.seqid = entry.chrom;
    read.strand = entry.get_strand();
    read.interval = entry.interval;
    read.mapq = entry.mapq;

    // Extract splice junctions from CIGAR N (REF_SKIP) operations
    size_t ref_pos = entry.interval.get_start();

    for (const auto& op : entry.cigar) {
        if (op.op == gio::cigar_op::REF_SKIP) {
            // N operation = intron/splice junction
            size_t donor = ref_pos;
            size_t acceptor = ref_pos + op.length;
            read.junctions.emplace_back(donor, acceptor);
        }

        // Advance reference position for reference-consuming operations
        if (op.consumes_reference()) {
            ref_pos += op.length;
        }
    }

    return read;
}

std::string processed_read::get_splice_signature() const {
    std::ostringstream ss;
    ss << seqid << ":" << strand << ":";

    if (junctions.empty()) {
        // Single-exon read - use position as signature
        ss << "SINGLE:" << interval.get_start() << "-" << interval.get_end();
    } else {
        for (size_t i = 0; i < junctions.size(); ++i) {
            if (i > 0) ss << ",";
            ss << junctions[i].donor << "-" << junctions[i].acceptor;
        }
    }

    return ss.str();
}

std::string processed_read::get_binned_signature(int bin_size) const {
    std::ostringstream ss;
    ss << seqid << ":" << strand << ":";

    if (junctions.empty()) {
        // Single-exon reads use larger bins
        int se_bin = bin_size * 10;  // 100bp default for single-exon
        size_t binned_start = (interval.get_start() / se_bin) * se_bin;
        size_t binned_end = (interval.get_end() / se_bin) * se_bin;
        ss << "SINGLE:" << binned_start << "-" << binned_end;
    } else {
        for (size_t i = 0; i < junctions.size(); ++i) {
            if (i > 0) ss << ",";
            size_t binned_donor = (junctions[i].donor / bin_size) * bin_size;
            size_t binned_acceptor = (junctions[i].acceptor / bin_size) * bin_size;
            ss << binned_donor << "-" << binned_acceptor;
        }
    }

    return ss.str();
}

// ============================================================================
// read_cluster implementation
// ============================================================================

void read_cluster::add_read(const processed_read* read) {
    if (members.empty()) {
        seqid = read->seqid;
        strand = read->strand;
    }

    members.push_back(read);

    // Update bounds
    if (read->interval.get_start() < start) {
        start = read->interval.get_start();
    }
    if (read->interval.get_end() > end) {
        end = read->interval.get_end();
    }
}

void read_cluster::finalize() {
    if (members.empty()) return;

    // Compute consensus junctions by averaging positions
    // First, collect all junction sets
    if (!members[0]->junctions.empty()) {
        size_t n_junctions = members[0]->junctions.size();
        consensus_junctions.resize(n_junctions);

        for (size_t j = 0; j < n_junctions; ++j) {
            // Average donor and acceptor positions across all members
            double sum_donor = 0.0;
            double sum_acceptor = 0.0;
            size_t count = 0;

            for (const auto* read : members) {
                if (j < read->junctions.size()) {
                    sum_donor += static_cast<double>(read->junctions[j].donor);
                    sum_acceptor += static_cast<double>(read->junctions[j].acceptor);
                    count++;
                }
            }

            if (count > 0) {
                consensus_junctions[j].donor = static_cast<size_t>(sum_donor / count + 0.5);
                consensus_junctions[j].acceptor = static_cast<size_t>(sum_acceptor / count + 0.5);
            }
        }
    }

    // Generate cluster ID
    std::ostringstream ss;
    ss << seqid << "_" << strand << "_" << start << "_" << end << "_n" << members.size();
    cluster_id = ss.str();
}

double read_cluster::mean_mapq() const {
    if (members.empty()) return 0.0;

    double sum = 0.0;
    for (const auto* read : members) {
        sum += read->mapq;
    }
    return sum / static_cast<double>(members.size());
}

const processed_read* read_cluster::get_representative() const {
    if (members.empty()) return nullptr;

    // Select by highest mapq, then by longest read
    const processed_read* best = members[0];
    for (const auto* read : members) {
        if (read->mapq > best->mapq ||
            (read->mapq == best->mapq &&
             read->interval.get_end() - read->interval.get_start() >
             best->interval.get_end() - best->interval.get_start())) {
            best = read;
        }
    }
    return best;
}

// ============================================================================
// read_clusterer implementation
// ============================================================================

read_clusterer::read_clusterer(const config& cfg)
    : cfg_(cfg) {}

void read_clusterer::add_read(const gio::sam_entry& entry) {
    stats_.total_reads++;

    // Filter by mapping quality
    if (entry.mapq < cfg_.min_mapq) {
        stats_.filtered_reads++;
        return;
    }

    // Process the read
    processed_read read = processed_read::from_sam_entry(entry);

    // Track single vs multi-exon
    if (read.is_single_exon()) {
        stats_.single_exon_reads++;
    } else {
        stats_.multi_exon_reads++;
    }

    // Store read and index by signature
    size_t idx = read_storage_.size();
    read_storage_.push_back(std::move(read));

    std::string sig = read_storage_.back().get_binned_signature(cfg_.junction_bin_size);
    signature_to_reads_[sig].push_back(idx);
}

std::vector<read_cluster> read_clusterer::finalize_chromosome() {
    auto clusters = refine_clusters();

    // Update stats
    for (const auto& cluster : clusters) {
        stats_.total_clusters++;
        if (cluster.read_count() > stats_.max_cluster_size) {
            stats_.max_cluster_size = cluster.read_count();
        }
        if (cluster.read_count() == 1) {
            stats_.single_read_clusters++;
        } else {
            stats_.multi_read_clusters++;
        }
    }

    // Clear state for next chromosome
    current_seqid_.clear();
    read_storage_.clear();
    signature_to_reads_.clear();

    return clusters;
}

std::vector<read_cluster> read_clusterer::refine_clusters() {
    std::vector<read_cluster> clusters;

    for (auto& [signature, read_indices] : signature_to_reads_) {
        if (read_indices.empty()) continue;

        // For each bin, we need to verify junction compatibility
        // and potentially split into sub-clusters

        // Simple approach: treat each bin as one cluster
        // More sophisticated: check junction tolerance and split if needed

        // Get first read to determine junction count
        const processed_read& first_read = read_storage_[read_indices[0]];

        if (first_read.is_single_exon()) {
            // Single-exon reads: just group by bin
            read_cluster cluster;
            for (size_t idx : read_indices) {
                cluster.add_read(&read_storage_[idx]);
            }
            cluster.finalize();
            clusters.push_back(std::move(cluster));
        } else {
            // Multi-exon reads: verify junction compatibility
            // Group reads with matching junctions within tolerance

            std::vector<bool> assigned(read_indices.size(), false);

            for (size_t i = 0; i < read_indices.size(); ++i) {
                if (assigned[i]) continue;

                const processed_read& seed = read_storage_[read_indices[i]];
                read_cluster cluster;
                cluster.add_read(&seed);
                assigned[i] = true;

                // Find compatible reads
                for (size_t j = i + 1; j < read_indices.size(); ++j) {
                    if (assigned[j]) continue;

                    const processed_read& candidate = read_storage_[read_indices[j]];

                    // Check junction compatibility
                    bool compatible = true;
                    if (seed.junctions.size() == candidate.junctions.size()) {
                        for (size_t k = 0; k < seed.junctions.size(); ++k) {
                            if (!seed.junctions[k].matches(candidate.junctions[k],
                                                           cfg_.junction_tolerance)) {
                                compatible = false;
                                break;
                            }
                        }
                    } else {
                        compatible = false;
                    }

                    if (compatible) {
                        cluster.add_read(&candidate);
                        assigned[j] = true;
                    }
                }

                cluster.finalize();
                clusters.push_back(std::move(cluster));
            }
        }
    }

    // Sort clusters by genomic position
    std::sort(clusters.begin(), clusters.end(),
              [](const read_cluster& a, const read_cluster& b) {
                  if (a.seqid != b.seqid) return a.seqid < b.seqid;
                  if (a.start != b.start) return a.start < b.start;
                  return a.end < b.end;
              });

    return clusters;
}

std::vector<read_cluster> read_clusterer::cluster_reads(gio::bam_reader& reader) {
    std::vector<read_cluster> all_clusters;

    for (const auto& entry : reader) {
        // Skip unmapped reads (should already be filtered by bam_reader)
        if (!entry.is_mapped()) continue;

        // Check if we've moved to a new chromosome
        if (!current_seqid_.empty() && entry.chrom != current_seqid_) {
            // Finalize previous chromosome
            auto chr_clusters = finalize_chromosome();
            all_clusters.insert(all_clusters.end(),
                               std::make_move_iterator(chr_clusters.begin()),
                               std::make_move_iterator(chr_clusters.end()));
        }

        current_seqid_ = entry.chrom;
        add_read(entry);
    }

    // Finalize last chromosome
    if (!current_seqid_.empty()) {
        auto chr_clusters = finalize_chromosome();
        all_clusters.insert(all_clusters.end(),
                           std::make_move_iterator(chr_clusters.begin()),
                           std::make_move_iterator(chr_clusters.end()));
    }

    // Calculate mean cluster size
    if (stats_.total_clusters > 0) {
        size_t total_reads_in_clusters = stats_.total_reads - stats_.filtered_reads;
        stats_.mean_cluster_size = static_cast<double>(total_reads_in_clusters) /
                                   static_cast<double>(stats_.total_clusters);
    }

    return all_clusters;
}
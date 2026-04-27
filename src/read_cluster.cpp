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
    read.interval = gdt::interval(entry.start, entry.end);
    read.mapq = entry.mapq;

    // Extract splice junctions from CIGAR N (REF_SKIP) operations
    size_t ref_pos = entry.start;

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
    std::string sig;
    sig.reserve(seqid.size() + 3 + junctions.size() * 20);
    sig += seqid;
    sig += ':';
    sig += strand;
    sig += ':';

    if (junctions.empty()) {
        sig += "SINGLE:";
        sig += std::to_string(interval.get_start());
        sig += '-';
        sig += std::to_string(interval.get_end());
    } else {
        for (size_t i = 0; i < junctions.size(); ++i) {
            if (i > 0) sig += ',';
            sig += std::to_string(junctions[i].donor);
            sig += '-';
            sig += std::to_string(junctions[i].acceptor);
        }
    }
    return sig;
}

std::string processed_read::get_binned_signature(int bin_size) const {
    std::string sig;
    sig.reserve(seqid.size() + 3 + junctions.size() * 20);
    sig += seqid;
    sig += ':';
    sig += strand;
    sig += ':';

    if (junctions.empty()) {
        int se_bin = bin_size * 10;
        size_t binned_start = (interval.get_start() / se_bin) * se_bin;
        size_t binned_end = (interval.get_end() / se_bin) * se_bin;
        sig += "SINGLE:";
        sig += std::to_string(binned_start);
        sig += '-';
        sig += std::to_string(binned_end);
    } else {
        for (size_t i = 0; i < junctions.size(); ++i) {
            if (i > 0) sig += ',';
            size_t binned_donor = (junctions[i].donor / bin_size) * bin_size;
            size_t binned_acceptor = (junctions[i].acceptor / bin_size) * bin_size;
            sig += std::to_string(binned_donor);
            sig += '-';
            sig += std::to_string(binned_acceptor);
        }
    }
    return sig;
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

    read_storage_.push_back(std::move(read));
}

std::vector<read_cluster> read_clusterer::finalize_chromosome() {
    auto clusters = build_clusters();

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

    return clusters;
}

std::vector<read_cluster> read_clusterer::build_clusters() {
    if (read_storage_.empty()) return {};

    // Build index sorted by (strand, junction_count, first_donor) so that
    // compatible reads are adjacent. Sweep forward collecting reads whose
    // first junction is within tolerance — no binning boundary artifacts.
    std::vector<size_t> indices(read_storage_.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        auto& ra = read_storage_[a];
        auto& rb = read_storage_[b];
        if (ra.strand != rb.strand) return ra.strand < rb.strand;
        if (ra.junctions.size() != rb.junctions.size())
            return ra.junctions.size() < rb.junctions.size();
        if (!ra.junctions.empty() && !rb.junctions.empty()) {
            if (ra.junctions[0].donor != rb.junctions[0].donor)
                return ra.junctions[0].donor < rb.junctions[0].donor;
        }
        return ra.interval.get_start() < rb.interval.get_start();
    });

    std::vector<read_cluster> clusters;
    std::vector<bool> assigned(indices.size(), false);

    for (size_t i = 0; i < indices.size(); ++i) {
        if (assigned[i]) continue;

        auto& seed = read_storage_[indices[i]];

        read_cluster cluster;
        cluster.add_read(&seed);
        assigned[i] = true;

        if (seed.is_single_exon()) {
            // Single-exon: sweep by start position proximity
            for (size_t j = i + 1; j < indices.size(); ++j) {
                if (assigned[j]) continue;
                auto& cand = read_storage_[indices[j]];
                if (cand.strand != seed.strand) break;
                if (!cand.is_single_exon()) break;
                if (cand.interval.get_start() - seed.interval.get_start()
                    > cfg_.single_exon_distance) break;
                cluster.add_read(&cand);
                assigned[j] = true;
            }
        } else {
            // Multi-exon: sweep while first donor is within tolerance
            for (size_t j = i + 1; j < indices.size(); ++j) {
                if (assigned[j]) continue;
                auto& cand = read_storage_[indices[j]];
                if (cand.strand != seed.strand) break;
                if (cand.junctions.size() != seed.junctions.size()) break;
                // First donor beyond tolerance → no more candidates
                if (cand.junctions[0].donor > seed.junctions[0].donor
                    + static_cast<size_t>(cfg_.junction_tolerance)) break;

                // Full junction compatibility check
                bool compatible = true;
                for (size_t k = 0; k < seed.junctions.size(); ++k) {
                    if (!seed.junctions[k].matches(cand.junctions[k],
                                                   cfg_.junction_tolerance)) {
                        compatible = false;
                        break;
                    }
                }
                if (compatible) {
                    cluster.add_read(&cand);
                    assigned[j] = true;
                }
            }
        }

        cluster.finalize();
        clusters.push_back(std::move(cluster));
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
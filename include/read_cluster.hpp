/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_READ_CLUSTER_HPP
#define ATROPLEX_READ_CLUSTER_HPP

// standard
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>

// genogrove
#include <genogrove/data_type/interval.hpp>
#include <genogrove/data_type/genomic_coordinate.hpp>
#include <genogrove/io/bam_reader.hpp>

namespace gdt = genogrove::data_type;
namespace gio = genogrove::io;

/**
 * Splice junction derived from CIGAR N (REF_SKIP) operations
 * Represents the boundary between two exons (intron coordinates)
 */
struct splice_junction {
    size_t donor;      // 5' end of intron (last exonic base + 1)
    size_t acceptor;   // 3' start of next exon

    splice_junction() : donor(0), acceptor(0) {}
    splice_junction(size_t d, size_t a) : donor(d), acceptor(a) {}

    size_t intron_length() const { return acceptor - donor; }

    bool operator==(const splice_junction& other) const {
        return donor == other.donor && acceptor == other.acceptor;
    }

    bool operator<(const splice_junction& other) const {
        if (donor != other.donor) return donor < other.donor;
        return acceptor < other.acceptor;
    }

    /**
     * Check if this junction matches another within a tolerance
     * @param other Junction to compare against
     * @param tolerance Maximum bp difference allowed for match
     * @return true if both donor and acceptor are within tolerance
     */
    bool matches(const splice_junction& other, int tolerance = 5) const {
        return std::abs(static_cast<int64_t>(donor) - static_cast<int64_t>(other.donor)) <= tolerance &&
               std::abs(static_cast<int64_t>(acceptor) - static_cast<int64_t>(other.acceptor)) <= tolerance;
    }
};

/**
 * Processed read with extracted splice junctions
 * Created from genogrove sam_entry by extracting N operations from CIGAR
 */
struct processed_read {
    std::string read_id;
    std::string seqid;
    char strand;
    gdt::interval interval;
    std::vector<splice_junction> junctions;
    uint8_t mapq;

    processed_read() : strand('.'), mapq(0) {}

    /**
     * Create processed_read from genogrove sam_entry
     * Extracts splice junctions from CIGAR N (REF_SKIP) operations
     */
    static processed_read from_sam_entry(const gio::sam_entry& entry);

    /**
     * Generate splice junction signature for clustering
     * Format: seqid:strand:d1-a1,d2-a2,...
     * For single-exon reads: seqid:strand:SINGLE:start-end
     */
    std::string get_splice_signature() const;

    /**
     * Generate binned splice signature for fuzzy clustering
     * Junction positions are binned to nearest bin_size
     * @param bin_size Bin size for position rounding (default 10bp)
     */
    std::string get_binned_signature(int bin_size = 10) const;

    /**
     * Check if this is a single-exon read (no splice junctions)
     */
    bool is_single_exon() const { return junctions.empty(); }

    /**
     * Get number of exons (junctions + 1)
     */
    size_t exon_count() const { return junctions.size() + 1; }
};

/**
 * Cluster of reads with similar splice junction signatures
 * Used for deduplication and consensus building
 */
struct read_cluster {
    std::string cluster_id;
    std::string seqid;
    char strand;
    size_t start;                              // Min start across all members
    size_t end;                                // Max end across all members
    std::vector<splice_junction> consensus_junctions;
    std::vector<const processed_read*> members; // Non-owning pointers

    read_cluster() : strand('.'), start(SIZE_MAX), end(0) {}

    /**
     * Add a read to this cluster
     * Updates bounds and member list
     */
    void add_read(const processed_read* read);

    /**
     * Finalize cluster after all reads added
     * Computes consensus junctions and selects representative
     */
    void finalize();

    /**
     * Get genomic coordinate for grove queries
     */
    gdt::genomic_coordinate get_coordinate() const {
        return gdt::genomic_coordinate(strand, start, end);
    }

    /**
     * Get read count in cluster
     */
    size_t read_count() const { return members.size(); }

    /**
     * Get mean mapping quality across members
     */
    double mean_mapq() const;

    /**
     * Get the representative read (highest mapq, longest)
     */
    const processed_read* get_representative() const;
};

/**
 * Read clustering engine
 * Groups reads by splice junction signature for efficient processing
 */
class read_clusterer {
public:
    struct config {
        int junction_bin_size = 10;    // Binning for initial clustering (bp)
        int junction_tolerance = 5;    // Max bp difference within cluster
        uint8_t min_mapq = 20;         // Minimum mapping quality
        int single_exon_bin_size = 100; // Larger bins for single-exon reads
    };

    read_clusterer() : read_clusterer(config{}) {}
    explicit read_clusterer(const config& cfg);

    /**
     * Cluster reads from a BAM file
     * Streams through file, processing chromosome-by-chromosome
     * @param reader Open bam_reader
     * @return Vector of clusters sorted by genomic position
     */
    std::vector<read_cluster> cluster_reads(gio::bam_reader& reader);

    /**
     * Cluster statistics
     */
    struct stats {
        size_t total_reads = 0;
        size_t filtered_reads = 0;      // Filtered due to mapq
        size_t single_exon_reads = 0;
        size_t multi_exon_reads = 0;
        size_t total_clusters = 0;
        size_t single_read_clusters = 0;
        size_t multi_read_clusters = 0;
        double mean_cluster_size = 0.0;
        double max_cluster_size = 0;
    };

    const stats& get_stats() const { return stats_; }

private:
    config cfg_;
    stats stats_;

    // Per-chromosome state
    std::string current_seqid_;
    std::vector<processed_read> read_storage_;
    std::unordered_map<std::string, std::vector<size_t>> signature_to_reads_;

    /**
     * Add a single read to current chromosome state
     */
    void add_read(const gio::sam_entry& entry);

    /**
     * Finalize current chromosome and return its clusters
     * Clears internal state for next chromosome
     */
    std::vector<read_cluster> finalize_chromosome();

    /**
     * Refine coarse clusters into final clusters
     * Verifies junction compatibility within bins
     */
    std::vector<read_cluster> refine_clusters();
};

#endif // ATROPLEX_READ_CLUSTER_HPP
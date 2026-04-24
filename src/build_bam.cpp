/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "build_bam.hpp"
#include "segment_builder.hpp"
#include "utility.hpp"

#include <sstream>

#include <genogrove/io/bam_reader.hpp>

namespace gio = genogrove::io;

void build_bam::build(grove_type& grove,
    const std::filesystem::path& filepath,
    std::optional<uint32_t> sample_id,
    chromosome_exon_caches& exon_caches,
    chromosome_segment_caches& segment_caches,
    size_t& segment_count,
    const expression_filters& filters,
    bool absorb,
    size_t fuzzy_tolerance,
    bool include_scaffolds,
    build_counters& counters,
    quant_sidecar::SampleStreamWriter* sidecar_writer,
    bool annotated_loci_only,
    const std::unordered_set<std::string>& chromosomes_filter) {

    const size_t segment_count_before = segment_count;

    // Step 1: Cluster reads by splice junction signature
    gio::bam_reader reader(filepath.string(), gio::bam_reader_options::primary_only());

    read_clusterer::config cfg;
    read_clusterer clusterer(cfg);
    auto clusters = clusterer.cluster_reads(reader);

    auto& stats = clusterer.get_stats();
    logging::info("Clustered " + std::to_string(stats.total_reads) + " reads into " +
                  std::to_string(stats.total_clusters) + " clusters from " +
                  filepath.filename().string());

    std::mutex grove_mutex;

    // Step 2: Convert each cluster to exon chain and create segments
    for (const auto& cluster : clusters) {
        // Skip single-exon clusters (no splice junctions)
        if (cluster.consensus_junctions.empty()) continue;

        // Scaffold filter (mirrors the GFF ingest path). Clusters on
        // unplaced contigs, alt haplotypes, fix patches, and decoys
        // are skipped unless --include-scaffolds is set. Each skipped
        // cluster represents one transcript-level unit, so we count
        // it against scaffold_filtered_transcripts and bail before
        // bumping input_transcripts.
        if (!is_main_chromosome(cluster.seqid, include_scaffolds)) {
            counters.scaffold_filtered_transcripts++;
            continue;
        }

        if (!chromosomes_filter.empty() &&
            chromosomes_filter.find(normalize_chromosome(cluster.seqid)) == chromosomes_filter.end()) {
            continue;
        }

        counters.input_transcripts++;

        // Filter by minimum read count (expression proxy). BAM clusters
        // only have a read-count signal, so only `--min-counts` applies.
        // Other attribute thresholds (TPM/FPKM/cov) are a no-op for BAM.
        float read_count = static_cast<float>(cluster.read_count());
        if (filters.min_counts >= 0 && read_count < filters.min_counts) {
            counters.discarded_transcripts++;
            continue;
        }

        std::string seqid = normalize_chromosome(cluster.seqid);

        // Derive exon coordinates from junctions
        auto exon_intervals = cluster_to_exon_coords(cluster);
        if (exon_intervals.size() < 2) {
            counters.discarded_transcripts++;
            continue;
        }

        // Resolve gene assignment by spatial query
        auto [gene_id, gene_name, gene_biotype] = resolve_gene(grove, cluster);

        // Generate a transcript ID for this cluster
        std::string tx_id = cluster.cluster_id.empty()
            ? seqid + "_" + std::to_string(cluster.start) + "_" + std::to_string(cluster.end)
            : cluster.cluster_id;

        // Insert exons and build chain
        std::vector<key_ptr> exon_chain;
        std::vector<gdt::genomic_coordinate> exon_coords;

        for (const auto& [ex_start, ex_end] : exon_intervals) {
            key_ptr exon_key = insert_exon(
                grove, grove_mutex, seqid, cluster.strand,
                ex_start, ex_end,
                gene_id, gene_name, gene_biotype,
                tx_id, exon_caches[seqid], sample_id, "BAM"
            );
            exon_chain.push_back(exon_key);
            exon_coords.emplace_back(cluster.strand, ex_start, ex_end);
        }

        // Create segment via segment_builder (same absorption rules as GTF)
        uint32_t seg_gene_idx = gene_registry::instance().intern(gene_id, gene_name, gene_biotype);
        segment_builder::create_segment(
            grove, grove_mutex, tx_id, seqid, cluster.strand,
            cluster.start, cluster.end,
            static_cast<int>(exon_intervals.size()),
            exon_coords, exon_chain,
            segment_caches[seqid], seg_gene_idx,
            sample_id, "BAM", segment_count,
            read_count, "",
            absorb, fuzzy_tolerance,
            counters, sidecar_writer, annotated_loci_only
        );
    }

    logging::progress_done(segment_count, segment_count - segment_count_before,
                           "Processed " + filepath.filename().string());
}

// ========================================================================
// Exon coordinate derivation
// ========================================================================

std::vector<std::pair<size_t, size_t>> build_bam::cluster_to_exon_coords(
    const read_cluster& cluster
) {
    std::vector<std::pair<size_t, size_t>> exons;

    if (cluster.consensus_junctions.empty()) {
        // Single exon
        exons.emplace_back(cluster.start, cluster.end);
        return exons;
    }

    // First exon: cluster start → first junction donor
    exons.emplace_back(cluster.start, cluster.consensus_junctions.front().donor);

    // Internal exons: between consecutive junctions
    for (size_t i = 0; i + 1 < cluster.consensus_junctions.size(); ++i) {
        exons.emplace_back(
            cluster.consensus_junctions[i].acceptor,
            cluster.consensus_junctions[i + 1].donor
        );
    }

    // Last exon: last junction acceptor → cluster end
    exons.emplace_back(cluster.consensus_junctions.back().acceptor, cluster.end);

    return exons;
}

// ========================================================================
// Gene resolution via spatial query
// ========================================================================

std::tuple<std::string, std::string, std::string> build_bam::resolve_gene(
    grove_type& grove,
    const read_cluster& cluster
) {
    std::string seqid = normalize_chromosome(cluster.seqid);
    gdt::genomic_coordinate query_coord(cluster.strand, cluster.start, cluster.end);

    // Spatial query: find overlapping segments in the grove
    auto result = grove.intersect(query_coord, seqid);

    // Find the best overlapping gene (by overlap length)
    std::string best_gene_id;
    std::string best_gene_name;
    std::string best_gene_biotype;
    size_t best_overlap = 0;

    for (auto* hit : result.get_keys()) {
        auto& feature = hit->get_data();
        if (!is_segment(feature)) continue;

        auto& seg = get_segment(feature);
        if (seg.absorbed) continue;

        auto& hit_coord = hit->get_value();
        size_t overlap_start = std::max(cluster.start, hit_coord.get_start());
        size_t overlap_end = std::min(cluster.end, hit_coord.get_end());

        if (overlap_end > overlap_start) {
            size_t overlap = overlap_end - overlap_start;
            if (overlap > best_overlap) {
                best_overlap = overlap;
                best_gene_id = seg.gene_id();
                best_gene_name = seg.gene_name();
                best_gene_biotype = seg.gene_biotype();
            }
        }
    }

    // If no overlap, create synthetic gene ID
    if (best_gene_id.empty()) {
        best_gene_id = "novel_" + seqid + "_" + std::to_string(cluster.start);
        best_gene_name = best_gene_id;
        best_gene_biotype = "unknown";
    }

    return {best_gene_id, best_gene_name, best_gene_biotype};
}

// ========================================================================
// Coordinate-based exon insertion (no GFF entry needed)
// ========================================================================

key_ptr build_bam::insert_exon(
    grove_type& grove,
    std::mutex& grove_mutex,
    const std::string& seqid,
    char strand,
    size_t start,
    size_t end,
    const std::string& gene_id,
    const std::string& gene_name,
    const std::string& gene_biotype,
    const std::string& transcript_id,
    exon_cache_type& exon_cache,
    std::optional<uint32_t> sample_id,
    const std::string& source
) {
    gdt::genomic_coordinate coord(strand, start, end);

    // Reuse existing exon if coordinates match
    auto cached = exon_cache.find(coord);
    if (cached != exon_cache.end()) {
        key_ptr exon_key = cached->second;
        auto& exon = get_exon(exon_key->get_data());
        exon.transcript_ids.insert(transcript_registry::instance().intern(transcript_id));
        if (sample_id.has_value()) {
            exon.add_sample(*sample_id);
        }
        if (!source.empty()) {
            exon.add_source(source);
        }
        return exon_key;
    }

    // Create new exon (gene_idx lives on segments, not exons)
    exon_feature new_exon;
    new_exon.transcript_ids.insert(transcript_registry::instance().intern(transcript_id));
    new_exon.exon_number = -1;  // Unknown from BAM

    if (sample_id.has_value()) {
        new_exon.add_sample(*sample_id);
    }
    if (!source.empty()) {
        new_exon.add_source(source);
    }

    genomic_feature feature = new_exon;
    key_ptr exon_key;
    {
        std::lock_guard<std::mutex> lock(grove_mutex);
        exon_key = grove.add_external_key(coord, feature);
    }
    exon_cache[coord] = exon_key;

    return exon_key;
}

// ========================================================================
// BAM header parsing for sample metadata
// ========================================================================

sample_info build_bam::parse_header(const std::filesystem::path& filepath) {
    sample_info info;
    info.source_file = filepath;
    info.id = filepath.stem().string();

    // Try to extract metadata from BAM header @RG and @PG lines
    try {
        gio::bam_reader reader(filepath.string());
        auto header_text = reader.get_header();

        // Parse @RG lines for sample/platform info
        std::istringstream stream(header_text);
        std::string line;
        while (std::getline(stream, line)) {
            if (line.substr(0, 3) == "@RG") {
                // Extract SM (sample), PL (platform), LB (library)
                auto extract_tag = [&](const std::string& tag) -> std::string {
                    size_t pos = line.find(tag + ":");
                    if (pos == std::string::npos) return "";
                    pos += tag.size() + 1;
                    size_t end = line.find('\t', pos);
                    return (end != std::string::npos) ? line.substr(pos, end - pos) : line.substr(pos);
                };

                std::string sm = extract_tag("SM");
                if (!sm.empty() && info.id == filepath.stem().string()) {
                    info.id = sm;
                }

                std::string pl = extract_tag("PL");
                if (!pl.empty()) info.platform = pl;
            }
            else if (line.substr(0, 3) == "@PG") {
                // Extract PN (program name), VN (version)
                auto extract_tag = [&](const std::string& tag) -> std::string {
                    size_t pos = line.find(tag + ":");
                    if (pos == std::string::npos) return "";
                    pos += tag.size() + 1;
                    size_t end = line.find('\t', pos);
                    return (end != std::string::npos) ? line.substr(pos, end - pos) : line.substr(pos);
                };

                std::string pn = extract_tag("PN");
                if (!pn.empty() && info.pipeline.empty()) info.pipeline = pn;

                std::string vn = extract_tag("VN");
                if (!vn.empty() && info.pipeline_version.empty()) info.pipeline_version = vn;
            }
        }
    } catch (const std::exception& e) {
        logging::warning("BAM header parse failed for " + filepath.string() + ": " + e.what());
    } catch (...) {
        logging::warning("BAM header parse failed for " + filepath.string() + ": unknown exception");
    }

    // BAM expression is always raw counts
    info.expr_type = sample_info::expression_type::COUNTS;

    return info;
}
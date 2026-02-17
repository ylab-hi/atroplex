/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "transcript_matcher.hpp"
#include "utility.hpp"

#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>
#include <functional>

// ============================================================================
// Helper functions for category string conversion
// ============================================================================

std::string to_string(structural_category cat) {
    switch (cat) {
        case structural_category::FSM: return "FSM";
        case structural_category::ISM: return "ISM";
        case structural_category::NIC: return "NIC";
        case structural_category::NNC: return "NNC";
        case structural_category::GENIC_INTRON: return "genic_intron";
        case structural_category::GENIC_GENOMIC: return "genic_genomic";
        case structural_category::ANTISENSE: return "antisense";
        case structural_category::INTERGENIC: return "intergenic";
        case structural_category::AMBIGUOUS: return "ambiguous";
        default: return "unknown";
    }
}

std::string subcategory::to_string() const {
    if (ism != ism_subcategory::NONE) {
        switch (ism) {
            case ism_subcategory::PREFIX: return "5prime_fragment";
            case ism_subcategory::SUFFIX: return "3prime_fragment";
            case ism_subcategory::INTERNAL: return "internal_fragment";
            case ism_subcategory::MONO_EXON: return "mono-exon";
            default: return "";
        }
    }
    if (nic != nic_subcategory::NONE) {
        switch (nic) {
            case nic_subcategory::COMBINATION: return "combination";
            case nic_subcategory::INTRON_RETENTION: return "intron_retention";
            case nic_subcategory::ALT_3PRIME: return "alternative_3end";
            case nic_subcategory::ALT_5PRIME: return "alternative_5end";
            case nic_subcategory::EXON_SKIPPING: return "exon_skipping";
            default: return "";
        }
    }
    if (nnc != nnc_subcategory::NONE) {
        switch (nnc) {
            case nnc_subcategory::NOVEL_DONOR: return "novel_donor";
            case nnc_subcategory::NOVEL_ACCEPTOR: return "novel_acceptor";
            case nnc_subcategory::NOVEL_BOTH: return "novel_donor_and_acceptor";
            case nnc_subcategory::NOVEL_EXON: return "novel_exon";
            default: return "";
        }
    }
    return "";
}

std::string match_result::category_string() const {
    std::string cat_str = ::to_string(category);
    std::string sub_str = subcat.to_string();
    if (!sub_str.empty()) {
        return cat_str + "_" + sub_str;
    }
    return cat_str;
}

// ============================================================================
// transcript_matcher implementation
// ============================================================================

transcript_matcher::transcript_matcher(grove_type& grove, const config& cfg)
    : grove_(grove), cfg_(cfg) {
    // Index known splice sites from reference
    index_splice_sites();
}

void transcript_matcher::index_splice_sites() {
    // TODO: Iterate through all exons in grove and record splice sites
    // For now this is a placeholder - would need grove iteration support
    // This would be populated from exon boundaries:
    // - donor = exon end (5' splice site)
    // - acceptor = exon start (3' splice site)
}

bool transcript_matcher::is_known_donor(const std::string& seqid, size_t position) const {
    // Check if position (with tolerance) matches any known donor site
    // For now, simplified - would use spatial index for efficiency
    for (int offset = -cfg_.splice_site_window; offset <= cfg_.splice_site_window; ++offset) {
        size_t check_pos = position + offset;
        // Hash: combine seqid hash with position
        size_t key = std::hash<std::string>{}(seqid) ^ (check_pos << 1);
        if (known_donor_sites_.count(key)) {
            return true;
        }
    }
    return false;
}

bool transcript_matcher::is_known_acceptor(const std::string& seqid, size_t position) const {
    for (int offset = -cfg_.splice_site_window; offset <= cfg_.splice_site_window; ++offset) {
        size_t check_pos = position + offset;
        size_t key = std::hash<std::string>{}(seqid) ^ (check_pos << 1);
        if (known_acceptor_sites_.count(key)) {
            return true;
        }
    }
    return false;
}

match_result transcript_matcher::match(const read_cluster& cluster) {
    match_result result;
    stats_.total_matches++;

    // Step 1: Find candidate segments via spatial query
    auto candidates = find_candidate_segments(cluster);

    if (candidates.empty()) {
        result.category = structural_category::INTERGENIC;
        update_stats(result);
        return result;
    }

    // Step 2-3: For each candidate, get exon chain and score junction match
    double best_score = 0.0;
    key_ptr best_segment = nullptr;
    std::string best_transcript_id;
    std::string best_gene_id;
    std::vector<key_ptr> best_exon_chain;
    int best_ref_junctions = 0;

    for (key_ptr seg_key : candidates) {
        if (!is_segment(seg_key->get_data())) continue;

        const auto& seg = get_segment(seg_key->get_data());
        if (seg.transcript_ids.empty()) continue;

        // All transcripts of a segment share the same exon chain — get it once
        auto exon_chain = get_exon_chain(seg_key, "");
        if (exon_chain.empty()) continue;

        double score = score_junction_match(cluster, exon_chain);

        if (score > best_score) {
            best_score = score;
            best_segment = seg_key;
            best_transcript_id = transcript_registry::instance().resolve(*seg.transcript_ids.begin());
            best_gene_id = seg.gene_id;
            best_exon_chain = exon_chain;
            best_ref_junctions = static_cast<int>(exon_chain.size()) - 1;
        }
    }

    // No matching segment found
    if (best_segment == nullptr) {
        // Check if we're within a gene but not matching any transcript
        // This could be GENIC_INTRON or GENIC_GENOMIC
        result.category = structural_category::INTERGENIC;
        update_stats(result);
        return result;
    }

    // Populate result with match info
    result.matched_segments.push_back(best_segment);
    result.matched_transcript_ids.push_back(best_transcript_id);
    result.matched_gene_ids.push_back(best_gene_id);
    result.reference_transcript = best_transcript_id;
    result.reference_gene = best_gene_id;
    result.junction_match_score = best_score;
    result.total_query_junctions = static_cast<int>(cluster.consensus_junctions.size());
    result.total_ref_junctions = best_ref_junctions;
    result.matching_junctions = static_cast<int>(best_score * result.total_query_junctions + 0.5);

    // Analyze splice sites for NIC/NNC classification
    analyze_splice_sites(result, cluster);

    // Find novel junctions
    result.novel_junctions = find_novel_junctions(cluster, best_exon_chain);

    // Classify using SQANTI categories
    classify_match(result, cluster, best_exon_chain);

    // Check for ambiguous matches (multiple segments with same best score)
    int equal_matches = 0;
    for (key_ptr seg_key : candidates) {
        if (!is_segment(seg_key->get_data())) continue;
        if (seg_key == best_segment) continue;

        const auto& seg = get_segment(seg_key->get_data());
        if (seg.transcript_ids.empty()) continue;

        auto exon_chain = get_exon_chain(seg_key, "");
        if (exon_chain.empty()) continue;

        double score = score_junction_match(cluster, exon_chain);
        if (std::abs(score - best_score) < 0.001 && score >= cfg_.min_junction_score) {
            equal_matches++;
            result.matched_segments.push_back(seg_key);
            result.matched_transcript_ids.push_back(
                transcript_registry::instance().resolve(*seg.transcript_ids.begin()));
            result.matched_gene_ids.push_back(seg.gene_id);
        }
    }

    if (equal_matches > 0 && result.is_fsm()) {
        result.category = structural_category::AMBIGUOUS;
    }

    update_stats(result);
    return result;
}

void transcript_matcher::classify_match(match_result& result,
                                         const read_cluster& cluster,
                                         const std::vector<key_ptr>& best_chain) {
    // Single-exon reads (no junctions)
    if (cluster.consensus_junctions.empty()) {
        // Check if entirely within intron or overlapping exon
        if (!best_chain.empty()) {
            // Check if read overlaps any exon
            bool overlaps_exon = false;
            for (const auto& exon_key : best_chain) {
                const auto& exon_coord = exon_key->get_value();
                if (cluster.end > exon_coord.get_start() &&
                    cluster.start < exon_coord.get_end()) {
                    overlaps_exon = true;
                    break;
                }
            }

            if (overlaps_exon) {
                // Mono-exon within gene region
                result.category = structural_category::GENIC_GENOMIC;
            } else {
                // Entirely within intron
                result.category = structural_category::GENIC_INTRON;
            }
        } else {
            result.category = structural_category::INTERGENIC;
        }
        return;
    }

    // Multi-exon reads
    auto ref_junctions = extract_reference_junctions(best_chain);

    // FSM: All query junctions match reference, same number of junctions
    if (result.junction_match_score >= 0.999 &&
        result.total_query_junctions == result.total_ref_junctions) {
        result.category = structural_category::FSM;
        return;
    }

    // ISM: All query junctions match reference, but fewer junctions (subset)
    if (result.junction_match_score >= 0.999 &&
        result.total_query_junctions < result.total_ref_junctions) {
        result.category = structural_category::ISM;
        result.subcat.ism = classify_ism(cluster, best_chain);
        return;
    }

    // Check if all splice sites are known (for NIC vs NNC)
    bool all_sites_known = (result.novel_donors == 0 && result.novel_acceptors == 0);

    if (all_sites_known) {
        // NIC: Novel combination of known splice sites
        result.category = structural_category::NIC;
        result.subcat.nic = classify_nic(cluster, best_chain);
    } else {
        // NNC: At least one novel splice site
        result.category = structural_category::NNC;
        result.subcat.nnc = classify_nnc(result);
    }
}

ism_subcategory transcript_matcher::classify_ism(const read_cluster& cluster,
                                                   const std::vector<key_ptr>& ref_chain) {
    if (cluster.consensus_junctions.empty()) {
        return ism_subcategory::MONO_EXON;
    }

    auto ref_junctions = extract_reference_junctions(ref_chain);
    if (ref_junctions.empty()) {
        return ism_subcategory::NONE;
    }

    // Find which reference junctions are matched
    bool matches_first = false;
    bool matches_last = false;

    for (const auto& read_junc : cluster.consensus_junctions) {
        if (read_junc.matches(ref_junctions.front(), cfg_.junction_tolerance)) {
            matches_first = true;
        }
        if (read_junc.matches(ref_junctions.back(), cfg_.junction_tolerance)) {
            matches_last = true;
        }
    }

    if (matches_first && !matches_last) {
        return ism_subcategory::PREFIX;  // 5' fragment (missing 3' end)
    } else if (!matches_first && matches_last) {
        return ism_subcategory::SUFFIX;  // 3' fragment (missing 5' end)
    } else if (!matches_first && !matches_last) {
        return ism_subcategory::INTERNAL;
    }

    return ism_subcategory::NONE;
}

nic_subcategory transcript_matcher::classify_nic(const read_cluster& cluster,
                                                   const std::vector<key_ptr>& ref_chain) {
    // Check for intron retention
    if (has_intron_retention(cluster, ref_chain)) {
        return nic_subcategory::INTRON_RETENTION;
    }

    auto ref_junctions = extract_reference_junctions(ref_chain);
    int query_junctions = static_cast<int>(cluster.consensus_junctions.size());
    int ref_junction_count = static_cast<int>(ref_junctions.size());

    // Check for exon skipping (query has fewer junctions but spans same region)
    if (query_junctions < ref_junction_count && query_junctions > 0) {
        // Could be exon skipping - check if query junctions skip reference junctions
        return nic_subcategory::EXON_SKIPPING;
    }

    // Default: novel combination of known splice sites
    return nic_subcategory::COMBINATION;
}

nnc_subcategory transcript_matcher::classify_nnc(const match_result& result) {
    if (result.novel_donors > 0 && result.novel_acceptors > 0) {
        return nnc_subcategory::NOVEL_BOTH;
    } else if (result.novel_donors > 0) {
        return nnc_subcategory::NOVEL_DONOR;
    } else if (result.novel_acceptors > 0) {
        return nnc_subcategory::NOVEL_ACCEPTOR;
    }
    return nnc_subcategory::NONE;
}

void transcript_matcher::analyze_splice_sites(match_result& result,
                                               const read_cluster& cluster) {
    result.known_donors = 0;
    result.known_acceptors = 0;
    result.novel_donors = 0;
    result.novel_acceptors = 0;

    for (const auto& junction : cluster.consensus_junctions) {
        // Donor = 5' end of intron (end of exon)
        if (is_known_donor(cluster.seqid, junction.donor)) {
            result.known_donors++;
        } else {
            result.novel_donors++;
        }

        // Acceptor = 3' end of intron (start of next exon)
        if (is_known_acceptor(cluster.seqid, junction.acceptor)) {
            result.known_acceptors++;
        } else {
            result.novel_acceptors++;
        }
    }
}

bool transcript_matcher::has_intron_retention(const read_cluster& cluster,
                                               const std::vector<key_ptr>& ref_chain) {
    // Intron retention: query spans an intron without the junction
    auto ref_junctions = extract_reference_junctions(ref_chain);

    for (const auto& ref_junc : ref_junctions) {
        // Check if query spans this intron
        if (cluster.start < ref_junc.donor && cluster.end > ref_junc.acceptor) {
            // Query spans the intron region - check if junction is present
            bool has_junction = false;
            for (const auto& query_junc : cluster.consensus_junctions) {
                if (query_junc.matches(ref_junc, cfg_.junction_tolerance)) {
                    has_junction = true;
                    break;
                }
            }
            if (!has_junction) {
                return true;  // Intron retention detected
            }
        }
    }
    return false;
}

std::vector<match_result> transcript_matcher::match_batch(
    const std::vector<read_cluster>& clusters) {

    std::vector<match_result> results;
    results.reserve(clusters.size());

    for (const auto& cluster : clusters) {
        results.push_back(match(cluster));
    }

    return results;
}

void transcript_matcher::update_grove(const read_cluster& cluster,
                                       const match_result& result,
                                       bool add_novel) {
    if (result.has_match()) {
        // Add read support to matched segments
        for (key_ptr seg_key : result.matched_segments) {
            if (is_segment(seg_key->get_data())) {
                auto& seg = get_segment(seg_key->get_data());
                if (const auto* rep = cluster.get_representative()) {
                    seg.add_read_support(rep->read_id);
                }
                stats_.segments_updated++;
            }
        }
    }

    // Create discovered segment for novel patterns if requested
    if (add_novel && result.is_novel()) {
        key_ptr new_seg = create_discovered_segment(cluster);
        if (new_seg != nullptr) {
            stats_.segments_created++;
        }
    }
}

void transcript_matcher::update_stats(const match_result& result) {
    switch (result.category) {
        case structural_category::FSM:
            stats_.fsm_matches++;
            break;
        case structural_category::ISM:
            stats_.ism_matches++;
            switch (result.subcat.ism) {
                case ism_subcategory::PREFIX: stats_.ism_prefix++; break;
                case ism_subcategory::SUFFIX: stats_.ism_suffix++; break;
                case ism_subcategory::INTERNAL: stats_.ism_internal++; break;
                case ism_subcategory::MONO_EXON: stats_.ism_mono_exon++; break;
                default: break;
            }
            break;
        case structural_category::NIC:
            stats_.nic_matches++;
            switch (result.subcat.nic) {
                case nic_subcategory::COMBINATION: stats_.nic_combination++; break;
                case nic_subcategory::INTRON_RETENTION: stats_.nic_intron_retention++; break;
                case nic_subcategory::ALT_3PRIME: stats_.nic_alt_3prime++; break;
                case nic_subcategory::ALT_5PRIME: stats_.nic_alt_5prime++; break;
                case nic_subcategory::EXON_SKIPPING: stats_.nic_exon_skipping++; break;
                default: break;
            }
            break;
        case structural_category::NNC:
            stats_.nnc_matches++;
            switch (result.subcat.nnc) {
                case nnc_subcategory::NOVEL_DONOR: stats_.nnc_novel_donor++; break;
                case nnc_subcategory::NOVEL_ACCEPTOR: stats_.nnc_novel_acceptor++; break;
                case nnc_subcategory::NOVEL_BOTH: stats_.nnc_novel_both++; break;
                case nnc_subcategory::NOVEL_EXON: stats_.nnc_novel_exon++; break;
                default: break;
            }
            break;
        case structural_category::GENIC_INTRON:
            stats_.genic_intron_matches++;
            break;
        case structural_category::GENIC_GENOMIC:
            stats_.genic_genomic_matches++;
            break;
        case structural_category::ANTISENSE:
            stats_.antisense_matches++;
            break;
        case structural_category::INTERGENIC:
            stats_.intergenic_matches++;
            break;
        case structural_category::AMBIGUOUS:
            stats_.ambiguous_matches++;
            break;
    }
}

std::vector<key_ptr> transcript_matcher::find_candidate_segments(
    const read_cluster& cluster) {

    std::vector<key_ptr> candidates;

    gdt::genomic_coordinate query = cluster.get_coordinate();
    auto result = grove_.intersect(query, cluster.seqid);

    for (auto* key : result.get_keys()) {
        if (!is_segment(key->get_data())) continue;
        if (get_segment(key->get_data()).absorbed) continue;

        // Check minimum overlap
        const auto& seg_coord = key->get_value();
        size_t overlap_start = std::max(query.get_start(), seg_coord.get_start());
        size_t overlap_end = std::min(query.get_end(), seg_coord.get_end());

        if (overlap_end > overlap_start) {
            size_t overlap = overlap_end - overlap_start;
            if (static_cast<int>(overlap) >= cfg_.min_overlap_bp) {
                candidates.push_back(key);
            }
        }
    }

    return candidates;
}

std::vector<key_ptr> transcript_matcher::get_exon_chain(
    key_ptr segment, const std::string& /*transcript_id*/) {

    std::vector<key_ptr> chain;

    const auto& seg = get_segment(segment->get_data());
    size_t edge_id = seg.segment_index;

    // Get SEGMENT_TO_EXON edge to find first exon (filter by segment index)
    auto first_exons = grove_.get_neighbors_if(segment,
        [edge_id](const edge_metadata& e) {
            return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON &&
                   e.id == edge_id;
        });

    if (first_exons.empty()) return chain;

    // Traverse EXON_TO_EXON edges (filter by segment index)
    key_ptr current = first_exons[0];
    while (current != nullptr) {
        chain.push_back(current);

        auto next_exons = grove_.get_neighbors_if(current,
            [edge_id](const edge_metadata& e) {
                return e.type == edge_metadata::edge_type::EXON_TO_EXON &&
                       e.id == edge_id;
            });

        current = next_exons.empty() ? nullptr : next_exons[0];
    }

    return chain;
}

double transcript_matcher::score_junction_match(
    const read_cluster& cluster,
    const std::vector<key_ptr>& exon_chain) {

    if (cluster.consensus_junctions.empty()) {
        // Single-exon read - check if it fits within one exon
        return exon_chain.size() >= 1 ? 1.0 : 0.0;
    }

    // Extract reference junctions from exon chain
    auto ref_junctions = extract_reference_junctions(exon_chain);

    if (ref_junctions.empty()) {
        return 0.0;
    }

    // Count matching junctions
    int matches = 0;
    for (const auto& read_junc : cluster.consensus_junctions) {
        for (const auto& ref_junc : ref_junctions) {
            if (read_junc.matches(ref_junc, cfg_.junction_tolerance)) {
                matches++;
                break;
            }
        }
    }

    return static_cast<double>(matches) /
           static_cast<double>(cluster.consensus_junctions.size());
}

std::vector<splice_junction> transcript_matcher::find_novel_junctions(
    const read_cluster& cluster,
    const std::vector<key_ptr>& exon_chain) {

    std::vector<splice_junction> novel;
    auto ref_junctions = extract_reference_junctions(exon_chain);

    for (const auto& read_junc : cluster.consensus_junctions) {
        bool found = false;
        for (const auto& ref_junc : ref_junctions) {
            if (read_junc.matches(ref_junc, cfg_.junction_tolerance)) {
                found = true;
                break;
            }
        }
        if (!found) {
            novel.push_back(read_junc);
        }
    }

    return novel;
}

key_ptr transcript_matcher::create_discovered_segment(const read_cluster& cluster) {
    // Build segment feature for discovered pattern
    segment_feature new_segment;
    new_segment.add_source("atroplex");  // Mark as discovered by atroplex

    std::ostringstream coord_ss;
    coord_ss << cluster.seqid << ":" << cluster.strand << ":"
             << cluster.start << "-" << cluster.end;
    new_segment.coordinate = coord_ss.str();

    new_segment.exon_count = static_cast<int>(cluster.consensus_junctions.size() + 1);

    // Add read support
    if (const auto* rep = cluster.get_representative()) {
        new_segment.add_read_support(rep->read_id);
    }
    new_segment.read_coverage = cluster.read_count();

    // Insert into grove
    gdt::genomic_coordinate coord(cluster.strand, cluster.start, cluster.end);
    genomic_feature feature = new_segment;

    key_ptr seg_key = grove_.insert_data(cluster.seqid, coord, feature);

    return seg_key;
}

std::vector<splice_junction> transcript_matcher::extract_reference_junctions(
    const std::vector<key_ptr>& exon_chain) {

    std::vector<splice_junction> junctions;

    for (size_t i = 0; i + 1 < exon_chain.size(); ++i) {
        const auto& exon_coord = exon_chain[i]->get_value();
        const auto& next_coord = exon_chain[i + 1]->get_value();

        // Junction is at end of current exon to start of next exon
        junctions.emplace_back(exon_coord.get_end(), next_coord.get_start());
    }

    return junctions;
}

// ============================================================================
// Output methods
// ============================================================================

namespace {
    std::string join_strings(const std::vector<std::string>& strs, const std::string& sep) {
        if (strs.empty()) return ".";
        std::ostringstream ss;
        for (size_t i = 0; i < strs.size(); ++i) {
            if (i > 0) ss << sep;
            ss << strs[i];
        }
        return ss.str();
    }

    std::string junctions_to_string(const std::vector<splice_junction>& junctions) {
        if (junctions.empty()) return ".";
        std::ostringstream ss;
        for (size_t i = 0; i < junctions.size(); ++i) {
            if (i > 0) ss << ",";
            ss << junctions[i].donor << "-" << junctions[i].acceptor;
        }
        return ss.str();
    }
}

void transcript_matcher::write_results(const std::string& filepath,
                                        const std::vector<read_cluster>& clusters,
                                        const std::vector<match_result>& results) {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        logging::error("Could not open output file: " + filepath);
        return;
    }

    // SQANTI-like header
    out << "isoform\t"
        << "chrom\t"
        << "strand\t"
        << "length\t"
        << "exons\t"
        << "structural_category\t"
        << "subcategory\t"
        << "associated_gene\t"
        << "associated_transcript\t"
        << "ref_length\t"
        << "ref_exons\t"
        << "junctions\t"
        << "junction_match_score\t"
        << "matching_junctions\t"
        << "query_junctions\t"
        << "ref_junctions\t"
        << "known_donors\t"
        << "known_acceptors\t"
        << "novel_donors\t"
        << "novel_acceptors\t"
        << "novel_junctions\t"
        << "read_count\t"
        << "representative_read\n";

    for (size_t i = 0; i < clusters.size(); ++i) {
        const auto& cluster = clusters[i];
        const auto& result = results[i];

        const auto* rep = cluster.get_representative();
        std::string rep_id = rep ? rep->read_id : ".";
        size_t length = cluster.end - cluster.start;

        out << cluster.cluster_id << "\t"
            << cluster.seqid << "\t"
            << cluster.strand << "\t"
            << length << "\t"
            << (cluster.consensus_junctions.size() + 1) << "\t"
            << to_string(result.category) << "\t"
            << result.subcat.to_string() << "\t"
            << result.reference_gene.value_or(".") << "\t"
            << result.reference_transcript.value_or(".") << "\t"
            << ".\t"  // ref_length - would need lookup
            << (result.total_ref_junctions + 1) << "\t"
            << junctions_to_string(cluster.consensus_junctions) << "\t"
            << result.junction_match_score << "\t"
            << result.matching_junctions << "\t"
            << result.total_query_junctions << "\t"
            << result.total_ref_junctions << "\t"
            << result.known_donors << "\t"
            << result.known_acceptors << "\t"
            << result.novel_donors << "\t"
            << result.novel_acceptors << "\t"
            << junctions_to_string(result.novel_junctions) << "\t"
            << cluster.read_count() << "\t"
            << rep_id << "\n";
    }

    out.close();
    logging::info("Wrote results to: " + filepath);
}

void transcript_matcher::write_summary(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        logging::error("Could not open summary file: " + filepath);
        return;
    }

    out << "# Atroplex Transcript Classification Summary (SQANTI-like)\n\n";

    out << "## Structural Categories\n";
    out << "Total clusters: " << stats_.total_matches << "\n\n";

    auto pct = [this](size_t count) {
        if (stats_.total_matches == 0) return 0.0;
        return 100.0 * static_cast<double>(count) / static_cast<double>(stats_.total_matches);
    };

    out << "FSM (Full Splice Match): " << stats_.fsm_matches
        << " (" << pct(stats_.fsm_matches) << "%)\n";
    out << "ISM (Incomplete Splice Match): " << stats_.ism_matches
        << " (" << pct(stats_.ism_matches) << "%)\n";
    out << "NIC (Novel In Catalog): " << stats_.nic_matches
        << " (" << pct(stats_.nic_matches) << "%)\n";
    out << "NNC (Novel Not in Catalog): " << stats_.nnc_matches
        << " (" << pct(stats_.nnc_matches) << "%)\n";
    out << "Genic Intron: " << stats_.genic_intron_matches
        << " (" << pct(stats_.genic_intron_matches) << "%)\n";
    out << "Genic Genomic: " << stats_.genic_genomic_matches
        << " (" << pct(stats_.genic_genomic_matches) << "%)\n";
    out << "Antisense: " << stats_.antisense_matches
        << " (" << pct(stats_.antisense_matches) << "%)\n";
    out << "Intergenic: " << stats_.intergenic_matches
        << " (" << pct(stats_.intergenic_matches) << "%)\n";
    out << "Ambiguous: " << stats_.ambiguous_matches
        << " (" << pct(stats_.ambiguous_matches) << "%)\n";

    out << "\n## ISM Subcategories\n";
    out << "5' fragment (prefix): " << stats_.ism_prefix << "\n";
    out << "3' fragment (suffix): " << stats_.ism_suffix << "\n";
    out << "Internal fragment: " << stats_.ism_internal << "\n";
    out << "Mono-exon: " << stats_.ism_mono_exon << "\n";

    out << "\n## NIC Subcategories\n";
    out << "Combination: " << stats_.nic_combination << "\n";
    out << "Intron retention: " << stats_.nic_intron_retention << "\n";
    out << "Alternative 3' end: " << stats_.nic_alt_3prime << "\n";
    out << "Alternative 5' end: " << stats_.nic_alt_5prime << "\n";
    out << "Exon skipping: " << stats_.nic_exon_skipping << "\n";

    out << "\n## NNC Subcategories\n";
    out << "Novel donor: " << stats_.nnc_novel_donor << "\n";
    out << "Novel acceptor: " << stats_.nnc_novel_acceptor << "\n";
    out << "Novel donor + acceptor: " << stats_.nnc_novel_both << "\n";
    out << "Novel exon: " << stats_.nnc_novel_exon << "\n";

    out << "\n## Grove Updates\n";
    out << "Segments updated with read support: " << stats_.segments_updated << "\n";
    out << "Novel segments created: " << stats_.segments_created << "\n";

    out.close();
    logging::info("Wrote summary to: " + filepath);
}
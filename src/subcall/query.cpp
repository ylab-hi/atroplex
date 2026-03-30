/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/query.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>

#include <genogrove/io/gff_reader.hpp>

#include "sample_info.hpp"
#include "utility.hpp"

namespace gio = genogrove::io;

namespace {

std::string expr_type_label(sample_info::expression_type t) {
    switch (t) {
        case sample_info::expression_type::COUNTS: return "counts";
        case sample_info::expression_type::TPM:    return "TPM";
        case sample_info::expression_type::FPKM:   return "FPKM";
        case sample_info::expression_type::RPKM:   return "RPKM";
        case sample_info::expression_type::CPM:    return "CPM";
        default:                                   return "expression";
    }
}

} // anonymous namespace

namespace subcall {

// ============================================================================
// Contrast parsing
// ============================================================================

query_contrast query_contrast::parse(const std::string& spec) {
    query_contrast c;
    size_t pos = spec.find(':');
    if (pos == std::string::npos) {
        throw std::invalid_argument("Invalid contrast format: " + spec +
                                    " (expected group_a:group_b)");
    }
    c.group_a = spec.substr(0, pos);
    c.group_b = spec.substr(pos + 1);

    if (c.group_a.empty() || c.group_b.empty()) {
        throw std::invalid_argument("Empty group in contrast: " + spec);
    }
    return c;
}

// ============================================================================
// CLI options
// ============================================================================

cxxopts::Options query::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex query",
        "Classify transcripts against the pan-transcriptome index");

    options.add_options("Input/Output")
        ("i,input", "Input GTF/GFF with transcripts to classify",
            cxxopts::value<std::string>())
        ;

    options.add_options("Classification")
        ("junction-tolerance", "Max bp difference for junction match",
            cxxopts::value<int>()->default_value("5"))
        ("min-junction-score", "Minimum junction match score (0-1)",
            cxxopts::value<double>()->default_value("0.8"))
        ("min-overlap", "Minimum overlap with reference segment (bp)",
            cxxopts::value<int>()->default_value("50"))
        ;

    options.add_options("Differential usage")
        ("c,contrast", "Group contrast for DTU (format: group_a:group_b, repeatable)",
            cxxopts::value<std::vector<std::string>>())
        ("group-by", "Manifest field to group samples by (default: condition)",
            cxxopts::value<std::string>()->default_value("condition"))
        ("fdr", "FDR threshold for significance",
            cxxopts::value<double>()->default_value("0.05"))
        ("min-expression", "Minimum expression to include a transcript in DTU",
            cxxopts::value<float>()->default_value("0"))
        ;

    add_common_options(options);

    return options;
}

void query::validate(const cxxopts::ParseResult& args) {
    // Must have a grove source
    if (!args.count("genogrove") && !args.count("build-from") && !args.count("manifest")) {
        throw std::runtime_error("Must provide -g/--genogrove, -m/--manifest, or -b/--build-from");
    }

    // Must have input transcripts
    if (!args.count("input")) {
        throw std::runtime_error("Must provide -i/--input with transcripts to classify");
    }

    std::string input = args["input"].as<std::string>();
    if (!std::filesystem::exists(input)) {
        throw std::runtime_error("Input file not found: " + input);
    }
}

// ============================================================================
// Main execution
// ============================================================================

void query::execute(const cxxopts::ParseResult& args) {
    std::string input_path = args["input"].as<std::string>();
    logging::info("Classifying transcripts from: " + input_path);

    // Step 1: Classify input transcripts against the grove
    auto results = classify_transcripts(input_path);

    // Log classification summary
    std::map<structural_category, size_t> category_counts;
    for (const auto& r : results) {
        category_counts[r.category]++;
    }

    logging::info("Classification complete: " + std::to_string(results.size()) + " transcripts");
    for (const auto& [cat, count] : category_counts) {
        logging::info("  " + to_string(cat) + ": " + std::to_string(count));
    }

    // Step 2: Output classification
    auto out_dir = resolve_output_dir(args, input_path);
    std::string basename = std::filesystem::path(input_path).stem().string();

    std::string class_path = (out_dir / (basename + ".query.tsv")).string();
    write_classification(class_path, results);

    // Step 3: Optional DTU analysis
    std::vector<std::pair<query_contrast, std::vector<dtu_result>>> all_dtu;

    if (args.count("contrast")) {
        double fdr_threshold = args["fdr"].as<double>();
        auto contrast_strs = args["contrast"].as<std::vector<std::string>>();

        for (const auto& spec : contrast_strs) {
            auto contrast = query_contrast::parse(spec);
            logging::info("Running DTU: " + contrast.group_a + " vs " + contrast.group_b);

            auto dtu = run_dtu(contrast, results, args);
            apply_fdr_correction(dtu, fdr_threshold);

            size_t sig = std::count_if(dtu.begin(), dtu.end(),
                [](const dtu_result& r) { return r.significant; });
            logging::info("  Significant: " + std::to_string(sig) +
                          "/" + std::to_string(dtu.size()));

            std::string dtu_path = (out_dir / (basename + "." +
                contrast.group_a + "_vs_" + contrast.group_b + ".dtu.tsv")).string();
            write_dtu_results(dtu_path, contrast, dtu);

            all_dtu.emplace_back(contrast, std::move(dtu));
        }
    }

    // Step 4: Summary
    std::string summary_path = (out_dir / (basename + ".query.summary.txt")).string();
    write_summary(summary_path, results, all_dtu);

    logging::info("Results written to: " + out_dir.string());
}

// ============================================================================
// Transcript classification
// ============================================================================

std::vector<query_result> query::classify_transcripts(const std::string& input_path) {
    std::vector<query_result> results;

    // Configure matcher
    transcript_matcher::config match_cfg;
    match_cfg.junction_tolerance = 5;
    match_cfg.min_junction_score = 0.8;
    match_cfg.min_overlap_bp = 50;

    transcript_matcher matcher(*grove, match_cfg, exon_caches_);

    // Read input GTF/GFF gene-by-gene and classify each transcript
    gio::gff_reader reader(input_path);

    // Group entries by transcript
    std::unordered_map<std::string, std::vector<gio::gff_entry>> transcripts;
    for (const auto& entry : reader) {
        if (entry.type != "exon") continue;
        auto tx_id = entry.get_transcript_id();
        if (tx_id.has_value()) {
            transcripts[tx_id.value()].push_back(entry);
        }
    }

    logging::info("Read " + std::to_string(transcripts.size()) +
                  " transcripts from input");

    for (auto& [tx_id, exon_entries] : transcripts) {
        // Sort exons in genomic order
        std::sort(exon_entries.begin(), exon_entries.end(),
            [](const gio::gff_entry& a, const gio::gff_entry& b) {
                return a.start < b.start;
            });

        // Skip mono-exon transcripts (consistent with build pipeline)
        if (exon_entries.size() < 2) continue;

        // Build a read_cluster from the exon chain (reuse matcher infrastructure)
        read_cluster cluster;
        cluster.seqid = normalize_chromosome(exon_entries.front().seqid);
        cluster.strand = exon_entries.front().strand.value_or('+');
        cluster.start = exon_entries.front().start;
        cluster.end = exon_entries.back().end;

        // Extract splice junctions from consecutive exons
        for (size_t i = 0; i + 1 < exon_entries.size(); ++i) {
            splice_junction junc(exon_entries[i].end, exon_entries[i + 1].start);
            cluster.consensus_junctions.push_back(junc);
        }

        cluster.cluster_id = tx_id;

        // Classify against the grove
        match_result match = matcher.match(cluster);

        // Build query result with per-sample context
        query_result qr;
        qr.transcript_id = tx_id;
        qr.category = match.category;
        qr.subcat = match.subcat;
        qr.junction_match_score = match.junction_match_score;
        qr.matching_junctions = match.matching_junctions;
        qr.total_query_junctions = match.total_query_junctions;
        qr.total_ref_junctions = match.total_ref_junctions;
        qr.known_donors = match.known_donors;
        qr.known_acceptors = match.known_acceptors;
        qr.novel_donors = match.novel_donors;
        qr.novel_acceptors = match.novel_acceptors;

        if (match.has_match()) {
            qr.gene_id = match.reference_gene.value_or(".");
            qr.gene_name = ".";

            // Collect per-sample presence and expression from matched segment
            auto* best_seg_key = match.matched_segments.front();
            const auto& seg = get_segment(best_seg_key->get_data());

            qr.gene_name = seg.gene_name();

            seg.sample_idx.for_each([&](uint32_t sid) {
                qr.sample_ids.push_back(sid);
                if (seg.expression.has(sid)) {
                    qr.sample_expression[sid] = seg.expression.get(sid);
                }
            });
        }

        results.push_back(std::move(qr));
    }

    return results;
}

// ============================================================================
// DTU analysis
// ============================================================================

std::vector<dtu_result> query::run_dtu(
    const query_contrast& contrast,
    const std::vector<query_result>& results,
    const cxxopts::ParseResult& /*args*/) {

    // TODO: Implement proper DTU statistical testing
    // For now, compute per-group proportions using expression from matched segments
    std::vector<dtu_result> dtu_results;

    logging::warning("DTU statistical testing not yet implemented — "
                     "proportions computed but p-values are placeholders");

    return dtu_results;
}

void query::apply_fdr_correction(std::vector<dtu_result>& results,
                                  double fdr_threshold) {
    if (results.empty()) return;

    std::vector<size_t> indices(results.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
              [&results](size_t a, size_t b) {
                  return results[a].p_value < results[b].p_value;
              });

    size_t n = results.size();
    double prev_fdr = 1.0;

    for (size_t i = n; i > 0; --i) {
        size_t idx = indices[i - 1];
        double fdr = results[idx].p_value * static_cast<double>(n) / static_cast<double>(i);
        fdr = std::min(fdr, prev_fdr);
        fdr = std::min(fdr, 1.0);
        results[idx].fdr = fdr;
        results[idx].significant = (fdr <= fdr_threshold);
        prev_fdr = fdr;
    }
}

// ============================================================================
// Output
// ============================================================================

void query::write_classification(const std::string& path,
                                  const std::vector<query_result>& results) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open output file: " + path);
        return;
    }

    // Collect all sample entries for column headers
    auto& registry = sample_registry::instance();
    std::vector<uint32_t> sample_ids;
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto& info = registry.get(static_cast<uint32_t>(i));
        if (info.type != "replicate") {
            sample_ids.push_back(static_cast<uint32_t>(i));
        }
    }

    // Header
    out << "transcript_id\tgene_id\tgene_name\t"
        << "structural_category\tsubcategory\t"
        << "junction_match_score\tmatching_junctions\t"
        << "query_junctions\tref_junctions\t"
        << "known_donors\tknown_acceptors\t"
        << "novel_donors\tnovel_acceptors\t"
        << "n_samples";

    // Per-sample presence columns
    for (uint32_t sid : sample_ids) {
        const auto& info = registry.get(sid);
        out << "\t" << info.id << ".present";
    }

    // Per-sample expression columns
    for (uint32_t sid : sample_ids) {
        const auto& info = registry.get(sid);
        if (info.type == "sample") {
            std::string expr_label = info.has_expression_type()
                ? expr_type_label(info.expr_type)
                : "expression";
            out << "\t" << info.id << "." << expr_label;
        }
    }
    out << "\n";

    // Rows
    for (const auto& r : results) {
        out << r.transcript_id << "\t"
            << r.gene_id << "\t"
            << r.gene_name << "\t"
            << to_string(r.category) << "\t"
            << r.subcat.to_string() << "\t"
            << r.junction_match_score << "\t"
            << r.matching_junctions << "\t"
            << r.total_query_junctions << "\t"
            << r.total_ref_junctions << "\t"
            << r.known_donors << "\t"
            << r.known_acceptors << "\t"
            << r.novel_donors << "\t"
            << r.novel_acceptors << "\t"
            << r.sample_ids.size();

        // Per-sample presence
        std::set<uint32_t> present_set(r.sample_ids.begin(), r.sample_ids.end());
        for (uint32_t sid : sample_ids) {
            out << "\t" << (present_set.count(sid) ? "1" : "0");
        }

        // Per-sample expression
        for (uint32_t sid : sample_ids) {
            const auto& info = registry.get(sid);
            if (info.type == "sample") {
                auto it = r.sample_expression.find(sid);
                if (it != r.sample_expression.end()) {
                    out << "\t" << it->second;
                } else {
                    out << "\t.";
                }
            }
        }
        out << "\n";
    }

    logging::info("Classification written to: " + path);
}

void query::write_dtu_results(const std::string& path,
                               const query_contrast& contrast,
                               const std::vector<dtu_result>& results) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open output file: " + path);
        return;
    }

    out << "gene_id\tgene_name\ttranscript_id\t"
        << "prop_" << contrast.group_a << "\t"
        << "prop_" << contrast.group_b << "\t"
        << "delta_proportion\tp_value\tfdr\tsignificant\n";

    for (const auto& r : results) {
        out << r.gene_id << "\t"
            << r.gene_name << "\t"
            << r.transcript_id << "\t"
            << r.prop_group_a << "\t"
            << r.prop_group_b << "\t"
            << r.delta_proportion << "\t"
            << r.p_value << "\t"
            << r.fdr << "\t"
            << (r.significant ? "yes" : "no") << "\n";
    }

    logging::info("DTU results written to: " + path);
}

void query::write_summary(const std::string& path,
                           const std::vector<query_result>& classification,
                           const std::vector<std::pair<query_contrast,
                               std::vector<dtu_result>>>& dtu_results) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open summary file: " + path);
        return;
    }

    out << "# Atroplex Query Summary\n\n";

    // Classification summary
    out << "## Structural Classification\n";
    out << "Total transcripts: " << classification.size() << "\n\n";

    std::map<structural_category, size_t> counts;
    for (const auto& r : classification) {
        counts[r.category]++;
    }

    auto pct = [&](size_t count) {
        if (classification.empty()) return 0.0;
        return 100.0 * static_cast<double>(count) /
               static_cast<double>(classification.size());
    };

    for (const auto& [cat, count] : counts) {
        out << to_string(cat) << ": " << count
            << " (" << pct(count) << "%)\n";
    }

    // DTU summary
    if (!dtu_results.empty()) {
        out << "\n## Differential Transcript Usage\n";
        for (const auto& [contrast, results] : dtu_results) {
            size_t sig = std::count_if(results.begin(), results.end(),
                [](const dtu_result& r) { return r.significant; });
            out << contrast.group_a << " vs " << contrast.group_b << ": "
                << sig << "/" << results.size() << " significant\n";
        }
    }

    logging::info("Summary written to: " + path);
}

} // namespace subcall
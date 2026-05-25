/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/subcall.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "utility.hpp"
#include "build/builder.hpp"
#include "build/build_gff.hpp"
#include "sample_manifest.hpp"

namespace subcall {

void subcall::add_common_options(cxxopts::Options& options) {
    options.add_options("Common")
        ("o,output-dir", "Output directory for results (default: input file directory)",
            cxxopts::value<std::string>())
        ("p,prefix", "Output file prefix (default: derived from manifest or first input file)",
            cxxopts::value<std::string>())
        ("t,threads", "Number of threads (0 = auto-detect)",
            cxxopts::value<uint32_t>()->default_value("1"))
        ("progress", "Show progress output")
        ("h,help", "Show help message")
        ;

    options.add_options("Grove")
        ("g,genogrove", "Directory containing a pre-built genogrove index (.ggx + .qtx)",
            cxxopts::value<std::string>())
        ("m,manifest", "Sample manifest file (TSV with metadata)",
            cxxopts::value<std::string>())
        ("b,build-from", "Build grove from file(s) without metadata (GFF/GTF)",
            cxxopts::value<std::vector<std::string>>())
        ("k,order", "Genogrove tree order",
            cxxopts::value<int>()->default_value("3"))
        ("min-counts", "Minimum `counts` value to keep a transcript that carries the `counts` attribute. Disabled by default.",
            cxxopts::value<float>()->default_value("-1"))
        ("min-TPM", "Minimum `TPM` value to keep a transcript that carries the `TPM` attribute. Disabled by default.",
            cxxopts::value<float>()->default_value("-1"))
        ("min-FPKM", "Minimum `FPKM` value to keep a transcript that carries the `FPKM` attribute. Disabled by default.",
            cxxopts::value<float>()->default_value("-1"))
        ("min-RPKM", "Minimum `RPKM` value to keep a transcript that carries the `RPKM` attribute. Disabled by default.",
            cxxopts::value<float>()->default_value("-1"))
        ("min-cov", "Minimum `cov` value to keep a transcript that carries the `cov` attribute. Disabled by default.",
            cxxopts::value<float>()->default_value("-1"))
        ("min-CPM", "Minimum `CPM` value to keep a transcript that carries the `CPM` attribute. Disabled by default.",
            cxxopts::value<float>()->default_value("-1"))
        ("filter-precedence", "Switch expression filtering to precedence mode. Comma-separated "
            "list of expression types (e.g. TPM,FPKM,counts). For each transcript, walk this list "
            "and evaluate ONLY the threshold for the first type the transcript carries; ignore "
            "other --min-X thresholds. Default (omitted): AND mode — every --min-X threshold the "
            "user set is evaluated against the transcript's matching attribute; drop if any fails.",
            cxxopts::value<std::string>())
        ("no-absorb", "Disable ISM (Incomplete Splice Match) segment absorption into longer parents")
        ("fuzzy-tolerance", "Max bp difference for fuzzy exon boundary matching during absorption (0 = exact only)",
            cxxopts::value<size_t>()->default_value("5"))
        ("prune-tombstones", "Physically remove absorbed (tombstoned) segments from the grove after build. "
            "Produces a smaller/cleaner .ggx at the cost of a slow post-build sweep "
            "(grove.remove_key is O(E) per call today). Default: off — tombstones stay in the tree and are filtered at query time.")
        ("include-scaffolds", "Keep transcripts on unplaced scaffolds, alt contigs, fix patches, "
            "and decoy sequences in addition to the main chromosomes (chr1..chr22, chrX, chrY, chrM). "
            "Default: off — scaffold features are filtered at ingest to keep the pan-transcriptome index "
            "focused on main-chromosome biology. Enable for non-human/non-mouse species or when you "
            "specifically need scaffold contributions.")
        ("annotated-loci-only", "Only keep sample transcripts that spatially overlap an annotation "
            "segment in the grove. Novel intergenic and antisense gene loci (which are predominantly "
            "long-read artifacts) are discarded. Requires at least one annotation entry in the manifest. "
            "Default: off — all transcripts are indexed.")
        ("chromosomes", "Restrict index to specific chromosomes (comma-separated, e.g. chr1,chr22,chrX). "
            "Features on unlisted chromosomes are silently skipped at ingest. Accepts both prefixed "
            "(chr1) and bare (1) names. Default: all chromosomes.",
            cxxopts::value<std::string>())
        ;
}

void subcall::apply_common_options(const cxxopts::ParseResult& args) {
    if (args.count("progress")) {
        logging::set_progress_enabled(true);
    }
}

build_options subcall::apply_grove_options(const cxxopts::ParseResult& args) {
    build_options opts;
    opts.threads = args["threads"].as<uint32_t>();

    // Expression thresholds (-1 = disabled).
    opts.filters.min_counts = args["min-counts"].as<float>();
    opts.filters.min_TPM    = args["min-TPM"].as<float>();
    opts.filters.min_FPKM   = args["min-FPKM"].as<float>();
    opts.filters.min_RPKM   = args["min-RPKM"].as<float>();
    opts.filters.min_cov    = args["min-cov"].as<float>();
    opts.filters.min_CPM    = args["min-CPM"].as<float>();

    // Optional precedence list. Empty = AND mode; non-empty = walk
    // user-given order, evaluate only the first carried type's threshold.
    // Throws on unknown type tokens so CLI typos fail fast.
    if (args.count("filter-precedence")) {
        std::string prec_arg = args["filter-precedence"].as<std::string>();
        std::istringstream ss(prec_arg);
        std::string token;
        while (std::getline(ss, token, ',')) {
            token.erase(0, token.find_first_not_of(' '));
            token.erase(token.find_last_not_of(' ') + 1);
            if (token.empty()) continue;
            auto t = sample_info::expression_type_from_attribute(token);
            if (t == sample_info::expression_type::UNKNOWN) {
                throw std::runtime_error(
                    "--filter-precedence: unrecognized expression type '" + token +
                    "' (expected one of: counts, TPM, FPKM, RPKM, cov, CPM)");
            }
            // Dedup: later occurrences of the same type are dead (the
            // first occurrence always wins in the build_gff precedence
            // walk). Silently drop them instead of letting the vector
            // grow with no-op entries.
            if (std::find(opts.filter_precedence.begin(),
                          opts.filter_precedence.end(), t)
                == opts.filter_precedence.end()) {
                opts.filter_precedence.push_back(t);
            }
        }
    }

    opts.absorb              = !args.count("no-absorb");
    opts.fuzzy_tolerance     = args["fuzzy-tolerance"].as<size_t>();
    opts.prune_tombstones    = args.count("prune-tombstones") > 0;
    opts.include_scaffolds   = args.count("include-scaffolds") > 0;
    opts.annotated_loci_only = args.count("annotated-loci-only") > 0;

    // Comma-separated chromosome allowlist; normalized so chr-prefixed
    // and bare names both work.
    if (args.count("chromosomes")) {
        std::string chr_arg = args["chromosomes"].as<std::string>();
        std::istringstream ss(chr_arg);
        std::string token;
        while (std::getline(ss, token, ',')) {
            token.erase(0, token.find_first_not_of(' '));
            token.erase(token.find_last_not_of(' ') + 1);
            if (!token.empty()) {
                opts.chromosomes_filter.insert(normalize_chromosome(token));
            }
        }
    }

    return opts;
}

std::filesystem::path subcall::resolve_output_dir(const cxxopts::ParseResult& args,
                                                   const std::string& fallback_input_path) const {
    std::filesystem::path dir;

    if (args.count("output-dir")) {
        dir = args["output-dir"].as<std::string>();
    } else if (!fallback_input_path.empty()) {
        dir = std::filesystem::path(fallback_input_path).parent_path();
    }

    if (dir.empty()) {
        dir = std::filesystem::current_path();
    }

    std::filesystem::create_directories(dir);
    return dir;
}

void subcall::run(const cxxopts::ParseResult& args) {
    validate(args);
    apply_common_options(args);
    setup_grove(args);
    write_build_summary(args);
    execute(args);
}

/// Try to open the quantification sidecar at the path implied by the
/// supplied grove path: replace the `.ggx` extension with `.qtx`. Sets
/// the subcall's `qtx_reader` member when successful; leaves it empty
/// (with an INFO log) when the file is absent or unreadable. Never
/// throws — sidecar absence is a legitimate state, not an error.
static void try_open_qtx_for(const std::filesystem::path& grove_path,
                             std::optional<quant_sidecar::Reader>& reader_out) {
    std::filesystem::path qtx_path = grove_path;
    qtx_path.replace_extension(".qtx");
    if (!std::filesystem::exists(qtx_path)) {
        logging::info("No .qtx sidecar found alongside grove (looked for "
                      + qtx_path.string()
                      + "); expression columns will be omitted");
        return;
    }
    try {
        reader_out.emplace(qtx_path);
        logging::info("Loaded quantification sidecar: " + qtx_path.string()
                      + " (" + std::to_string(reader_out->samples().size())
                      + " samples, "
                      + std::to_string(reader_out->segment_block_count())
                      + " segment blocks)");
    } catch (const std::exception& e) {
        logging::warning(std::string("Failed to open .qtx sidecar at ")
                         + qtx_path.string() + ": " + e.what()
                         + " — expression columns will be omitted");
        reader_out.reset();
    }
}

void subcall::setup_grove(const cxxopts::ParseResult& args) {
    if (args.count("genogrove")) {
        std::filesystem::path grove_dir = args["genogrove"].as<std::string>();

        if (!std::filesystem::is_directory(grove_dir)) {
            throw std::runtime_error(
                "-g/--genogrove expects a directory, not a file: " +
                grove_dir.string());
        }

        // Scan for exactly one .ggx in the directory.
        std::filesystem::path gg_path;
        for (const auto& entry : std::filesystem::directory_iterator(grove_dir)) {
            if (entry.path().extension() == ".ggx") {
                if (!gg_path.empty()) {
                    throw std::runtime_error(
                        "Multiple .ggx files in directory: " + grove_dir.string());
                }
                gg_path = entry.path();
            }
        }
        if (gg_path.empty()) {
            throw std::runtime_error(
                "No .ggx file found in directory: " + grove_dir.string());
        }

        logging::info("Loading grove from: " + gg_path.string());
        load_grove(gg_path.string());
        try_open_qtx_for(gg_path, qtx_reader);
        return;  // Grove loaded from file, skip building
    }

    int order = args["order"].as<int>();

    // Collect samples from both sources
    std::vector<sample_info> all_samples;

    // From manifest (metadata in TSV)
    if (args.count("manifest")) {
        std::string manifest_path = args["manifest"].as<std::string>();
        logging::info("Loading manifest: " + manifest_path);
        sample_manifest manifest(manifest_path);
        logging::info("Found " + std::to_string(manifest.size()) + " sample(s) in manifest");

        for (const auto& info : manifest) {
            all_samples.push_back(info);
        }
    }

    // From build-from (metadata parsed from GFF headers)
    if (args.count("build-from")) {
        auto build_files = args["build-from"].as<std::vector<std::string>>();
        logging::info("Parsing headers from " + std::to_string(build_files.size()) + " file(s)");

        for (const auto& filepath : build_files) {
            sample_info info = build_gff::parse_header(filepath);
            all_samples.push_back(std::move(info));
        }
    }

    // Build grove if we have samples
    if (!all_samples.empty()) {
        build_options opts = apply_grove_options(args);

        logging::info("Creating grove with order: " + std::to_string(order));

        if (opts.filters.any_active()) {
            std::string note = "Expression filters:";
            if (opts.filters.min_counts >= 0) note += " counts>=" + std::to_string(opts.filters.min_counts);
            if (opts.filters.min_TPM    >= 0) note += " TPM>="    + std::to_string(opts.filters.min_TPM);
            if (opts.filters.min_FPKM   >= 0) note += " FPKM>="   + std::to_string(opts.filters.min_FPKM);
            if (opts.filters.min_RPKM   >= 0) note += " RPKM>="   + std::to_string(opts.filters.min_RPKM);
            if (opts.filters.min_cov    >= 0) note += " cov>="    + std::to_string(opts.filters.min_cov);
            if (opts.filters.min_CPM    >= 0) note += " CPM>="    + std::to_string(opts.filters.min_CPM);
            logging::info(note);
        }
        if (!opts.filter_precedence.empty()) {
            std::string note = "Filter mode: precedence (";
            for (size_t i = 0; i < opts.filter_precedence.size(); ++i) {
                if (i) note += ",";
                note += sample_info::expression_type_attribute(opts.filter_precedence[i]);
            }
            note += ") — evaluate only the first type each transcript carries";
            logging::info(note);
            if (!opts.filters.any_active()) {
                logging::warning("--filter-precedence set but no --min-X thresholds active — "
                                 "filter will be a no-op");
            }
        } else if (opts.filters.any_active()) {
            logging::info("Filter mode: AND — drop transcript if any present-and-thresholded type fails");
        }

        if (!opts.absorb) {
            logging::info("ISM segment absorption disabled");
        } else if (opts.fuzzy_tolerance > 0) {
            logging::info("Fuzzy absorption tolerance: " + std::to_string(opts.fuzzy_tolerance) + "bp");
        }
        if (opts.prune_tombstones) {
            logging::info("Physical tombstone removal enabled (--prune-tombstones)");
        }
        if (opts.include_scaffolds) {
            logging::info("Scaffold filter disabled (--include-scaffolds): all seqids retained");
        } else {
            logging::info("Scaffold filter active: main chromosomes only (chr1..chr22, chrX, chrY, chrM)");
        }
        if (!opts.chromosomes_filter.empty()) {
            std::string chr_list;
            for (const auto& c : opts.chromosomes_filter) {
                if (!chr_list.empty()) chr_list += ", ";
                chr_list += c;
            }
            logging::info("Chromosome filter active: " + chr_list);
        }
        if (opts.annotated_loci_only) {
            logging::info("Annotated-loci-only mode: sample transcripts at novel loci will be discarded");
        }
        grove = std::make_unique<grove_type>(order);

        // Derive the final single-file quantification sidecar path.
        {
            std::string fallback_for_outdir;
            if (args.count("manifest")) {
                fallback_for_outdir = args["manifest"].as<std::string>();
            } else if (args.count("build-from")) {
                fallback_for_outdir = args["build-from"].as<std::vector<std::string>>().front();
            }
            auto out_dir = resolve_output_dir(args, fallback_for_outdir);
            std::string prefix = resolve_prefix(args);
            if (!out_dir.empty()) {
                opts.qtx_path = (out_dir / (prefix + ".qtx")).string();
            }
        }

        auto build_start = std::chrono::steady_clock::now();
        build_stats = builder::build_from_samples(*grove, all_samples, opts, nullptr);
        auto build_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - build_start).count();
        build_stats->build_time_seconds = build_elapsed;

        // Record build parameters for summary provenance
        build_stats->build_parameters["order"] = std::to_string(order);
        build_stats->build_parameters["absorb"] = opts.absorb ? "yes" : "no";
        if (opts.absorb && opts.fuzzy_tolerance > 0)
            build_stats->build_parameters["fuzzy_tolerance"] = std::to_string(opts.fuzzy_tolerance) + "bp";
        build_stats->build_parameters["include_scaffolds"] = opts.include_scaffolds ? "yes" : "no";
        build_stats->build_parameters["annotated_loci_only"] = opts.annotated_loci_only ? "yes" : "no";
        if (!opts.chromosomes_filter.empty()) {
            std::string chr_list;
            for (const auto& c : opts.chromosomes_filter) {
                if (!chr_list.empty()) chr_list += ",";
                chr_list += c;
            }
            build_stats->build_parameters["chromosomes"] = chr_list;
        }
        if (opts.filters.min_counts >= 0)
            build_stats->build_parameters["min_counts"] = std::to_string(static_cast<int>(opts.filters.min_counts));
        if (opts.filters.min_TPM >= 0)
            build_stats->build_parameters["min_TPM"] = std::to_string(static_cast<int>(opts.filters.min_TPM));
        if (opts.filters.min_FPKM >= 0)
            build_stats->build_parameters["min_FPKM"] = std::to_string(static_cast<int>(opts.filters.min_FPKM));
        if (opts.filters.min_RPKM >= 0)
            build_stats->build_parameters["min_RPKM"] = std::to_string(static_cast<int>(opts.filters.min_RPKM));
        if (opts.filters.min_cov >= 0)
            build_stats->build_parameters["min_cov"] = std::to_string(static_cast<int>(opts.filters.min_cov));
        if (opts.filters.min_CPM >= 0)
            build_stats->build_parameters["min_CPM"] = std::to_string(static_cast<int>(opts.filters.min_CPM));
        if (!opts.filter_precedence.empty()) {
            std::string prec_list;
            for (auto t : opts.filter_precedence) {
                if (!prec_list.empty()) prec_list += ",";
                prec_list += sample_info::expression_type_attribute(t);
            }
            build_stats->build_parameters["filter_precedence"] = prec_list;
        }

        logging::info("Grove ready with spatial index and graph structure");

        // Open the qtx reader we just produced so the subclass execute()
        // can use it directly without reloading from disk.
        if (!opts.qtx_path.empty()) {
            try_open_qtx_for(opts.qtx_path, qtx_reader);
        }
    }
}

std::string subcall::resolve_prefix(const cxxopts::ParseResult& args) const {
    if (args.count("prefix")) {
        return args["prefix"].as<std::string>();
    }

    // Derive from manifest filename
    if (args.count("manifest")) {
        return std::filesystem::path(args["manifest"].as<std::string>()).stem().string();
    }

    // Derive from first build-from file
    if (args.count("build-from")) {
        return std::filesystem::path(
            args["build-from"].as<std::vector<std::string>>().front()).stem().string();
    }

    // Derive from genogrove index directory (find the .ggx, use its stem).
    if (args.count("genogrove")) {
        std::filesystem::path grove_dir = args["genogrove"].as<std::string>();
        if (std::filesystem::is_directory(grove_dir)) {
            for (const auto& entry : std::filesystem::directory_iterator(grove_dir)) {
                if (entry.path().extension() == ".ggx")
                    return entry.path().stem().string();
            }
        }
        return grove_dir.stem().string();
    }

    return "atroplex";
}

void subcall::load_grove(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open()) {
        throw std::runtime_error("Cannot open index file: " + path);
    }

    // Read magic bytes + version
    char magic[4];
    ifs.read(magic, 4);
    if (!ifs.good()) {
        throw std::runtime_error("Failed to read .ggx header (file too short?): " + path);
    }
    if (std::string(magic, 4) != "AGRX") {
        throw std::runtime_error("Invalid .ggx file (bad magic): " + path);
    }
    uint16_t version;
    ifs.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (!ifs.good()) {
        throw std::runtime_error("Failed to read .ggx version (file truncated?): " + path);
    }
    if (version != 1) {
        throw std::runtime_error("Unsupported .ggx version: " + std::to_string(version));
    }

    // Deserialize registries (order must match save_grove)
    (void)gene_registry::deserialize(ifs);
    source_registry::deserialize(ifs);
    (void)transcript_registry::deserialize(ifs);
    (void)sample_registry::deserialize(ifs);

    // Deserialize grove (handles its own zlib decompression)
    grove = std::make_unique<grove_type>(grove_type::deserialize(ifs));

    logging::info("Loaded index from: " + path);
    logging::info("Chromosomes: " + std::to_string(grove->get_root_nodes().size()) +
                  ", Genes: " + std::to_string(gene_registry::instance().size()) +
                  ", Sources: " + std::to_string(source_registry::instance().size()) +
                  ", Transcripts: " + std::to_string(transcript_registry::instance().size()) +
                  ", Samples: " + std::to_string(sample_registry::instance().size()));
}

void subcall::build_grove(const std::vector<std::string>& files, int order, uint32_t threads) {
    logging::info("Creating grove with order: " + std::to_string(order));
    grove = std::make_unique<grove_type>(order);
    auto build_start = std::chrono::steady_clock::now();
    build_stats = builder::build_from_files(*grove, files, threads);
    auto build_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - build_start).count();
    build_stats->build_time_seconds = build_elapsed;
    logging::info("Grove ready with spatial index and graph structure");
}

void subcall::save_grove(const std::string& path) {
    if (!grove) {
        logging::warning("No grove to save");
        return;
    }

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.is_open()) {
        throw std::runtime_error("Cannot create index file: " + path);
    }

    // Write magic bytes + version
    ofs.write("AGRX", 4);  // Atroplex GRove indeX
    uint16_t version = 1;
    ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));

    // Serialize registries (order must match load_grove)
    gene_registry::instance().serialize(ofs);
    source_registry::instance().serialize(ofs);
    transcript_registry::instance().serialize(ofs);
    sample_registry::instance().serialize(ofs);

    // Serialize grove (handles its own zlib compression)
    grove->serialize(ofs);

    logging::info("Index saved to: " + path + " (" +
                  std::to_string(std::filesystem::file_size(path) / 1024) + " KB)");
}

void subcall::write_build_summary(const cxxopts::ParseResult& args) {
    if (!build_stats) return;

    std::string prefix = resolve_prefix(args);
    std::string fallback;
    if (args.count("manifest")) {
        fallback = args["manifest"].as<std::string>();
    } else if (args.count("build-from")) {
        fallback = args["build-from"].as<std::vector<std::string>>().front();
    }

    auto out_dir = resolve_output_dir(args, fallback);
    std::string summary_path = (out_dir / (prefix + ".ggx.summary")).string();
    build_stats->write_summary(summary_path);
}

} // namespace subcall
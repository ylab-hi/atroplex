/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of genogrove and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_UTILITY_HPP
#define ATROPLEX_UTILITY_HPP

// standard
#include <chrono>
#include <string>

/**
 * Normalize chromosome name to UCSC/GENCODE style
 * - "1" → "chr1", "X" → "chrX", "Y" → "chrY"
 * - "MT" → "chrM" (Ensembl mitochondria)
 * - "chr1" → "chr1" (already normalized)
 */
std::string normalize_chromosome(const std::string& seqid);

/**
 * Return true if `seqid` names a canonical main chromosome on the
 * UCSC/GENCODE human assembly — chr1..chr22, chrX, chrY, chrM (or the
 * unprefixed forms 1..22, X, Y, M, MT). Everything else — unplaced
 * scaffolds (`chr1_KI270706v1_random`), unlocalized contigs
 * (`chrUn_*`), fix patches (`*_fix`), alt haplotypes (`*_alt`), decoy
 * sequences (`chrEBV`), and so on — returns false.
 *
 * The main-chromosome set is hard-coded for human; mouse and other
 * species with different naming conventions should pass
 * `include_scaffolds=true` to disable the filter entirely.
 *
 * Intended to be called with a raw GTF seqid **before** it is run
 * through `normalize_chromosome()`, so that both the prefixed and
 * unprefixed forms match. Scaffold names with extra suffixes (the
 * long tail that blows up pan-transcriptome catalog size at cohort
 * scale) all fail the strict match and are filtered out.
 *
 * @param seqid           The chromosome name as read from the GTF
 * @param include_scaffolds  If true, always return true (no filtering)
 */
bool is_main_chromosome(const std::string& seqid, bool include_scaffolds = false);

namespace logging {
    void info(const std::string& message);
    void warning(const std::string& message);
    void error(const std::string& message);

    /**
     * Enable or disable progress output (carriage return updates)
     * Progress is disabled by default for compatibility with non-interactive
     * environments (SLURM, CI, etc.). Use --progress flag to enable.
     * @param enabled Whether to show progress updates
     */
    void set_progress_enabled(bool enabled);

    /**
     * Display progress with carriage return (overwrites current line)
     * Shows line count with comma formatting and lines/sec rate
     * Does nothing if progress is disabled via set_progress_enabled(false)
     * @param lines Number of lines read
     * @param prefix Prefix message (e.g., "Processing gencode.gtf")
     */
    void progress(size_t lines, const std::string& prefix = "Processing");

    /**
     * Clear the progress line and print final summary
     * @param segments Total segments created
     * @param prefix Prefix message
     */
    void progress_done(size_t segments, const std::string& prefix = "Processed");

    /**
     * Reset the progress timer (call before starting a new progress sequence)
     */
    void progress_start();
}

#endif //ATROPLEX_UTILITY_HPP

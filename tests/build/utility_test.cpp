/*
 * Unit tests for `utility.hpp` helpers used by the build path.
 *
 * Currently covers `is_main_chromosome`, which gates whether GFF/BAM
 * entries enter the grove. Edge cases matter here because a false
 * positive silently drops biology, and a false negative blows up the
 * cohort-scale catalog with scaffold noise.
 */

#include <gtest/gtest.h>

#include "utility.hpp"

// ── Prefixed autosomes ──────────────────────────────────────────────

TEST(IsMainChromosomeTest, PrefixedAutosomes) {
    EXPECT_TRUE(is_main_chromosome("chr1"));
    EXPECT_TRUE(is_main_chromosome("chr9"));
    EXPECT_TRUE(is_main_chromosome("chr10"));
    EXPECT_TRUE(is_main_chromosome("chr22"));
}

TEST(IsMainChromosomeTest, UnprefixedAutosomes) {
    EXPECT_TRUE(is_main_chromosome("1"));
    EXPECT_TRUE(is_main_chromosome("9"));
    EXPECT_TRUE(is_main_chromosome("10"));
    EXPECT_TRUE(is_main_chromosome("22"));
}

// ── Sex chromosomes + mitochondria ──────────────────────────────────

TEST(IsMainChromosomeTest, SexAndMito) {
    EXPECT_TRUE(is_main_chromosome("chrX"));
    EXPECT_TRUE(is_main_chromosome("chrY"));
    EXPECT_TRUE(is_main_chromosome("chrM"));
    EXPECT_TRUE(is_main_chromosome("X"));
    EXPECT_TRUE(is_main_chromosome("Y"));
    EXPECT_TRUE(is_main_chromosome("M"));
    EXPECT_TRUE(is_main_chromosome("MT"));  // Ensembl mitochondria
}

// ── Out-of-range autosomes ──────────────────────────────────────────

TEST(IsMainChromosomeTest, RejectsOutOfRangeAutosomes) {
    EXPECT_FALSE(is_main_chromosome("chr0"));
    EXPECT_FALSE(is_main_chromosome("0"));
    EXPECT_FALSE(is_main_chromosome("chr23"));
    EXPECT_FALSE(is_main_chromosome("23"));
    EXPECT_FALSE(is_main_chromosome("chr100"));
    EXPECT_FALSE(is_main_chromosome("chr999"));
}

// ── Scaffolds / alt haplotypes / fix patches / decoys ──────────────

TEST(IsMainChromosomeTest, RejectsUnplacedContigs) {
    EXPECT_FALSE(is_main_chromosome("chr1_KI270706v1_random"));
    EXPECT_FALSE(is_main_chromosome("chr14_GL000009v2_random"));
    EXPECT_FALSE(is_main_chromosome("chrUn_KI270302v1"));
    EXPECT_FALSE(is_main_chromosome("chr11_KI270721v1_random"));
}

TEST(IsMainChromosomeTest, RejectsAltAndFixPatches) {
    EXPECT_FALSE(is_main_chromosome("chr1_KN196474v1_fix"));
    EXPECT_FALSE(is_main_chromosome("chr19_KI270938v1_alt"));
    EXPECT_FALSE(is_main_chromosome("HG2232_PATCH"));
}

TEST(IsMainChromosomeTest, RejectsDecoys) {
    EXPECT_FALSE(is_main_chromosome("chrEBV"));
    EXPECT_FALSE(is_main_chromosome("hs37d5"));
    EXPECT_FALSE(is_main_chromosome("NC_007605"));
}

// ── Lexical/format edge cases ───────────────────────────────────────

TEST(IsMainChromosomeTest, RejectsEmpty) {
    EXPECT_FALSE(is_main_chromosome(""));
}

TEST(IsMainChromosomeTest, RejectsBarePrefix) {
    EXPECT_FALSE(is_main_chromosome("chr"));
}

TEST(IsMainChromosomeTest, RejectsLeadingZeroAutosomes) {
    // Canonical GTFs never use zero-padded autosome names; rejecting
    // them keeps the filter strict and avoids accidentally matching
    // unusual fixture names.
    EXPECT_FALSE(is_main_chromosome("01"));
    EXPECT_FALSE(is_main_chromosome("chr01"));
    EXPECT_FALSE(is_main_chromosome("022"));
    EXPECT_FALSE(is_main_chromosome("chr022"));
}

TEST(IsMainChromosomeTest, RejectsMixedAlphanumericSuffix) {
    // First non-digit kicks the autosome check out — these all fail.
    EXPECT_FALSE(is_main_chromosome("1X"));
    EXPECT_FALSE(is_main_chromosome("chr1a"));
    EXPECT_FALSE(is_main_chromosome("22-extra"));
}

TEST(IsMainChromosomeTest, RejectsLowercaseSex) {
    // Intentionally strict — GFFs use uppercase X/Y/M. Lowercase is
    // almost always a pipeline bug that should not silently match.
    EXPECT_FALSE(is_main_chromosome("chrx"));
    EXPECT_FALSE(is_main_chromosome("y"));
    EXPECT_FALSE(is_main_chromosome("chrm"));
}

// ── include_scaffolds override ──────────────────────────────────────

TEST(IsMainChromosomeTest, IncludeScaffoldsOverridesEverything) {
    EXPECT_TRUE(is_main_chromosome("chr1_KI270706v1_random", /*include_scaffolds=*/true));
    EXPECT_TRUE(is_main_chromosome("chrEBV",                 /*include_scaffolds=*/true));
    EXPECT_TRUE(is_main_chromosome("",                       /*include_scaffolds=*/true));
    EXPECT_TRUE(is_main_chromosome("chr100",                 /*include_scaffolds=*/true));
    EXPECT_TRUE(is_main_chromosome("anything_literally",     /*include_scaffolds=*/true));
}

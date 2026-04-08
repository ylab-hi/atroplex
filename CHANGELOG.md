# Changelog

All notable changes to atroplex will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **Splice site hash collisions**: replaced weak XOR-based hash (`seqid_hash ^ position<<1`) with proper `(seqid, position)` composite keys using boost::hash_combine pattern — eliminates false positive "known" splice site hits across chromosomes that caused NIC/NNC misclassification
- Fixed potential `size_t` underflow in splice site tolerance window when position is near zero
- Fixed infinite loop in `sample_bitset` iterator: `advance()` did not mark iterator as exhausted when all bits were consumed, causing range-based for loops over `sample_idx` to hang
- Fixed missing `absorb` parameter in final `process_gene()` call per GTF file (last gene always skipped absorption)

### Added
- Discover test suite: 11 unit tests constructing `read_cluster` objects directly for each SQANTI category (FSM, ISM prefix/suffix, NIC, NNC, GENIC_INTRON, GENIC_GENOMIC, INTERGENIC), plus 6 integration tests exercising the full SAM → cluster → classify pipeline
- Query classification test suite: 8 tests verifying per-sample presence tracking (shared, annotation-only, sample-only, novel), expression propagation, and gene assignment
- Test fixtures: minimal chr22 reference GTF (2 genes, 3 transcripts), sample GTF with expression, query GTF, and 23-read SAM file covering all classification categories
- Updated genogrove dependency to v0.20.0 (fixes B+ tree `split_node` corruption, adds grove iteration API and serialization)
- Updated cxxopts dependency to v3.3.1
- Adapted to genogrove v0.17.0+ API: `gff_entry`/`sam_entry` use direct `start`/`end` fields; `gff_entry.attributes` uses `std::less<>` transparent comparator

### C++20 Modernization
- Replaced GCC `__builtin_popcountll`/`__builtin_ctzll`/`__builtin_ctz`/`__builtin_popcount` with portable `std::popcount` and `std::countr_zero` from `<bit>`
- Used `std::string_view` + `std::from_chars` in `chromosome_compare` to eliminate heap allocations in sort comparator
- Used `std::string_view` in `build_gff::parse_header` trim lambda
- Added `[[nodiscard]]` to `segment_builder`, `transcript_matcher`, and `read_clusterer` return-value functions
- Fixed `read_clusterer::stats::max_cluster_size` type from `double` to `size_t`

### Changed
- **Read clustering**: replaced binned clustering with sort+sweep algorithm — eliminates boundary artifacts where reads 1bp apart could land in different bins and never be compared; sorts by `(strand, junction_count, first_donor)` then sweeps within `junction_tolerance`
- Removed `--junction-bin` CLI option (binning is no longer used internally)
- Single-pass candidate scoring in `transcript_matcher::match()` — avoids redundant graph traversals for ambiguity detection
- CI workflow now runs `ctest` after build (tests were compiled but never executed)
- **Absorption rules v2**: rewrote ISM absorption as 9 structured rules (0-8) with ref/sample distinction, fuzzy subsequence matching, and proper mono-exon classification (see `absorption_rules.txt`)
- Annotations are now always processed before samples during grove construction (sorted in `builder::build_from_samples`)
- Transcripts within a gene are sorted by exon count (descending) before processing, ensuring multi-exon segments exist before mono-exon fragments are checked
- `classify_subsequence()` now distinguishes 3' ISM (1-2 missing exons, absorb) from 3' degradation (3+ missing, drop vs ref / keep vs sample)

### Refactored
- **Extracted `segment_builder`** from `build_gff`: all segment creation and absorption logic moved to format-agnostic `segment_builder.hpp/cpp`, enabling future BAM input support without duplicating absorption rules
- `build_gff::process_transcript()` now pre-computes seqid/strand/span and delegates to `segment_builder::create_segment()`
- Named constants for absorption thresholds (`TERMINAL_TOLERANCE_BP`, `FUZZY_TOLERANCE_BP`, `MAX_ISM_MISSING_EXONS`)
- Hoisted `gene_index.find()` lookup in `create_segment()` — single lookup shared across Steps 2/3/4
- Eliminated unnecessary vector copy in `process_gene()` transcript sorting (sort keys only, not full entry vectors)
- **Absorption rule 0 (FSM)**: identical exon coordinates merged (pointer identity, then fuzzy ≤5bp fallback)
- **Absorption rule 1 (5' ISM)**: contiguous subset at 5' end kept as separate segment (potential APA)
- **Absorption rule 2 (3' ISM)**: contiguous subset at 3' end, missing 1-2 exons, absorbed (RT dropout)
- **Absorption rule 3 (3' degradation)**: contiguous subset at 3' end, missing 3+ exons, dropped vs reference / kept vs sample
- **Absorption rule 4 (internal fragment)**: contiguous subset missing exons from both ends, dropped vs reference / kept vs sample
- **Absorption rule 5 (terminal variant)**: same intron chain with TSS/TES differing <50bp, absorbed
- **Absorption rule 6 (mono-exon gene overlap)**: mono-exon overlapping a multi-exon gene without crossing an intron boundary, dropped
- **Absorption rule 7 (mono-exon intron retention)**: mono-exon spanning from one exon across an intron into the next exon, kept
- **Absorption rule 8 (mono-exon intergenic)**: mono-exon with no gene overlap, dropped
- **Fuzzy subsequence matching**: coordinate-proximity fallback (≤5bp tolerance) integrated into Rules 0-4 when pointer identity fails, handling cross-caller coordinate variation
- `reset()` methods on `transcript_registry`, `gene_registry`, and `source_registry` for test isolation
- GoogleTest-based test suite for absorption rules (`tests/build/absorption_test.cpp`)
- Test fixtures for each absorption rule (`tests/build/fixtures/`)
- `absorption_rules.txt` documenting all rules with biological rationale, ref/sample distinction, and execution order
- Biological replicate merging via `--min-replicates N` and manifest `group` column
- Auto-inference of replicate groups from sample IDs (strips `_repNN` suffix)
- Subcommand architecture: `build`, `analyze`, `discover`, `query`
- Variant-based genomic feature system (`exon_feature`, `segment_feature`)
- Edge metadata and graph structure with segment-level edge deduplication
- Gene-by-gene GFF/GTF processing with cross-file feature deduplication
- Sample manifest (TSV) with ENCODE-aligned metadata fields
- GFF header parsing (`##property: value` format) for `--build-from` workflows
- Pan-transcriptome sample tracking via `sample_registry`
- Entry type distinction (`type`: "sample" vs "annotation") affecting stats and output
- Source tracking (GFF column 2: HAVANA, ENSEMBL, TALON, etc.) on all features
- Expression auto-detection from GFF transcript attributes (counts/TPM/FPKM/RPKM/cov)
- Expression type tracking per sample in output column headers
- `--min-expression` filter for excluding low-count transcripts during build
- `--no-absorb` flag to disable ISM segment absorption
- Two-tier statistics: quick build summary from builder caches, full `atroplex analyze` via tree traversal
- Full analysis output: per-sample/per-source CSVs, exon/segment sharing TSVs
- Conserved exon detail TSV (per-exon rows with per-sample transcripts and expression)
- Splicing hub analysis: exons with >10 downstream branches, per-sample entropy, PSI, expression
- Branch detail TSV: per-(hub x target) fraction and expression
- Isoform diversity via Jaccard distance (analyze subcommand)
- Read clustering by splice junction signature from BAM input
- Transcript matching with coarse-to-fine spatial + graph queries
- Splice site indexing from exon caches for NIC/NNC classification
- Novel segment creation during discovery phase
- Query subcommand: classify input transcripts against the index with per-sample presence/expression; optional DTU via `--contrast`
- Index serialization (.ggx format): AGRX magic + version, registries + zlib-compressed grove
- Per-file memory usage logging during grove construction
- Build time reporting in index summary
- R visualization script for statistics (`scripts/visualize_stats.R`)
- Expression injection script (`scripts/inject_expression.py`)
- ENCODE test data with expression-injected GTFs
- Dockerfile for automated builds

### Optimized
- Gene string interning via `gene_registry` (gene_id, gene_name, gene_biotype stored once, features carry `uint32_t`)
- Source tracking via `source_registry` with `uint16_t` bitfield (replaces `unordered_set<string>`)
- Sample tracking via `sample_bitset` dynamic bitset (replaces `unordered_set<uint32_t>`, supports 100+ samples)
- Expression storage via lazy `expression_store` (null pointer for features without expression, flat vector with sentinel)
- Removed stored coordinate strings from features (computed on demand from `genomic_coordinate`)
- Removed unused `overlapping_features` field from `exon_feature`
- `sorted_vec` replaces `unordered_set<uint32_t>` for `transcript_ids` on exons and segments (~13% memory reduction)
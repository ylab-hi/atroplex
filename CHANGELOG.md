# Changelog

All notable changes to atroplex will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Subcommand architecture: `build`, `analyze`, `discover`, `dtu`
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
- ISM segment absorption: contiguous exon-chain subsets merged into parent segments (forward + reverse)
- Mono-exon transcript filtering during build
- Two-tier statistics: quick `--stats` from builder caches, full `atroplex analyze` via tree traversal
- Full analysis output: per-sample/per-source CSVs, exon/segment sharing TSVs
- Conserved exon detail TSV (per-exon rows with per-sample transcripts and expression)
- Splicing hub analysis: exons with >10 downstream branches, per-sample entropy, PSI, expression
- Branch detail TSV: per-(hub x target) fraction and expression
- Isoform diversity via Jaccard distance (analyze subcommand)
- Read clustering by splice junction signature from BAM input
- Transcript matching with coarse-to-fine spatial + graph queries
- Match classification (EXACT, COMPATIBLE, NOVEL_JUNCTION, NOVEL_EXON, INTERGENIC, AMBIGUOUS)
- Novel segment creation during discovery phase
- Per-file memory usage logging during grove construction
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
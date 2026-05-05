# Atroplex

A pan-transcriptome indexing and analysis toolkit for long-read sequencing data. Atroplex builds a spatial index with graph structure from multiple annotation sources, discovers novel transcripts from BAM input, classifies transcripts against the index, exports per-sample GTF files from the index, and provides detailed analysis of exon/segment sharing, splicing complexity, and isoform diversity.

## Installation

### Docker (recommended)

```bash
docker build -t atroplex .
docker run --rm -v $(pwd):/data atroplex --help
```

### From source

Requires CMake 3.14+, a C++20 compiler, and htslib installed via pkg-config.

```bash
cmake -B build -S .
cmake --build build
./build/atroplex --help
```

### Dependencies

- **cxxopts** (v3.3.1): Command-line argument parsing (fetched automatically)
- **genogrove** (v0.21.0): Genomic interval data structures and graph structure (fetched automatically)
- **htslib**: Reading BAM/SAM files (system dependency via pkg-config)
- **zlib**: Compression support

## Quick Start

```bash
# Build index from a single annotation file
atroplex build -b annotation.gtf

# Build pan-transcriptome from multiple sources via manifest
atroplex build -m manifest.tsv

# Build with expression filter
# (manifest's `expression_attribute` column declares which GFF attribute
# to read on each sample — here: `counts` for TALON samples)
atroplex build -m manifest.tsv --min-counts 3

# Build only specific chromosomes (for targeted analysis)
atroplex build -m manifest.tsv --chromosomes chr1,chr22,chrX

# Run full per-sample inspection (sharing, splicing hubs, diversity)
atroplex inspect -m manifest.tsv -o results/

# Classify transcripts against the index
atroplex query -i transcripts.gtf -m manifest.tsv -o results/

# Classify with differential transcript usage
atroplex query -i transcripts.gtf -m manifest.tsv --contrast treated:control

# Discover novel transcripts from long-read data
atroplex discover -i reads.bam -m manifest.tsv -o results/

# Export per-sample GTF files from a built index
atroplex export -g index.ggx -o export/

# Compact a built index by physically removing absorbed segments
atroplex compact -g index_dir/ -o compacted/
```

## Subcommands

### `atroplex build` — Build pan-transcriptome index

Builds a genogrove index from one or more annotation files or BAM files. Each input file is treated as a separate entry in the pan-transcriptome, with features deduplicated across files and sample/source provenance tracked on every exon and segment.

```bash
# From manifest (full metadata per sample)
atroplex build -m ENCODE/manifest.tsv -o results/

# From annotation files directly (metadata parsed from GFF headers)
atroplex build -b gencode.gtf -b sample1.gtf -b sample2.gtf

# With expression filtering — thresholds are per-attribute and each sample
# declares which attributes to read via the manifest's `expression_attribute`
# column. A TALON manifest row with `expression_attribute = counts` is
# filtered by --min-counts; a StringTie row with `expression_attribute = cov,TPM`
# is filtered by both --min-cov and --min-TPM (AND semantics, missing
# attributes on a given transcript are pass-through).
atroplex build -m manifest.tsv --min-counts 3
atroplex build -m manifest.tsv --min-cov 1 --min-TPM 0.5

# Build only specific chromosomes
atroplex build -m manifest.tsv --chromosomes chr1,chr22

# Drop sample transcripts at novel loci (keep only annotated gene regions)
atroplex build -m manifest.tsv --annotated-loci-only

# Disable ISM absorption
atroplex build -m manifest.tsv --no-absorb

# Physical tombstone removal (slower, produces smaller .ggx for distribution)
atroplex build -m manifest.tsv --prune-tombstones
```

Build always produces a serialized index (`.ggx`) and a build summary (`.ggx.summary`).

### `atroplex inspect` — Full pan-transcriptome inspection

Performs detailed per-sample inspection: overview / per-source / biotype
breakdowns, exon & segment sharing, conserved-exon detail, and splicing hubs
(per-sample PSI + entropy + branch fan-out). Splicing events (cassette / alt-5′ / alt-3′ / IR / alt-terminal / mutually-exclusive) are opt-in via `--events`.

```bash
# Full inspection from manifest
atroplex inspect -m ENCODE/manifest.tsv -o results/

# Full inspection from annotation files
atroplex inspect -b gencode.gtf -b sample1.gtf -o results/

# Filter to segments in >= 5 samples (annotations always kept)
atroplex inspect -m manifest.tsv -o results/ --min-samples 5

# Relax the conservation threshold: conserved == present in >= 95% of samples
atroplex inspect -m manifest.tsv -o results/ --conserved-fraction 0.95

# Include splicing event catalog (cassette, alt-5'/3', IR, etc.)
atroplex inspect -m manifest.tsv -o results/ --events
```

Inspect-specific options: `--min-samples` (skip segments in < N samples, annotations always kept), `--conserved-fraction <(0,1]>` (fraction of sample-typed entries a feature must appear in to be classified as conserved; default `1.0` = strict "in every sample"; relax for a dropout-tolerant conserved core), `--events` (write per-gene splicing event catalog, off by default)

### `atroplex query` — Classify transcripts against the index

Classifies input transcripts (GTF/GFF) against the pan-transcriptome index using SQANTI-like structural categories (FSM, ISM, NIC, NNC, etc.). Optionally performs differential transcript usage (DTU) analysis between sample groups.

```bash
# Classify transcripts
atroplex query -i transcripts.gtf -m manifest.tsv -o results/

# With differential transcript usage between groups (chi-squared + BH-FDR)
# group_a and group_b must be values of the manifest's `group` column
# (or auto-inferred from the `_repNN` suffix on sample IDs). DTU requires
# a .qtx sidecar to have been written during build.
atroplex query -i transcripts.gtf -m manifest.tsv --contrast treated:control --fdr 0.05
```

### `atroplex discover` — Discover novel transcripts

Clusters aligned long reads by splice junction signature and matches them against the index to classify transcripts as known, compatible, or novel.

```bash
atroplex discover -i reads.bam -m manifest.tsv -o results/
```

### `atroplex export` — Reconstruct per-sample GTFs from the index

Walks a built index and writes one GTF file per sample with gene, transcript, and exon lines. Expression values (when available) are emitted as GTF attributes. Supports filters to restrict output by sample, gene, region, biotype, source, or sample frequency.

```bash
# Export all samples from a pre-built index (pass the directory containing the .ggx)
atroplex export -g index_dir/ -o export/

# Export a specific sample, protein-coding genes only
atroplex export -g index_dir/ --sample HL60_M1_rep1 --biotype protein_coding

# Export conserved features in a genomic region
atroplex export -g index_dir/ --region chr22:10000000-15000000 --conserved-only

# Export the dropout-tolerant conserved core (segments in >= 95% of samples)
atroplex export -g index_dir/ --conserved-only --conserved-fraction 0.95

# Export features present in at least 3 samples, from HAVANA
atroplex export -g index_dir/ --min-samples 3 --source HAVANA
```

Export-specific options: `--sample`, `--gene`, `--region chr:start-end`, `--min-samples`, `--conserved-only`, `--conserved-fraction <(0,1]>` (sample-typed-entry fraction the segment must hit for `--conserved-only`; default `1.0` = strict; mirrors the inspect convention with annotations excluded from the denominator), `--biotype`, `--source` (all filters are AND'd).

### `atroplex compact` — Compact a built index

Physically removes absorbed (tombstoned) segments from an existing `.ggx`, producing a smaller index without changing query semantics. Use this when an index was built without `--prune-tombstones` and you want to reclaim space later. The companion `.qtx` is already remapped against live segments at build time, so it is copied through unchanged alongside `.ggx.summary`.

```bash
# Compact an existing index into a separate output directory
atroplex compact -g index_dir/ -o compacted/

# Compact a structure-only index (no .qtx alongside the .ggx)
atroplex compact -g index_dir/ -o compacted/ --no-qtx
```

By default, `compact` fails fast if no `.qtx` is found alongside the input `.ggx` — pass `--no-qtx` to opt out. It also refuses to write to the input directory, so a partial write can never overwrite the original.

Compact-specific options: `--no-qtx`.

### Common Options

| Option | Description |
|--------|-------------|
| `-o, --output-dir` | Output directory (default: input file directory) |
| `-p, --prefix` | Output file prefix (default: derived from manifest or first input) |
| `-t, --threads` | Number of threads (default: 1) |
| `--progress` | Show progress output |
| `-m, --manifest` | Sample manifest file (TSV) |
| `-b, --build-from` | Build from GFF/GTF file(s) |
| `-g, --genogrove` | Directory containing a pre-built genogrove index (`.ggx` + optional `.qtx` sidecar) |
| `-k, --order` | Genogrove tree order (default: 3) |
| `--min-counts` | Minimum `counts` value for transcripts whose sample declares `counts` in its manifest `expression_attribute` column (default: -1, disabled) |
| `--min-TPM` | Minimum `TPM` value for transcripts whose sample declares `TPM` (default: -1, disabled) |
| `--min-FPKM` | Minimum `FPKM` value for transcripts whose sample declares `FPKM` (default: -1, disabled) |
| `--min-RPKM` | Minimum `RPKM` value for transcripts whose sample declares `RPKM` (default: -1, disabled) |
| `--min-cov` | Minimum `cov` value for transcripts whose sample declares `cov` (default: -1, disabled) |
| `--no-absorb` | Disable ISM segment absorption into longer parent segments |
| `--fuzzy-tolerance` | Max bp difference for fuzzy exon boundary matching (default: 5) |
| `--prune-tombstones` | Physically remove absorbed segments from the grove post-build (slower, smaller .ggx) |
| `--include-scaffolds` | Keep transcripts on unplaced scaffolds, alt contigs, fix patches, and decoy sequences. Default: off — GFF/BAM ingest is filtered to canonical main chromosomes only (`chr1..chr22`, `chrX`, `chrY`, `chrM`). Enable for non-human/non-mouse species or when you specifically need scaffold contributions. |
| `--chromosomes` | Restrict index to specific chromosomes (comma-separated, e.g. `chr1,chr22,chrX`). Accepts both prefixed and bare names. Default: all chromosomes. |
| `--annotated-loci-only` | Only keep sample transcripts that overlap an annotation segment. Novel intergenic loci are discarded; novel isoforms at annotated loci inherit the annotation gene identity. |

## Input Files

### Sample Manifest (TSV)

A tab-separated file specifying input files with metadata. This is the recommended way to build a pan-transcriptome with full sample tracking.

```tsv
file	id	type	assay	biosample	condition	species	platform	pipeline	expression_attribute	description
gencode.v49.gtf	GENCODE_v49	annotation	.	.	.	Homo sapiens	.	GENCODE	.	Reference annotation
sample1.gtf	HL60_M1_rep1	sample	RNA-seq	HL-60	M1 macrophage	Homo sapiens	PacBio Sequel II	TALON	counts	HL-60 M1 replicate 1
sample2.gtf	HL60_M1_rep2	sample	RNA-seq	HL-60	M1 macrophage	Homo sapiens	PacBio Sequel II	TALON	counts	HL-60 M1 replicate 2
strtie1.gtf	STRTIE_01	sample	RNA-seq	brain	healthy	Homo sapiens	Illumina NovaSeq	StringTie	cov,TPM	StringTie sample with both cov and TPM
```

- Tab-separated, `"."` for empty values (VCF convention)
- `file` column is required; all others optional
- `type`: `"annotation"` for reference annotations, `"sample"` for experimental data (default: `"sample"`)
- Column names are case-insensitive
- Relative paths resolved from manifest directory
- `id` auto-generated from filename if not provided
- Optional `group` column for replicate grouping; if absent, groups are auto-inferred by stripping `_repNN` suffix from sample IDs
- Optional `expression_attribute` column declares which GFF attributes carry quantitative expression for this sample: `counts`, `TPM`, `FPKM`, `RPKM`, `cov`, or a comma-separated list like `cov,TPM` for StringTie samples that carry multiple quantifications. Empty or `.` means no expression filtering for the sample. The FIRST declared attribute is what gets stored per-feature and appears in per-sample output column headers (e.g., `SAMPLE1.counts`, `STRTIE_01.cov`)

The `type` field matters: only entries with `type = "sample"` count toward "conserved" thresholds. By default a feature is conserved when it is present in every sample-typed entry; relax with `--conserved-fraction` (e.g. `0.95` for a dropout-tolerant core). Reference annotations participate in structural analysis but don't inflate sample counts.

### GFF/GTF Files

When using `--build-from` without a manifest, metadata is parsed from GFF/GTF headers:

```
##id: GENCODE_v49_GRCh38
##description: Evidence-based annotation of the human genome
##species: Homo sapiens
##annotation_source: GENCODE
##annotation_version: v49
##type: annotation
```

#### Required Attributes

Each exon feature must have `gene_id` and `transcript_id` in column 9.

#### Expression Attributes

Each sample declares which GFF attribute(s) carry quantitative expression via the manifest's `expression_attribute` column. Supported values: `counts`, `TPM`, `FPKM`, `RPKM`, `cov`. Multiple may be listed comma-separated (e.g., `cov,TPM` for StringTie samples). Empty or `.` disables expression filtering and storage for that sample.

- The **first** declared attribute is what gets stored per-feature via `expression_store` and appears as the column header suffix in per-sample outputs (e.g., `ENCSR_01.counts`, `STRTIE_01.cov`).
- All declared attributes are **evaluated against their matching CLI threshold** (`--min-counts`, `--min-TPM`, etc.) with **AND semantics** — a transcript is kept only if every (declared attribute that has an active threshold) meets it. Missing attributes on a given transcript are pass-through.
- Samples that declare no `expression_attribute` (or use `.`) are never filtered on expression — same pass-through semantics as annotations.
- A CLI threshold that no manifest sample declares emits a warning at build time so you know the filter had no effect.

Example: a mixed ENCODE+StringTie build where TALON samples are filtered on raw read counts and StringTie samples are filtered on StringTie coverage, while GENCODE passes through unfiltered:

```tsv
file             id          type        expression_attribute
gencode.gtf      GENCODE     annotation  .
encsr_01.gtf     ENCSR_01    sample      counts
encsr_02.gtf     ENCSR_02    sample      counts
strtie_01.gtf    STRTIE_01   sample      cov,TPM
```

```bash
atroplex build -m manifest.tsv --min-counts 3 --min-cov 1 --min-TPM 0.5
```

Use `scripts/inject_expression.py` to add TALON quantification counts into GTF files before building.

### BAM/SAM Files

BAM files can be used as input to `build` (via manifest or `--build-from`). Reads are clustered by splice junction signature, and clusters are converted to exon chains and segments using the same absorption rules as GFF input. Read counts serve as expression values.

#### Chromosome Name Normalization

Atroplex normalizes chromosome names to UCSC/GENCODE style automatically:

| Input | Normalized |
|-------|------------|
| `1`, `2`, ... `22` | `chr1`, `chr2`, ... `chr22` |
| `X`, `Y` | `chrX`, `chrY` |
| `MT` | `chrM` |
| `chr1` (already prefixed) | `chr1` (unchanged) |

This allows mixing annotations from different sources (e.g., GENCODE + Ensembl).

## ISM Absorption

During index construction, Incomplete Splice Match (ISM) transcripts — truncated versions of full-length transcripts — are absorbed into their parent segments rather than creating separate entries. This reduces noise from degradation artifacts and technical truncation in long-read data.

Absorption rules (in execution order):

| Rule | Pattern | Action |
|------|---------|--------|
| 0 | Identical exon structure (FSM) | Merge metadata |
| 5 | Same intron chain, TSS/TES within 50bp | Absorb |
| 6 | Mono-exon overlapping gene, no intron crossing | Drop |
| 7 | Mono-exon spanning exon-intron-exon (intron retention) | Keep |
| 8 | Mono-exon intergenic | Drop |
| 1 | 5' ISM (contiguous subset at 5' end) | Keep |
| 2 | 3' ISM (1-2 missing exons from 5' end) | Absorb |
| 3 | 3' degradation (3+ missing from 5' end) | Drop vs ref, Keep vs sample |
| 4 | Internal fragment (both ends missing) | Drop vs ref, Keep vs sample |

After creating a new segment, reverse absorption applies the same rules to existing shorter segments. Matching uses pointer identity first, then fuzzy coordinate matching within `--fuzzy-tolerance` bp. Absorption can be disabled with `--no-absorb`.

## Output Files

### Build output (all subcommands that build a grove)

| File | Description |
|------|-------------|
| `{prefix}.ggx` | Serialized grove index |
| `{prefix}.ggx.summary` | Build summary (genes, transcripts, segments, exons, biotypes, per-chromosome) |

### `atroplex inspect` output

Outputs land directly under `{output-dir}/`, grouped by category
subfolders (no wrapper folder — matches the convention of the other
subcommands):

```
{output-dir}/
  overview/
    {basename}.overview.tsv          Global counts (genes, transcripts, segments, exons)
    {basename}.per_sample.tsv        Per-sample metrics (one row per sample)
    {basename}.per_source.tsv        Per-GFF-source metrics (HAVANA / ENSEMBL / TALON / ...)
    {basename}.biotype.tsv           Gene + transcript biotype breakdown (long-form)
  sharing/
    {basename}.exon_sharing.tsv      Exon sharing summary (metrics × samples)
    {basename}.segment_sharing.tsv   Segment sharing summary (metrics × samples)
    {basename}.conserved_exons.tsv   Per-exon detail for exons meeting the conservation threshold (`--conserved-fraction`, default = all samples)
  splicing_hubs/
    {basename}.splicing_hubs.tsv     Hub exons (>10 downstream branches) with per-sample PSI + entropy
    {basename}.branch_details.tsv    Per-(hub × target) branch fraction + expression
  splicing_events/                   (only with --events)
    {basename}.splicing_events.tsv   Classified events: cassette / alt-5′ / alt-3′ / IR / alt-terminal / mutex
```

#### Exon Sharing Summary (`.exon_sharing.tsv`)

Metrics as rows, samples as columns, plus a `total` column:

| Metric | Description |
|--------|-------------|
| `total` | Total exons in this sample |
| `exclusive` | Exons only in this sample |
| `shared` | Exons in 2+ but not all samples |
| `conserved` | Exons meeting the conservation threshold (`--conserved-fraction`, default = all samples) |
| `constitutive` | Exons in all transcripts of their gene |
| `alternative` | Exons in only some transcripts of their gene |

#### Segment Sharing Summary (`.segment_sharing.tsv`)

Same format: `total`, `exclusive`, `shared`, `conserved`.

#### Conserved Exons Detail (`.conserved_exons.tsv`)

One row per exon present in **all samples**. Per-sample columns include transcript counts and expression values (expression columns only for sample-type entries).

#### Splicing Hubs (`.splicing_hubs.tsv`)

Exons with more than 10 unique downstream exon targets, indicating complex alternative splicing decision points. Per-sample columns include branch counts, shared/unique classification, transcript counts, Shannon entropy, PSI, and expression.

#### Branch Details (`.branch_details.tsv`)

One row per (hub exon, downstream target) pair with per-sample branch usage fractions and expression.

### `atroplex query` output

| File | Description |
|------|-------------|
| `{basename}.query.tsv` | Per-transcript classification with per-sample presence/expression |
| `{basename}.query.summary.txt` | Classification summary (counts per category) |
| `{basename}.{group_a}_vs_{group_b}.dtu.tsv` | DTU results (when `--contrast` is provided) |

Only matched transcripts are emitted in `query.tsv` — unmatched transcripts (intergenic, antisense) are counted in the summary but omitted from the TSV since they carry no per-sample or structural information.

#### Classification columns (`query.tsv`)

| Column | Description |
|--------|-------------|
| `transcript_id` | Query transcript identifier (from input GTF `transcript_id` attribute) |
| `gene_id` | Gene ID of the best-matching segment in the index (may be an annotated gene like ENSG… or a tool-generated ID like MSTRG.* for novel loci) |
| `gene_name` | Gene name of the best-matching segment (`.` when no gene name is available) |
| `structural_category` | SQANTI-like classification: FSM, ISM, NIC, NNC, genic_intron, genic_genomic |
| `subcategory` | Refinement: ISM → 5prime_fragment / 3prime_fragment / internal_fragment; NIC → combination / intron_retention / exon_skipping / alternative_3end / alternative_5end; NNC → novel_donor / novel_acceptor / novel_both / novel_exon |
| `junction_match_score` | Fraction of query splice junctions matching the best reference segment (0.0–1.0) |
| `matching_junctions` | Number of query junctions that match the reference |
| `query_junctions` | Total splice junctions in the query transcript |
| `ref_junctions` | Total splice junctions in the best-matching reference segment |
| `known_donors` | Query donor sites (exon 3' ends) matching any known donor in the index |
| `known_acceptors` | Query acceptor sites (exon 5' starts) matching any known acceptor in the index |
| `novel_donors` | Query donor sites not found in the index |
| `novel_acceptors` | Query acceptor sites not found in the index |
| `n_samples` | Number of samples in the index that contain the matched segment |
| `{sample_id}.present` | Per-sample presence (1/0) — one column per entry in the index |
| `{sample_id}.{expr_type}` | Per-sample expression at the matched segment (only when `.qtx` sidecar is available; `.` = no data) |

#### DTU columns (`*.dtu.tsv`)

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier |
| `gene_name` | Gene name |
| `transcript_id` | Transcript being tested |
| `prop_{group_a}` | Proportion of gene expression from this transcript in group A |
| `prop_{group_b}` | Proportion of gene expression from this transcript in group B |
| `delta_proportion` | Difference in proportions (group A − group B) |
| `p_value` | Chi-squared test p-value (Wilson-Hilferty approximation) |
| `fdr` | Benjamini-Hochberg adjusted p-value |
| `significant` | `yes` / `no` based on `--fdr` threshold |

### `atroplex discover` output

| File | Description |
|------|-------------|
| `{input}.atroplex.tsv` | Per-cluster match results (SQANTI-like) |
| `{input}.atroplex.summary.txt` | Match statistics |

### `atroplex export` output

One GTF file per exported sample in the output directory:

| File | Description |
|------|-------------|
| `{sample_id}.gtf` | Reconstructed GTF with gene/transcript/exon lines, expression as attributes |

## Key Concepts

### Pan-Transcriptome Index

Atroplex builds a combined index from multiple annotation sources. Each input file (reference annotation or sample assembly) is registered with metadata, and every feature (exon, segment) tracks which samples and sources it came from. Features at identical coordinates are deduplicated:

- **Exons**: Deduplicated by genomic coordinates (chromosome + strand + start + end)
- **Segments**: Deduplicated by exon structure (ordered list of exon coordinates within a transcript)

### Two-Level Feature System

1. **Segments** are spatially indexed transcript paths used for coarse queries. A segment represents one or more transcripts that share the same exon structure. Segments participate in interval tree queries.

2. **Exons** are graph-only external keys used for fine-grained verification. They are linked into chains via edges and are not spatially indexed.

The graph connects segments to their first exon (SEGMENT_TO_EXON edge), and exons to subsequent exons (EXON_TO_EXON edges). Each edge carries a numeric segment index for ID-based traversal at branching exons.

### Structural Categories (SQANTI-like)

The `query` and `discover` subcommands classify transcripts into structural categories:

| Category | Description |
|----------|-------------|
| FSM | Full Splice Match — all junctions match a reference transcript |
| ISM | Incomplete Splice Match — subset of reference junctions |
| NIC | Novel In Catalog — novel combination of known splice sites |
| NNC | Novel Not in Catalog — at least one novel splice site |
| GENIC_INTRON | Mono-exon entirely within an intron |
| GENIC_GENOMIC | Overlaps both intron and exon regions |
| ANTISENSE | Overlaps gene on opposite strand |
| INTERGENIC | No gene overlap |

### Sample Types

Entries are classified as `"annotation"` (reference catalogs like GENCODE) or `"sample"` (experimental assemblies). This distinction affects:

- **Conserved/sharing statistics**: Only sample-type entries count toward "all samples" thresholds
- **Expression output**: Expression columns only appear for sample-type entries in TSV output files
- **Absorption rules**: Rules 3 and 4 treat annotation parents differently from sample parents

## Scripts

### inject_expression.py

Injects expression counts from TALON quantification TSV into GTF files:

```bash
python scripts/inject_expression.py \
  input.gtf counts.tsv output.expr.gtf \
  --header id=SAMPLE_001 \
  --header species="Homo sapiens"
```

Matches on `talon_transcript` attribute in GTF to `transcript_ID` in TSV. Adds `counts "N"` attribute to transcript lines.

### visualize_stats.R

Generates multi-panel visualization from analysis output:

```bash
Rscript scripts/visualize_stats.R <stats_dir> [output_prefix]
```

Requires: ggplot2, tidyr, dplyr, scales, patchwork. Produces PDF and PNG.

## License

GPLv3. See [LICENSE](LICENSE) for details.
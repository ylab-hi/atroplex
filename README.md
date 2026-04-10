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
- **genogrove** (v0.20.2): Genomic interval data structures and graph structure (fetched automatically)
- **htslib**: Reading BAM/SAM files (system dependency via pkg-config)
- **zlib**: Compression support

## Quick Start

```bash
# Build index from a single annotation file
atroplex build -b annotation.gtf

# Build pan-transcriptome from multiple sources via manifest
atroplex build -m manifest.tsv

# Build with replicate merging and expression filter
atroplex build -m manifest.tsv --min-replicates 2 --min-expression 3

# Run full per-sample analysis (sharing, splicing hubs, diversity)
atroplex analyze -m manifest.tsv -o results/

# Classify transcripts against the index
atroplex query -i transcripts.gtf -m manifest.tsv -o results/

# Classify with differential transcript usage
atroplex query -i transcripts.gtf -m manifest.tsv --contrast treated:control

# Discover novel transcripts from long-read data
atroplex discover -i reads.bam -m manifest.tsv -o results/

# Export per-sample GTF files from a built index
atroplex export -g index.ggx -o export/
```

## Subcommands

### `atroplex build` — Build pan-transcriptome index

Builds a genogrove index from one or more annotation files or BAM files. Each input file is treated as a separate entry in the pan-transcriptome, with features deduplicated across files and sample/source provenance tracked on every exon and segment.

```bash
# From manifest (full metadata per sample)
atroplex build -m ENCODE/manifest.tsv -o results/

# From annotation files directly (metadata parsed from GFF headers)
atroplex build -b gencode.gtf -b sample1.gtf -b sample2.gtf

# With replicate merging (require feature in >= 2 replicates)
atroplex build -m manifest.tsv --min-replicates 2

# With expression filtering (skip transcripts below threshold)
atroplex build -m manifest.tsv --min-expression 3

# Disable ISM absorption
atroplex build -m manifest.tsv --no-absorb
```

Build always produces a serialized index (`.ggx`) and a build summary (`.ggx.summary`).

### `atroplex analyze` — Full pan-transcriptome analysis

Performs detailed per-sample analysis including Jaccard isoform diversity, exon/segment sharing statistics, and splicing hub detection.

```bash
# Full analysis from manifest
atroplex analyze -m ENCODE/manifest.tsv -o results/

# Full analysis from annotation files
atroplex analyze -b gencode.gtf -b sample1.gtf -o results/
```

### `atroplex query` — Classify transcripts against the index

Classifies input transcripts (GTF/GFF) against the pan-transcriptome index using SQANTI-like structural categories (FSM, ISM, NIC, NNC, etc.). Optionally performs differential transcript usage (DTU) analysis between sample groups.

```bash
# Classify transcripts
atroplex query -i transcripts.gtf -m manifest.tsv -o results/

# With differential transcript usage between conditions
atroplex query -i transcripts.gtf -m manifest.tsv --contrast treated:control --fdr 0.05

# Group by a different manifest field
atroplex query -i transcripts.gtf -m manifest.tsv --contrast HeLa:K562 --group-by biosample
```

### `atroplex discover` — Discover novel transcripts

Clusters aligned long reads by splice junction signature and matches them against the index to classify transcripts as known, compatible, or novel.

```bash
atroplex discover -i reads.bam -m manifest.tsv -o results/
```

### `atroplex export` — Reconstruct per-sample GTFs from the index

Walks a built index and writes one GTF file per sample with gene, transcript, and exon lines. Expression values (when available) are emitted as GTF attributes. Supports filters to restrict output by sample, gene, region, biotype, source, or sample frequency.

```bash
# Export all samples from a pre-built index
atroplex export -g index.ggx -o export/

# Export a specific sample, protein-coding genes only
atroplex export -g index.ggx --sample HL60_M1_rep1 --biotype protein_coding

# Export conserved features in a genomic region
atroplex export -g index.ggx --region chr22:10000000-15000000 --conserved-only

# Export features present in at least 3 samples, from HAVANA
atroplex export -g index.ggx --min-samples 3 --source HAVANA
```

Export-specific options: `--sample`, `--gene`, `--region chr:start-end`, `--min-samples`, `--conserved-only`, `--biotype`, `--source` (all filters are AND'd).

### Common Options

| Option | Description |
|--------|-------------|
| `-o, --output-dir` | Output directory (default: input file directory) |
| `-p, --prefix` | Output file prefix (default: derived from manifest or first input) |
| `-t, --threads` | Number of threads (default: 1) |
| `--progress` | Show progress output |
| `-m, --manifest` | Sample manifest file (TSV) |
| `-b, --build-from` | Build from GFF/GTF file(s) |
| `-g, --genogrove` | Load pre-built genogrove index (.ggx) |
| `-k, --order` | Genogrove tree order (default: 3) |
| `--min-expression` | Minimum expression to include a transcript (default: -1, disabled) |
| `--no-absorb` | Disable ISM segment absorption into longer parent segments |
| `--fuzzy-tolerance` | Max bp difference for fuzzy exon boundary matching (default: 5) |
| `--min-replicates` | Merge biological replicates; require features in >= N replicates (default: 0, no merge) |

## Input Files

### Sample Manifest (TSV)

A tab-separated file specifying input files with metadata. This is the recommended way to build a pan-transcriptome with full sample tracking.

```tsv
file	id	type	assay	biosample	condition	species	platform	pipeline	description
gencode.v49.gtf	GENCODE_v49	annotation	.	.	.	Homo sapiens	.	GENCODE	Reference annotation
sample1.gtf	HL60_M1_rep1	sample	RNA-seq	HL-60	M1 macrophage	Homo sapiens	PacBio Sequel II	TALON	HL-60 M1 replicate 1
sample2.gtf	HL60_M1_rep2	sample	RNA-seq	HL-60	M1 macrophage	Homo sapiens	PacBio Sequel II	TALON	HL-60 M1 replicate 2
```

- Tab-separated, `"."` for empty values (VCF convention)
- `file` column is required; all others optional
- `type`: `"annotation"` for reference annotations, `"sample"` for experimental data (default: `"sample"`)
- Column names are case-insensitive
- Relative paths resolved from manifest directory
- `id` auto-generated from filename if not provided
- Optional `group` column for replicate grouping; if absent, groups are auto-inferred by stripping `_repNN` suffix from sample IDs

The `type` field matters: only entries with `type = "sample"` count toward "conserved" thresholds (features present in ALL samples). Reference annotations participate in structural analysis but don't inflate sample counts.

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

#### Expression Auto-Detection

Expression values are automatically parsed from transcript-level GFF entries. Priority order: `counts` > `TPM` > `FPKM` > `RPKM` > `cov`. Values are stored per-sample on segment features. The expression type (counts, TPM, etc.) is recorded on the sample and appears in output column headers.

Use `scripts/inject_expression.py` to add expression counts from TALON quantification TSV files into GTF files before building.

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

## Biological Replicate Merging

When `--min-replicates N` is provided, biological replicates within the same experiment group are merged after grove construction:

- Groups are determined by the `group` column in the manifest, or auto-inferred by stripping `_repNN` suffix from sample IDs
- Features must be present in >= N replicates to survive (threshold capped at group size)
- Expression values are averaged across replicates
- Original replicate entries are excluded from downstream statistics

## Output Files

### Build output (all subcommands that build a grove)

| File | Description |
|------|-------------|
| `{prefix}.ggx` | Serialized grove index |
| `{prefix}.ggx.summary` | Build summary (genes, transcripts, segments, exons, biotypes, per-chromosome) |

### `atroplex analyze` output

All output goes to `{output-dir}/analysis/` with organized subfolders:

```
analysis/
  {basename}.analysis.txt            Full text report
  {basename}.sample_stats.csv        Per-sample metrics (samples as columns, metrics as rows)
  {basename}.source_stats.csv        Per-source metrics (sources as columns, metrics as rows)
  sharing/
    {basename}.exon_sharing.tsv      Exon sharing summary
    {basename}.segment_sharing.tsv   Segment sharing summary
    {basename}.conserved_exons.tsv   Conserved exon detail
  splicing_hubs/
    {basename}.splicing_hubs.tsv     Splicing hub exons
    {basename}.branch_details.tsv    Per-branch detail
```

#### Exon Sharing Summary (`.exon_sharing.tsv`)

Metrics as rows, samples as columns, plus a `total` column:

| Metric | Description |
|--------|-------------|
| `total` | Total exons in this sample |
| `exclusive` | Exons only in this sample |
| `shared` | Exons in 2+ but not all samples |
| `conserved` | Exons in ALL samples |
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
| `{basename}.query.summary.txt` | Classification summary |
| `{basename}.{group_a}_vs_{group_b}.dtu.tsv` | DTU results (when `--contrast` is provided) |

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
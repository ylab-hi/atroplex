# Atroplex

A pan-transcriptome indexing and analysis toolkit for long-read sequencing data. Atroplex builds a spatial index with graph structure from multiple annotation sources, tracks features across samples, and provides detailed analysis of exon/segment sharing, splicing complexity, and isoform diversity.

## Installation

### Docker (recommended)

```bash
docker pull <dockerhub-username>/atroplex:latest
docker run --rm -v $(pwd):/data atroplex build -b /data/annotation.gtf --stats -o /data/output/
```

Or build the image locally:

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

- **cxxopts** (v3.2.1): Command-line argument parsing (fetched automatically)
- **genogrove** (v0.14.0): Genomic interval data structures and graph structure (fetched automatically)
- **htslib**: Reading BAM/SAM/FASTQ files (system dependency via pkg-config)
- **zlib**: Compression support

## Quick Start

```bash
# Build index from a single annotation file with summary statistics
atroplex build -b annotation.gtf --stats

# Build pan-transcriptome from multiple sources via manifest
atroplex build -m manifest.tsv --stats

# Run full per-sample analysis (sharing, splicing hubs, diversity)
atroplex analyze -m manifest.tsv -o results/

# Discover novel transcripts from long-read data
atroplex discover -i reads.bam -m manifest.tsv
```

## Subcommands

### `atroplex build` - Build pan-transcriptome index

Builds a genogrove index from one or more annotation files. Each input file is treated as a separate entry in the pan-transcriptome, with features deduplicated across files and sample/source provenance tracked on every exon and segment.

```bash
# From manifest (full metadata per sample)
atroplex build -m ENCODE/manifest.tsv --stats -o results/

# From annotation files directly (metadata parsed from GFF headers)
atroplex build -b gencode.gtf -b sample1.gtf -b sample2.gtf --stats

# Both combined
atroplex build -m manifest.tsv -b extra_annotation.gtf --stats
```

Use `--stats` to write a compact index summary to the output directory.

### `atroplex analyze` - Full pan-transcriptome analysis

Performs detailed per-sample analysis including Jaccard isoform diversity, exon/segment sharing statistics, and splicing hub detection. This is more expensive than `--stats` and produces multiple output files organized in subfolders.

```bash
# Full analysis from manifest
atroplex analyze -m ENCODE/manifest.tsv -o results/

# Full analysis from annotation files
atroplex analyze -b gencode.gtf -b sample1.gtf -o results/
```

### `atroplex discover` - Discover novel transcripts

Clusters aligned long reads by splice junction signature and matches them against the index to classify transcripts as known, compatible, or novel.

```bash
atroplex discover -i reads.bam -m manifest.tsv -o results/
```

### `atroplex dtu` - Differential transcript usage

*(Under development)*

### Common Options

| Option | Description |
|--------|-------------|
| `-o, --output-dir` | Output directory (default: input file directory) |
| `-t, --threads` | Number of threads (default: 1) |
| `-s, --stats` | Write compact index summary |
| `--progress` | Show progress output |
| `-m, --manifest` | Sample manifest file (TSV) |
| `-b, --build-from` | Build from GFF/GTF file(s) |
| `-g, --genogrove` | Load pre-built genogrove index (.gg) |
| `-k, --order` | Genogrove tree order (default: 3) |

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

Expression values are automatically parsed from transcript-level GFF entries. Priority order: `counts` > `TPM` > `FPKM` > `RPKM` > `cov`. Detected values are:

1. Stored on segments (transcript-level)
2. Propagated to exons (accumulated across transcripts sharing the exon)
3. The expression type (counts, TPM, etc.) is recorded on the sample and appears in output column headers

Use `scripts/inject_expression.py` to add expression counts from TALON quantification TSV files into GTF files before building.

#### Chromosome Name Normalization

Atroplex normalizes chromosome names to UCSC/GENCODE style automatically:

| Input | Normalized |
|-------|------------|
| `1`, `2`, ... `22` | `chr1`, `chr2`, ... `chr22` |
| `X`, `Y` | `chrX`, `chrY` |
| `MT` | `chrM` |
| `chr1` (already prefixed) | `chr1` (unchanged) |

This allows mixing annotations from different sources (e.g., GENCODE + Ensembl).

## Output Files

### `--stats` (compact summary)

Written directly to the output directory:

```
{basename}.index_stats.txt
```

A compact text summary with:
- Overview (chromosomes, genes, transcripts, segments, exons, edges)
- Input entries listed by name and type
- Genes by biotype
- Transcripts per gene / exons per segment distributions
- Per-chromosome breakdown

### `atroplex analyze` (full analysis)

All output goes to `{output-dir}/analysis/` with organized subfolders:

```
analysis/
  {basename}.analysis.txt            Full text report
  {basename}.sample_stats.csv        Per-sample metrics (CSV)
  {basename}.source_stats.csv        Per-source metrics (CSV)
  sharing/
    {basename}.exon_sharing.tsv      Exon sharing summary
    {basename}.segment_sharing.tsv   Segment sharing summary
    {basename}.conserved_exons.tsv   Conserved exon detail
  splicing_hubs/
    {basename}.splicing_hubs.tsv     Splicing hub exons
    {basename}.branch_details.tsv    Per-branch detail
```

#### Full Text Report (`.analysis.txt`)

Human-readable report with all sections: overview, biotype breakdown, distributions, exon sharing (constitutive/alternative), segment sharing (conserved/shared/sample-specific), graph structure, isoform diversity (Jaccard), per-sample comparison table, per-source comparison table, and top splicing hubs.

#### Per-Sample CSV (`.sample_stats.csv`)

Samples as columns, metrics as rows. Includes metadata rows (source_file, assay, biosample, etc.) followed by statistics: segments, exclusive segments, exons, exclusive exons, genes, transcripts, isoform diversity, deduplication ratio, mean expression, expressed segments.

#### Per-Source CSV (`.source_stats.csv`)

GFF sources (column 2: HAVANA, ENSEMBL, TALON, etc.) as columns, metrics as rows: segments, exclusive segments, exons, exclusive exons, genes.

#### Exon Sharing Summary (`.exon_sharing.tsv`)

Metrics as rows, samples as columns, plus a `total` column with global values:

| Metric | Description |
|--------|-------------|
| `total` | Total exons in this sample |
| `exclusive` | Exons only in this sample |
| `shared` | Exons in 2+ but not all samples |
| `conserved` | Exons in ALL samples |
| `constitutive` | Exons in all transcripts of their gene |
| `alternative` | Exons in only some transcripts of their gene |

#### Segment Sharing Summary (`.segment_sharing.tsv`)

Same format as exon sharing:

| Metric | Description |
|--------|-------------|
| `total` | Total segments in this sample |
| `exclusive` | Segments only in this sample |
| `shared` | Segments in 2+ but not all samples |
| `conserved` | Segments in ALL samples |

#### Conserved Exons Detail (`.conserved_exons.tsv`)

One row per exon that is present in **all samples**. Provides per-sample transcript counts and expression values for cross-sample comparison of core exons.

| Column | Description |
|--------|-------------|
| `exon_id` | Exon identifier |
| `gene_name` | Gene symbol |
| `gene_id` | Gene identifier |
| `chromosome` | Chromosome |
| `coordinate` | Genomic coordinate (chr:strand:start-end) |
| `n_transcripts` | Total transcripts using this exon |
| `constitutive` | `yes` if present in all transcripts of the gene, `no` if alternative |
| `{sample}.transcripts` | Number of transcripts using this exon in the sample |
| `{sample}.{expr_type}` | Expression value (e.g., `.counts`, `.TPM`); only for sample-type entries |

- Expression columns only appear for entries with `type = "sample"` (not annotations)
- Missing expression shown as `.` (present but not quantified)
- Sorted by chromosome then coordinate

#### Splicing Hubs (`.splicing_hubs.tsv`)

Exons with more than 10 unique downstream exon targets, indicating complex alternative splicing decision points. One row per hub exon.

| Column | Description |
|--------|-------------|
| `gene_name` | Gene symbol |
| `gene_id` | Gene identifier |
| `exon_id` | Hub exon identifier |
| `coordinate` | Genomic coordinate |
| `exon_number` | 1-based position in a representative transcript chain |
| `total_exons` | Total exons in the representative segment |
| `total_branches` | Total unique downstream targets across all samples |
| `total_transcripts` | Total transcripts using this hub exon |

Per-sample columns (one set per entry):

| Column | Description |
|--------|-------------|
| `{sample}.branches` | Unique downstream targets in this sample |
| `{sample}.shared` | Branches also present in at least one other sample |
| `{sample}.unique` | Branches only in this sample |
| `{sample}.transcripts` | Transcripts using this hub exon in this sample |
| `{sample}.entropy` | Shannon entropy of branch usage: H = -Sigma(f_i * log2(f_i)) |
| `{sample}.psi` | Traditional PSI (Percent Spliced In): hub transcripts / gene transcripts |
| `{sample}.{expr_type}` | Expression at hub exon (sample-type entries only) |

- **Entropy** measures splicing complexity at the junction. Higher entropy means more evenly distributed branch usage; 0 means all transcripts take the same path.
- **PSI** measures exon inclusion rate. PSI = 1.0 means all transcripts of the gene include this exon.
- Missing values shown as `.` (hub exon not present in that sample).

#### Branch Details (`.branch_details.tsv`)

One row per (hub exon, downstream target) pair. Shows how transcript flow is distributed across individual branch targets.

| Column | Description |
|--------|-------------|
| `hub_gene_name` | Gene symbol |
| `hub_gene_id` | Gene identifier |
| `hub_exon_id` | Hub exon identifier |
| `hub_coordinate` | Hub genomic coordinate |
| `target_exon_id` | Downstream target exon identifier |
| `target_coordinate` | Target genomic coordinate |
| `{sample}.fraction` | Branch usage fraction: target transcripts / hub transcripts |
| `{sample}.{expr_type}` | Expression at target exon (sample-type entries only) |

- Fraction values sum to ~1.0 across all targets for a given hub and sample
- `.` indicates the target is not present in that sample

### `atroplex discover` output

```
{input}.atroplex.tsv            Per-cluster match results
{input}.atroplex.summary.txt    Match statistics
```

## Key Concepts

### Pan-Transcriptome Index

Atroplex builds a combined index from multiple annotation sources. Each input file (reference annotation or sample assembly) is registered with metadata, and every feature (exon, segment) tracks which samples and sources it came from. Features at identical coordinates are deduplicated:

- **Exons**: Deduplicated by genomic coordinates (chromosome + strand + start + end)
- **Segments**: Deduplicated by exon structure (ordered list of exon coordinates within a transcript)

### Two-Level Feature System

1. **Segments** are spatially indexed transcript paths used for coarse queries. A segment represents one or more transcripts that share the same exon structure. Segments participate in interval tree queries.

2. **Exons** are graph-only external keys used for fine-grained verification. They are linked into chains via edges and are not spatially indexed.

The graph connects segments to their first exon (SEGMENT_TO_EXON edge), and exons to subsequent exons (EXON_TO_EXON edges).

### Expression Handling

Expression values can come from two sources:

1. **Transcript-level** (current): Auto-detected from GFF transcript entries (counts, TPM, FPKM, RPKM, cov). Stored on segments and propagated to exons by accumulation across transcripts.
2. **Exon-level** (future): Direct per-exon quantification (e.g., from short-read counting). Stored directly on exon features.

When transcript-level expression is propagated to exons, it is **accumulated** — if an exon appears in 3 transcripts with counts 100, 50, and 200, the exon's expression for that sample will be 350. For constitutive exons (present in all transcripts of a gene), this equals total gene expression.

### Sample Types

Entries are classified as `"annotation"` (reference catalogs like GENCODE) or `"sample"` (experimental assemblies). This distinction affects:

- **Conserved/sharing statistics**: Only sample-type entries count toward "all samples" thresholds
- **Expression output**: Expression columns only appear for sample-type entries in TSV output files
- **Column headers**: Expression columns include the detected type (e.g., `SAMPLE1.counts`, `SAMPLE2.TPM`)

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
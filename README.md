# atroplex

A tool for breakpoint-aware transcript discovery from long-read sequencing data.


## Quick Start

```bash
# Build from annotation file
atroplex --build-from annotation.gtf --order 10

# Build from multiple files (pan-transcriptome)
atroplex --build-from gencode.gtf --build-from ensembl.gtf --order 10
```


## Input File Format (GFF/GTF)

Atroplex accepts GFF3 and GTF annotation files. The following requirements must be met:

### Required Attributes

Each feature must have these attributes in column 9:

| Attribute | Description | Example |
|-----------|-------------|---------|
| `gene_id` | Unique gene identifier | `ENSG00000139618` |
| `transcript_id` | Unique transcript identifier | `ENST00000380152` |

### Required Feature Types

| Feature Type | Description |
|--------------|-------------|
| `exon` | Exonic regions (required for segment construction) |

### Feature Types (Currently Skipped)

The following feature types are present in GFF/GTF files but are **not yet processed**:

| Feature Type | Description | Status |
|--------------|-------------|--------|
| `CDS` | Coding sequence regions | TODO |
| `five_prime_UTR` / `5UTR` | 5' untranslated regions | TODO |
| `three_prime_UTR` / `3UTR` | 3' untranslated regions | TODO |
| `start_codon` | Start codon position | TODO |
| `stop_codon` | Stop codon position | TODO |

These will be used in future versions to annotate exons with overlapping CDS/UTR regions.

### Supported Annotation Sources

Atroplex has been tested with:

- **GENCODE** (recommended): `gencode.v49.annotation.gtf`
- **Ensembl**: `Homo_sapiens.GRCh38.110.gtf`

### Chromosome Name Normalization

Atroplex automatically normalizes chromosome names to UCSC/GENCODE style:

| Input | Normalized |
|-------|------------|
| `1`, `2`, ... `22` | `chr1`, `chr2`, ... `chr22` |
| `X`, `Y` | `chrX`, `chrY` |
| `MT` | `chrM` |
| `chr1` (already prefixed) | `chr1` (unchanged) |

This allows mixing annotations from different sources (e.g., GENCODE + Ensembl) in pan-transcriptome builds.

### Example GTF Entry

```
chr1    GENCODE exon    11869   12227   .   +   .   gene_id "ENSG00000290825"; transcript_id "ENST00000456328"; gene_name "DDX11L2";
chr1    GENCODE exon    12613   12721   .   +   .   gene_id "ENSG00000290825"; transcript_id "ENST00000456328"; gene_name "DDX11L2";
chr1    GENCODE CDS     12613   12721   .   +   0   gene_id "ENSG00000290825"; transcript_id "ENST00000456328"; gene_name "DDX11L2";
```

### Deduplication Behavior

When building from multiple files:

- **Exons**: Deduplicated by genomic coordinates (chr + strand + start + end)
- **Segments**: Deduplicated by exon structure (ordered list of exon coordinates)

If the same exon or segment appears in multiple files, it is stored once with sample tracking metadata.



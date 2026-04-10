# Atroplex: Differentiation Roadmap

Goal: Position atroplex as a **graph-based structural splicing analysis framework** — distinct
from isopedia (read-level genotyper/lookup tool). Lead narrative: "isopedia tells you what exists,
atroplex tells you what's happening and why."

---

## Phase 0: Genogrove Prerequisites ✅ DONE

### 0.1 Fix `split_node` bug ✅
Fixed in genogrove v0.19.0.

### 0.2 Grove serialization/deserialization (.ggx files) ✅
Implemented: AGRX magic + version, registries + zlib-compressed grove. Save via
`subcall::save_grove()`, load via `subcall::load_grove()`.

### 0.3 Grove iteration API ✅
Using option (c): builder caches passed downstream. `chromosome_exon_caches` available
via `exon_caches_` member on subcall base class.

---

## Phase 1: Discovery Pipeline ✅ DONE

### 1.1 Implement `index_splice_sites()` ✅
Implemented in `src/transcript_matcher.cpp`. Walks `chromosome_exon_caches`,
stores `(seqid, position)` composite keys with proper hash combining (fixed
cross-chromosome collision bug). Called from `transcript_matcher` constructor.

### 1.2 Verify discovery pipeline end-to-end ✅
Verified via discover test suite (`tests/discover/`):
- 11 unit tests: direct `read_cluster` construction for each SQANTI category
  (FSM, ISM prefix/suffix, NIC, NNC, GENIC_INTRON, GENIC_GENOMIC, INTERGENIC)
- 6 integration tests: full SAM → `read_clusterer` → `transcript_matcher` pipeline
- Fixtures: minimal chr22 reference GTF (2 genes, 3 transcripts), 23-read SAM
- Read clustering uses sort+sweep (no binning boundary artifacts)

---

## Phase 1b: Refined Absorption Rules ✅ DONE

Clean structural rules for ISM absorption during index creation. These prevent ISM
inflation (isopedia's weakness) while preserving biologically meaningful isoforms.

### Absorption rules (implemented in `segment_builder.cpp`)

| Rule | Pattern | Action | Status |
|------|---------|--------|--------|
| 0 | FSM — identical exon coordinates | Merge metadata | ✅ |
| 1 | 5' ISM — contiguous subset at 5' end | Keep (potential APA) | ✅ |
| 2 | 3' ISM — 1-2 missing exons from 5' | Absorb (RT dropout) | ✅ |
| 3 | 3' degradation — 3+ missing from 5' | Drop vs ref / Keep vs sample | ✅ |
| 4 | Internal fragment — both ends missing | Drop vs ref / Keep vs sample | ✅ |
| 5 | Terminal variant — same introns, TSS/TES <50bp | Absorb | ✅ |
| 6 | Mono-exon gene overlap — no intron crossing | Drop | ✅ |
| 7 | Mono-exon intron retention — spans exon-intron-exon | Keep | ✅ |
| 8 | Mono-exon intergenic — no gene overlap | Drop | ✅ |

- Fuzzy subsequence matching (≤5bp tolerance) for Rules 0-4
- Reverse absorption when parent arrives after ISM
- Ref/sample distinction: annotations processed before samples
- Tests: `tests/build/absorption_test.cpp` with per-rule fixtures
- Documented in `absorption_rules.txt`

### Deferred
- **NMD absorption veto** — requires `--genome` + CDS info (Phase 2c)
- **Splice site dinucleotide storage** — requires `--genome` FASTA (Phase 2b.2)

---

## Phase 2: Alternative Splicing Event Catalog ✅ DONE

### 2.1 Define splicing event types ✅
Implemented in `include/splicing_catalog.hpp`:
- CASSETTE_EXON, MUTUALLY_EXCLUSIVE, ALT_5PRIME, ALT_3PRIME
- INTRON_RETENTION, ALT_FIRST_EXON, ALT_LAST_EXON

### 2.2 Implement graph walker ✅
Implemented in `src/splicing_catalog.cpp`. Walks `chromosome_gene_segment_indices`,
compares exon chains across segments within each gene, classifies events by topology.

### 2.3 Integrate into `analyze` subcommand ✅
Splicing events collected during analysis, event counts in summary output.

### 2.4 Extend hub analysis with event classification ✅
`event_type` column added to `splicing_hubs.tsv`.

---

## Phase 2b: BAM Input for `atroplex build` ✅ PARTIALLY DONE

### 2b.1 Extract shared build helpers ✅
`segment_builder` extracted from `build_gff` — format-agnostic segment creation
and absorption logic shared by GFF and BAM builders.

### 2b.2 Splice site validator ❌ NOT STARTED
- Requires `--genome` FASTA (genogrove FASTA reader being added upstream)
- Canonical dinucleotide validation (GT-AG, GC-AG, AT-AC)
- Short-read junction support index

### 2b.3 BAM builder module ✅
Implemented in `include/build_bam.hpp`, `src/build_bam.cpp`:
- Clusters BAM reads via `read_clusterer`
- Derives exon coordinates from consensus junctions
- Resolves gene assignment via spatial query
- Same absorption rules as GTF path via `segment_builder`
- Expression = read count, source = "BAM"

### 2b.4 CLI integration ✅
`--build-from` auto-detects BAM via `filetype_detector`. GTF processed before BAM
in mixed builds. BAM header parsed for sample metadata (@RG/@PG).

### 2b.5 Short-read quality layer ❌ NOT STARTED
- Junction support filtering, novel splice site evidence
- Requires `--genome`/`--short-reads`/`--junctions` CLI options

---

## Phase 2c: NMD Prediction on Segments ❌ NOT STARTED

Blocked on `--genome` FASTA support (Phase 2b.2 prerequisite).

### 2c.1 Reference CDS integration
- Parse CDS features from reference annotation GTF
- Store per-gene CDS start/stop coordinates
- Map reference CDS onto segment exon chains to determine reading frame

### 2c.2 ORF prediction per segment
- Requires `--genome` (FASTA) from Phase 2b
- Extract spliced sequence, translate, find PTC

### 2c.3 NMD classification (50nt rule)
- PTC > 50nt upstream of last exon-exon junction → NMD candidate
- Add `nmd_candidate: bool` + `ptc_position: uint32_t` to `segment_feature`

### 2c.4 NMD catalog output
- Per-gene NMD summary, poison exon identification
- Cross-sample NMD tracking, connect to splicing event catalog

### 2c.5 Unproductive splicing analysis
- Fraction of NMD targets per gene per sample
- Differential unproductive splicing between conditions

---

## Phase 3: Sequence-Aware Discovery ❌ NOT STARTED

Neither isopedia nor SQANTI3 reports sequence content of discovered isoforms.
Long reads carry the full sequence — use it.

### 3.1 Extract consensus sequences per read cluster
- Extend `read_cluster::finalize()` with consensus sequence from member reads
- Simple majority-vote per position (or BAM MD/CS tags)

### 3.2 Report sequences in discovery output
- Add `consensus_sequence` and `sequence_variants` columns to `.atroplex.tsv`

### 3.3 Allele-specific isoform detection
- Extract SNVs from consensus vs reference
- Optional VCF input for haplotype phasing

### 3.4 RNA editing detection
- Compare consensus to reference at known editing sites
- Optional known editing sites BED file

---

## Phase 4: Differential Transcript Usage ✅ PARTIALLY DONE

### 4.1 Basic DTU from graph + expression ✅
Implemented in `src/subcall/query.cpp`:
- Chi-squared test on transcript proportions between condition groups
- Benjamini-Hochberg FDR correction
- `--contrast group_a:group_b` CLI, `--group-by` manifest field, `--fdr` threshold
- Output: `.dtu.tsv` with proportions, delta, p-value, FDR, significance
- Query test suite verifies per-sample presence/expression tracking

### 4.2 Isoform switching detection ❌ NOT STARTED
- Find dominant segment per gene per sample
- Flag genes where dominant segment changes between conditions

### 4.3 Exon-level differential inclusion ❌ NOT STARTED
- Per-exon PSI per sample, differential PSI between conditions
- Connects to splicing event catalog (Phase 2)

---

## Phase 5: Benchmarking ❌ NOT STARTED

### 5.1 LRGASP simulated data benchmark
- Download LRGASP simulated PacBio + ONT reads
- Build grove from GENCODE V38 (the LRGASP reference)
- Run `atroplex discover` on simulated reads
- Compare: atroplex classification accuracy vs StringTie+SQANTI3 pipeline
- Show: atroplex does classification in one step, others need two tools
- Compare against isopedia: show it can only re-identify, not classify

### 5.2 ENCODE cohort analysis
- Build grove from GENCODE + ENCODE TALON assemblies (use existing manifest)
- Run full analysis: splicing events, hubs, sharing
- Show biological findings: tissue-specific splicing patterns, hub complexity
  differences, condition-specific events
- Target: 50-100 samples minimum for credibility

### 5.3 Discovery + classification showcase
- Take ENCODE BAM files → discover against GENCODE-only grove
- Show: atroplex finds and classifies novel transcripts that isopedia
  would need a pre-assembled GTF to even look for
- Quantify: "X% of NIC events are intron retentions, Y% are exon skipping" —
  structural insight isopedia cannot provide

---

## Phase 6: Nice-to-have (strengthens paper, not essential)

### 6.1 Exon essentiality scoring
- Per exon: fraction of gene's transcript paths that include it
- Score 1.0 = constitutive (all paths), 0.1 = rarely included
- Already have the data — just need to compute and output

### 6.2 Exon co-occurrence matrix
- Per gene: pairwise exon inclusion/exclusion across all transcript paths
- Identify obligate pairs (always co-occur) and exclusive pairs (never co-occur)

### 6.3 Fix genogrove `split_node` bug ✅
Fixed in genogrove v0.19.0.

### 6.4 Fusion detection as graph edges
- Model fusions as cross-gene SEGMENT_TO_EXON edges
- Track fusion breakpoint at exon resolution
- Per-sample fusion tracking via existing sample_bitset

---

## Phase 7: LLM Query Interface [HIGH PRIORITY — UNIQUE CAPABILITY]

Natural language queries over the pan-transcriptome index. No other tool offers this.
"Which genes have tissue-specific splicing hubs?", "Show me all NMD candidates in brain",
"What isoforms of BRCA1 are shared across cancer samples?"

### 7.1 Structured query API (genogrove level)
- Programmatic query interface returning structured results from the grove
- High-level queries beyond raw `intersect()` and edge traversal
- Gene-centric queries: all segments/exons for a gene, per-sample presence
- Feature filtering: by category, sample count, expression, biotype
- **Where:** genogrove — generic interval/graph query layer, any grove user benefits

### 7.2 Domain-specific query layer (atroplex level)
- Biological semantics on top of structured queries
- Query types:
  - Gene queries: isoforms, splicing events, DTU, expression per sample
  - Splicing queries: hubs, branch targets, PSI, entropy by condition
  - Sharing queries: conserved/exclusive features, sample overlap
  - NMD queries: unproductive isoforms, poison exons (Phase 2c)
  - Comparison queries: "how does gene X differ between tissue A and B"
- Returns structured JSON for LLM consumption

### 7.3 MCP server (`atroplex serve`)
- New subcommand: starts an MCP (Model Context Protocol) server
- Exposes query functions as MCP tools with typed inputs/outputs
- LLM (Claude, etc.) calls tools via function calling, gets JSON back
- Workflow: natural language → LLM translates to tool calls → execute →
  LLM interprets structured results → natural language response
- Requires: loaded grove (.ggx), optionally genome FASTA

### 7.4 Interactive mode (`atroplex chat`)
- Optional: built-in chat interface using Claude API
- Loads grove + genome, starts conversation
- User asks questions in natural language, atroplex answers using the index
- Could also work as a Claude Code MCP tool (atroplex as MCP server,
  Claude Code as client)

---

## Phase 8: Streaming Analyze Architecture [CRITICAL — OOM at 1000 samples]

Current `index_stats::collect()` accumulates massive per-sample maps that blow up memory.
Killed by OOM on 128GB with ~1000 samples — fails in Phase 1 before any analysis begins.

### 8.1 Single-pass streaming traversal
Rewrite `collect()` to traverse the grove once, processing each segment + its exon chain
inline. No `segment_keys` vector, no `tx_to_samples` map, no `gene_sample_tx` map.

```
for each chromosome:
  for each segment in B+ tree leaves:
    increment global + per-sample counters (flat arrays)
    walk exon chain via edges
      increment exon counters
      detect hubs inline
    → done, nothing stored
```

### 8.2 Replace maps with flat counters
Per-sample stats currently use `unordered_map<uint32_t, unordered_set<string>>` (~40KB per
sample per map). Replace with flat `vector<counter_struct>` indexed by sample_id. Counter
struct holds only `size_t` fields (segments, exons, exclusive, shared, conserved, etc.).

### 8.3 Per-gene diversity inline
Entropy and effective isoform metrics can be computed per-gene during traversal: accumulate
exon usage counts for the current gene, compute entropy when the gene boundary is crossed
(segments are sorted by coordinate within each chromosome), add to running mean.

### 8.4 Hub detection inline
Already partially streaming (branch_details.tsv written inline). Extend to compute PSI
on the fly using the current gene's segment data rather than precomputed `gene_sample_tx`.

### 8.5 Split index_stats (issue #24)
Separate `build_summary` (lightweight, always) from `analysis_report` (heavyweight, analyze
only). Remove the `detailed` flag. Each has its own `collect()` with appropriate streaming.

### Memory targets
| Samples | Current peak | Target peak |
|---------|-------------|-------------|
| 100 | ~10-20 GB | ~1-2 GB |
| 1000 | OOM (>128 GB) | ~5-10 GB |
| 20000 | impossible | ~20-30 GB |

---

## Current Priority Order (for paper submission)

| Priority | Phase | Status | Blocker |
|----------|-------|--------|---------|
| ~~1~~ | ~~0.1-0.3 Genogrove prerequisites~~ | ✅ Done | — |
| ~~2~~ | ~~1.1 Splice site indexing~~ | ✅ Done | — |
| ~~3~~ | ~~1.2 Discovery e2e~~ | ✅ Done | — |
| ~~4~~ | ~~1b Absorption rules v2~~ | ✅ Done | — |
| ~~5~~ | ~~2.1-2.4 Splicing event catalog~~ | ✅ Done | — |
| ~~6~~ | ~~2b.1/2b.3/2b.4 BAM input~~ | ✅ Done | — |
| 7 | 8.1-8.5 Streaming analyze (OOM fix) | ❌ | None — critical for 1000+ samples |
| 8 | 2b.2/2b.5 Splice site validation | ❌ | `--genome` FASTA in genogrove |
| 9 | 2c NMD prediction | ❌ | `--genome` FASTA |
| 10 | 7.1-7.3 LLM/MCP query interface | ❌ | None — high priority |
| 11 | 5.1 LRGASP benchmark | ❌ | None — needed for paper |
| 12 | 4.2-4.3 Isoform switching + exon PSI | ❌ | None — DTU foundation exists |
| 13 | 3.1-3.2 Sequence extraction | ❌ | None |
| 14 | 5.2-5.3 ENCODE cohort analysis | ❌ | None — biological story |
| 15 | Everything else | ❌ | As time allows |
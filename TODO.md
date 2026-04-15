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

Legacy `index_stats::collect()` accumulated sample×feature maps that blew up memory and
OOM'd on 128GB with ~1000 samples before any analysis began. Phase 8 replaces that with a
single-pass streaming `analysis_report` whose memory scales with feature count, not sample
count. Golden rule: **flat per-sample counters or streamed detail rows, never a matrix
indexed by (sample_id, feature_id)**. Per-gene accumulators are cleared at every
chromosome boundary.

### 8.1 Basic overview stats (global + per-sample + per-chromosome) ✅
New `analysis_report` struct in `include/analysis_report.hpp` + `src/analysis_report.cpp`.
Single-pass B+ tree traversal: segment counting, inline exon chain walk with pointer
dedup via `unordered_set<const void*> visited_exons`, lightweight per-gene accumulators
finalized at chromosome boundaries. `sample_counters` is a flat vector indexed by
sample_id — no maps keyed on features.

Output:
- `analysis/overview/{basename}.overview.tsv` — global metric/value pairs
  (samples, genes, transcripts, segments, exons, edges, absorbed, distributions)
- `analysis/overview/{basename}.per_sample.tsv` — rows per sample with gene/segment/
  exon counts, exclusive/shared/conserved, expression

Wired into the `atroplex analyze` subcommand (`src/subcall/analyze.cpp`). The legacy
`index_stats::collect` call is commented out behind a `TODO(Phase 8)` block in that
file.

### 8.1b Builder path refactor + tombstone physical removal ✅
Not originally in the Phase 8 plan; added during the refactor because the legacy
build-summary output was entangled with `index_stats` and needed to move before the
streaming analyze work could cleanly delete the old struct.

- **Bumped genogrove to v0.21.0** (CMakeLists.txt:28) for the new `grove.remove_key()`
  API from upstream PR #305. Enables physical tombstone cleanup.
- **`builder::remove_tombstones()`** (src/builder.cpp) runs after `merge_replicates`,
  before summary collection. Walks `gene_indices` (authoritative source — unlike
  `segment_caches` which `try_reverse_absorption` eagerly erases), calls
  `grove.remove_key(seqid, key)` for each absorbed segment, then a single
  `grove.remove_edges_if()` pass strips orphan `EXON_TO_EXON` edges carrying the
  tombstone's `segment_index`. Prunes `gene_indices` of the absorbed entries and
  defensively scans `segment_caches` for any stragglers. Returns the removed count.
- **`include/build_summary.hpp` + `src/build_summary.cpp`** (new) — extracted the
  ~150-line inline stats collection from `builder::build_from_samples`. Owns the
  `build_summary` struct, the `.ggx.summary` text writer, and the new `build_counters`
  struct (see below). Replaces `index_stats` in the build path entirely;
  `subcall::build_stats` is now `std::optional<build_summary>`.
- **`build_counters`** struct threaded explicitly (no default args) through
  `build_gff::build / process_gene / process_transcript`, `build_bam::build`, and
  `segment_builder::create_segment`. Tracks five categories:
  - `input_transcripts` — raw count seen
  - `merged_transcripts` — folded into an existing segment at FSM/terminal/ISM
    forward-absorption (Rules 0, 2, 5, fuzzy-FSM)
  - `absorbed_segments` — tombstoned by reverse absorption and swept
  - `discarded_transcripts` — dropped entirely (`--min-expression`, mono-exon
    Rules 6/8, Rules 3/4 fragment drops vs reference)
  - `replicates_merged` — collapsed by `--min-replicates`
  All five surface in the `.ggx.summary` "Transcript processing:" section.
- **`tests/build/builder_pipeline_test.cpp`** (new, 4 tests) — the first tests to
  invoke `builder::build_from_samples` end-to-end (every earlier test bypassed the
  orchestrator and called `build_gff::build` directly):
  1. `BuilderFullPipeline_CountersPopulated` — asserts each counter field
  2. `RemoveTombstones_PhysicallyRemoved` — walks grove leaves, asserts no
     `absorbed == true` survives
  3. `RemoveTombstones_OrphanEdgesPruned` — swept-vs-unswept `edge_count()`
     comparison proves the `remove_edges_if` pass ran
  4. `BuildSummary_WrittenFileContainsCounters` — re-reads the written
     `.ggx.summary` file and greps for the processing section
  Generated fixture files on-the-fly in a per-test tmp dir; updated pre-existing
  tests (absorption/serialization/discover/query) to pass `fuzzy_tolerance` + a
  `build_counters&` explicitly now that those are no longer defaulted.

Landed as commit `d81d1f0`.

### 8.2 Per-sample sharing stats ✅
Pure TSV transposition of the flat `per_sample` counters already populated in 8.1.
Zero new traversal, zero new memory beyond the small sample-id + label vectors the
writers build (same size as `write_per_sample`).

- **`analysis_report::write_exon_sharing(path)`** — metric × (total + samples) TSV
  with rows `total`, `exclusive`, `shared`, `conserved`, `constitutive`, `alternative`.
  The last two rows are populated by 8.4's per-gene pass (below); the row structure
  was added in 8.2 and the counter sources in 8.4.
- **`analysis_report::write_segment_sharing(path)`** — same layout, rows `total`,
  `exclusive`, `shared`, `conserved`.

Output:
- `analysis/sharing/{basename}.exon_sharing.tsv`
- `analysis/sharing/{basename}.segment_sharing.tsv`

Wired via `subcall::analyze` after `collect()` returns.

### 8.3 Splicing hub detection inline ✅
Branching exons with `> MIN_HUB_BRANCHES` (10) unique downstream targets are
detected inline during the exon chain walk, buffered per-gene, and written to
disk at gene finalization. Rows are emitted from inside `finalize_gene` and the
gene's `pending_hubs` vector is destroyed at the chromosome boundary via
`active_genes.clear()` — nothing ever accumulates across genes or chromosomes.

**Detection:** on each exon's first global visit (guarded by the existing
`visited_exons` dedup set), a second `grove.get_neighbors_if()` call with an
unfiltered predicate returns every outgoing `EXON_TO_EXON` edge. Targets are
deduplicated into a local `unordered_set<key_ptr>`; if the count exceeds the
threshold, a `pending_hub { key_ptr exon, vector<key_ptr> targets, chain_pos,
chain_total }` is registered on the current gene's accumulator. `chain_pos` /
`chain_total` are tracked by counters on the segment's exon chain loop.

**Emission at finalize_gene:** for each pending hub, the three outer-scope
per-sample buffers (`hub_branches`, `hub_shared`, `hub_unique`, each sized
`num_samples`) are `std::fill`-reset and populated by iterating the target
bitsets once — each target contributes to every sample in its `sample_idx`,
with shared/unique classified by the target's own `sample_count()`. The row is
written immediately and the next hub reuses the same buffers.

**Memory discipline:**
- `pending_hubs` is scoped inside `gene_acc`; cleared at every chromosome
  boundary. Peak ≈ one chromosome's worth of hub pointer lists, typically
  ~250 KB.
- `hub_branches` / `hub_shared` / `hub_unique` are allocated **once** in the
  outer `collect()` scope and reused across every hub in every gene. Peak =
  `O(num_samples)` × 3 vectors, *never* `O(num_hubs × num_samples)`.

**Skipped per-sample columns** (require transcript→sample attribution the flat
data model doesn't support without a `sample × feature` map):
- `.transcripts` — per-sample hub transcript count
- `.entropy` — per-sample Shannon entropy of branch usage
- `.psi` — traditional PSI (hub_transcripts / gene_transcripts per sample)

These can be reintroduced later via a bounded per-gene segment→sample tally if
needed, but for now the R viz scripts can derive fractions by grouping
`branch_details.tsv`.

**Skipped global columns (deferred to 8.6):**
- `event_type` — cassette / alt_5' / alt_3' classification overlaps with the
  splicing event catalog work.

**Outputs:**
- `analysis/splicing_hubs/{basename}.splicing_hubs.tsv` — one row per hub.
  Global columns: `gene_name`, `gene_id`, `exon_id`, `coordinate`,
  `exon_number`, `total_exons`, `total_branches`, `total_transcripts`.
  Per-sample columns (non-replicate entries only): `.branches`, `.shared`,
  `.unique`, plus `.{expr_type}` for sample-type entries.
- `analysis/splicing_hubs/{basename}.branch_details.tsv` — one row per
  (hub × target) pair. Global columns: `hub_gene_name`, `hub_gene_id`,
  `hub_exon_id`, `hub_coordinate`, `target_exon_id`, `target_coordinate`.
  Per-sample columns: `.present` (1 or `.`) plus `.{expr_type}` for sample-type
  entries.

Wired via `analysis_report::begin_splicing_hub_streams(hubs_path, branches_path)`
called from `src/subcall/analyze.cpp` before `collect()`.

### 8.4 Per-gene diversity metrics + constitutive/alternative ✅
Runs inside the same single traversal as 8.1. Per-gene accumulator extended with:
- `segment_tx_counts: vector<size_t>` — per-segment transcript counts used for the
  effective-isoforms calculation
- `exon_seg_counts: unordered_map<key_ptr, size_t>` — segments-in-this-gene that
  contain each exon (incremented on every exon visit regardless of first-visit
  dedup, so accurate)
- `first_visit_exons: vector<key_ptr>` — exons whose global first visit fell in
  this gene; used at finalization for constitutive/alternative classification
  without double counting

At `finalize_gene` for multi-segment genes:
- **Effective isoforms** = `2^H(segment→tx-count distribution)`
- **Gene exon entropy** = `H(exon usage fraction distribution)`
- Both sums accumulated into `sample_counters.effective_isoforms_sum` /
  `gene_exon_entropy_sum` with a `multi_segment_genes` denominator; `write_per_sample`
  emits the mean as `mean_effective_isoforms` / `mean_gene_exon_entropy`.
- **Constitutive vs alternative**: for each exon in `first_visit_exons`, if
  `exon_seg_counts[e] == gene.segment_count` it's constitutive, else alternative.
  Per-sample `constitutive_exons` / `alternative_exons` counters are incremented
  for every sample in `exon.sample_idx`.

Writer updates:
- `write_per_sample` adds four new columns: `constitutive_exons`, `alternative_exons`,
  `mean_gene_exon_entropy`, `mean_effective_isoforms`.
- `write_exon_sharing` adds the `constitutive` and `alternative` rows introduced
  in 8.2.

Memory cost: per active gene, `exon_seg_counts` (~20 entries) + `segment_tx_counts`
(~10 entries) + `first_visit_exons` (~20 entries). Cleared at every chromosome
boundary. Scales with feature count, never with sample count.

### 8.5 Conserved exon detail streaming ✅
Per-exon rows streamed inline during the exon chain walk in `collect()`. No
accumulator — rows are emitted the moment an exon is first visited if its
`sample_count() == total_samples_for_conserved`.

- **`analysis_report::begin_conserved_exon_stream(path)`** — called from `analyze.cpp`
  *before* `collect()` so the file handle and sample-label header are ready when
  traversal begins. Enumerates non-replicate registry entries in registry order,
  writes header `exon_id\tgene_name\tgene_id\tchromosome\tcoordinate\tn_transcripts\t{label}.{expr_type}…`.
- The `unique_ptr<ofstream>` lives on the `analysis_report` instance and is closed
  automatically on destruction. `conserved_stream_sample_ids` + `conserved_stream_is_sample`
  are bounded by sample count.
- Inside the existing exon chain walk: if `exon_conserved && stream_open`, write
  one row with gene info, coordinate, total `transcript_ids.size()`, and per-sample
  expression columns (`.` for missing) — **sample-type entries only** (annotations
  get no expression columns).

Output: `analysis/sharing/{basename}.conserved_exons.tsv`

**Scoped-out (deferred):**
- Per-sample transcript count columns — require a `transcript_id → {sample_ids}`
  map, which is exactly the OOM-prone sample×feature structure the rewrite exists
  to avoid. Not reintroducing.
- `constitutive` flag column — 8.4 computes constitutive/alternative at *gene
  finalization* while 8.5 emits rows at *first exon visit*, so the flag isn't
  known at emission time. Adding it would require buffering per-gene conserved
  exon candidates until finalization. Small bounded cost; acceptable for a small
  follow-up if biologists want it, but not in this pass. The sample-level
  `constitutive_exons` / `alternative_exons` counts are already in exon_sharing.tsv
  and per_sample.tsv.

### 8.6 Splicing event catalog inline ❌
At gene finalization: compare segment exon chains (still in gene accumulator) to
detect cassette exons, alt 5'/3', intron retention, MXE. No separate pass — reuse
the per-gene data already built for 8.4.

Output: `analysis/splicing_events/{basename}.splicing_events.tsv`

Implementation notes:
- Requires per-gene segment exon chains at finalization. Currently `gene_acc` only
  has `exon_seg_counts` (usage counts, not the chains themselves). To detect
  cassettes and alt-splice events we need the ordered exon sequence per segment.
  Option A: store `vector<vector<key_ptr>>` in gene_acc (per-gene bound by
  `#segments × avg_exons_per_segment` — a few KB per gene). Option B: re-walk the
  chains at finalization by iterating the grove's graph again — doubles traversal
  cost. Prefer A; memory stays O(features-in-gene).
- Legacy `splicing_catalog::collect_from_grove` already classifies events against
  a pre-built grove — can be adapted or reused inline.

### 8.7 Remove legacy `index_stats::collect()` ❌
Once 8.3 + 8.6 land and `analysis_report` covers all output the old system produced,
delete `src/index_stats.cpp::collect()` and the `detailed` flag. The `index_stats`
struct itself can go too — `build_summary` (from 8.1b) already replaced its
build-path usage, and analyze no longer references it. Relates to issue
[#24](https://github.com/ylab-hi/atroplex/issues/24).

Cleanup checklist at that point:
- Delete `src/index_stats.cpp` entirely
- Delete `include/index_stats.hpp`
- Remove commented-out legacy block in `src/subcall/analyze.cpp`
- Remove `#include "index_stats.hpp"` from `src/subcall/analyze.cpp`
- Remove any lingering `index_stats::` references (grep first)
- Close issue #24

### Memory targets
| Samples | Legacy peak | Streaming target | Status |
|---------|------------|-----------------|--------|
| 100 | ~10-20 GB | ~1-2 GB | Unverified (needs benchmark run after 8.6) |
| 1000 | OOM (>128 GB) | ~5-10 GB | Unverified |
| 20000 | impossible | ~20-30 GB | Unverified |

Benchmarking deferred until 8.6 + 8.7 land so we measure the full streaming pipeline
at once, not intermediate states.

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
| ~~7a~~ | ~~8.1 + 8.1b + 8.2 + 8.3 + 8.4 + 8.5 Streaming analyze (core + hubs)~~ | ✅ Done | — |
| 7b | 8.6 splicing event catalog + 8.7 legacy cleanup | ❌ | None — completes streaming analyze |
| 8 | 2b.2/2b.5 Splice site validation | ❌ | `--genome` FASTA in genogrove |
| 9 | 2c NMD prediction | ❌ | `--genome` FASTA |
| 10 | 7.1-7.3 LLM/MCP query interface | ❌ | None — high priority |
| 11 | 5.1 LRGASP benchmark | ❌ | None — needed for paper |
| 12 | 4.2-4.3 Isoform switching + exon PSI | ❌ | None — DTU foundation exists |
| 13 | 3.1-3.2 Sequence extraction | ❌ | None |
| 14 | 5.2-5.3 ENCODE cohort analysis | ❌ | None — biological story |
| 15 | Everything else | ❌ | As time allows |
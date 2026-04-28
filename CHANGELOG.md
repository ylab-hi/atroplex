# Changelog

All notable changes to atroplex will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Removed
- **Build-time replicate merging** (#47): `--min-replicates` and `--min-replicate-fraction` removed from the build path. The `merge_replicates()` function destructively cleared individual sample bits and replaced them with group bits, losing per-sample information needed for within-group analysis (e.g., boxplots). It also didn't save memory (bitset size unchanged) or remove dead segments. Replaced by inspect-time `--min-samples` (#48).

### Added
- **`sharing/conserved_segments.tsv`**: per-segment counterpart to the existing `conserved_exons.tsv`. Streamed inline during `analysis_report::collect()` for any segment passing `is_conserved(min_required_for_conserved_)`. Columns: `segment_id`, `gene_id`, `gene_name`, `gene_biotype`, `chromosome`, `coordinate`, `exon_count`, `n_transcripts`, `n_samples`, `sources` (semicolon-joined from the source bitfield). Per-sample expression columns appear when a `.qtx` sidecar is loaded, mirroring the conserved-exons gating. Lets the conserved transcript-structure core of the pan-transcriptome be enumerated and downstream-annotated (biotype, function, expression) rather than just summarized as a count.
- **`--conserved-fraction` for inspect**: configurable conservation threshold. Fraction in `(0, 1]` of sample-typed registry entries (annotations excluded) that a feature must appear in to be classified as conserved. Default `1.0` preserves the strict "in every sample" definition (back-compatible). Relaxing (e.g. `--conserved-fraction 0.95`) yields a dropout-tolerant "conserved core" — at 21K samples, strict equality drops to ~6 segments because a single missed sample disqualifies an otherwise universal structure. The threshold applies to both `is_conserved()` checks (segments and exons), so per-sample `conserved_segments` / `conserved_exons` counts and the conserved-detail TSVs all use the same definition. `is_conserved(N)` semantics changed from `count == N` to `count >= N`; at default fraction=1.0 the two are equivalent.
- **Query: skip unmatched rows, document output format** ([#54](https://github.com/ylab-hi/atroplex/pull/54)): `query.tsv` now omits intergenic/antisense rows (zero per-sample data); counts remain in the summary. README documents all `query.tsv` and `dtu.tsv` columns. New test verifies query against unannotated loci returns FSM with tool-generated gene_id.
- **Performance: candidate-based splice site classification, ostringstream removal** ([#71](https://github.com/ylab-hi/atroplex/pull/71), closes [#15](https://github.com/ylab-hi/atroplex/issues/15), partial [#12](https://github.com/ylab-hi/atroplex/issues/12)): splice sites for NIC/NNC classification are now derived from spatial candidates already found during `match()`, eliminating the O(segments × exons) grove traversal at matcher construction. More correct than the global index (only considers overlapping loci). Signature methods use `string::reserve` + `+=` instead of `ostringstream`.
- **C++ modernization: string_view, [[nodiscard]], named constants** ([#66](https://github.com/ylab-hi/atroplex/pull/66), closes [#60](https://github.com/ylab-hi/atroplex/issues/60)): `string_view` for 8 utility/logging functions; `[[nodiscard]]` on registry `intern()` and `Reader::lookup()`; named constants `BITS_PER_WORD`, `MAX_SOURCES`, `PROGRESS_INTERVAL`; `walk_exon_chain()` takes `const grove_type&`.
- **Splicing event catalog opt-in** ([#65](https://github.com/ylab-hi/atroplex/pull/65)): `--events` flag gates the per-gene splicing event TSV (cassette, alt-5'/3', IR, alt-terminal, mutually exclusive). Off by default — the file can reach tens of GB on cohort-scale builds. Hub analysis always runs regardless. When disabled, both exon chain capture and event detection are skipped.
- **Always-on per-chromosome logging in inspect** ([#64](https://github.com/ylab-hi/atroplex/pull/64)): `Inspecting chrN (X/Y)...` now prints via `logging::info` regardless of `--progress`. Fixes garbled output when `--progress` was enabled (missing newline after final progress line).
- **Per-leaf progress for inspect** ([#51](https://github.com/ylab-hi/atroplex/pull/51)): `--progress` now shows per-chromosome leaf-node and segment counts during the grove traversal, so long-running inspections show visible progress.
- **`--chromosomes` build filter** ([#50](https://github.com/ylab-hi/atroplex/pull/50)): restrict index building to specific chromosomes (e.g. `--chromosomes chr1,chr22,chrX`). Features on unlisted chromosomes are silently skipped at ingest. Accepts both prefixed (`chr1`) and bare (`1`) names. Recorded in `.ggx.summary` build parameters. Applied in both GFF and BAM build paths.
- **Build parameters in `.ggx.summary`** ([#49](https://github.com/ylab-hi/atroplex/pull/49)): the summary now records order, absorb, fuzzy_tolerance, scaffold/loci filters, expression thresholds, and chromosome filter for reproducibility.
- **`--min-samples` for inspect** (#48): non-destructive sample-count filter applied during the grove traversal. Segments present in fewer than N samples are skipped — all downstream outputs (overview, sharing, hubs, events) reflect the threshold. The grove is not modified, so different thresholds can be explored without rebuilding.
- **Integration tests for .qtx remap and transitive chain** ([#55](https://github.com/ylab-hi/atroplex/pull/55), closes [#37](https://github.com/ylab-hi/atroplex/issues/37), [#38](https://github.com/ylab-hi/atroplex/issues/38)): end-to-end tests exercising reverse absorption → `absorbed_into_idx` → `remove_tombstones` remap → `merge_to_qtx` heap re-injection → `Reader::lookup()`. Tier 2.1 verifies two-sample remap; Tier 2.2 verifies A→B→C transitive path compression.
- **Absorption test coverage** ([#52](https://github.com/ylab-hi/atroplex/pull/52), closes [#10](https://github.com/ylab-hi/atroplex/issues/10)): exact assertions for Rules 0/3/4/8 with dedicated fixtures. Rules 3 and 4 now verify ref/sample distinction (drop vs keep). Rule 8 has an explicit intergenic mono-exon test.

- **`--annotated-loci-only` build filter** (#43): drops sample transcripts that don't spatially overlap any annotation segment on the same strand. Novel isoforms at annotated loci inherit the annotation's `gene_idx`, so sample-specific gene_ids (MSTRG.*, ENCLB*) no longer inflate the gene count. Annotation transcripts always pass. Docker CI now pushes PR images tagged `pr-<N>` for HPC testing; build log shows per-file segment delta; overview renames `transcripts` to `source_transcript_ids`.

### Changed
- **Decompose `collect()` and `finalize_gene` into focused methods** ([#70](https://github.com/ylab-hi/atroplex/pull/70), closes [#59](https://github.com/ylab-hi/atroplex/issues/59)): splits the 740-line `collect()` and 332-line `finalize_gene` lambda into 8 extracted methods. Internal types (`pending_hub`, `gene_acc`) and hub buffers promoted to class members. `collect()` is now a ~60-line orchestrator.
- **Decompose `create_segment` into rule helpers** ([#68](https://github.com/ylab-hi/atroplex/pull/68), partial [#59](https://github.com/ylab-hi/atroplex/issues/59)): extracts `try_terminal_variant_merge`, `handle_mono_exon`, `try_subsequence_absorption`, `apply_gene_idx_inheritance` from the 233-line function. `create_segment` is now a ~40-line orchestrator.
- **Consolidate build parameters into `build_options` struct** ([#67](https://github.com/ylab-hi/atroplex/pull/67), closes [#57](https://github.com/ylab-hi/atroplex/issues/57)): replaces 13-15 individual parameters on `build_from_samples`, `build_gff::build`, `build_bam::build`, `process_gene`, and `process_transcript` with a single struct. Net -63 lines.
- **Gene_idx inheritance always applied** ([#63](https://github.com/ylab-hi/atroplex/pull/63)): sample transcripts that spatially overlap an annotation segment now always inherit its gene_idx, regardless of `--annotated-loci-only`. Previously this was gated behind the flag, causing 19M inflated gene counts on 21K-sample builds without it (each StringTie gene_id created a separate gene entry). `--annotated-loci-only` now only controls whether novel-locus transcripts are kept or dropped.

### Fixed
- **Cross-sample gene_idx inheritance via spatial overlap**: `apply_gene_idx_inheritance` now falls back to the first overlapping sample-derived segment when no annotation parent is present. Previously the inheritance loop was gated on `is_parent_annotation(...)`, so two samples both contributing novel transcripts at the same locus (e.g. StringTie's `MSTRG.1` in two TCGA files spanning the same coordinates) created two separate `gene_idx` entries. On the 21K-sample comprehensive build this produced 2.8M `gene_id` strings instead of cross-sample-collapsed loci. Affects every per-gene downstream metric (`total_genes`, `single/multi_isoform_genes`, `mean/median/max_segments_per_gene`, `mean_gene_exon_entropy`, `mean_effective_isoforms`, per-source `genes`). Annotation-source parents still take priority for inheritance — the fallback only fires when no overlapping segment is annotation-typed. Tested via new `CrossSampleInheritance_NovelLocus` (two non-annotation fixtures with overlapping spans but disjoint exon chains).
- **`--min-samples` counted annotations toward threshold** ([#69](https://github.com/ylab-hi/atroplex/pull/69)): `--min-samples` now counts only `type="sample"` entries. Annotation-only segments (GENCODE reference isoforms) always pass regardless of the threshold. Previously `sample_count()` included annotations, so reference isoforms were filtered at `--min-samples 2`. Overview TSV gains `exon_links` metric and a `note` column.
- **Duplicate .qtx records at merge time** ([#62](https://github.com/ylab-hi/atroplex/pull/62), fixes [#35](https://github.com/ylab-hi/atroplex/issues/35), [#58](https://github.com/ylab-hi/atroplex/issues/58)): `flush_block` in `run_final_merge` now collapses duplicate `(segment, sample)` pairs by summing values before writing. Removes inconsistent reader-side workarounds. Also standardizes file-open error logging in `analysis_report` and adds stream read checks to `.ggx` deserialization.
- **`is_annotation_sample` checked `annotation_source` field, not just `type`**: samples with `annotation_source` filled in the manifest (common in StringTie GTFs built against a reference) were misclassified as annotations, bypassing the `--annotated-loci-only` filter and gene_idx inheritance. Fix: check `info.type == "annotation"` only.
- **Annotation mono-exon genes dropped by Rules 6/7/8**: GENCODE mono-exon genes (miRNAs, snoRNAs, ~4K genes) were classified as intergenic/gene-overlap and discarded because no multi-exon segment existed at their locus yet. Fix: annotation transcripts skip mono-exon Rules 6/7/8 — they always create segments as the reference.
- **Gene_idx inheritance ran for annotations**: overlapping GENCODE genes on the same strand would steal each other's gene_idx via the inheritance block. Fix: inheritance guarded by `!is_annotation_sample`.
- **Build summary overcounted genes**: `build_summary::collect()` did not skip absorbed or zero-attribution segments. Fix: added skip matching inspect's logic.

### Changed
- **Spatial absorption: gene_id removed as structural key** (#41, fixes #40): absorption candidate lookup in `segment_builder::create_segment` and `try_reverse_absorption` now uses `grove.intersect()` (spatial) instead of the `gene_segment_index` (gene_id hash). The `gene_segment_index` data structure (`chromosome_gene_segment_indices`) and all threading through build_gff / build_bam / builder / subcall are removed. `gene_idx` dropped from `exon_feature` (gene info lives on segments only; exons derive it from the parent segment during traversal). Breaking: `.ggx` serialization format changes (exon no longer carries gene_idx) — rebuild required. Fixes 3.8M inflated gene count on mixed-source cohorts (StringTie MSTRG IDs no longer create separate gene entries when their exon structure matches an existing annotation segment). Enables cross-annotation builds (GENCODE + RefSeq + custom GTFs merge spatially). Net -36 lines across 21 files.
- **`-g/--genogrove` now takes a directory** (#39): instead of passing a `.ggx` file path, pass the directory containing the index. The command scans for exactly one `.ggx` (errors on 0 or >1) and auto-discovers the sibling `.qtx` sidecar. Matches the write-side convention where `build` writes to `-o` with `--prefix`. Passing a file directly errors with a clear message. One-way break — update scripts that pass `.ggx` paths.
- **Rename `analyze` → `inspect` + flatten output layout** (#33): the pan-transcriptome analysis subcommand is now `atroplex inspect` to better describe what it does (structural inspection: overview / sharing / hubs / events — not "analysis" in the statistical sense). Outputs land directly under `-o/--output-dir` in category subfolders (`overview/`, `sharing/`, `splicing_hubs/`, `splicing_events/`) — the previous `analysis/` wrapper directory is gone, matching the convention of `build` / `query` / `discover` / `export`. The old subcommand name is not aliased (one-way break — update scripts).
- **Sidecar-backed expression restored in `inspect`** (#33, batch 3b): per-sample `.{expr_type}` columns (e.g. `SAMPLE1.counts`) are back on `splicing_hubs.tsv`, `branch_details.tsv`, and `conserved_exons.tsv`; `expressed_segments` + `mean_expression` return to `per_sample.tsv`. Values are read from the `{prefix}.qtx` sidecar at inspect time via a new `quant_sidecar::Reader*` parameter on `analysis_report::collect()`, with per-exon aggregation via a per-chromosome `exon_expr_sum` accumulator (keyed by `key_ptr`, cleared at every chromosome boundary — bounded by exons-per-chromosome, no sample×feature matrix). Columns are only emitted when the sidecar loaded successfully, so no-sidecar runs produce clean TSVs without all-`.` columns.
- **Quantitative DTU restored in `query --contrast`** (#33, batch 3b): chi-squared test on the 2×k gene contingency table with Wilson-Hilferty p-value approximation + Benjamini-Hochberg FDR. Replaces the error-throwing stub added in 3a. Per-sample expression for each matched transcript is looked up from the `.qtx` sidecar in `classify_transcripts` and stored on `query_result::sample_expression`; DTU consumes it directly. `--contrast` without a sibling `.qtx` now fails fast with a clear error pointing at the rebuild requirement.
- **Drop `query --group-by`** (#33): `sample_info::group` is now the only grouping axis for DTU. `--contrast A:B` uses group-column values directly (manifest `group` column, or auto-inferred from `_repNN` suffix). Samples whose `group` is empty or `.` are skipped in DTU. Simplifies the CLI and matches the product decision to use `group` everywhere (including the future Phase 8.8 group-aware `inspect` aggregation).

### Fixed
- **Reverse absorption preserves .qtx expression records** (#36, fixes #34): when `try_reverse_absorption::absorb_into_parent` tombstones a candidate segment, any `.qtx` records already written under the candidate's `segment_index` were dropped by `merge_to_qtx`'s tombstone filter. New `absorbed_into_idx: optional<size_t>` on `segment_feature` records the parent link at absorption time. `builder::remove_tombstones` derives a `tombstone_remap` map with transitive path compression and hands it to `merge_to_qtx`. `run_final_merge` re-injects remapped records into the merge heap at the parent's index — they sort naturally into the parent's block since parents always have a higher `segment_index` (monotonic creation order). Tombstones without a parent (Rule 3/4 drops via `tombstone_candidate`) remain in the drop set as before. Breaking: `.ggx` serialization format gains the new field (rebuild required). 4 new tests cover basic remap, remap-over-drop priority, merge-into-existing-parent, and empty-remap backward compatibility.
- **Expression dropped on every sample-side merge** (#33, build-path bug from 3a): `segment_builder::merge_into_segment` received `expression_value` but never wrote it to the sidecar — only the "Step 5: Create new segment" branch of `create_segment` wrote. Because the builder sorts annotations first and samples after, every segment was created from an annotation transcript (no expression → nothing written) and every sample transcript then merged via Rule 0 exact / Rule 5 terminal / fuzzy-FSM / Rule 2 ISM_3PRIME, silently dropping its `counts`. Result: `.qtx` files stayed at the 84-byte empty-header baseline even on fully populated manifests with valid `expression_attribute` declarations. Fix threads `sidecar_writer` through `merge_into_segment` and appends the record there when `expression_value >= 0`. Known narrower edge case still open: `try_reverse_absorption` loses expression for the tombstoned-candidate case (sample→sample only; annotation-first ordering avoids it in practice).

### Changed
- **Breaking: expression no longer stored in the grove** (batch 3a): `expression_store` removed from `segment_feature` and `exon_feature`. Per-sample expression values (counts / TPM / FPKM / RPKM / cov, after the `--min-*` filters pass at ingest) are now written to per-sample `.qtx` sidecar files at `{output_dir}/{prefix}.quants/{sample_id}.qtx` during build, instead of being stored in the in-memory grove. This is the primary memory fix that unblocks cohort-scale TCGA/GTEx builds — the projected ~150 GB of per-segment flat `vector<float>` storage at 21K samples is gone. The `.ggx` format no longer carries expression fields on segments or exons (one-way break — rebuild required). Scaffolding for filter logic (`--min-counts`, `--min-TPM`, etc.) and `sample_info::expression_attributes` / `expr_type` is preserved; the value is still extracted from GTF attributes for the filter gate, just no longer stored on the segment after passing. `query --contrast` is temporarily disabled with a clear error pointing at the follow-up PR that adds sidecar-read quantitative DTU. `analyze` outputs drop the per-sample `.{expr_type}` columns from `splicing_hubs.tsv`, `branch_details.tsv`, `conserved_exons.tsv`, and the `mean_expression` / `expressed_segments` columns from `sample_stats.csv` (all restored in the follow-up PR via sidecar reads). The `--min-expression` legacy flag is removed.

### Added
- **Per-hub PSI and entropy in `splicing_hubs.tsv`**: each hub row now includes two new per-sample columns, `{label}.psi` and `{label}.entropy`, alongside the existing `.branches`/`.shared`/`.unique`. PSI is the fraction of the sample's segments in the hub's gene that pass through the hub exon (segment-level, not read-weighted — sparse-friendly at cohort scale). Entropy is the Shannon entropy over the sample's usage distribution across the hub's downstream branch targets: `H = −Σ p_i · log₂(p_i)`, where `p_i` is the fraction of the sample's hub-touching segments that continue to target `i`. The computation runs in two passes over `gene_acc.segments_in_gene`: Pass 1 builds the per-sample PSI denominator once per gene (shared across every hub in the gene); Pass 2 walks each hub's targets per sample, accumulating numerator + branch counts. Three new outer-scope buffers (`hub_psi_num`, `hub_psi_den`, `branch_counts`) are allocated once in `collect()` and reset with `std::fill` between hubs, so peak memory stays at `O(num_samples × global_max_targets)` independent of hub count. Analyze output chain capture is now also triggered by `hub_stream` being armed (previously only `splicing_events_stream`), because the PSI/entropy passes need the per-segment exon chain
- **Scaffold filter at ingest** (#30): new `--include-scaffolds` CLI flag (default: off) gates GFF/BAM ingest to canonical main chromosomes only — `chr1..chr22`, `chrX`, `chrY`, `chrM` (and the bare/`MT` variants). Everything else (unplaced contigs, alt haplotypes, fix patches, decoys) is skipped at the entry loop and counted as `scaffold_filtered_transcripts` in `build_counters`. The filter runs before any grove insertion, so scaffold features never touch the spatial index or the registries. Cohort-scale TCGA/GTEx builds typically see ~10–15% segment-count reduction with no loss of main-chromosome biology. New `is_main_chromosome()` helper in `utility.hpp` handles both prefixed and unprefixed forms with strict autosome validation (rejects any suffix like `chr1_KI270706v1_random`, `chrUn_*`, `*_alt`, `*_fix`, `chrEBV`). `.ggx.summary` now emits a `Scaffold-filtered: N` line when any were filtered
- **`quant_sidecar` module** (#30): new `include/quant_sidecar.hpp` + `src/quant_sidecar.cpp` defining the `.qtx` per-sample quantification sidecar format (`AQTX` magic, 32-byte header, 12-byte `(segment_index, value)` records sorted by `segment_index`). `Writer` supports buffered append + sort-and-flush with atomic tmp/rename; `Reader` supports O(log N) `lookup` via binary search and templated `for_each` linear scan. Ships as an unused primitive in this PR — the follow-up PR wires sidecars into the build path and removes `expression_store` from the in-memory grove to address cohort-scale OOM

### Fixed
- **CI runner OOM during test binary compilation** (#27): every test binary previously re-compiled the full `src/` glob via `ATROPLEX_LIB_SOURCES`, so `cmake --build --parallel` kicked off `N × cores` parallel g++ processes. After adding a fourth test binary (`builder_pipeline_tests`), GCC 13/14 Debug and Release jobs silently OOM'd on the 16GB ubuntu-24.04 runners — the log went quiet for 47 minutes before "runner lost communication" death. Clang passed because its per-TU memory is materially lower. Fixed by promoting the source glob to a new `atroplex_core` STATIC library target that the main executable and every test binary link against — library sources compile exactly once instead of 5×
- **Splice site hash collisions**: replaced weak XOR-based hash (`seqid_hash ^ position<<1`) with proper `(seqid, position)` composite keys using boost::hash_combine pattern — eliminates false positive "known" splice site hits across chromosomes that caused NIC/NNC misclassification
- Fixed potential `size_t` underflow in splice site tolerance window when position is near zero
- Fixed infinite loop in `sample_bitset` iterator: `advance()` did not mark iterator as exhausted when all bits were consumed, causing range-based for loops over `sample_idx` to hang
- Fixed missing `absorb` parameter in final `process_gene()` call per GTF file (last gene always skipped absorption)
- Fixed silent error swallowing (#18): replaced bare `catch (...) {}` with `catch (const std::exception& e)` + `logging::warning` in GFF expression parsing and BAM header parsing — malformed values now produce visible diagnostics instead of silent pass-through
- Fixed discover read coverage (#20): matched segments now record full cluster read count instead of just 1 (the representative); removed `supporting_reads` vector to avoid memory blowup on large BAMs

### Added
- **Streaming `analysis_report` for `atroplex analyze`** (#27, Phase 8.1): new `include/analysis_report.hpp` + `src/analysis_report.cpp` replaces legacy `index_stats::collect()` which OOM'd on ~1000 samples. Single-pass B+ tree traversal with flat `vector<sample_counters>` indexed by sample_id, per-gene accumulators cleared at each chromosome boundary, and `transcripts_per_gene` / `exons_per_segment` distribution vectors bounded by feature count. Zero sample×feature maps. Output: `analysis/overview/{basename}.overview.tsv` + `{basename}.per_sample.tsv`. Legacy `index_stats::collect()` commented out in `subcall::analyze` pending phases 8.2–8.7
- **`build_summary` class extracted from builder** (#27): `include/build_summary.hpp` + `src/build_summary.cpp`. Owns the `.ggx.summary` text writer and collects its data from builder caches post-sweep. Replaces `index_stats` in the build path; `subcall::build_stats` is now `std::optional<build_summary>`. ~150-line inline stats block moved out of `builder::build_from_samples`
- **`build_counters` transcript processing tracker** (#27): new struct threaded explicitly through `build_gff`/`build_bam`/`segment_builder` (no default args — the five categories are always populated). Tracks `input_transcripts`, `merged_transcripts` (Rules 0/2/5/fuzzy-FSM folded into existing segments), `absorbed_segments` (tombstoned and swept), `discarded_transcripts` (min-expression, mono-exon, fragment drops), and `replicates_merged`. Surfaces in the `.ggx.summary` "Transcript processing:" section
- **Physical tombstone removal** (#27): `builder::remove_tombstones()` runs after `merge_replicates`, before summary collection. Walks `gene_indices` (authoritative tombstone source — `segment_caches` are eagerly erased by `try_reverse_absorption`), calls `grove.remove_key()` on each absorbed segment, then a single `grove.remove_edges_if()` pass strips orphan `EXON_TO_EXON` edges carrying the tombstone's `segment_index`. Requires genogrove v0.21.0
- **Builder pipeline test suite** (#27): `tests/build/builder_pipeline_test.cpp` with 4 tests covering the first-ever end-to-end exercise of `builder::build_from_samples`: counter aggregation (`BuilderFullPipeline_CountersPopulated`), physical tombstone removal from the B+ tree (`RemoveTombstones_PhysicallyRemoved`), orphan edge pruning via swept-vs-unswept `edge_count()` comparison (`RemoveTombstones_OrphanEdgesPruned`), and `.ggx.summary` text file contents (`BuildSummary_WrittenFileContainsCounters`). Generates temporary GTF fixtures in a per-test tmp directory
- Updated genogrove dependency to v0.21.0 (adds `grove.remove_key()` B+ tree key removal API with rebalancing from upstream #305)
- `atroplex export` subcommand (#19): reconstruct per-sample GTF files from a `.ggx` index with filters (`--sample`, `--gene`, `--region`, `--min-samples`, `--conserved-only`, `--biotype`, `--source`)
- Reusable `gtf_writer` utility for GTF line formatting (gene/transcript/exon)
- CI now triggers only once per PR (push to main + pull_request, no duplicate runs)
- Serialization roundtrip test suite (#17): 8 tests verifying .ggx index write/read roundtrip preserves registries, spatial queries, graph edges, expression values, and sample membership
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

### Performance
- Replaced Jaccard diversity with O(N) entropy metrics (#25): eliminates O(N² × E) pairwise computation in Phase 4b; new `gene_exon_entropy` and `effective_isoforms` metrics scale to large groves without the 200-segment cap
- Added progress logging to analyze phases (#25): `collect()` now emits per-phase status messages for long-running analyses
- Added early-exit filters to absorption matching (#23): O(1) exon count and span guards skip candidates that cannot possibly be subsequence matches, avoiding expensive O(N_parent × N_sub) scans
- Added explicit `subsequence_type::FSM` for fuzzy full-length matches (#23): replaces the old `ISM_3PRIME` shortcut with a semantically correct type
- Cached exon counts before transcript sort (#22): eliminates O(T·E·log T) rescans in `process_gene()`, speeding up per-gene transcript ordering

### Refactored
- Removed `exon_caches` dependency from `transcript_matcher` (#21): splice site indexing now walks the grove directly, enabling `--genogrove` for discover/query
- Added `splicing_catalog::collect_from_grove()` (#21): reconstructs per-gene segment chains from grove traversal, enabling `--genogrove` for analyze

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
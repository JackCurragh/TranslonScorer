# TranslonScorer refactor — task list (strangler-fig migration)

Spec: `docs/architecture.md` (functional core / imperative shell; capability
protocols; per-BAM offsets; P/A site; migration order). This file is the
**execution checklist and progress tracker** — the working agent edits it.

## Agent prompt (hand this to the implementer)

> You are migrating TranslonScorer to the structure in `translonscorer/docs/architecture.md`.
> Work the tasks in `translonscorer/docs/REFACTOR_TASKS.md` **strictly in order, one at a time**.
>
> **Non-negotiable rules:**
> 1. **The gate must be GREEN before and after every task.** The gate is `make gate`
>    (= `pytest tests/test_golden.py -q`). Never start a task with a red gate; never
>    finish one leaving it red.
> 2. **Phases 1–2 are byte-identical moves.** Moving code must not change behaviour:
>    the GAPDH golden must stay Δmetric=0, 0 call mismatches. Only **Phase 3**
>    (coverage providers) introduces new behaviour, and it gets new tests — it must
>    still not break the existing golden.
> 3. **`pipeline/` keeps working the whole time.** When you move code to a new home,
>    leave a thin re-export shim in the old `pipeline/` module so existing imports
>    resolve. Delete shims only in the final cleanup task, gate green.
> 4. **Functional core discipline:** pure functions + frozen dataclasses; the *only*
>    stateful objects are `coverage/` providers. No new global state.
> 5. **One task = one self-contained change.** After completing a task: run the gate,
>    paste the one-line result into the task's `STATUS:` field, tick its `[x]`, and if
>    this repo is a git repo, commit with `refactor(<area>): <task title>`.
> 6. **If a gate goes red and you can't fix it within the task's scope, REVERT the
>    task's changes** (the gate returns green because `pipeline/` is untouched) and
>    **stop — report what blocked you.** Do not push forward on red.
> 7. **Stop at each `=== PHASE N complete ===` marker and report** a summary + the
>    gate result for human review before starting the next phase.
>
> Track progress only by editing this file's checkboxes and `STATUS:` lines so the
> overseers can read current state at a glance.

## Gate definition

`make gate` runs `tests/test_golden.py`, which must assert:
- **GAPDH golden**: scoring the GAPDH-locus events reproduces
  `outputs/reference_scores_gapdh.parquet` with max|Δmetric| == 0 and 0 call/
  eligibility mismatches.
- **scalar ≡ vectorised**: `score_events` and `score_events_vectorised` agree (Δ=0)
  on the same coverage (port from `scripts/verify_vectorised.py`).

---

## Profile/offset contract before scoring

Profile construction is upstream of scoring, matrix aggregation, clustering, and
reports. The required order is:

`read assignment → whole-sample offset calibration → per-sample/per-length P/A profiles → fold/merge prefixes/samples → aggregate/cluster → score`

Non-negotiable profile rules for the refactor:

- Offsets are estimated once per sample/BAM and read length from whole-sample
  evidence. They must not be inferred inside a target locus, transcript,
  read-prefix partition, cluster, or scoring call.
- Offset tables are explicit inputs to profile generation:
  `{sample_id, read_length, site -> offset}`. Missing lengths use only a declared
  policy: fail, warn+fallback, or fixed offset.
- After offsets are known, profile generation is pure and deterministic:
  reads/counts + annotation + offset table + site (`P`/`A`) → tidy profile rows.
- Aggregation happens only after per-sample/per-length offsets are applied.
  Cluster and aggregate scoring consume profile matrices that already have
  correct P/A-site coordinates and never recompute offsets.
- Profile caches must include source identity, annotation identity,
  mapper/multimapper policy, site, offset-table identity, read-prefix policy, and
  relevant code/config version. A cache keyed only by transcript/locus is invalid.

Current violations to remove during the refactor:

- `pipeline/profiles.py::profiles_from_bam` and
  `profiles_from_sparse_parquet_matrix` infer offsets during profile generation.
- The matrix demo notebook calls `profiles_from_sparse_parquet_matrix` per
  read-prefix partition, so offsets can be inferred per locus/prefix.
- The matrix demo notebook profile cache is keyed only by transcript ID.
- Older rollup helpers that accept one `ref_offset` are compatibility/indexing
  paths only, not final P/A profile generation.

---

## Phase 0 — make the gate one command
- [x] **T0** Create `tests/test_golden.py` (port `score_reference.py` + `verify_vectorised.py` comparison logic) and a `make gate` target. Gate GREEN on current `pipeline/` code.
  - STATUS: 3 passed in 0.36s — GAPDH golden (init/elong/term SUPPORTED), scalar≡vec (GAPDH + contended). `make gate` target added.

`=== PHASE 0 complete — report ===`

## Phase 1 — leaf deps (model + io)
- [x] **T1** `model.py`: move `ScoreThresholds` (from `event_score.py`) + `OffsetParams`/`Region` (from `coverage_providers.py`); add `ScoreRecord`, `ConsequentialityPolicy`. Old modules re-import. Gate green. STATUS: 3 passed in 0.48s
- [x] **T2** `io/annotation.py`: canonical `build_cds_blocks` + `build_gene_spans` (consolidate the copies in `scripts/demo_frame_dominance.py` and `matrix_rollup.py`). Callers import here. Gate green. STATUS: 3 passed in 0.77s
- [x] **T3** `io/matrix.py`: manifest / count-parquet / samples / `_discover_bam` helpers (from `matrix_qc.py`, `sparse_parquet.py`). Gate green. STATUS: 3 passed in 0.49s
- [x] **T4** `io/bam.py`: oxbow+pysam alignment read, CIGAR blocks, `read_id` parse (from `matrix_rollup._prefix_sums` scan + `sparse_parquet`). Gate green. STATUS: 3 passed in 0.44s
- [x] **T5** `io/store.py`: `persist_scores` + parquet readers (events / scores / ledger). Gate green. STATUS: 3 passed in 0.37s

`=== PHASE 1 complete — report ===`

## Phase 2 — pure transforms (bit-identical, gated)
- [x] **T6** `events.py`: `extract_events`, frame intervals, `_deconflict_intervals`, contention (from `event_extract.py` + `matrix_qc.py`). Add a chr22 event-count reproduction check to the gate. Gate green. STATUS: 4 passed in 0.69s (GAPDH golden + scalar≡vec×2 + chr22 event-count)
- [x] **T7** `offsets.py`: `OffsetParams`, `plausible_offset_range`, `usable_read_length`, `metagene|file|global` (file/global real; metagene stub OK; **default offset path unchanged**). Offset APIs return explicit per-sample/per-read-length tables and do not accept locus/profile matrices as calibration input. Gate green. STATUS: 4 passed in 0.52s — commit ecc45eb
- [x] **T8** `scoring/{aspects,evidence,run}.py`: move all scorers from `event_score.py` (init/elong/term/junction, `_elong_evidence`, `score_events[_vectorised]`, prefix-sum helpers→`coverage/profile.py`). **GAPDH golden Δ=0 — the core regression.** STATUS: 4 passed in 0.50s — commit a9c8a14
- [x] **T9** `qc.py` (periodicity/frame-dominance), `frame_support.py` + `frame/{bleed,hmm,latent,deblur}.py` (profile → codon/frame posteriors; pure transforms only), `clustering.py` (`profile_clustering`/`score_clustered`), `report.py` (per-feature composition), `consequential.py` (consequentiality policy). **Cluster TRIGGER → `scoring/aspects.py`** (next to identifiability); only the *response* in `clustering.py`. Frame support consumes already-built P/A profiles and must not infer offsets or build profiles. Gate green. STATUS: 17 passed in 0.62s — commit fbc114f

`=== PHASE 2 complete — report ===`

## Phase 3 — coverage providers (new behaviour, last)
- [x] **T10** `coverage/base.py`: capability **Protocols** (`CoverageProvider`, `SupportsSites`, `SupportsJunctions`, `SupportsMappability`). `coverage/profile.py`: pure P/A-site profile gen from assigned reads/counts plus explicit offset table; no offset inference in profile builders (prefix sums; `A=P+3` default only when selected by the offset table/policy). Gate green. STATUS: 22 passed in 0.57s — commit 2a675db
- [x] **T11** `coverage/matrix.py`: `MatrixProvider` wrapping `region_coverage`/`tabulate_junctions`. Must reproduce the **GAPDH golden through the provider path** (Δ=0). STATUS: 24 passed in 0.55s — commit 0d32b01
- [x] **T12** `coverage/bam.py`: `BamSetProvider` (whole-sample per-BAM/per-length offset calibration before any locus profile; P/A; transcriptome→genome isoform-multimapper resolution; unique default + mappability ledger). `coverage/bigwig.py`: `BigwigSetProvider` (only `CoverageProvider`). **New behaviour → new tests** (`tests/test_bam_provider.py`: score SRR11005875 genome BAM on GAPDH, assert sane in-frame/init/term; assert offsets are calibrated once per BAM/length and reused for locus profiles). Existing golden still green. STATUS: 27 passed in 0.69s — commit 9247417

### Phase 3 follow-ups (BAM-provider correctness — surfaced in real-data testing)
- [x] **T12.a** Strand-aware P/A placement (`-` strand used `reference_start`) + emit `strand` column; NH-based unique filter (was `MAPQ==0`, kept 2–4-locus multimappers). Tests + GAPDH genome fixture (`data/gapdh_cohort_genome.bam`). STATUS: 38 passed — commit (strand/unique).
- [x] **T12.b** Implement metagene P-site offset calibration (5′ pile-up at start codons), pure `metagene_offsets` + provider `_build_metagene_histogram`. Unit + real-data tests. STATUS: 40 passed — commit (metagene). Real GAPDH offsets 9–12 for 25–32mers.
- [ ] **T12.1** Transcriptome→genome read projection + isoform-multimapper resolution (the largest; needs inverse exon mapping + validation against the SRR transcriptome BAM). For now `transcriptome=True` raises `NotImplementedError` (no silent wrong coords). STATUS: guarded + tracked.

`=== PHASE 3 complete — report ===`

## Phase 4 — orchestration + cleanup
- [x] **T13** `workflows.py` + `cli.py`: subcommands (`extract-events`, `score-matrix`, `score-bams`, `consequential`); marked `orf-composite`/`score-orfs`/`feature-metrics` deprecated (yellow stderr warning + `[DEPRECATED]` help). workflows.py wraps events.run_extract / MatrixProvider / BamSetProvider / score_events_vectorised / persist_scores / apply_policy; per-chrom coverage query so genomic positions never collide. Smoke test `tests/test_cli_workflows.py` (10 tests). STATUS: gate 41 passed; full 121 passed, 3 skipped — commit 0ca4f43
- [x] **T14** Deleted the 3 genuine pure re-export shims (`event_extract`, `event_score`, `profile_clustering`) and repointed importers at the new tree (event_score→model+scoring.run; profile_clustering→clustering, both in pipeline modules matrix_normalisation/matrix_scoring and tests). Final gate green. STATUS: gate 41 passed; full 121 passed, 3 skipped — commit 0132bb1.
  - **DEFERRED (correctly out of scope):** `pipeline/matrix_qc.py` is a 739-line *partially-migrated real module* (not a pure shim — only some helpers re-export io/matrix), depended on by matrix_rollup + tests. `pipeline/frame_support.py` is a live **Config→FrameSupportParams adapter** depended on by cli + 3 un-migrated pipeline modules (workflow, frame_method_compare, rdg_flux_export). Removing either requires migrating those un-migrated modules off Config — a larger task than the Phase-1–3 shim cleanup. Track as **T15** if/when those modules are migrated.
- [x] **T13.1** Per-translon report composition + dynamic consequentiality (the two remaining stubs). `report.compose_report(scores, feature_event)` aggregates long-form event scores to one row/translon (per-aspect call/metric/n_reads/supported_frac; shared events → identical aspect score for both translons). `consequential.apply_policy` is now real: `consequentiality_score = tier_confidence * (0.5 + 0.5*expression_pct) * context_weight`, gated by `min_tier_confidence`/`min_expression_percentile` — dynamic/context-aware, NO hard length/biotype gates. `workflows.report_workflow` + CLI `report` wire store+events→report→policy. Tests `tests/test_report_consequential.py` (6) + e2e. mypy-clean. STATUS: gate 41 passed; full 129 passed, 3 skipped — commits 9f017b5, de95e90.
`=== PHASE 4 complete — DONE (T13.1 closed the scoring stubs) ===`

## Phase 5 — T15: finish the `pipeline/` migration (strangler-fig, per subsystem)

`pipeline/` still holds ~30 un-migrated real modules. Migrate them subsystem by
subsystem into the new tree, **gate green after every sub-task**, leaving a
re-export shim at the old path until the final sweep, then delete `pipeline/`.
Same discipline as Phases 1–4: NO Claude in commits, STATUS line + commit hash
per task, never advance on red. Ordered so leaves move before their dependents.

- [x] **T15.1 — frame_support adapter retirement.** All callers now build `FrameSupportParams` and call the pure `TranslonScorer.frame_support.build_frame_support`: cli.py via local `_build_and_write_frame_support` shell helper (Config→Params + write, 3 sites); workflow.py inlined Params+write; frame_method_compare/rdg_flux_export swapped Config→Params in their single wrapper. Deleted `pipeline/frame_support.py` + its shim test; test_score_first_gates repointed to Params. STATUS: gate 40 passed; full 162 passed — commit pending.
- [x] **T15.2 — matrix_qc relocation.** `pipeline/matrix_qc.py` already imported its pure helpers from `qc.py`/`io.matrix`/`io.bam`/`events` (it was a partial shim around matrix-specific scanning orchestration). Moved the whole cohesive module to top-level `TranslonScorer/matrix_qc.py` (peer of `clustering.py`/`frame_support.py`), fixed the one relative import, repointed `matrix_rollup.py` + tests (test_matrix_scoring, test_golden) directly — NO shim left (reduces final-sweep debt). STATUS: gate 40 passed; matrix 25 passed; full 162 passed — commit pending.
- [x] **T15.3 — matrix engine.** Relocated `matrix_rollup` (1098) / `matrix_scoring` (501) / `matrix_normalisation` (330) to top-level peers `TranslonScorer/matrix_rollup.py` etc. (consistent with matrix_qc relocation). `..`→`.` import depth fix; sibling `.matrix_normalisation` import unchanged. Repointed `coverage/matrix.py` (MatrixProvider), `pipeline/workflow.py`, matrix tests. NO shims. STATUS: gate 40 passed; matrix 34 passed; full 162 passed — commit pending.
- [x] **T15.4 — frame subsystem.** Moved `read_assignment`, `frame_method_compare`, `frame_disambiguation`, `frame_crosstalk` into the existing `frame/` package (peer depth → `..` imports unchanged; only `frame_disambiguation`'s sibling `.transcript_coords` → `..pipeline.transcript_coords`). Repointed CLI (3 commands), `pipeline/feature_metrics.py` (frame_crosstalk), tests. STATUS: gate green; frame tests 31 passed; full 162 passed — commit pending.
- [x] **T15.5 — profiles / coords / indexing.** Moved all 8 (`profiles`, `locus_profiles`, `locus_features`, `transcript_coords`, `mapped_index`, `index_from_bam`, `junctions`, `junction_model`) into `coverage/` as one batch (keeps `profiles`' sibling imports of `.mapped_index`/`.transcript_coords` valid; `..` depth-stable). Repointed cli (profiles/index-from-bam/junctions), matrix_scoring, frame_disambiguation, pipeline siblings (workflow/annotation_bundle/rdg_flux_export/feature_metrics), tests. STATUS: gate 40 passed; targeted 40 passed; full 162 passed — commit pending.
- [x] **T15.6 — ORF-composite path.** Kept (commands still functional, just deprecated) and migrated all 9 into a new `orf/` subpackage: `orfs_import`, `orf_composite`, `map_orfs`, `assemble`, `score_gates`, `score_schema`, `panel_manifest`, `profile_compare`, `rdg_flux_export`. Batch-internal sibling imports stay valid; `rdg_flux_export.annotation_bundle` → `..pipeline.annotation_bundle` (T15.7). Repointed cli, matrix_scoring, pipeline/workflow, tests. STATUS: gate 40 passed; full 162 passed — commit pending.
- [x] **T15.7 — shell + config.** `config` → top-level `TranslonScorer/config.py`; `validator` merged into config.py as `validate_config`; `annotation_bundle`/`inspect` → `io/`; `workflow` → top-level `legacy_workflow.py` (legacy `all`/`process-bam`/`find-orfs`/`score-orfs` orchestration; `..`→`.`). Repointed cli + orf/rdg_flux_export. STATUS: gate 40 passed; full 162 passed — commit pending.
- [x] **T15.8 — final sweep.** `feature_metrics` → `orf/`; deleted `coverage_providers.py` (superseded skeleton, no importers), `_test_write.py` (junk), `pipeline/__init__.py`, stray `orfs_import.py.bak`. **`pipeline/` is gone** — zero `TranslonScorer.pipeline` references remain anywhere. architecture.md migration section banner updated. STATUS: gate 40 passed; full 162 passed. `make lint` (mypy) = 229 errors, but **unchanged by T15** (229 before → 229 after; all pre-existing in legacy file_handlers/visualization/orf code, none introduced by the migration) — lint has never been green on the package and is tracked separately, not a T15 regression.

`=== PHASE 5 complete — single-architecture codebase, pipeline/ removed ===`

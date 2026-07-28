# Changelog

All notable changes to TranslonScorer are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Removed
- **Flat-offset matrix coverage.** `score-matrix --ref-offset` is gone and
  `--psite-index` is now **required** (likewise in `pipeline` matrix mode).
  A single offset applied across all read lengths smears the P-site across
  frames and pinned `elong_in_frame` to the ~0.33 random floor, so the mode's
  only behaviour was to produce noise that looked like a score.
- **The FrameRollup *scoring* path.** `score-matrix --gtf`, `pipeline
  --cds-gtf/--cds-bigbed`, `workflows.score_matrix_rollup_workflow` and
  `scoring.run.score_elongation_from_rollup` are removed. It got offsets right
  but ran as a separate transcript-coordinate pipeline that never went through
  `_score_events_over_provider`: elongation only, `--context-gtf` and
  `--mappability-bigwig` ignored, `identifiability`/`breadth` hardcoded. The
  P-site index delivers the same offset accuracy for every event type through
  the shared scorer.
  The FrameRollup **utilities** (`build_frame_rollup`, `calibrate_offsets`,
  `score_frame_rollup`) are unchanged and still used by `build-psite-index`
  and the QC/figure scripts.

- **The second scorer.** `scoring.run.score_events_vectorised` is renamed to
  `score_events` and is now the only scorer in product code. The scalar
  implementation — which carried a "must reproduce this" docstring and so was a
  standing hand-sync drift risk — moved to `tests/reference_scorer.py` as
  `score_events_scalar`, where it stays useful as an independent oracle the
  shipped scorer is diffed against on every run.

### Changed
- **CLI surface: 27 top-level commands → 15 + 2 groups.** Nothing is deleted or
  deprecated; twelve commands moved under two groups so the top level reads as
  the product:
  - `research` — `score-compare-frame`, `score-compare-existing`,
    `compare-profiles`, `compare-frame-methods`, `frame-disambiguation`,
    `compare-read-assignment`, `validate-panel`, `export-rdg-flux`
  - `legacy` — `plot`, `orfs-import`, `assemble`, `map-orfs`

  Invoke as `translonscorer research compare-profiles …`. Scripts calling these
  at the top level need the group name inserted.
- `MatrixProvider` now requires `psite_index_dir` and raises if it is missing,
  rather than silently defaulting to the flat offset.
- Matrix and BAM/bigWig remain two first-class coverage strategies that differ
  only at `provider.coverage()`; everything downstream is one shared path.
- **`score_events` signature**: coverage is now a `cov_df` DataFrame
  (`pos`, `count`), matching what providers return. The old scalar took a
  `{pos: count}` dict; that shape survives only in the test reference.
- The GAPDH golden test now runs the **shipped** scorer. It previously ran the
  scalar reference, so the product was only validated against the golden
  transitively via the scalar≡vectorised comparison.

## [0.2.0] - 2026-06-26

### Added
- **FrameRollup scoring path** (`matrix_rollup.py`): reads projected to
  transcriptome coordinates, `phase0 = tx_pos % 3` stored offset-free.
  `build_frame_rollup` replaces the flat-offset genomic scan for elongation
  scoring.
- **`calibrate_offsets(target_frame=0)`**: per-(sample, length) P-site offset
  calibration that correctly finds offsets mapping in-frame reads to frame 0.
  Fixes the root cause of the ~33% random-floor periodicity on mixed-protocol
  cohorts. Optimised: evaluates 3 mod-3 classes instead of 9 raw offsets.
- **`score_frame_rollup`**: fully vectorised Polars implementation (join +
  pivot) replacing the Python dict accumulation loop.
- **`prevalence_from_rollup`**: per-feature sample-axis prevalence (FR5) — 
  fraction of samples showing in-frame periodicity above threshold. Pure Polars.
- **`score_elongation_from_rollup`** (`scoring/run.py`): maps feature-level
  FrameRollup scores to per-event elongation evidence in the standard schema.
- **`score_matrix_rollup_workflow`** (`workflows.py`): end-to-end FrameRollup
  scoring pipeline wired for use from the CLI.
- **`build_coverage_index`**: full `(sample, feature_id, tx_pos, length) →
  count` primary substrate (FR6) for per-translon read profiles and clustering.
  Shares `_fr_worker_init` / `_ci_worker` parallelism with `build_frame_rollup`.
- **`profile_from_index`**: aggregate CoverageIndex to `(feature_id, tx_pos) →
  count` with optional per-sample and per-length filters and normalisation.
- **`score-matrix --gtf`** CLI option: pass a GTF to activate the FrameRollup
  path with calibrated offsets; legacy flat-offset path retained as fallback.
- FR3 positive-control gate test: canonical GAPDH CDS scores SUPPORTED via
  FrameRollup. 25 tests in `tests/test_frame_rollup.py` (all green).

### Fixed
- Matrix elongation frame scoring was at the random 33% floor for mixed-protocol
  cohorts. Flat `ref_offset=15` applied to samples with calibrated offsets in
  different mod-3 classes (12, 13, 14) averages to 1/3. Fixed by
  `build_frame_rollup` + `calibrate_offsets(target_frame=0)`.

## [0.1.1]

- Initial alpha: CLI for scoring translational events from Ribo-seq data.

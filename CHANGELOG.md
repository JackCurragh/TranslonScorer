# Changelog

All notable changes to TranslonScorer are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2026-08-09

**Contains breaking changes.** `score-matrix --ref-offset` and the FrameRollup
scoring path are gone, and twelve CLI commands moved under `research`/`legacy`
groups. Pre-1.0, so these land in a minor bump; see Removed/Changed below.

### The output contract

Because the score store feeds downstream annotation work, the output schema is
treated as a public contract from this release on. One row per
`(event_id, aspect)`, with exactly **one headline metric per row** carried in
`metric` and named by `metric_name`. Secondary statistics live in the
`evidence` JSON blob. The four headline metrics are unchanged from 0.2.0:

| aspect | `metric_name` | meaning | scale |
| --- | --- | --- | --- |
| init | `init_rise` | consensus log2 fold-change, body / outer flank, over flanks (9, 18, 30, 60) | log2, unbounded |
| term | `term_drop` | consensus log2 fold-change, body / UTR; positive = drop | log2, unbounded |
| elongation | `elong_in_frame` | frame-0 share over *uncontended* positions, falling back to all positions | [0, 1] |
| junction | `junc_confident_spanning` | confident spanning-read support | count |

Metrics are on **different scales**; a consumer must switch on `metric_name`
before comparing or thresholding. This was already true and is now written down.

### Added
- **Translation signature scores** from Chothani et al., *Molecular Cell* (2022),
  as `TranslonScorer/scoring/signature.py`: PIF, CIF and ribosome drop-off, as
  pure transforms over a signal vector.
- **`cif`** in the `evidence` blob of every elongation row. Null when the span
  is not a whole number of codons.
- **`dropoff`** in the `evidence` blob of every termination row. It rides on
  termination rather than elongation because it needs positions past the ORF
  end, which is where the splice-aware flank machinery already lives. Null when
  the 33-nt window runs off the contig.
- `tests/reference_signature.py`, a deliberately unrefactored transliteration of
  the published R, used as an independent oracle. Never imported by product code.
- Evidence-blob agreement tests. `_assert_identical` compares only
  `metric`/`call`/`eligibility`, so everything in `evidence` was previously
  pinned by no test at all; the batched and scalar scorers are now diffed field
  by field across both fixtures.

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

### Fixed
- `_flank_positions` returned ascending genomic order on the minus strand,
  making the non-splice-aware fallback the exact reverse of the splice-aware
  path for every minus-strand flank. **No score changed** — the only consumer
  was `_codon_levels`, which bins in 3s and takes median/max/total, and every
  flank length in use is a multiple of 3, so the bin set was identical either
  way. Fixed because the next order-sensitive consumer would have read every
  minus-strand window backwards; `dropoff` is that consumer.

### Notes for consumers
- **PIF was already being emitted, under another name.** The existing
  `overall_in_frame` evidence field is exactly Chothani PIF
  (`sum(psites[seq(1,n,3)]) / sum(psites)`). No new field was added for it. The
  headline `elong_in_frame` is *not* PIF: it is the same ratio restricted to
  positions not contended by an overlapping event in a different frame. The two
  differ exactly where events overlap.
- **`term_drop` and `dropoff` both mean "signal falls after the stop" and are
  not interchangeable.** `term_drop` is an unbounded multi-flank consensus log2
  fold-change over all positions, flanks out to 60 nt. `dropoff` is a bounded
  [0,1] ratio over a fixed 33-nt window using frame-0 positions only. Neither
  supersedes the other; report which one you used.
- **CIF reproduces two quirks of the published R deliberately.** Its threshold
  is the literal `33.33`, not `100/3`, so a codon with exactly equal signal in
  all three frames counts as frame-0 dominant. Its denominator is *every* codon
  including zero-signal ones, so uncovered codons lower the score rather than
  being excluded. Both are pinned by tests; changing either would silently move
  every CIF value away from the published definition.
- **Drop-off assumes the translon span includes its stop codon**, since the
  window is frame-locked so index 17 is the terminal stop nucleotide. If an
  upstream annotation excludes the stop, every value is 3 nt out of register —
  worth checking per source. Verified against
  `translon_db/translons.sqlite` (8,852,481 rows): 99.9% have
  `terminal_codon_class = 'stop'` and 100.0% have `length_mod3 = 0`, so it
  holds there. `reference_cds` carries both conventions explicitly (`bed_*`
  stop-included, `stop_excluded_*`) plus a `reference_stop_policy` column.
- **Not yet validated against published Chothani values.** These agree with the
  R *formulas* (fuzz-tested against the transliteration); they have not been
  compared to numbers from the original pipeline on real data.

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

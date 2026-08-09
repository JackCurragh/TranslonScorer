# Changelog

All notable changes to TranslonScorer are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Because the score store is consumed by downstream annotation work, **the output
schema is treated as a public contract**. Every release states what a consumer
can rely on, and any change to it appears here first.

## [Unreleased]

## [0.3.0] — 2026-08-09

### The output contract

One row per `(event_id, aspect)`, with exactly **one headline metric per row**
carried in `metric` and named by `metric_name`. Secondary statistics live in
the `evidence` JSON blob. The four headline metrics are unchanged from 0.2.0:

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
- **CIF is emitted** in the `evidence` blob of every elongation row, as `cif`.
  Null when the span is not a whole number of codons.
- **Drop-off is emitted** in the `evidence` blob of every termination row, as
  `dropoff`. Null when the 33-nt window runs off the contig.
- `tests/reference_signature.py`, a deliberately unrefactored transliteration of
  the published R, used as an independent oracle in the test suite. Never
  imported by product code.
- Evidence-blob agreement tests. `_assert_identical` compares only
  `metric`/`call`/`eligibility`, so everything in `evidence` was previously
  unpinned by any test; the batched and scalar scorers are now diffed field by
  field across both fixtures.

### Notes for consumers

- **PIF was already being emitted, under another name.** The existing
  `overall_in_frame` evidence field is exactly Chothani PIF
  (`sum(psites[seq(1,n,3)]) / sum(psites)`). No new field was added for it. The
  headline `elong_in_frame` is *not* PIF: it is the same ratio restricted to
  positions not contended by an overlapping event in a different frame. The two
  differ exactly where events overlap.
- **`term_drop` and `dropoff` are both "signal falls after the stop" and are not
  interchangeable.** `term_drop` is a multi-flank consensus log2 fold-change
  over all positions, flanks out to 60 nt. `dropoff` is a bounded [0,1] ratio
  over a fixed 33-nt window using frame-0 positions only. Neither supersedes the
  other; report which one you used.
- **CIF reproduces two quirks of the published R deliberately.** Its threshold
  is the literal `33.33`, not `100/3`, so a codon with exactly equal signal in
  all three frames counts as frame-0 dominant. Its denominator is *every* codon
  including zero-signal ones, so uncovered codons lower the score rather than
  being excluded. Both are pinned by tests; changing either would silently move
  every CIF value away from the published definition.
- **Drop-off assumes the translon span includes its stop codon.** The window is
  frame-locked so index 17 is the terminal stop nucleotide, taken from the
  event's `term_pos`. If an upstream annotation *excludes* the stop codon, every
  `dropoff` value is 3 nt out of register — so this is worth checking per
  source. Verified against `translon_db/translons.sqlite` (8,852,481 rows):
  99.9% have `terminal_codon_class = 'stop'` and 100.0% have `length_mod3 = 0`,
  so the assumption holds there. `reference_cds` carries both conventions
  explicitly (`bed_*` stop-included, `stop_excluded_*`) alongside a
  `reference_stop_policy` column; use the stop-included pair.
- CIF returns null for spans that are not a whole number of codons. On the
  translon database above that is 21 rows out of 8.85M.
- Not yet validated against published Chothani values. The scores agree with the
  R *formulas* (fuzz-tested against the transliteration); they have not been
  compared to numbers from the original pipeline on real data.

### Fixed

- `_flank_positions` returned ascending genomic order on the minus strand, making
  the non-splice-aware fallback the exact reverse of the splice-aware path for
  every minus-strand flank. **No score changed** — the only consumer was
  `_codon_levels`, which bins in 3s and takes median/max/total, and every flank
  length in use is a multiple of 3, so the bin set was identical either way. It
  is fixed because the next order-sensitive consumer would have read every
  minus-strand window backwards.

## [0.2.0]

Prior releases predate this changelog.

[Unreleased]: https://github.com/JackCurragh/TranslonScorer/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/JackCurragh/TranslonScorer/compare/v0.2.0...v0.3.0

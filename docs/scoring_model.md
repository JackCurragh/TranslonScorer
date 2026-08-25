# Scoring model

## Unit

The scored unit is the **event**: a deduplicated piece of genomic geometry that
many features can claim. Not the translon.

| event | geometry | identity key |
|---|---|---|
| `elongation` | one translation block, with phase | `E\|chrom\|strand\|start\|end\|phase` |
| `junction` | one (donor, acceptor) intron | `J\|chrom\|strand\|donor\|acceptor` |
| `init` | one genomic position (translon 5' end) | `I\|chrom\|strand\|pos` |
| `term` | one genomic position (translon 3' end) | `T\|chrom\|strand\|pos` |

Identity is a deterministic hash of the key, so the same geometry from any feature
source produces the same `event_id`. `feature_event` maps features to the events
they claim; `event_overlap` records elongation events that overlap in a different
frame.

Scoring per event rather than per translon is the reason the engine is tractable
at annotation scale, and it is not negotiable for performance reasons alone.

## What a block-level score is

A per-block score is a statement **about that block**. It is not an estimator of
any translon-level quantity, and must not be presented as one.

This has a specific consequence for the Chothani/Sonia scores. Block-level CIF is
a valid measurement of a block. Sonia CIF is defined over the complete
transcript-oriented ORF vector. These are different objects; averaging the first
does not produce the second, and the result must not be labelled Sonia CIF.

Belief in a translon is an assessment over the **set of blocks** it claims,
derived from per-block evidence — not a re-measurement at translon scale.

## Measurement and decision are separate

Every event-aspect produces one long-form row:

```
event_id, aspect, group, tier, n_reads, metric, metric_name,
eligibility, call, evidence(json), thresholds_version,
map_track_mean, map_track_low
```

- `eligibility ∈ {NA, INSUFFICIENT, ELIGIBLE}` — could this be tested.
- `call ∈ {SUPPORTED, UNSUPPORTED, AMBIGUOUS}` — did it pass.

These are never conflated: "couldn't test" is not "failed". `NA` is missing
evidence — never zero, never a negative score.

The call is a thin, re-derivable function of the continuous evidence plus an
explicit versioned `ScoreThresholds`. Change a threshold and re-derive calls;
never re-measure.

`metric` is the natural per-aspect quantity. It is **not** normalised across
aspects, and there is no composite score.

## Current metrics

| aspect | metric | evidence carried alongside |
|---|---|---|
| `init` | `init_rise` — log2((body+α)/(leader+α)), median codon-bin levels, median over flanks 9/18/30/60 | `rise_by_flank`, `stability`, `flank_peakiness`, `flank_spliced`, `boundary_axes_by_flank` |
| `elongation` | `elong_in_frame` — in-frame fraction where unconfounded | `clean_in_frame`, `overall_in_frame` (= Sonia PIF), `cif` (first-class column), `n_codons` (first-class column), `breadth`, `identifiability`, `competitor_share`, `noise_share` |
| `term` | `term_drop` — log2((body+α)/(UTR+α)) | same step machinery as init, plus `dropoff` (Sonia, bounded), `boundary_axes_by_flank` |
| `junction` | `junc_confident_spanning` — reads with ≥6 nt aligned both sides | `short_spanning`, `unspliced`, `psi` |

Review flags — `stability`, `flank_peakiness`, `map_track_mean`, `map_track_low` —
are recorded as evidence and never gate a call.

## Boundary evidence axes (init/term)

The depth-only step (`consensus_rise`) can't distinguish a real start/stop
from anything else that changes density there. `boundary_axes_by_flank` adds,
per flank length, both sides of the boundary:

- **periodicity** — frame-0 share of the flank's signal (PIF's definition,
  applied to a flank instead of a whole ORF span).
- **uniformity** — peakiness (max/median codon-bin) and Gini coefficient
  (`signature.gini`) over codon bins. Gini reacts to the whole distribution,
  not just the tallest bin.
- **breadth** — fraction of codon bins with any signal.
- **periodicity significance** — one-sided Mann-Whitney U p-value comparing
  per-codon frame-0 share between the two flanks (`aspects._periodicity_significance`);
  `None` below a 5-codon floor per side, or if scipy is unavailable — scipy
  stays optional, this axis just goes missing rather than crashing scoring.
  Inspired by RiboCode's use of a nonparametric periodicity test, not a
  reproduction of it.

All of this is evidence, with one narrow exception: in the existing
borderline band (`0 < consensus_rise < rise_thr`), a significant periodicity
result at enough flank lengths (`periodicity_min_agree_frac` of the available
`periodicity_p` values below `periodicity_significance_alpha`) resolves the
call from AMBIGUOUS to SUPPORTED — auditable via
`evidence["periodicity_resolved_ambiguous"]`. This can only move
AMBIGUOUS → SUPPORTED; an event already SUPPORTED or UNSUPPORTED by depth
alone is untouched.

## Frame and strand

`abs_frame(phase, strand)` gives the genomic `pos % 3` residue carrying in-frame
signal. Codons are three consecutive *transcript* positions, so on the minus
strand a codon's first base sits at the highest genomic coordinate; vectors are
reversed into transcript order before frame indexing.

Initiation flanks are anchored at the start codon's **first** base. Termination
flanks are anchored at the terminal nucleotide — the stop codon's **last** base.
The two anchors are offset by two positions relative to codon structure, so
frame-0 does not sit at the same index in both.

Leader and UTR flanks are projected across introns via `splice_context` rather
than read as raw flanking genomic bases.

## Composition

`feature_event ⨝ scores` → per-translon, per-aspect evidence. Elongation's
`metric` is read-count-weighted across its segment events (not length —
`compose_report`/`_compose_per_translon`); a feature's `{aspect}_call` is
SUPPORTED when at least `supported_frac_min` (default 0.5) of its events are
individually SUPPORTED. Composition reports evidence; it does not make
annotation decisions.

CIF composes differently, deliberately: `elongation_cif_approx` is
codon-count-weighted, not read-weighted, because CIF's denominator is
codons — read-weighting would let a short, deep block outvote a long,
modest-coverage one. It is still only an approximation (labelled as such):
the real translon-level CIF needs codons concatenated across blocks in
transcript order, crossing the intron, which is the deferred isoform work.
`report.compose_block_detail` returns the un-collapsed per-block CIF values
this approximates, one row per block, nothing lost.

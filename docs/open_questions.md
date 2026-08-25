# Open questions

Everything raised and not yet settled. Nothing here is agreed. Items move out of
this file into one of the three settled docs when decided.

---

## 1. Established observations, not yet acted on

Facts checked against code or data. They constrain the decisions below.

- **Codons straddling block boundaries are dropped.** `orf_signal_vector` trims
  partial codons at each block end, losing up to one codon per boundary. Currently
  invisible. Still open — this is exactly what the real (non-approximate)
  translon-level CIF in §5 would need to fix, by concatenating codons across
  blocks instead of trimming at each one.
- **Drop-off denominator dies on 14.22% of rows.** `zero_dropoff_denominator` in
  the pancreas run (83,354 / 586,105). The bounded [0,1] ratio is degenerate at
  that rate.
- **Two scorers, different estimators.** Pancreas start-rise uses flank 45, mean,
  α=0.1. TranslonScorer uses median-of-codon-bins, flanks 9/18/30/60, α=1.0.
  Definitions are reconciled; empirics are not transferable.
- **Provider incomparability has a calibration consequence.** Since offsets enter
  at three different points, any threshold or reference distribution must be built
  per provider, on the same track.

## 2. Boundary metrics

- **Frame-aware start rise / stop drop — done, in a different shape than
  proposed.** Not a frame-0-only log2FC step; `boundary_axes_by_flank` adds
  periodicity (frame-0 share), uniformity (peakiness + Gini), breadth, and a
  Mann-Whitney significance test per flank length, both sides — see
  `scoring_model.md`. `init_rise`/`term_drop` themselves are unchanged
  (still total-read); the frame-resolved read is additional evidence, and
  can resolve a borderline call (see `_decide_step`). Not yet established
  against reviewed examples beyond the one synthetic case in
  `test_splice_aware_flank.py` that motivated the AMBIGUOUS-resolver in the
  first place — real-data calibration is still open.
- **Short-window noise floor — done.** `_periodicity_significance` returns
  `None` (not a number) below a 5-codon-per-side floor
  (`ScoreThresholds.periodicity_min_codons`).
- **Multiscale windows.** Still 9/18/30/60 only — the 2/3/5/10-codon variant
  was not built.
- **What to derive from the multiscale vector.** Partially addressed: every
  new axis is computed and kept per flank length (nothing collapsed early),
  so agreement/disagreement across scales is visible. A specific
  short-range/long-range/sign-consistency summary statistic is still not
  derived — the raw per-flank-length picture is what's available today.

## 3. Mappability-adjusted CIF

- Usable codon = all three bases map to the exon model **and** meet a mappability
  requirement. Zero-signal-but-mappable codons stay in the denominator — removing
  them deletes the negative evidence.
- Must always travel with `raw CIF`, `usable_codon_count`, `usable_codon_fraction`.
  `CIF_mappable = 1.0` over 3 codons is not `0.8` over 80.
- Needs a minimum denominator (proposed: ≥10 usable codons, or ≥50% of ORF codons)
  and stratification by ORF length before it could be a headline metric.
- Requires promoting fields from evidence JSON to schema columns.
- Thresholds must be re-derived under the adjusted definition; current CIF
  thresholds do not carry over.

## 4. Competing starts sharing a stop

- Competition is real and holds under hypothesis testing. Supplying the candidate
  set is what makes competitors **enumerable rather than latent** — it bounds the
  problem rather than removing it.
- `identifiability` / `competitor_share` / `noise_share` already do this for
  overlapping frames. Candidates sharing a stop is the same mechanism at wider
  scope.
- **Grouping key is not decided.** `(chrom, strand, stop, frame)` is insufficient.
  Also relevant: transcript/exon compatibility, splice path, exact stop geometry,
  whether candidate bodies actually overlap, alternative stops, and duplicate
  geometry submitted by different tools.
- Measuring a start's downstream body truncated before the next competing start is
  the practical intermediate.
- **Piecewise step model** (background + step per start − step at stop, with sparse
  regularisation) is not the next production metric. The cheap precursor is a
  diagnostic: per shared-stop group, plot codon-level frame-0 signal, overlay all
  candidate starts, show each start-rise and the distance to its neighbours. That
  tests whether the 45-nt start-rise is being contaminated without committing to a
  model.
- Proposed evidence ordering: raw local rise < frame-specific local rise <
  frame-specific rise sustained across codons < sustained body plus consistent stop
  drop with competing starts considered.

## 5. Sonia scores at feature level

- Sonia PIF/CIF/drop-off are defined over the complete transcript-oriented ORF
  vector. Block-level values are valid block measurements but are not those.
- A feature-level layer would need explicit fields rather than evidence JSON, and
  requires transcript projection — which TranslonScorer does not do, and which
  touches the deferred isoform question. The pancreas recompute enumerates all
  compatible transcripts without selecting one as truth. **Partial progress for
  CIF**: `cif`/`n_codons` are now first-class columns and
  `elongation_cif_approx` gives a codon-weighted translon-level number — but
  it's explicitly an approximation (codons aren't concatenated across the
  intron), not the real transcript-projected value this bullet is about.
  `report.compose_block_detail` keeps the un-collapsed per-block values.
- **Absence taxonomy.** `sonia_unique_status` distinguishes six outcomes
  (`ok`, `no_transcript_exon_mapping`, `zero_total_signal`,
  `zero_dropoff_denominator`, `insufficient_dropoff_flank`,
  `stop_not_overlapping_transcript`). TranslonScorer's `eligibility` has three and
  collapses all of these into `INSUFFICIENT`. Port direction is pancreas →
  TranslonScorer.

## 6. CDS and non-CDS assessment

- Run the same scorer over canonical CDS translons and non-CDS candidates, same
  track, one row per translon, with `feature_class` provenance and source tool.
- Do not calibrate CDS and non-CDS separately.
- Keep native TranslonScorer metrics alongside Sonia, to test whether the
  event-centric scores add information or reproduce the Sonia signal.
- Partly built already: CDS-control pif/cif/dropoff distributions
  (`transcode_pooled_baseline.py`) and Sonia-vs-native correlations
  (`build_pancreas_score_correlations.py`).

## 7. Thresholds

- Position each facet against a **reference distribution** rather than an absolute
  cutoff: CDS first, then high-confidence non-CDS as a second reference population
  so short and non-canonical candidates have a reference that is not structurally
  unlike them.
- Reference distributions are per provider and per track.

## 8. Undecided mechanisms

- **Block → translon.** Belief in a translon is an assessment over the set of
  blocks it claims. The rule for that is not decided. It cannot reduce to variance:
  real ORFs are not uniform (5' ramp, pausing, per-block depth).
- **Scorable vs scored.** Which facets a candidate must be scorable on to be
  assessed at all, versus merely scored where available. This determines what a
  reference distribution can be built from, and how features with few scorable
  facets are handled.
- **Output schema.** Which quantities are first-class columns and which stay in
  evidence JSON.

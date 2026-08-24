# Significance testing plan

Not yet implemented. This is a forward-looking plan, written the same way the
boundary-axes work was planned before it was built — captured here so it
survives past one conversation, not because anything below is agreed in the
same sense as `annotation_decisions.md`.

## Context

The init/term boundary axes added this session include a periodicity
significance test (`aspects._periodicity_significance`, one-sided
Mann-Whitney U, evidence-only except as an AMBIGUOUS-resolver). Elongation,
CIF, and termination dropoff still rely on fixed-ratio thresholds computed
straight from raw counts, which are unreliable exactly where depth is
lowest — a 2-read codon can trivially clear CIF's 33.33% cutoff by chance.
This doc proposes extending the same "replace a bare threshold with a
data-size-aware test" idea to the remaining aspects, plus feeding the result
into the feature-level believability hierarchy in `consequential.py`.

## 1. Elongation: frame test

- **Current**: `elong_in_frame >= thr.elong_in_frame` (fixed 0.5) and
  `breadth >= thr.elong_breadth`.
- **Proposal**: per elongation event, test the observed frame-0 read count
  against the null of a uniform 1/3 split, using the `cont_by_frame`/frame
  tallies `_elong_evidence` already receives — no new coverage walk needed.
  - Chi-square goodness-of-fit across all 3 frame counts vs uniform expected
    (`scipy.stats.chisquare`) is the better fit over a plain binomial test:
    it uses the full 3-way split already computed, and catches signal that's
    unevenly split between the two *non*-frame-0 positions, which a
    frame-0-vs-rest binomial test can't see.
  - Evidence-only field (`elong_frame_p`) in v1, not a call-gate — same
    "resolver, not a blind switch" pattern as the boundary axes. Natural fit:
    resolve the existing `identifiability`-driven AMBIGUOUS band, not replace
    the `elong_in_frame` threshold outright.
- **Open decision**: resolver-only, or eventually replaces the fixed
  threshold? Flagged, not answered here.

## 2. CIF: per-codon significance

- **Current**: `cif()` in `signature.py` is a per-codon deterministic
  threshold (frame-0 share > 33.33%), deliberately reproducing the Chothani
  R reference exactly — including its low-depth noise problem.
- **Proposal**: a **separate** new metric, not a change to `cif()` (which
  must stay reference-faithful) — e.g. `cif_significant(signal,
  min_codon_reads) -> Optional[float]`: for each codon with enough reads, a
  binomial test (frame-0 count vs total, p=1/3); the metric is the fraction
  of *testable* codons where that test is significant. This is CIF's
  low-depth problem addressed without touching the function whose whole job
  is matching a published reference.
- **Open decision**: minimum reads-per-codon floor (the per-codon analogue
  of `periodicity_min_codons`, which is per-flank).

## 3. Termination dropoff: significance

- **Current**: `dropoff()` — bounded [0,1] ratio, before/(before+after),
  `None` at a zero denominator (14.22% of rows in the pancreas run, per
  `open_questions.md`).
- **Proposal**: reuse `_periodicity_significance`'s shape (or the function
  itself) on the same 33-nt window `_dropoff_at` already builds — before-stop
  vs after-stop frame-0 counts, one test instead of a new implementation.
  Cheapest item on this list.

## 4. Junction: deferred, explicitly

- No natural uniform null the way frame has 1/3 — a meaningful test needs a
  reference/background spanning-rate distribution, which is the same
  deferred work as `open_questions.md` §7 (reference-distribution
  thresholds). Named here so it isn't silently dropped, not attempted
  half-built against a null that doesn't mean anything yet.

## 5. Feeding the believability hierarchy (`consequential.py`)

- `apply_policy()`'s `tier_confidence` currently only sees binary
  `{aspect}_call` — a translon `SUPPORTED` by the periodicity-resolver at one
  flank length ranks identically to one `SUPPORTED` cleanly at every scale.
- **Proposal**: a per-event `confidence` field (e.g. count of significance
  tests that passed, out of how many were computable for that event),
  promoted through composition the same way `cif`/`n_codons` were promoted
  this session, then blended into `tier_confidence` as a weight instead of a
  flat mean-of-SUPPORTED.
- This is the piece that actually answers "is this event believable" beyond
  a binary call. Deliberately sequenced last — the composition-weighting
  question here is exactly the kind of thing the CIF read-vs-codon-weighting
  discussion this session got right by not rushing, and it depends on 1–3
  actually landing first so there's real per-event confidence data to weight.

## Sequencing

1. Termination dropoff significance — reuses `_periodicity_significance`
   almost verbatim, cheapest.
2. Elongation frame test — reuses existing frame-tally machinery.
3. CIF per-codon significance — new function, needs its own spec for the
   reads-per-codon floor.
4. Hierarchy integration (§5) — depends on 1–3 landing.
5. Junction — revisit only once the reference-distribution work
   (`open_questions.md` §7) is scoped; not before.

## Open decisions needing sign-off before implementation

- Elongation: resolver-only vs. eventually replacing the fixed threshold.
- CIF: the per-codon reads floor value.
- Whether `cif_significant` is a new metric alongside `cif()` (recommended)
  or changes CIF's definition (not recommended — breaks reference fidelity).

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

## 0. Design principle: one question, several convergent lenses

A plain p-value only says "this doesn't look random." It says nothing about
*whose* signal it is. Nearby real biology — an overlapping ORF, a uORF,
genuine stop-codon readthrough, a duplicated locus — can make a window look
non-random for reasons that have nothing to do with the event under test. No
choice of null distribution fixes this: the ambiguity isn't in the null, it's
in the alternative hypothesis, which a single test was never built to
resolve.

The fix isn't a smarter single test — it's asking the *same underlying
question* through several measurements that would diverge if the signal
came from somewhere else, and requiring them to agree. Each aspect below is
organized the same way:

- **Core question** — the one thing that phase of translation actually
  asks, stated in plain terms.
- **Lenses** — 2–4 measurements that each answer that same question from a
  different angle (magnitude, spatial distribution/uniformity, read-length
  consistency). A real event moves all of them together; an artifact tends
  to satisfy only the one lens its own quirk happens to match.
- **Attribution, not exclusion** — where a neighboring frame or region also
  shows signal, the goal is *not* to decide a single winner. Overlapping
  translation is real biology (dual-coding regions, uORFs with independent
  starts, genuine readthrough) and is exactly the kind of thing this tool
  exists to capture, not launder into a single-frame answer. Each candidate
  event's battery is evaluated on its own terms; two events can both come
  back positive, and both get reported.

This generalizes something already in the codebase in miniature: the
boundary periodicity test is already run at several flank lengths and
requires agreement across most of them before it resolves an AMBIGUOUS call
(`_decide_step`, `evidence.py`). That's a one-axis version (flank length) of
the same idea — the plan below adds more axes and applies it per aspect.

## 1. Initiation

**Core question**: does in-frame density increase downstream of this
candidate start, and hold?

Three lenses, each already represented by a piece of the codebase rather
than invented fresh:

- **Level** — the existing coverage step (`_step_score`: body level vs
  outer level, log2 fold-change). Did density go up at all.
- **Phase** — the existing boundary periodicity test
  (`_periodicity_significance`): is the increase specifically *in-frame*,
  not just more reads generally.
- **Uniformity** — not yet wired to the boundary, but `signature.gini()`
  already exists. A genuine start should settle into steady, roughly
  uniform in-frame elongation downstream; a stray pileup tends to stay
  peaky/patchy. Reusing `gini` on the downstream in-frame codon series adds
  this as a third leg at no new-metric cost.

  **Implementation note (added after landing, see
  docs/significance_testing_results.md):** "uniformity" turned out to bundle
  two genuinely different questions, not one. Gini measures MAGNITUDE
  CONCENTRATION — is the mass piled into one bin, regardless of where. It
  says nothing about whether the in-frame SHARE itself holds steady across
  the window: broadly-spread signal that is nonetheless much more in-frame
  near the boundary than further out has unremarkable Gini but is not
  spatially uniform in the sense this lens is meant to catch. So the
  initiation battery carries both: `gini_body_inframe` (magnitude) and
  `body_share_consistency_p` (spatial frame-0-SHARE consistency between the
  near and far halves of the body flank, via the same `_share_consistency`
  two-sided test elongation and termination use — see below). Elongation's
  §2 "Uniformity" (`body_uniformity_p`) and termination's §4 "Continuity"
  (`dropoff_continuity_p`) were ALREADY the share-consistency kind, not the
  Gini kind — they were built first and reused for initiation once the
  distinction became visible by comparing the three lenses side by side.
  A future aspect's "uniformity" lens should pick deliberately between the
  two rather than assuming they're interchangeable.

These aren't three independent confound checks — they're three
measurements of one claim ("density increased and stabilized here"), which
is why they're expected to move together for a real start.

**Attribution — neighboring/overlapping starts**: when a candidate or
annotated start sits nearby in a different frame, don't pick a winner.
Run the same three-lens battery independently in each frame. Three
distinguishable outcomes fall out:

1. Only this frame shows the increase → ordinary single-frame initiation;
   the neighboring signal was noise or unrelated.
2. Both frames independently pass their own battery → genuine overlapping
   translation; report both events.
3. This frame's apparent increase disappears once measured on its own
   terms (i.e. it was really the tail of the neighbor's signal) → leakage,
   not a real start here.

This only works if each frame's battery is computed on its own terms rather
than as a head-to-head "frame A vs frame B" comparison — a head-to-head
framing structurally can't return outcome 2.

## 2. Elongation

**Core question**: given translation has started, does in-frame density
stay elevated and consistent across the whole body?

- **Level** — chi-square goodness-of-fit of the 3-way frame tally
  (`cont_by_frame`) against uniform 1/3, reusing `_elong_evidence`'s
  existing tallies. Preferred over a frame-0-vs-rest binomial test because
  it also catches signal unevenly split between the two non-frame-0
  positions.
- **Uniformity** — split the body into halves (or a sliding window) and
  test whether frame-0 dominance is consistent along its length. A
  confound that only overlaps part of the ORF shows up as a local patch,
  not a body-wide signature; genuine elongation should look similar
  throughout. (This is the SPATIAL SHARE-CONSISTENCY kind of "uniformity",
  not the Gini/magnitude-concentration kind — see the implementation note
  under §1 Initiation.)
- **Cross-check against existing region flags** — the codebase already
  tracks a low-mappability flag per locus (`map_track_low`,
  `coverage/base.py`). Wire it into the elongation confidence explicitly:
  a "significant" frame call sitting on a flagged low-mappability region
  should count for less, since multi-mapped reads are a known way to
  manufacture apparent periodicity that doesn't belong to this locus.

**Attribution — overlapping ORFs in the body**: same per-frame-on-its-own-
terms approach as initiation. If a different frame also shows a
significant, uniform signature over some or all of the body, that's
consistent with real dual-coding/overlapping translation, not evidence
against this frame — report both if both pass their own battery.

**Open decision**: resolver-only (feeds the AMBIGUOUS band, same as the
boundary axes) vs. eventually replacing the fixed `elong_in_frame`
threshold outright. Flagged, not answered here.

## 3. CIF: per-codon significance

**Core question**: how much of the ORF's codon-by-codon signal is
individually frame-0-dominant beyond what depth alone would produce by
chance?

- **Level** — per codon with enough reads, a binomial test (frame-0 count
  vs total, p=1/3); the metric is the fraction of *testable* codons where
  that test is significant. Kept as a **separate** metric from `cif()`,
  which must stay reference-faithful to the Chothani R implementation —
  this doesn't change that function, it adds a low-depth-aware companion.
- **Uniformity** — are the significant codons spread across the ORF, or
  clustered in one stretch? Genuine elongation should give broad, scattered
  significance; a local confound (e.g. a short overlapping element)
  clusters it. Testable with a simple contiguity/run-length statistic over
  the per-codon significance calls.
- **Cross-aspect agreement** — does this line up with the elongation
  chi-square test over the same window? Two independently-computed aspects
  agreeing is stronger evidence than either alone, and cheap since both
  already read the same underlying per-codon frame tallies.

**Attribution**: if a stretch of codons is significantly dominant in a
*different* frame, don't treat that as noise to explain away — check
whether that stretch's own frame passes the same per-codon battery on its
own terms. A consistent block of off-frame significance is a candidate for
a real local overlapping element, not necessarily contamination of this
ORF's CIF.

**Open decision**: the minimum reads-per-codon floor (the per-codon
analogue of `periodicity_min_codons`, which is per-flank, not per-codon).

## 4. Termination dropoff

**Core question**: does in-frame density fall after the stop, relative to
before — the mirror image of initiation's question.

- **Level** — reuse `_periodicity_significance`'s shape (or the function
  itself) on the same 33-nt window `_dropoff_at` already builds: before-stop
  vs after-stop frame-0 counts. Cheapest item on this list, and the most
  direct extension of an existing test.
- **Continuity (uniformity, mirrored)** — does the "before" window's
  profile look like a continuation of ordinary upstream elongation, or is
  it an isolated pileup right at the stop with no elongation signature
  behind it? A real termination event should have the elongation
  signature (§2) holding just upstream; a stop-proximal queuing/stalling
  artifact often doesn't. (Also the share-consistency kind, not Gini — see
  §1's implementation note.)
- **Downstream attribution** — does the "after" window overlap a separate
  annotated/candidate ORF (a dORF)?

**Attribution — why this one has three outcomes, not two**: "signal after
the stop" is not automatically an artifact. Genuine stop-codon readthrough
is real, documented translational biology, not noise. So the after-window
result should be read as:

1. Clean drop, no sustained after-signal → ordinary termination.
2. Sustained *in-frame* signal after the stop, continuous with before →
   readthrough — a real event in its own right, worth flagging/reporting,
   not just a low dropoff score.
3. After-signal exists but doesn't connect to this event (different frame,
   or attributable to a separate downstream ORF's own start) → a distinct
   translation event nearby, not readthrough of this one.

This is the same three-way shape as initiation's neighbor-start check,
applied to the downstream side of termination instead of the upstream side
of initiation.

## 5. Junction

**Core question**: does the reading frame stay continuous across the
splice junction, consistent with one continuous translated ORF — or does
each side behave independently?

Previously deferred outright for lacking a natural uniform null (no
equivalent of frame's 1/3 for spanning rate). That's still true, and a
reference-distribution approach is still blocked on `open_questions.md`
§7. But the same convergent-lenses idea doesn't need an external null — it
needs **internal consistency**, which is checkable now:

- **Frame match** — is frame-0 share comparable immediately upstream and
  immediately downstream of the junction? A single continuous ORF should
  look similar on both sides.
- **Sensitivity to the junction reads themselves** — does the frame call
  change materially if junction-spanning reads are excluded and only
  exonic reads on either side are used? If the call only holds up *because*
  of the spanning reads, that's weaker evidence than a call that's robust
  either way.

**Attribution**: if the two sides don't match, that isn't automatically a
bad junction call — it can mean each exon is independently translated
(e.g. alternative frame usage after splicing, or an exon translated on its
own). Same principle as everywhere else: test each side on its own terms
before concluding the junction is spurious rather than reporting two
distinct segments.

This is a smaller, nearer-term step than the full reference-distribution
work, and doesn't need to wait on it.

## 6. Feeding the believability hierarchy (`consequential.py`)

- `apply_policy()`'s `tier_confidence` currently only sees binary
  `{aspect}_call` — a translon `SUPPORTED` by the periodicity-resolver at
  one flank length ranks identically to one `SUPPORTED` cleanly at every
  scale.
- **Proposal**: a per-event `confidence` field built from the lens
  batteries above — not "count of one test type passed," but "how many of
  the *distinct, convergent* lenses for this aspect agreed," promoted
  through composition the same way `cif`/`n_codons` were promoted this
  session, then blended into `tier_confidence` as a weight instead of a
  flat mean-of-SUPPORTED.
- Where attribution produced more than one reportable event at a locus
  (overlapping initiation, readthrough, an independently-translated exon),
  each gets its own confidence from its own battery — this hierarchy
  should not collapse them back into one answer.
- Deliberately sequenced last — depends on §1–5 actually landing first so
  there's real per-event, per-lens data to weight.

## Sequencing

1. Termination dropoff level test — reuses `_periodicity_significance`
   almost verbatim, cheapest.
2. Elongation level test (chi-square) — reuses existing frame-tally
   machinery.
3. Initiation uniformity lens (`gini` at the boundary) — reuses an
   existing function, no new statistics.
4. CIF per-codon significance — new function, needs its own spec for the
   reads-per-codon floor.
5. Uniformity/continuity lenses for elongation and termination, and the
   attribution logic (per-frame-on-its-own-terms batteries) for
   initiation, elongation, CIF, and termination — this is the bulk of the
   new work, and is what actually delivers "capture overlapping
   translation" rather than just "reject noise."
6. Hierarchy integration (§6) — depends on 1–5 landing.
7. Junction internal-consistency checks — independent of the
   reference-distribution work, can proceed in parallel once §5's
   attribution pattern exists to reuse.

## Open decisions needing sign-off before implementation

- Elongation: resolver-only vs. eventually replacing the fixed threshold.
- CIF: the per-codon reads floor value.
- Whether `cif_significant` is a new metric alongside `cif()` (recommended)
  or changes CIF's definition (not recommended — breaks reference
  fidelity).
- Per-event, per-read-length frame tallies don't currently exist (only
  sample/cohort-level per-length QC does, in `matrix/qc.py`). A read-length-
  consistency lens for any aspect needs this instrumented first — bigger
  lift than the rest of this plan, not included in the batteries above
  until scoped separately.
- How many lenses need to agree, per aspect, before a battery counts as
  "passed" for attribution purposes (e.g. 2 of 3, all of them) — not fixed
  here, needs a value per aspect once real data exists to tune against.

## Addendum: what §1–§7 landing and real-data validation changed

Everything above is the original plan, left as written. This section
records what building it and running it against the real GAPDH cohort BAM
(see docs/significance_testing_results.md for full numbers) actually
surfaced — three things the plan didn't anticipate.

1. **Initiation's "uniformity" bundled two different questions** (see the
   implementation note under §1) — Gini (magnitude concentration) and
   spatial frame-0-share consistency are not the same test, and elongation/
   termination's "uniformity"/"continuity" lenses were already the
   share-consistency kind. Fixed by giving initiation both.

2. **3+-way overlapping candidate frames are not a rare edge case in real
   annotation.** On the GAPDH region fixture, 150 of 336 elongation events
   (45%) had two or more competing candidate frames once `event_overlap`
   was wired in (2 to 8 competitors, most commonly 2–3) — a gene-dense
   region, but not a contrived one. §0/§6's "attribution not exclusion"
   principle was written with the two-frame case foremost in mind; at this
   prevalence, the N>2 case needs its own real design pass at some point,
   not indefinite deferral. This implementation deliberately does not
   attempt one (per its own brief) — `pair_elongation_neighbors`/
   `pair_cif_neighbors`/`pair_boundary_neighbors` flag `multi_way_overlap`
   and report each frame's own battery result independently, nothing more.

3. **`event_overlap` (the table `_contention` already builds) is not wired
   into the production `score-bams` path.** `score_bams_workflow` /
   `_score_events_over_provider` (workflows.py) never read it and never
   pass an `overlaps_df` into `scoring.run.score_events`, so every
   elongation event scored via the standard CLI path today has
   `competitor_share={}` regardless of real overlaps on disk — contention
   scoring and this session's neighbor attribution only actually run when a
   caller builds `overlaps_df` by hand (as this session's validation script
   and `tests/test_golden.py`'s contended-pair tests do). This predates this
   session's work — the elongation contention machinery itself was already
   in place — but it means the two are silently disconnected in the shipped
   CLI, and finding #2 above (which needed `overlaps_df` to fire at all)
   would otherwise have been invisible. Wiring `event_overlap` (plus a
   `comp_phase` join back onto `events`, since the on-disk table doesn't
   carry it either -- see `run_gapdh_validation.py`'s manual join for
   the shape needed) into `score_bams_workflow`/`score_matrix_workflow` is
   flagged as follow-up, not attempted here.

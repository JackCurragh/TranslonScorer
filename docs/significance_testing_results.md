# Significance testing: implementation results

Companion to `docs/significance_testing_plan.md` (§0–§7, "Sequencing",
"Open decisions"). That doc is the spec; this one records what landed,
what real data showed, and what to revisit.

## What landed

Every lens and cross-check in the plan's §1–§7 is implemented and unit
tested, landed as one commit per aspect on `feat/significance-testing`:

| Aspect | Lenses | Where |
|---|---|---|
| Termination | Level (`dropoff_significance_p`), Continuity (`dropoff_continuity_p`), Downstream/readthrough (`dropoff_after_share` + `classify_termination_downstream`) | `scoring/aspects.py`, `scoring/attribution.py` |
| Elongation | Level (`frame_chisq_p`), Uniformity (`body_uniformity_p`), Mappability downweight | `scoring/aspects.py`, `scoring/attribution.py` |
| Initiation | Level (`consensus_rise`, pre-existing), Phase (`periodicity_p`, pre-existing), Uniformity — magnitude (`gini_body_inframe`) **and** spatial (`body_share_consistency_p`, added after a design-gap fix — see plan §1's implementation note) | `scoring/aspects.py` |
| CIF | Level (`cif_codon_significance`), Uniformity/contiguity (`cif_codon_contiguity`), Cross-aspect agreement with elongation's chi-square | `scoring/signature.py`, `scoring/attribution.py` |
| Junction | Frame-match, spanning-sensitivity proxy | `scoring/aspects.py::junction_internal_consistency` |
| Attribution | Per-aspect batteries (`init_battery`/`term_battery`/`elong_battery`/`cif_battery`), 3-outcome neighbor classification (`classify_neighbor_outcome`), pairing over elongation/CIF overlaps and a proximity heuristic for init/term (`pair_elongation_neighbors`/`pair_cif_neighbors`/`pair_boundary_neighbors`) | `scoring/attribution.py` |
| Hierarchy (§6) | `confidence` — first-class schema column, read-weighted through `compose_report`, blended into `consequential.py`'s `tier_confidence` as a per-aspect weight | `scoring/evidence.py`, `report.py`, `consequential.py` |

Every new lens is **evidence-only** — none of them changed `_decide_step`
or `_elong_evidence`'s eligibility/call logic. That was a deliberate
constraint (see CODE VOICE / scope), and it means one thing up front:

**No per-event `call` (SUPPORTED/AMBIGUOUS/UNSUPPORTED) changed on GAPDH
with the new layer ON vs OFF, by design.** The pre-existing
`periodicity_resolved_ambiguous` mechanism in `_decide_step` (init/term
AMBIGUOUS→SUPPORTED via the boundary periodicity test) predates this
session and is unaffected. What the new layer changes is `confidence`
(new field, feeds `tier_confidence`) and attribution labels (new,
previously unrepresentable). Those are covered below instead of a
call-count table that would just read "0 changed everywhere."

## Real-data validation: GAPDH cohort BAM

`data/gapdh_cohort_genome.bam` (24,197 alignments, chr12) scored against
`data/chr12_gapdh_region.gtf` (chrom renamed `12`→`chr12` to match the
BAM; CDS feature type) via `extract_events_workflow` →
`score_bams_workflow` → `report_workflow` → `consequential_workflow` (the
real, shipped workflow functions, not synthetic fixtures).

- 636 events extracted: 40 init, 336 elongation, 62 term, 198 junction,
  across 116 translons (this GTF window covers GAPDH plus several
  neighboring genes, most with zero reads in this cohort-subsetted BAM).
- 24,187 of the 24,197 total alignments overlap the GAPDH gene span itself
  (chr12:6,534,500–6,538,400); real depth in this fixture is
  GAPDH-specific, not spread across the window.
- The real GAPDH transcript, `ENST00000229239`, has 8,893 total scored
  reads and is the only translon in the fixture with substantial signal
  end to end.

### Confidence distribution

51 of 636 events got a non-null `confidence` (585 are null — see the
floor-firing findings below for why). Distribution: min 0.0, max 1.0,
elongation mean 0.47 (n=11 with any reads), init mean 0.0375 (n=40, but
39 of them are zero-read candidates).

**High-confidence example** — elongation event `8524976404572690219`
(part of `ENST00000229239`, span 6537583–6537996, 686 reads): confidence
1.0. Both lenses agree cleanly: `frame_chisq_p` is essentially zero and
`body_uniformity_p` is well above 0.05 (consistent across both halves of
the body). `metric` (0.448) is just *below* `elong_in_frame` (0.50), so
the call is still UNSUPPORTED under the fixed-ratio threshold — but the
lens battery shows this is not noise, it's periodicity that the ratio
threshold alone doesn't capture. (A second, nearby block in the same
translon, `14065136805926077183`, 3,741 reads, is the same story in more
detail: `frame_chisq_p ≈ 4×10⁻⁵²`, `cif_significance.frac_significant =
0.44`, `cif_contiguity.max_run_frac = 0.125` — broadly scattered
significant codons, not a local artifact — pooled confidence 0.6.)

**Moderate-confidence example** — init event `4298823927098781090`
(`ENST00000229239`'s own start, chr12:6,534,832, 99 reads): confidence
0.5, call AMBIGUOUS. `consensus_rise` is positive (0.23, weak) so Level
passes; `body_share_consistency_p` (0.37, only computable at the 60-nt
flank) says the near/far halves agree, so spatial Uniformity passes; but
Phase (`periodicity_p`, none of the four flank lengths clear α=0.05) and
magnitude Uniformity (`gini_body_inframe`, mostly >0.6) both disagree.
2 of 4 lenses agree → confidence 0.5, which is exactly the "some signal,
not fully corroborated" case the battery is meant to distinguish from
either extreme, and matches the event's own AMBIGUOUS call.

**Low-confidence / caveat** — a handful of events with `call=null`
(eligibility INSUFFICIENT, below `min_reads=20`) still get a non-null
confidence, including one at exactly 6 reads scoring confidence 1.0 (its
tiny per-codon vector happened to agree on both lenses by chance).
`confidence` is computed from raw evidence regardless of eligibility —
it does **not** imply enough data to trust. This is harmless in
`consequential.py`'s `tier_confidence` (which only weights aspects with a
non-null `call`, so INSUFFICIENT events are excluded from that average
regardless of their raw confidence), but a caller reading the per-event
`confidence` column directly should gate on `eligibility`/`call` first —
worth calling out since it wasn't obvious until real depth-starved events
surfaced it.

### Attribution

`pair_elongation_neighbors`/`pair_cif_neighbors` need `overlaps_df` (built
from `event_overlap` joined back to `events` for `comp_phase` — see the
plan doc's addendum for why the standard `score-bams` CLI path doesn't
supply this on its own). With it wired in manually for this validation:

| outcome | elongation | CIF |
|---|---|---|
| no_neighbor | 152 | 152 |
| multi_way_overlap | 150 | 150 |
| neither_supported | 34 | 34 |
| both_independent | 0 | 0 |

`multi_way_overlap` (2–8 competing frames, median 2–3) fired on 45% of
elongation events in this locus — **not** a rare edge case; see the STOP
AND FLAG note below. `both_independent` never fired, but that's a depth
artifact of this fixture (only 11 of 336 elongation events have any reads
at all, so two competing frames both having real signal is essentially
unobservable here) rather than evidence the outcome doesn't occur — the
unit tests in `test_attribution.py` exercise it directly with synthetic
evidence.

`pair_boundary_neighbors` (init/term, 60-nt proximity heuristic) found no
same-locus different-frame pairs in this fixture at either aspect —
ran cleanly, just nothing to report at this window size in this window's
annotation.

`classify_termination_downstream`: all 62 termination events classified
`clean_drop` — a direct consequence of the termination floor-firing
finding below (no after-window signal to call readthrough on).

### Floor-firing findings (STOP AND FLAG)

**Termination's `dropoff_significance_p`/`dropoff_continuity_p` never
fired: 0 of 62 termination events reached the 5-codon-per-side floor.**
Traced to the cause, not just observed: P-site coverage in the ±60-nt
window around GAPDH's own stop codon (chr12:6,538,166) is effectively
zero under this pipeline run's metagene-calibrated offsets, even though
(a) 1,792 raw alignments overlap chr12:6,538,000–6,538,300 and (b)
querying the same region directly with the offset calibration's *global*
fallback (bypassing metagene) recovers 63 units of P-site signal nearby.
This is upstream P-site offset-calibration behavior (likely thin
metagene calibration data — only 40 init events, most zero-read, to
calibrate per-read-length offsets against), not a bug in the new
significance lenses and not a reason to lower the floor, per the brief.
Flagging plainly: on this fixture, under the pipeline's default offset
settings, termination's new lenses are implemented and tested but never
exercised end-to-end on real data.

**CIF's per-codon floor (`cif_codon_min_reads=10`) fired sparsely but did
fire**: 7 of 336 elongation events had at least one testable codon (the
other events are the zero/near-zero-read flanking transcripts). Where it
fired, it fired richly — the 3,741-read example above has 72 testable
codons. Not a floor problem, a depth-coverage-of-the-fixture problem: most
of the 336 elongation blocks in this GTF window belong to genes this
particular cohort BAM subset simply doesn't cover.

**3+-way overlapping candidate frames are common, not rare** (STOP AND
FLAG per the brief): 150/336 elongation events, up to 8 competitors. Per
the brief, this was **not** resolved — `multi_way_overlap` is flagged and
each frame's own battery result is reported independently, no policy
invented for N>2. See the plan doc's addendum for the recommendation that
this needs a real design pass, not indefinite deferral, given the
prevalence.

### What wasn't exercised on real data

- **Junction internal consistency** (`junction_internal_consistency`):
  implemented and unit tested, but not wired into the production
  `score_junction_event` call site in `scoring/run.py` this session —
  needs acceptor-side phase plumbing through `feature_event`'s shared
  schema across all four event types, which is more schema churn than
  this pass covers (see the commit message on
  `93c16c3` and the plan doc addendum). Not exercised on GAPDH.
- **Mappability downweight** (`compose_confidence`'s `map_track_low`
  multiplier): unit tested (`test_compose_confidence_downweights_low_mappability`),
  not exercised on GAPDH since no `--mappability-bigwig` was supplied to
  this validation run. Didn't force it, per the brief.
- **bigwig cross-check**: not attempted. The larger real bigwigs mentioned
  in the brief (`/Users/jackt/SRD5A1/...`) are RiboCrypt/BodyMap tracks;
  reconciling their build/chrom-naming against this GTF fixture and
  confirming genuine coverage at a specific locus was not a small
  side-task given the time already spent getting the GAPDH BAM path fully
  wired and validated end to end, and the GAPDH BAM already gave a real,
  traceable result (including the offset-calibration finding above) that
  a bigwig run would not have surfaced (bigwig has no read-level
  information, so P-site offset calibration doesn't apply the same way).
  Relying on the GAPDH BAM as primary real-data evidence, as the brief
  allows.
- `data/translons.bb`/`translons.bed12`: checked for overlap with the
  GAPDH locus (chr12:6,534,512–6,538,374) — zero rows overlap. This call
  set looks to be lncRNA/uORF-focused (`LINC00115`, `LINC01128`, ... in
  the entries near the top of the file), not a GAPDH CDS call set, so
  there was nothing to compare against here. Not forced.

## Follow-up: pancreas cohort bigwig validation

Same locus as the GAPDH validation (chr12, GAPDH ± flanking genes), but a
different question: was the GAPDH BAM's termination-floor-never-fires
finding a thin-calibration artifact of that fixture, or something more
fundamental? Real pancreas-cohort data
(`/Users/jackt/projects/all-RiboSeq/hpc_pancreas_local/`, outside this
repo) makes that testable directly, since the bigwig path has **no**
offset-calibration step at all — coverage is whatever depth the bigwig
already encodes.

### Setup

- **Chrom naming / genome build**: verified, not assumed. The pancreas
  bigwigs (`merged_bigwigs/*.bw`) use `chr1`…`chr22`,`chrX`,`chrY`,`chrM`
  (`chr`-prefixed), matching the bed12 annotations directly — no rename
  needed this time (unlike the GAPDH GTF, which needed `12`→`chr12`). Gene
  coordinates (e.g. `ENSG00000111640.15` / GAPDH at
  chr12:6,534,832–6,538,170) line up with GRCh38, consistent with
  `inputs/gencode.v47.annotation.gtf` in the same tree and with the GRCh38
  coordinates used in the GAPDH validation — same build throughout.
- **bed12 needs no GTF conversion.** `TranslonScorer/io/feature_sources.py`
  already has `from_bed12()`, wired up via
  `extract_events_workflow(..., bed12_path=...)` — it reads columns 0/1/3/5/
  10/11 (chrom/start/name/strand/blockSizes/blockStarts) directly and
  ignores the rest, so the `itemRgb="0"` (vs a real RGB triplet) in these
  files doesn't matter. No new conversion code was needed or written.
- **Pairing chosen**: `merged_bigwigs/merged_good_unique_with_junction.{forward,reverse}.bw`
  paired with `iRibo.Pancreas_pooled.bed12`. Reasoning: "unique" (not
  multi-mapped) avoids exactly the manufactured-periodicity-via-multimapping
  concern the mappability cross-check (item 8) exists for; "good" (the
  broader of the two quality tiers, "good"/"great") maximizes depth for a
  floor-firing test, which is the whole point of this follow-up; "with
  junction" keeps the fullest depth picture even though bigwig can't use
  spanning reads for splice detection either way. `iRibo` was chosen over
  the other four orthogonal callers because `translonscorer_merged_bigwig/
  iRibo.good.unique.with_junction/` already existed as a partially-run
  pairing on disk (events extracted genome-wide, scoring aborted before
  completion — evidently the combo someone already considered primary) and
  its calls include GAPDH itself, preserving direct continuity with the
  GAPDH BAM validation. Region-restricted to chr12:6,400,000–6,700,000 (13
  real genes + 25 short "candidate_orfN" calls) rather than genome-wide,
  since a full genome-wide bigwig score run is exactly what the prior
  partial run had aborted out of.
- **Is the bigwig pre-P-site-shifted?** Checked directly, not assumed:
  raw per-position bigwig depth over GAPDH's biggest exon block
  (chr12:6,537,583–6,537,996) concentrates 86% of signal into a single
  genomic-residue-mod-3 bin. Un-shifted footprint coverage (raw 5′ ends or
  whole-footprint pileups) does not produce that kind of crisp single-frame
  concentration — so yes, this bigwig is pre-shifted at generation time.
  `ribometric_pancreas_offsets.csv` (28/29/30nt → offset 15) is reference/
  provenance for how, not something to feed into `score_bigwigs_workflow`
  (which has no offset-calibration parameter on the bigwig path at all —
  `site` is accepted but ignored, per its own docstring). **But** — see the
  frame-registration finding below: pre-shifted does not mean shifted to
  the SAME register this codebase expects.

### Finding 1: the depth/floor-firing question is answered cleanly — yes

| | GAPDH BAM (thin) | Pancreas bigwig (real depth) |
|---|---|---|
| Termination events with any dropoff lens value | 0 / 62 | 8 / 38 (21%) |
| Elongation events with ≥1 CIF-testable codon | 7 / 336 | 159 / 180 (88%) |
| Total CIF-testable codons (all elongation events) | not computed (too few events) | 5,153 |
| Elongation events with n_reads > 0 | 11 / 336 | 176 / 180 (98%) |

Termination's dropoff lenses (`dropoff_significance_p`/`dropoff_continuity_p`)
now fire on real data — the floor genuinely was sitting idle for lack of
depth, not because of anything wrong with the lens or its floor value. The
GAPDH BAM's specific finding (P-site coverage collapsing to ~0 right at the
stop under that fixture's thin metagene calibration) is confirmed to be a
calibration/depth artifact of that fixture: it doesn't reproduce here,
where there's no calibration step to go wrong and depth is orders of
magnitude higher. (The 8 values that did fire aren't individually
significant at α=0.05 — 0.54, 0.70, 0.94 — which just says these
particular stop codons don't show a dramatic drop in this narrow window;
plausible real biology, not evidence against the lens. Dropoff's window is
anchored directly at the true stop-codon position from raw exon
coordinates, not derived from `phase`/`abs_frame`, which turns out to
matter — see Finding 2.)

### Finding 2 (the headline result): a systematic, strand-dependent frame-registration mismatch — not a calibration artifact, something more fundamental, and NOT a bug in this session's lenses

Checked every elongation event with ≥500 reads in the region (117 events)
by comparing the SCORER's target frame (`abs_frame(phase, strand)`, `a_e`)
against the ACTUALLY dominant genomic-residue-mod-3 bin in the raw bigwig
coverage over that event's span:

**0 of 117 events had their dominant frame match `a_e`.** Every one of
them was strongly periodic (72–98% of signal in one frame — this is
excellent P-site precision, not noise) but the dominant frame was
systematically offset from `a_e` by exactly one nucleotide, and the
*direction* of the offset depends cleanly on strand:

- strand `-1` events: true dominant frame = `(a_e + 1) % 3` in 78 of 79
  checked events; the one exception (892 reads, the lowest-depth event
  checked, only 44% concentration in its own dominant frame vs 72–98%
  everywhere else) is also the one elongation event this validation's
  attribution pass flagged `both_independent` — plausibly genuine
  cross-frame contention rather than a clean single-frame signal to
  register-check in the first place.
- strand `+1` events: true dominant frame = `(a_e + 2) % 3` (i.e. `a_e −
  1`), 38 of 38 checked events, no exceptions.

Concretely, GAPDH's own biggest exon block (event `14065136805926077183`,
same coordinates and same hashed `event_id` as in the GAPDH BAM
validation — the bed12-derived phase and the GENCODE-CDS-derived phase
from the GTF agree exactly, phase 2, `a_e=1`): the pancreas bigwig gives
frame 0 as dominant (86% of 337,292 reads) — frame 1 (`a_e`) carries only
5.3%. Re-querying the **same exact locus from the GAPDH BAM** (built via
this codebase's own offset placement, not an externally pre-shifted
track) gives frame 1 (`a_e`) as dominant (44.1% of 3,834 reads, noisier
because depth there is ~90× lower, but the right register).

This pins the mismatch down cleanly: **TranslonScorer's own phase
computation and its own P-site placement (via `coverage/bam.py`) are
internally self-consistent** — verified independently, not assumed, by
recomputing the identical locus through the BAM path and finding it lands
on `a_e` correctly. The discrepancy is specific to this externally
pre-generated pancreas bigwig's P-site placement convention disagreeing
with this codebase's convention by one nucleotide, in a strand-dependent
way (the signature of, e.g., a coordinate off-by-one that isn't itself
strand-symmetric — a 0-based/1-based slip, or an offset computed from the
wrong read end for one strand during bigwig generation). Isolating the
exact external cause (something in the pipeline that produced
`merged_bigwigs/`, not in `BigwigSetProvider` or `abs_frame`, given the
BAM-path self-consistency check above) is out of scope for this
validation and was not attempted — flagging it plainly rather than
guessing further or patching around it, per the brief.

**What this does and doesn't invalidate:**

- `frame_chisq_p` (the elongation Level lens) tests "is the 3-way tally
  non-uniform at all", not "does `a_e` specifically dominate" — so its
  extremely small p-values on this data (essentially every high-depth
  event) remain a valid statement ("real periodicity exists here") even
  though the register is off. This is exactly the distinction that lens
  was designed around (see its docstring, "preferred over a frame-0-vs-
  rest binomial test") — a useful, unplanned confirmation that it measures
  something meaningfully different from the frame-0-specific metrics.
- Everything that specifically tests "is `a_e` dominant" — the elongation
  `metric` (`effective_in_frame`), `cif_codon_significance`'s one-sided
  binomial test, `body_uniformity_p` (which does still run correctly as a
  *within*-body-in-the-wrong-frame consistency check, just not meaningful
  for "is this really translated") — is measuring the wrong register on
  this specific bigwig+bed12 pairing and should not be read as "GAPDH and
  its neighbors aren't translated here" (they clearly are, per the raw
  periodicity). This shows up starkly in the call distribution: 1 of 180
  elongation events came back SUPPORTED despite 176/180 having real reads
  and the region being unambiguously, strongly periodic.
- Init/term boundary lenses (`periodicity_p`, `gini_body_inframe`,
  `body_share_consistency_p`) and both dropoff lenses are anchored
  directly at the true start/stop genomic position from raw exon
  coordinates, not derived from `phase`/`abs_frame` — so they are NOT
  subject to this mismatch. This lines up with the data: init calls went
  from 1/40 SUPPORTED (GAPDH BAM) to 14/38 SUPPORTED (pancreas bigwig),
  and termination went from 0/62 to 13/38 SUPPORTED + 1 AMBIGUOUS — real
  improvement from real depth, on lenses immune to Finding 2.

**STOP AND FLAG**: this is a real, structural mismatch between an
externally-generated coverage track and this codebase's frame convention,
discovered by this validation, not a bug introduced by this session's
lens work (every affected lens is computing exactly what it says it
computes; the input register was ambiguous). No code changes were made in
response to this — the cause is genuinely upstream/external and unclear
without dedicated investigation (which BAM/bigwig-generation pipeline
produced `merged_bigwigs/`, and what P-site convention it used, are
questions for whoever owns `hpc_pancreas_local/`, not something to guess
at and patch here). Anyone using `merged_bigwigs/` with TranslonScorer
should re-derive or verify P-site registration before trusting any
frame-0-specific metric (elongation `metric`, CIF) from it; chi-square-
style "any periodicity" and boundary-anchored lenses remain trustworthy.

### Finding 3: 3+-way overlap prevalence does not reproduce here — but likely reflects annotation density, not a contradicted result

`event_overlap` for this bed12-derived region has only 4 rows (vs 478 for
the GAPDH BAM run's GTF+ORF-finder-derived event set covering the same
genomic window). With `overlaps_df` wired in the same manual way as the
GAPDH validation: 2 `both_independent`, 1 `neither_supported`, 1
`only_this_frame`, 176 `no_neighbor` — **zero `multi_way_overlap`**,
unlike the GAPDH BAM run's 150/336 (45%).

Read this as a difference in annotation SOURCE density, not a reversal of
the GAPDH finding: `iRibo.Pancreas_pooled.bed12` is a curated, filtered
set of final ORF calls (13 real genes + a modest number of short
`candidate_orfN` entries in this window), while the GAPDH BAM validation's
events came from a GTF plus this codebase's own de-novo ORF-finder
candidate generation — inherently far more redundant/overlapping by
construction. The 3+-way overlap design gap (flagged in the plan doc
addendum) stands; this follow-up doesn't resolve or contradict it, it
just shows the prevalence is annotation-source-dependent, which is worth
knowing when interpreting either number.

### Finding 4: junction bigwig ceiling confirmed as expected, not a bug

All 142 junction events in this region came back `eligibility=
INSUFFICIENT, call=null` — exactly the documented, pre-existing bigwig
ceiling (no CIGAR, no way to count spanning reads from per-base depth
alone; see `_warn_bigwig_cannot_score_junctions` in `workflows.py` and the
project's own bigwig-mappability-gap history). Expected, not investigated
further, per the brief.

### Termination downstream classification, revisited

`classify_termination_downstream` returned `clean_drop` for all 62 events
in the GAPDH BAM run (nothing to classify — no after-window signal at
all). With real depth: 32 `clean_drop`, 6 `readthrough`. Since
`dropoff_after_share` is anchored at the true stop position (immune to
Finding 2), these 6 are a genuine first observation of the readthrough
branch actually firing on data, though verifying whether they represent
real stop-codon readthrough biology versus some other after-window signal
source would need follow-up beyond this validation's scope.

## Open decisions: defaults picked and why

All added to `ScoreThresholds` (`model.py`), each commented
first-pass/unvalidated:

| Field | Default | Reasoning |
|---|---|---|
| `cif_codon_min_reads` | 10 | Order-of-magnitude below the whole-event `min_reads=20` floor — a single codon should need much less than a whole ORF to be individually testable, but 10 still gives a one-sided binomial test real power at p=1/3 (roughly 8+/10 frame-0 hits needed for α=0.05). |
| `lens_agree_frac` | 0.5 | "2 of 3" / "half of the lenses" as the default battery-pass bar, kept uniform across aspects rather than tuned per aspect — real data (this validation) wasn't enough to justify aspect-specific values yet; CIF's 3-lens battery and elongation's 2-lens battery both use it as-is. Revisit once a larger cohort exists. |
| `elong_uniformity_min_codons` | 5 | Matches `periodicity_min_codons` (the pre-existing per-side codon floor for the boundary test) — same order of magnitude, same underlying Mann-Whitney U power requirement. |
| `init_uniformity_gini_max` | 0.6 | Below the Gini midpoint (0.5) but with headroom — a perfectly uniform in-frame series (Gini≈0) is unrealistic even for a clean start (finite codon count, Poisson-ish noise), so the cutoff sits above 0, but well below "half the mass in one bin" (~0.5–0.7 territory empirically, unvalidated). |
| `cif_level_min_frac` | 0.5 | Majority of testable codons significant — same "half" convention as `lens_agree_frac`, chosen for consistency rather than an independent derivation. |
| `cif_contiguity_max_run_frac` | 0.7 | A run covering more than 70% of all significant codons reads as "basically one cluster"; below that, calling it "scattered" seemed reasonable. Loosest of the new thresholds since contiguity's real distribution on a real cohort is unknown. |
| `mappability_confidence_penalty` | 0.7 | Downweight, not a veto (per the plan's explicit instruction) — 0.7 keeps most of an elongation event's confidence intact under low mappability while still moving it, rather than picking something closer to 0 (near-veto) or 1 (no-op). |
| new alpha/significance cutoffs | reuse `periodicity_significance_alpha` (0.05) everywhere | Per the brief's default-unless-real-reason-to-deviate instruction. Nothing in the GAPDH validation gave a concrete reason to pick a different α for chi-square, binomial, or two-proportion tests, so none of them introduced a second alpha field. |
| `confidence`: JSON-only vs first-class | **first-class** (`_RECORD_SCHEMA`) | Same rule cif/n_codons were promoted under: `consequential.py`'s `tier_confidence` needs to weight/aggregate by it across events within a translon, which an evidence-JSON-only field can't support without re-parsing JSON in the composition hot path. Documented in the schema comment and the `177b750` commit message. |

## Design gaps surfaced by real data, not resolved here

- **3+-way overlapping frames are common** (45% of elongation events on
  the GAPDH locus) — flagged per the brief, not resolved. See the plan
  doc's addendum for the recommendation.
- **`event_overlap` isn't wired into `score_bams_workflow`/
  `score_matrix_workflow`.** Pre-existing gap (the contention-scoring
  machinery itself predates this session), but it means elongation
  attribution is invisible via the standard CLI path today without a
  caller manually building `overlaps_df`. Flagged as follow-up in the plan
  doc's addendum.
- **Junction internal consistency's acceptor-side phase plumbing** is not
  wired into the production path (see above) — the lens itself works and
  is tested against synthetic donor/acceptor coverage.
- **`confidence` on INSUFFICIENT events** can be non-null and even 1.0 on
  very low read counts (see the low-confidence/caveat example above) —
  harmless in current composition (excluded via the null `call`), but a
  latent footgun for any future direct consumer of the raw column.

## Plan doc updates

`docs/significance_testing_plan.md` gained:

- An implementation note under §1 (Initiation) explaining the Gini vs.
  share-consistency distinction (the mid-session design-gap fix, commit
  `1db12eb`), with forward-references from §2/§4 so a future reader
  doesn't have to reverse-engineer it from diffs.
- An "Addendum: what §1–§7 landing and real-data validation changed"
  section at the end, covering the three findings above (uniformity
  design gap, 3+-way prevalence, `event_overlap` not wired into
  production) in the plan's own terms.

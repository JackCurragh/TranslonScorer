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

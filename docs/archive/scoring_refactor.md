# Scoring refactor — measure, report everything, decide late

Status: **agreed direction** (2026-07-16). Supersedes the "score-and-decide-in-one-pass, aggregate-only, three-matrix-paths" model. Companion docs: `event_scoring_model.md`, `matrix_scoring_redesign.md`, `matrix_normalisation_strategy.md` (the last one's normalisation contract still applies).

This doc records the target design and the staged plan to get there. It is written against a concrete snag list (below) found by reading the current scoring stack; each snag maps to a correction and a stage.

---

## 0. Principles (the whole thing in four lines)

1. **A scorer measures, it does not decide.** Each per-aspect scorer is a pure function `coverage + event → {named metrics}`. No thresholds, eligibility, or call inside it.
2. **Report every metric that is calculated.** All metrics are first-class, long-form `(event_id, aspect, metric_name, value)`. Nothing hides in a JSON blob; a new metric is new rows, never a schema change. Metrics stay useful regardless of any rule applied on top.
3. **Decide late, from data.** Eligibility, call, and consequentiality are separate passes over the metric table, driven entirely by threshold/policy data. Change a threshold → re-derive calls without re-measuring. No literal cutoff lives in code.
4. **One coverage shape, one scorer path.** Every evidence source (BAM, bigwig, matrix) answers exactly one question — "per-position coverage over these regions" — and flows through the same scorer. No source-specific scoring path.

A corollary the whole codebase should follow: **prefer plain functions in modules over classes/Protocols/dataclasses that don't earn their keep.** Structure should be developer-grown, not scaffolding.

---

## 1. Target pipeline

```
sources ─► events ─► per-position coverage (ONE shape, correctly offset)
             │
             ▼
   pure per-aspect measurement  (init / term / elongation / junction)
             │   emits ALL metrics, no decisions
             ▼
   long-form metric store  (event_id, aspect, metric_name, value, group)
             │
             ├─► calls        = late, data-driven transform (thresholds)
             ├─► feature report = late, aspect-agnostic composition
             └─► consequentiality = late, policy transform
```

`group ∈ {aggregate, sample, cluster}` is a parameter over the same path, not a code branch — see §4.

---

## 2. Snag list → corrections

IDs match the review of 2026-07-16. `file:line` anchors are current at time of writing.

### A. Extensibility (the "many scores, logic later" goal)

- **A1 — aspect set hardcoded at composition.** `report.py:23` fixes `("init","elongation","term")`; `consequential.py:35` fixes the call columns. A new aspect is scored and stored, then silently dropped.
  **Fix:** composition discovers aspects via `group_by("aspect")` and emits one column-block per aspect found. Chain order (init→elong→term) becomes optional display ordering only.

- **A2 — one headline `metric` + JSON `evidence` blob.** `evidence.py:32,174`. Secondary numbers (`clean_in_frame`, `breadth`, `identifiability`, `competitor_share`, `noise_share`, `psi`, …) are buried in a JSON string.
  **Fix:** each scorer returns a flat dict of every number it computes; store all as long-form metric rows. Every metric is queryable. Eligibility/call become their own derived rows/columns, not baked into the measurement.

- **A3 — magic cutoffs.** `report.py:126` (`supported_frac >= 0.5`), `consequential.py:87` (`0.5 + 0.5*expr`).
  **Fix:** no literal survives. `supported_frac` cutoff → `ScoreThresholds`; the consequentiality blend → `ConsequentialityPolicy`, applied in the late pass. Calls/consequentiality recomputable from stored metrics.

### B. Path divergence (make matrix analogous to bigwig/aggregate)

- ~~**B1 — flat-offset default → ~0.33 in-frame floor.**~~ ~~**B2 — FrameRollup scores elongation only.**~~ ~~**B3 — FrameRollup elongation hardcodes `identifiability=None, breadth=1.0, covered_nt=0`.**~~
  **ALL DONE 2026-07-28.** Both the flat-offset mode and the FrameRollup scoring path are deleted (`score_elongation_from_rollup`, `score_matrix_rollup_workflow`, `--gtf`, `--ref-offset`, `pipeline --cds-gtf/--cds-bigbed`). `MatrixProvider` requires `psite_index_dir` and raises without it. Matrix coverage is per-position genomic exactly like bigwig/BAM; the three sources now diverge at a single line (`provider.coverage()` in `_score_events_over_provider`) and are identical downstream.
  **Open decision RESOLVED:** precompute (psite-index) is the only matrix path; the `build-psite-index` step is required. BAM/bigwig still calibrate on the fly per file, which is the right trade at 1–20 files.
  The FrameRollup *utilities* (`build_frame_rollup`, `calibrate_offsets`, `score_frame_rollup`) are retained — `build-psite-index` is built on them and the QC/figure scripts use them.

- ~~**B4 — two scorers (scalar `score_events` + vectorised `score_events_vectorised`) hand-synced.**~~ **DONE 2026-07-28.** Product ships one `scoring.run.score_events` (elongation batched via prefix sums; init/term/junction a per-event loop). The scalar implementation moved to `tests/reference_scorer.py` as an independent oracle — kept deliberately, because the prefix-sum elongation kernel is not obviously correct by inspection and a second implementation to diff against is worth more than it costs. It must never be imported by product code. The GAPDH golden now runs the **shipped** scorer (it previously ran the scalar, reaching the product only transitively).

### C. Normalisation / cross-sample

- **C1 — `size_factors()` is a stub returning `{}`.** `matrix.py:148`. No library-size normalisation; deep samples dominate aggregate coverage.
  **Fix:** derive per-sample size factors from the per-read totals `build-matrix-cache` already computes. Report the size factor / normalised depth as metrics (per principle 2) so normalisation is visible.

- **C2 — expression percentile ranks raw `total_reads`.** `consequential.py:74`. Follows from C1.
  **Fix:** use depth-normalised expression once C1 lands; report both raw and normalised depth as metrics; policy picks one, both stay visible.

- **C3 — aggregate-only; per-sample and cluster unwired.** `workflows.py:292` always calls `coverage(..., by_sample=False)`. Provider supports `by_sample=True`.
  **Fix:** thread `group_level ∈ {aggregate, sample, cluster}` through `_score_events_over_provider → provider.coverage(by_sample=…)`, scoring once per group. Same wiring lands **clustering** (label → k aggregates → score each through the one path) and **per-sample prevalence** for free. Only the matrix provider does `by_sample` cheaply; BAM/bigwig are already per-file. See §4.

### D. Smaller

- **D1 — unstranded coverage feeds both strands.** `workflows.py:307`. Silent ± mixing. **Fix:** require strandedness metadata; raise, don't mix.
- **D2 — map-track gap reads as unmappable.** `workflows.py:211`. Diagnostic-only; low priority. Distinguish "no data" from "zero" only if the track carries it.
- **D3 — junction `psi` computed but unused in the call.** `aspects.py:377`. Dissolves into A2/A3: `psi` reported regardless; whether it gates is a threshold.

---

## 3. Metric store schema

Long-form, self-describing, extensible without migration:

```
event_id     UInt64
aspect       Utf8      # "init" | "term" | "elongation" | "junction" | <any new scorer>
group        Utf8      # "aggregate" | sample_name | cluster_id
metric_name  Utf8      # "clean_in_frame", "breadth", "identifiability", "init_rise", ...
value        Float64
```

Calls are derived rows/columns produced by the late decision pass, not stored by the scorer. `evidence` JSON blob is removed — its contents become metric rows. `map_track_mean`/`map_track_low` become ordinary metrics.

`ScoreRecord` (model.py) — the never-instantiated typed mirror of the old schema — is deleted.

---

## 4. `group_level` and clustering

One parameter, three behaviours, one scorer:

- `aggregate` — every sample collapsed to one pseudo-sample; score once (today's behaviour).
- `sample` — `group_id = sample_name`; one scored set per sample → **prevalence**.
- `cluster` — at a locus, label samples by profile *shape* (`clustering.cluster_locus_profiles` on the shape-normalised `X_cluster`), aggregate within cluster on the scale-preserving `X_score` (`aggregate_cluster_profiles`), then score each of the k cluster aggregates through the **same** per-position scorer.

Label and aggregate stay two steps (different matrices; labels have standalone QC value; lets a "not confusing → skip clustering" trigger short-circuit). `clustering.py` + `matrix_normalisation.py` are kept and become leaf utilities feeding this. The deleted `matrix_scoring.py` was a false start at this against the old bigwig engine.

---

## 5. De-ceremony census (grep-verified 2026-07-16)

- **Delete:** `ScoreRecord` (never instantiated); `DictCoverageProvider` (zero references).
- **Collapse — capability Protocols** (`CoverageProvider`/`SupportsSites`/`SupportsJunctions`/`SupportsMappability`): the whole hierarchy is load-bearing at exactly one `isinstance` (`workflows.py:170`), which *still* catches `NotImplementedError` because bigwig defines `junction_support` only to raise. Replace with a plain "does this source provide junctions? else skip" check. Delete all four Protocols.
- **Simplify providers:** each is config + `coverage()` (+ a stub `size_factors`, + `junction_support`, + `mappability_ledger` — verify the last is live before cutting). A coverage source is "config + a `coverage(regions, site)` function"; prefer a module function / closure over a class implementing an interface.
- **Keep:** `ScoreThresholds` (used; expand it to absorb the A3 cutoffs). `ConsequentialityPolicy` (built with real args 3× in cli; applied late/separately — already correct shape).

---

## 6. Staged plan

Each stage gate-green before the next; no literal cutoffs reintroduced.

1. **Measure/report/decide-late inversion + un-hardcode (A1, A2, A3, D3).** Scorers become pure measurement; long-form all-metrics store; aspect-agnostic composition; cutoffs into thresholds/policy. Mostly `scoring/` + `report.py` + `consequential.py`. Highest leverage; directly serves the goal. — **NOT STARTED** (A1 partly done: composition is already aspect-agnostic).
2. ~~**Path convergence (B1, B2, B3, D1).**~~ — **DONE 2026-07-28**, except D1 (unstranded coverage still mixes strands).
3. ~~**One scorer + de-ceremony (B4, §5).**~~ — **DONE.** B4 done 2026-07-28; `ScoreRecord`/`DictCoverageProvider`/the Protocols were removed 2026-07-16.
4. **Size factors + `group_level` wiring (C1, C2, C3, §4).** Implement size factors; thread `group_level`; land per-sample prevalence and clustering through the one path. — **NOT STARTED.**

---

## 7. Open decisions

- ~~**Offsets: precompute vs on-the-fly.**~~ **RESOLVED 2026-07-28:** precompute (psite-index) is the only matrix path; BAM/bigwig calibrate per file.
- **Clustering trigger:** when is a locus "confusing enough" to cluster rather than aggregate? Lives next to identifiability in the decision pass. Heuristic TBD.
- **Init/term vectorisation:** left as a per-event loop; revisit only if it shows up in profiles.

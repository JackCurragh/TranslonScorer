# Matrix scoring — capability design & performance audit

Status: **draft / working design** (2026-06-25, iteration 2). Supersedes the implicit "aggregate-and-score-once" model in the matrix path. Companion docs: `event_scoring_model.md` (event/aspect model), `matrix_query_engine.md` (query-engine vision), `matrix_normalisation_strategy.md`.

This document:
1. **Characterises the requirements** — functional (FR) and non-functional (NFR) — for making the most of a ~6k-sample sparse RiboSeq matrix.
2. **Audits how to satisfy them quickest and best** — what to pre-compute vs compute on demand, so that adding correctness (per-sample, per-length, offset-calibrated periodicity) does **not** reintroduce the slowness we just removed.

---

## 0. Triggering finding (why this exists)

On the real 6k matrix (chr22, TransCODE v45 comprehensive), `score-matrix` aggregate tier produced **canonical protein-coding ORFs (PT) at `elong_in_frame` ≈ 0.342 / 3% SUPPORTED — the random 3-frame floor**, *below* speculative uORFs. Known-translated CDS cannot show no periodicity ⇒ pipeline defect, not biology.

Structural root cause:
- Aggregate cache `read_totals = (read_id, total_count)` **discards read length**.
- `region_coverage`/`_cov_worker` applies a **flat `ref_offset=15`** to all lengths (5′→P-site).
- Counts summed over ~6k samples with **no per-sample offset calibration**.

Periodicity needs the P-site offset *per read length* (28→12, 29→12…), *per sample*. A flat offset smears the P-site across frames → in-frame → ~1/3. Shape metrics (`init_rise`, `term_drop`) survive because they don't depend on frame.

**The periodicity-aware, per-sample, offset-calibrated machinery already exists** (`matrix_rollup.py`: `build_matrix_rollup`, `tabulate_rollup`, `calibrate_offsets`, `rollup_to_periodicity`, `periodicity_qc_scalable`, `tabulate_profiles`, `tabulate_junctions`, `profile_matrix`). The scoring path bypasses it. Most of this work is **unification + the right cache**, not greenfield.

---

## 1. Current-state audit

**Works:** event model (`scoring/aspects.py`, `scoring/evidence.py`); shape aspects `init_rise`/`term_drop` (discriminate, offset-robust); `read_totals` abundance fast path (~flat in samples, ~90 s warm chr22); the periodicity rollup substrate in `matrix_rollup.py` (offset-independent `phase0`, per-(sample,length,strand,phase0) rollups, GBA, profiles, junctions) — built but **not wired to scoring**.

**Broken / missing:** frame/periodicity scoring unsound (§0); **two divergent coverage paths** (scoring's flat-offset `region_coverage` vs the calibrated rollup) that must converge; junctions dead in scoring (`n_reads=0`, `psi=null`); **no per-sample summarisation** (the matrix's whole point); per-chromosome full-matrix re-query (whole-genome ~50× chr22); report policy degenerate (582/582 consequential); no empirical calibration harness.

---

## 2. Requirements

Notation: each requirement names the **sufficient statistic** it needs — the minimal per-read reduction from which it can be computed. This is what §3 caches.

### 2.1 Functional requirements

- **FR1 — Event extraction.** Build init/elongation/termination/junction events from any feature source (GTF/`--feature-type`, BED12, bigBed, FASTA, sqlite), in **genomic or transcriptome** coordinates. *[exists, genomic]*
- **FR2 — Metric families.** Per event, quantify ribosome-occupancy evidence in these families:
  - **Abundance** — coverage / RPF density / expression percentile. *Stat: `(feature) → count`.*
  - **Initiation** — `init_rise`, flank peakiness, stability. *Stat: coverage by position near start.*
  - **Termination** — `term_drop`. *Stat: coverage by position near stop.*
  - **Elongation / frame** — in-frame fraction, periodicity, breadth, frame dominance. *Stat: `(length, strand, phase0) → count` where `phase0=(pos−cds_start)%3`.*
  - **Junction / splicing** — confident-spanning, ψ (PSI), unspliced. *Stat: CIGAR-N-spanning reads per junction.*
  - **Sub-codon profile** — metagene, A/P/E-site profiles. *Stat: `(length, position) → count`.*
- **FR3 — Calibrated frame.** Frame/periodicity metrics MUST use **per-(sample,length) calibrated P-site offsets**, validated such that **canonical CDS score as translated** (positive-control gate). Offset applied analytically, not baked (see NFR2).
- **FR4 — Group levels.** Every metric summarisable at three levels with `group_level` as a **parameter, not a separate code path**: `aggregate`, `sample` (→ prevalence), `cluster`.
- **FR5 — Prevalence outputs.** For per-sample summaries emit `n_eligible_samples`, `n_supported_samples`, `support_fraction`, and the **per-sample metric distribution** (quantiles, not just mean). Eligibility is per-sample (read-count gated); prevalence denominator = eligible samples.
- **FR6 — Profiles & clustering.** Produce per-locus **[sample × position] profile matrices**; support **sample clustering** (tissue/condition/quality strata) and **ORF/event clustering** (co-translation).
- **FR7 — Differential translation.** Contrast prevalence/metrics between sample groups (e.g. condition A vs B at a locus, or genome-wide).
- **FR8 — Multimapping policy.** `unique` (default) and `guilty_by_association` rescue, **preserving `phase0`** so frame survives redistribution. *[exists in `build_matrix_rollup`]*
- **FR9 — Normalisation.** Per-sample depth + composition normalisation before any cross-sample/cluster comparison. *(`matrix_normalisation_strategy.md`)*
- **FR10 — Calibration harness.** Positive (canonical CDS) and negative (intron / shuffled / UTR) control sets; thresholds set by **separation** (ROC/Youden), per-metric anchoring; phase-2 empirical calibration.
- **FR11 — Report composition.** Compose per-event evidence into per-translon reports with a **calibrated** consequential policy (fix the current 582/582).
- **FR12 — Serving.** CLI today (`extract-events`, build-index, `score-matrix`, `report`, `pipeline`); converge on query primitives `query/expression/profile/compare/qc` over one substrate, group-level parameterised; eventually MCP/public API.

### 2.2 Non-functional requirements

- **NFR1 — Flat in samples.** Scoring a cohort of any size touches **output-sized** data, not full nnz. Prevalence is a group-by on a sample key over a bounded rollup.
- **NFR2 — Offset is analytic, not baked.** Recalibrating the P-site offset MUST NOT force a full index rebuild. Achieved by storing `phase0`/length and applying `frame=(phase0+offset)%3` at query time.
- **NFR3 — One-time, parallel build.** Index build is once per (reference/annotation, matrix), parallel, bounded by **unique reads × small categorical fan-out**.
- **NFR4 — Bounded memory.** **No genome-wide `(sample × position)` materialisation.** Sample-resolved position data is scoped to features/regions.
- **NFR5 — Deterministic.** No order-dependent multimapping; reproducible outputs.
- **NFR6 — Migration safety.** Bit-identical (or characterised-delta) equivalence vs the reference path on golden + HAS_MATRIX tests during the cutover.
- **NFR7 — Cache provenance.** Annotation-bound caches are versioned/rebuildable and stamped with (matrix, annotation, offset-model) provenance; the abundance cache stays annotation-independent.

---

## 3. Performance audit & cache design — "will it become slow?" → no, with the right substrate

**The principle that gave us flat-in-samples scoring still holds; we cached the wrong reduction.** Every metric is an associative reduction over reads; speed = cache the **minimal sufficient statistic per family, keyed by exactly the dimensions it needs**, then query = `filter → combine`, flat in samples.

`(read_id, total_count)` is sufficient for abundance **only** — it was fast *because* it collapsed length and phase, which is *why* frame broke. Speed came from discarding signal.

### 3.1 Primary substrate — a per-reference, feature-scoped, length-resolved coverage index (transcriptome-first)

Build, **in concert with the BAM and per reference/annotation**, scoped to the feature set (CDS/ORFs — not the whole genome):

```
CoverageIndex:  (sample, feature_id, tx_pos, length)  →  count
```
where `tx_pos` is the read 5′ end in **transcriptome (feature-relative) coordinates**.

Why this is the right primary cache:
- **Not bigger than today.** It re-keys the existing `(read_id, sample, count)` store from read_id to position and collapses reads sharing a position (that *is* coverage). Scoped to features (~1–5 % of genome; tens of Mb in transcriptome space) it is **≤ current read-level nnz**, usually much smaller (NFR3, NFR4).
- **Length retained ⇒ offset analytic** (NFR2): never bakes a P-site choice; recalibration is free.
- **Transcriptome coordinates** make frame natural (`phase0=(tx_pos−cds_start)%3`, no CDS-phase bookkeeping, no splicing fragmentation), shrink the position space ~30×, and match RiboMetric's space for cross-tool calibration.
- **Build cost:** projecting genome-aligned reads → transcript(s) via the annotation is the build step (generalises `build_alignment_index`); multi-transcript reads follow the multimap policy (FR8). Annotation-bound by nature (NFR7) — correct, since frame/profile metrics are intrinsically annotation-relative.

### 3.2 Derived reductions (one build, all families)

From `CoverageIndex`, derive on build (or lazily, cached):
- **FrameRollup** — collapse `tx_pos → phase0`: `(sample, feature_id, length, strand, phase0) → count`. Tiny categorical fan-out; offset analytic; serves FR2-elongation at scale, flat in samples (NFR1).
- **Abundance** — collapse length & position: `(sample, feature_id) → count`; or keep the annotation-independent genome-level `read_totals` for pure expression.
- **Profile[sample × position]** — filter a feature, pivot `tx_pos` with the per-(sample,length) offset applied → A/P-site profile matrix; input to FR6 clustering. Scoped to analysis loci (NFR4).
- **Junctions** — built alongside from CIGAR-N-spanning reads (`tabulate_junctions`); serves FR2-junction/ψ.

The FrameRollup is a *strict reduction* of the CoverageIndex, so the two are consistent by construction — not competing caches.

### 3.3 Coordinate split (genomic vs transcriptome)
- **Genomic, annotation-independent:** `read_totals` abundance — keep for expression-only queries that must not depend on an annotation.
- **Transcriptome, annotation-bound (primary):** the feature-scoped `CoverageIndex` + derived FrameRollup/Profiles/Junctions — the substrate for everything frame/profile. Rebuilt per annotation; stamped with provenance (NFR7).

### 3.4 Offset handling
Per-(sample,length) offset, calibrated **centrally on the reduced FrameRollup** (it sees all partitions per sample, which per-partition workers cannot — `calibrate_offsets`), applied analytically `frame=(phase0±offset)%3` (`_frame_at`). No rebuild on recalibration (NFR2). Caveat: analytic shift is approximate within ~offset nt of a feature boundary — negligible for aggregate stats.

### 3.5 Query-time fixes
- Hoist the per-partition index read out of the per-chromosome / per-feature loop — index read is reference-independent, must happen once per partition per run, not once per chromosome (current `_score_events_over_provider` defect).
- `group_level` (aggregate/sample/cluster) is a group-by parameter over the same rollup pass — prevalence is the sample-keyed group-by (NFR1).

### 3.6 Target cost model
- **Build** (once per annotation): CoverageIndex + FrameRollup ≈ current `build-matrix-cache` order (minutes, parallel). Profiles built region-scoped on demand.
- **Score query:** `filter + combine` over FrameRollup → flat in samples, single pass per reference.
- **Profiles/clustering:** region-scoped fetch (`tabulate_profiles` already ~5.4 s for GAPDH across 257 partitions).

Net: correctness (per-sample, per-length, calibrated) is recovered **without** losing flat-in-samples — by caching the length-resolved, feature-scoped, transcriptome `CoverageIndex` and deriving the phase0 FrameRollup from it. The corner to avoid is a genome-wide, offset-baked, length-collapsed `(sample, position)` cache — maximally expensive *and* maximally lossy.

---

## 4. Target architecture

```
feature sources ──► events (init/elong/term/junction), genomic or transcriptome   [FR1]
matrix + annotation ──► build_index (per reference, in concert with BAM):
        CoverageIndex (sample, feature_id, tx_pos, length) → count   [primary, feature-scoped]
        ├─ FrameRollup    (sample, feature_id, length, strand, phase0) → count   [derived; FR2/FR3]
        ├─ Abundance      (sample, feature_id) → count   (+ genome read_totals)   [FR2]
        └─ Junctions      per-junction spanning support                           [FR2]
        offset calibration  per (sample,length)  [calibrate_offsets; FR3/NFR2]
score(events, index, group_level): ONE substrate
        frame ← FrameRollup + analytic offset ;  shape ← init/term ;  junction ← Junctions
        group_level ∈ {aggregate, sample→prevalence, cluster}   [FR4/FR5]
profiles(loci) ──► [sample×position] ──► cluster samples / ORFs   [FR6]; differential [FR7]
normalise [FR9] ──► calibrate vs controls [FR10] ──► report/policy [FR11] ──► serve [FR12]
```

Single structural change underpinning it all: **scoring consumes the length-resolved transcriptome CoverageIndex (and its FrameRollup), not the flat-offset aggregate coverage. One substrate, `group_level` as a parameter.**

---

## 5. Sequencing

0. **Prove the root cause, in transcriptome coordinates** (cheap; de-risks everything *and* validates the projection): take a well-expressed canonical chr22 CDS, project a few partition BAMs' reads to transcript coordinates, compute frame fraction **by read length** — show periodicity appears at the correct per-length offset and vanishes at flat-15.
1. **Build the primary substrate** — feature-scoped, length-resolved transcriptome `CoverageIndex` + derived FrameRollup (generalise `build_alignment_index`/`build_matrix_rollup`); per-(sample,length) `calibrate_offsets`.
2. **Unify scoring** onto it: frame + junction via the index; keep `init_rise`/`term_drop` and `read_totals`. **Positive-control gate**: canonical CDS must score SUPPORTED.
3. **Sample axis**: `group_level=sample` → prevalence (FR5); `cluster`.
4. **Profiles + clustering** (FR6) and differential translation (FR7).
5. **Calibration harness** (FR10); fix report policy (FR11) and per-chromosome re-query.
6. **Serving / query primitives** (FR12) → engine/MCP.

## 6. Open questions / risks
- **Read→transcript projection**: multi-transcript reads, overlapping ORFs, isoform ambiguity — define the policy (unique-transcript vs distribute; interaction with FR8 GBA).
- **Boundary approximation** of the analytic offset shift near feature ends (NFR2 caveat).
- **Per-sample eligibility** thresholds — most samples have too few reads per event; define the denominator and minimum-read gates (FR5).
- **Sample×position memory** ceiling — where we refuse genome-wide sample-resolved profiles (NFR4).
- **Genomic vs transcriptome for non-canonical / novel ORFs** without a transcript model — may need a genome-feature CoverageIndex variant.
- **Normalisation** model before cross-sample comparison (FR9) — depth + composition.
- **Provenance/versioning** of annotation-bound caches (NFR7).
```

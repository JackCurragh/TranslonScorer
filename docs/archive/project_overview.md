# TranslonScorer & RiboSeq Pipeline — Project Overview

*Last updated: 2026-06-27. This is a working snapshot; code and HPC state evolve faster than docs.*

---

## 1. The Ecosystem

Four tools live under `/Users/jackt/projects/all-RiboSeq/`:

| Project | Path | Role |
|---------|------|------|
| **ensembl-genes-nf** | `ensembl-genes-nf/` | Main Nextflow pipeline: acquisition → QC → alignment → matrix/BAM outputs |
| **TranslonScorer** | `translonscorer/` | ORF scoring engine — the primary focus of active development |
| **RiboMetric** | `RiboMetric/` | QC metrics (periodicity, offsets, read-length distribution, CDS enrichment) |
| **get-RPF** | `get-RPF/` | Read acquisition and collapsing |

TranslonScorer is NOT called from ensembl-genes-nf — it runs as a separate step consuming the pipeline's outputs.

---

## 2. Data: The 6k-Sample Matrix

The main dataset is a 6,000-sample Ribo-seq cohort on the EMBL-EBI HPC at:

```
/hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/
```

### global_partitioned — the scoring substrate

```
global_partitioned/<PREFIX>/          # 257 partitions, 4-nt prefix sharding (AAAA..TTTT)
  reads.parquet                       # read_id → genomic coords + sequence + length
  counts.parquet                      # read_id × sample_id → count (sparse)
  samples.parquet                     # sample_id → sample_name, study_id
```

- ~1 TB for 5 partitions; not copyable locally — all work runs on HPC
- reads are sharded by the first 4 nt of their sequence (so one read always lands in exactly one partition)
- `counts.parquet` is the nnz-sparse join table; the counts themselves are small but the join with reads is large
- rRNA reads were removed from the matrix (~56M reads across the cohort) — but this had zero effect on frame QC because rRNA reads don't overlap CDS

### studies_partitioned — NOT used for scoring

CSR npz format, no genomic coordinates, per-study. Cannot be scored by TranslonScorer.

### Key HPC output files (translonscorer_6k/)

```
psite_index_filtered/
  qc_per_sample.parquet               # 5228 samples, 30 QC columns
  offsets.parquet                     # 96,476 (sample_name, length, offset) rows
  per_sample_length_frames.parquet    # 76,751 (sample_id, length, n_reads, f0/f1/f2/f0_pct)
  usable_sample_lengths.parquet       # 11,088 (sample_id, length) pairs, 2270 samples
  usable_samples_pif.parquet          # 2270-row summary
```

---

## 3. P-site Index Pipeline (`build-psite-index`)

Produces calibrated P-site offsets and per-sample QC before scoring. Three phases:

**Phase 1a — matrix rollup** (`matrix_rollup.py`):
- Scans all 257 partition BAMs in parallel (spawn workers, not fork — polars threadpool deadlocks fork)
- Assigns each in-CDS read to `(sample, length, strand, phase0)` where `phase0 = (cds_phase + 5'pos) % 3` — this is offset-independent
- Groups into `FrameRollup`: per-(sample, length) frame distribution
- `calibrate_offsets()` finds the per-(sample, length) offset that maximises frame-0 reads
  - Only 3 possible offsets (10, 11, 12 — mod-3 classes in range 10-18)
  - Requires ≥50 CDS reads per (sample, length) to calibrate
  - Added cohort-modal fallback: for (sample, length) pairs below the threshold, use the most common offset at that length across all calibrated samples
  - Result: 96,476 (sample, length) pairs covered (was 79,970 before fallback)

**Phase 1b — fast_reads_qc** (`matrix_qc.py`):
- Derives 5' dinucleotide from partition directory name (partitions are sharded by 4-nt prefix so `dir.name[:2]` = 5' dinuc) — avoids reading the sequence column for all reads
- Only reads sequence for 3' dinucleotide
- Aggregates inside the lazy plan before collect (`streaming=True`) to avoid OOM
- Produces per-sample length distribution and ligation bias metrics

**Phase 2 — bake index**:
- Combines calibrated offsets with QC metrics into `qc_per_sample.parquet`

### Key bugs fixed in this pipeline

- `fast_reads_qc` originally OOM'd because it joined the sequence column for all 5228 samples per partition (~200GB RAM killed). Fix: derive 5' dinuc from dir name, only read sequence for 3' dinuc, stream
- `rfd` values in psite_index.py are lists `[f0, f1, f2]` not dicts — calling `.items()` on them crashed. Fix: `enumerate(fd)`
- Polars `group_by()` returns group keys as tuples, not scalars — `for sname, grp` failed silently (str of tuple `"('DRR428600',)"` never matched sample IDs). Fix: `for (sname,), grp`

---

## 4. Sample QC: Per-Length PIF Filter

**Why aggregate periodicity fails**: summing frame counts across all read lengths collapses to ~0.33 (random) because different lengths peak in different frames under a uniform offset. The per-length signal is far stronger.

**Per-length PIF** (Proportion In Frame, `f0_pct = frames[0]/total * 100`):
- Computed for each (sample, length) pair at lengths 25-37 nt, min 100 reads
- A (sample, length) pair is "usable" if its f0_pct > the cohort 75th percentile at that length — a data-driven noise floor
- A sample is "usable" if it has ≥10,000 total usable reads

### 6k cohort breakdown (5228 samples total)

| Category | Count | Description |
|----------|-------|-------------|
| No CDS reads | 901 | No frame data possible |
| Good Ribo-seq | 1438 | f0 ≥ 60%, rpf_28_32_prop ≥ 0.3 |
| Bad processing | 725 | RPF lengths OK but f0 ~36% flat — genuinely aperiodic, not fixable |
| Misclassified | 410 | Disome/long footprints at 33-37 nt, f0 ~66-68% — recoverable |
| Borderline | 1754 | Mixed signal |
| **Usable (PIF filter)** | **2270** | 11,088 (sample, length) pairs |

The PIF filter recovered 883 samples missed by the old periodicity ≥ 0.5 threshold.

---

## 5. Scoring: `score-matrix` Command

### One path: `--psite-index` (required since 2026-07-28)

`MatrixProvider` reads coverage from the P-site index
(`psite_index.query_genomic_coverage`), which stores per-(sample, length)
offset-corrected genomic positions. Scoring then runs through the ordinary
`_score_events_over_provider` pipeline, so init/term/junction/mappability are
scored alongside elongation. `(sample, length)` pairs in
`usable_sample_lengths.parquet` are honoured automatically if that file is
present in the index dir.

```bash
translonscorer score-matrix \
  --events-dir /path/to/events \
  --matrix-dir /hps/.../global_partitioned \
  --store-dir /path/to/scores \
  --data-version v1 \
  --psite-index /hps/.../translonscorer_6k/psite_index_filtered
```

Omitting `--psite-index` is now a usage error rather than a silent fall back.

### Why the two removed paths were removed

**Flat offset (`--ref-offset`, formerly the default).** It applied
`ref_offset=15` to every read length. True P-site offsets are per read length
(28 nt → 12, 29 nt → 12, …); a single flat offset smears the P-site across
frames, collapsing `elong_in_frame` to ~0.33 (random). Canonical protein-coding
CDS scored ~0.342 in-frame — indistinguishable from noise. `init_rise` and
`term_drop` were unaffected (coverage-shape metrics, not frame-dependent).
A mode whose only behaviour is to produce noise is a footgun, not an option.

**FrameRollup (`--gtf`).** It got the offsets right, but as a structurally
separate pipeline in transcript coordinates that never went through
`_score_events_over_provider` — so it scored elongation only, ignored
`--context-gtf`/`--mappability-bigwig`, and hardcoded `identifiability=None,
breadth=1.0`. The P-site index gives the same offset accuracy with every event
type through the shared scorer, so FrameRollup had nothing left to offer as a
*scoring* path.

The FrameRollup *utilities* (`build_frame_rollup`, `calibrate_offsets`,
`score_frame_rollup` in `matrix/rollup.py`) are unchanged and still in use:
`build-psite-index` is built on them, and the QC/figure scripts use
`score_frame_rollup` to demonstrate calibrated-vs-flat periodicity. Only the
scoring path was removed.

### The two coverage strategies (deliberately kept)

`MatrixProvider` (annotation scale, thousands of samples via a prebuilt index)
and `BamSetProvider`/`BigwigSetProvider` (1–20 files read directly) are both
first-class. They differ **only** in how per-locus coverage is obtained — one
`provider.coverage(regions, site=...)` call at `workflows.py:294`. Events,
scorer, report, consequentiality and store below that line are identical.

---

## 6. Event Scoring Model

The core architectural decision: **score deduplicated genomic events once, compose into per-translon reports**.

### Why

8.85M translons in the TransCODE v45 catalog → 54.5M blocks. After deduplication:
- Elongation: 54.5M → 681k (80× reduction)
- Initiations: 8.85M → 388k (23×)
- Terminations: → 176k (50×)
- ~78× total dedup; 1,518,455 unique events genome-wide

Shared evidence is computed once and is consistent across all translons that use it.

### Event types and keys

| Aspect | Key | Dedup |
|--------|-----|-------|
| Initiation | (chrom, start, strand) | 388k |
| Termination | (chrom, stop, strand) | 176k |
| Elongation | (chrom, seg_start, seg_end, strand, **phase**) | 687k |
| Junction | (chrom, donor, acceptor, strand) | 267k |

`event_id = xxh3_64(canonical_key)` — deterministic, stable across annotation rebuilds.

### Scoring per aspect

- **Initiation**: log2 fold-change `log2((body_level + α) / (out_level + α))` at multiple flank lengths (9/18/30/60 nt). Median codon-bin levels (spike-robust). SUPPORTED if log2FC ≥ 1.0 (≥2×).
- **Elongation**: `elong_in_frame` (in-frame A-site fraction) + `breadth` (covered_nt/span_nt). SUPPORTED if in_frame ≥ 0.5 AND breadth ≥ 0.2.
- **Termination**: mirror of initiation (step-down).
- **Junction**: confident_spanning reads (≥6 nt overhang both sides). SUPPORTED if ≥20 confident spanning.

### Contention / identifiability

~49% of elongation events are contended (uORF-over-CDS overlaps). Each contended elongation event carries:
- `clean_in_frame`: in-frame signal in the unconfounded portion
- `identifiability`: fraction of contended coverage attributable to this event's frame
- `status = AMBIGUOUS` if identifiability < 0.5

### Storage

- Structure: `translons.sqlite` extended with `events`, `feature_event`, `event_overlap` tables
- Scores: append-only partitioned Parquet `fact_event_score/` (grain: event_id, group_id, aspect, tier, annotation_version, data_version)
- Re-scoring = new partition rows, never mutating existing data

### Current results (genome-wide aggregate tier, chr12 validation)

Support rates (SUPPORTED/eligible):
- Junction: 91% (122k supported)
- Elongation: 68% (379k supported, 130k ambiguous from contention)
- Initiation: 53% (98k supported)
- Termination: 52% (32k supported)

### Consequential set

8.85M translons → 536,731 distinct features → **10,503 'most consequential'** (novel non-CDS, full-chain-supported). Refined to **6,329** after reproducibility filter (supported in ≥3 samples individually — the aggregate tier over-counts; single-sample signal artifacts are real).

---

## 7. Matrix Query Engine

Full design note: `translonscorer/docs/matrix_query_engine.md`

### Cost model (256-partition, 6k samples)

| Phase | Cost | Scales with |
|-------|------|-------------|
| Phase A (annotate CDS frames) | ~640 s once | N_unique reads, NOT sample count |
| Phase B old (Python loop) | ~83 s/sample | Sample count — was the bottleneck |
| Phase B new (Polars streaming) | ~0.1 s/sample | Effectively flat |

Full-scale warm query: **7.7 s vs old 3120 s = 405× faster**.

### Key components in `matrix_rollup.py`

- `build_frame_rollup()`: per-partition parallel BAM scan → `(sample_name, length, strand, phase0, feature_id, count)` rollup. spawn workers (not fork).
- `calibrate_offsets(agg)`: find per-(sample, length) offset from reduced rollup
- `score_frame_rollup(rollup, offsets)`: apply offsets analytically, aggregate to per-(feature, length) `elong_in_frame`
- `tabulate_profiles()`: per-(sample, gene, position) A-site coverage matrix for clustering
- `tabulate_junctions()`: per-junction spanning read counts with overhang classification

---

## 8. Periodicity QC Deconfliction

**Root cause of near-zero periodicity scores**: CDS interval conflict rate was 713%. All CDS exon intervals from all transcript isoforms were pooled into one per-(chrom, strand) list. Alternative splice isoforms and overlapping genes caused most positions to be covered by 7+ conflicting intervals with different phases — frame assignment was arbitrary, collapsing to uniform ~0.33.

**Fix** (`_deconflict_intervals()` in `matrix_qc.py`): event-sweep that only emits intervals where the phase is unanimous across all overlapping transcripts. 3.67M raw intervals → 335k deconflicted segments.

**Why aggregate periodicity still understates signal**: per-length frame counts partially cancel when summed (different lengths peak in different frames). The per-length signal is far stronger: 29 nt dominant-frame fraction = 0.94 cohort-wide.

---

## 9. RiboMetric Score Redesign

Active work in `RiboMetric/` (two checkouts: `RiboMetric/` and `RiboMetric-v120/`; work targets `RiboMetric/`).

**Problem**: three disagreeing scoring systems (cards / badges / gate), scores not anchored (calibration ranges inert), arbitrary transforms.

**Design**: every score anchored 0 (null model) → 1 (ideal), unified config-driven resolver, tiered metrics.

| Tier | Metrics | Gate? |
|------|---------|-------|
| Tier 1 (identity) | periodicity_dominance, periodicity_information, CDS enrichment | Yes |
| Tier 2 (usability) | recommended_read_proportion, coverage breadth, saturation | No |
| Tier 3 (caveats) | dup/multimap/softclip/terminal bias | No |

**Key scoring changes done (S0–S5)**:
- `periodicity_dominance`: raw fraction (identity method), pass ≥ 0.70, warn ≥ 0.50; 1/3 baseline marker on plot
- CDS enrichment ratio: `1 - 1/E` where E = observed_CDS_fraction / expected_CDS_fraction (corrects organism confound)
- Terminal bias: `inverse_linear(max_value=2.0)` KL divergence
- `sqrt` dropped from periodicity_information
- Report restructured into tier sections + context strip + diagnostics

**Parked**: empirical recalibration (Phase 2, needs corpus). `frame_dominance_rescaled` still in code but should be switched back to `identity` when work resumes.

---

## 10. Known Bugs and Gaps

### Active / high priority

| Bug | Location | Status |
|-----|----------|--------|
| `elong_in_frame` at ~0.33 random floor | matrix scoring flat offset | **Fixed 2026-07-14, path removed 2026-07-28** (still pending real-cohort validation): `MatrixProvider.coverage()` reads the P-site index (`psite_index.query_genomic_coverage`) and there is no longer a flat-offset alternative to fall into — `--psite-index` is required and `--ref-offset`/`--gtf` are gone. See §5. |
| Junction scoring dead (all n_reads=0) | `workflows._score_events_over_provider` | **Fixed 2026-07-14**: root cause was `_score_events_over_provider` never calling `provider.junction_support()` — both `MatrixProvider`/`BamSetProvider` implement it, it was just never invoked. See `_junction_support_for_chrom`. |
| init_rise/term_drop spuriously inflated near a splice site | `scoring/aspects.py` leader/UTR flank | **Fixed 2026-07-14**: flank was flat genomic, reading intronic sequence when the true leader/UTR was spliced. Now splice-aware via `--context-gtf` on score-matrix/score-bams/pipeline. See `_project_flank`. |
| Bigwig scoring unreachable from the event pipeline | `coverage/bigwig.py`, `workflows.py` | **Fixed 2026-07-14**: `BigwigSetProvider` existed but was never imported by `workflows.py`, and its `stranded` flag was a stub (always sum-merged fwd+rev). New `score_bigwigs_workflow` + `score-bigwig` CLI command; `coverage()` now emits a real `strand` column. Still lossy by design: no P/A site, no junction spanning, no mappability. |
| `--offset-method metagene` silently dead on score-bams/pipeline | `workflows.score_bams_workflow` | **Fixed 2026-07-14**: `BamSetProvider` needs genomic start-codon coords to calibrate; `score_bams_workflow` never passed them, so it always fell back to `global_offset`. Now derives them from extracted `init` events. Only works for genome-aligned BAMs (`transcriptome=False`) — transcriptome-aligned still falls back to global. |
| `hrf` metric computes raw max count, not frame periodicity | `core/scoring.py:115` | **Legacy-only** (2026-07-14 audit): only reachable via the deprecated `profiles`/`profiles_from_bigwig` command, not the current event-scoring pipeline (`score-matrix`/`score-bams`/`score-bigwig`) — lower priority than it reads. Unfixed. |
| Composite score = raw unbounded sum, depth-dominated | `score_orfs_workflow` | **Legacy-only**: `score-orfs` command no longer exists (removed; see `test_deprecated_commands_removed`). Only relevant if this legacy path is ever revived. Unfixed. |
| 582/582 consequential (100% pass rate) | `consequential.py` policy | Policy was a stub; `apply_policy` is now real but may need threshold tuning |
| Per-sample reproducibility not in policy | `consequential.py` | `ConsequentialityPolicy` has no reproducibility knobs yet |
| Mappability-based eligibility does not exist | `coverage/base.py` `SupportsMappability`, `scoring/evidence.py` | **Diagnostic annotation added 2026-07-14** (not the same thing as the `SupportsMappability`/`mappability_ledger()` BAM-NH-tag stub, which is still unimplemented and untouched): `--mappability-bigwig` on score-matrix/score-bams/score-bigwig/pipeline attaches a precomputed mappability track; every scored event gets `map_track_mean`/`map_track_low` (never affects eligibility/call — annotation only), surfaced in the composed report as `{aspect}_map_track_mean`/`{aspect}_map_track_low`. See `workflows._map_track_for_chrom`. |

### Infrastructure / pipeline

| Issue | Detail |
|-------|--------|
| `build-psite-index` teardown hang | spawn Pool / polars threadpool join hang after writing last partition; output is complete |
| Cache write not atomic | `sink_parquet` goes straight to final path — killed job leaves corrupt file that reruns silently skip |
| `build-matrix-cache` teardown hang | Same pattern |

---

## 11. Code Structure (key files)

```
translonscorer/
  TranslonScorer/
    matrix_rollup.py      # FrameRollup, calibrate_offsets, tabulate_profiles, tabulate_junctions
    matrix_qc.py          # fast_reads_qc, periodicity QC, _deconflict_intervals
    psite_index.py        # build_psite_index orchestration (Phases 1a/1b/2)
    workflows.py          # score_matrix_rollup_workflow (--psite-index wired here)
    cli.py                # CLI commands including score-matrix --psite-index
    pipeline/
      event_extract.py    # extract_events, run_extract
      event_score.py      # score_initiation_event, score_elongation_event, score_termination_event
      matrix_rollup.py    # older module — superseded by TranslonScorer/matrix_rollup.py on the new path
    scoring/
      evidence.py         # event record schema
      run.py              # score_elongation_from_rollup
    report.py             # compose_report (per-translon aspect rollup)
    consequential.py      # apply_policy, ConsequentialityPolicy
  docs/
    event_scoring_model.md
    matrix_query_engine.md
    matrix_scoring_redesign.md
    upstream_qc_from_matrix.md
    RELEASE.md
    REFACTOR_TASKS.md
  outputs/
    reference_scores_gapdh.parquet   # golden reference (22 event-aspect rows)
  scripts/
    score_reference.py
    verify_vectorised.py
    per_translon_report.py
    consequential_set.py
    refine_consequential.py
    expression_consequential.py
    finalize_consequential.py
```

---

## 12. Release Machinery

- `bump2version` bumps both `__init__.py` and `pyproject.toml` (dynamic version via attr is unsafe due to cli import chain)
- CI: py3.10/py3.12 + samtools + coverage → Codecov + mypy (report-only for now)
- Release: tag `v*` → GitHub Actions → twine → PyPI Trusted Publishing
- **Not yet wired**: register PyPI Trusted Publisher + create `pypi` GitHub environment
- Current version: 0.1.1
- Runbook: `translonscorer/docs/RELEASE.md`

Active development branch: `local/read-assignment-prototype` on `github.com/JackCurragh/TranslonScorer`

---

## 13. Immediate Next Steps

1. **Validate `--psite-index` scoring on HPC**: reinstall from branch, run `score-matrix` with `--psite-index psite_index_filtered/` on a chr subset, check that `elong_in_frame` for canonical CDS (PT) rises above ~0.33 baseline
2. **Diagnose junction scoring**: all 154 confident_spanning events show n_reads=0, psi=null — separate from the frame offset bug
3. **Fix atomic cache writes**: write to `.tmp` path then rename to avoid corrupt-file silent skip on reruns
4. **Fix teardown hang**: spawn Pool teardown after polars-heavy workers — likely a Rust threadpool join issue; try `pool.terminate()` instead of `pool.join()` after collecting results
5. **Per-sample reproducibility in consequential policy**: wire the 40%-cut refinement (elong supported in ≥3 samples) into `ConsequentialityPolicy`

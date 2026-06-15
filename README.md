# TranslonScorer

Score translation events (**translons**) from Ribo-seq data — for **individual
samples** (1–20 BAMs or bigWigs) and at **annotation scale** (a sparse
multi-sample matrix of thousands of samples).

TranslonScorer scores each *aspect* of a translation event chain separately —
**initiation** (P-site at the start codon), **elongation** (A-site, in-frame),
**termination** (A-site at the stop), and **splice junctions** — and reports
per-translon evidence plus a dynamic, context-aware *consequentiality* score. It
does **not** apply hard length/biotype gates: a start-stop ORF or a translated
small RNA is scored on its evidence like anything else.

---

## Installation

```sh
pip install "TranslonScorer @ git+https://github.com/JackCurragh/TranslonScorer"
# optional extras: [bigwig] [zarr] [viz] [fastbam] [dev] or [full]
pip install "TranslonScorer[full] @ git+https://github.com/JackCurragh/TranslonScorer"
```

Core needs `pysam`, `polars`, `numpy`. `scipy` is optional (hierarchical
clustering of confusing loci).

---

## Concepts (read this first)

The scoring model is **event-centric**. Translons are decomposed into
*deduplicated genomic events* (one start codon shared by many isoforms is one
init event, etc.), each event is scored **once**, and the scores are composed
back into per-translon reports. This is ~78× less work than scoring every
translon independently and guarantees shared evidence is consistent.

The pipeline is the same regardless of data source — only the **coverage
provider** changes:

| source | provider | command |
|---|---|---|
| 1–20 genome/transcriptome BAMs | `BamSetProvider` | `score-bams` |
| 1–20 bigWigs (coverage only) | `BigwigSetProvider` | (via API) |
| annotation-scale sparse matrix | `MatrixProvider` | `score-matrix` |

All three feed one scoring core and produce the same `fact_event_score` store
and per-translon report.

The end-to-end flow is four steps:

```
extract-events        annotation (sqlite) ─► genomic event store
   │
score-bams / score-matrix   events + coverage ─► append-only score store
   │
report                scores + events ─► per-translon report (+ consequentiality)
   │
consequential         re-gate an existing report under a different policy
```

`pipeline` runs **extract → score → report** in one call (see below).

---

## Quick start

### A. Individual files (your own BAMs)

Score a handful of genome-aligned Ribo-seq BAMs against an annotation:

```sh
# one-shot
translonscorer pipeline \
  --sqlite annotation.translons.sqlite \
  --bam sampleA.bam --bam sampleB.bam \
  --out-dir results/ \
  --offset-method global --global-offset 12

# results/report.parquet now has one row per translon with per-aspect calls
```

Transcriptome-aligned BAMs (reads positioned along the mRNA) are projected to
genome coordinates — UTRs included, so 5′UTR uORFs land correctly:

```sh
translonscorer pipeline \
  --sqlite annotation.translons.sqlite \
  --bam tx_aligned.bam --transcriptome --annotation annotation.gtf \
  --out-dir results/
```

### B. The matrix (annotation scale)

Score events against the sparse multi-sample matrix (partition directories):

```sh
translonscorer pipeline \
  --sqlite annotation.translons.sqlite \
  --matrix matrix/partition_000 --matrix matrix/partition_001 \
  --chrom chr12 \
  --out-dir results_chr12/ \
  --data-version matrix_v1
```

`--chrom` restricts extraction (handy for a tractable run); omit it for the
whole genome. See [`runs/run_matrix_scoring.sh`](runs/run_matrix_scoring.sh) for
a reproducible matrix run you can copy.

### C. Step by step (full control)

```sh
# 1. events (once per annotation; reused across all scoring runs)
translonscorer extract-events --sqlite anno.sqlite --out-dir run/events --chrom chr12

# 2a. score from BAMs ...
translonscorer score-bams  --events-dir run/events --bam a.bam --bam b.bam \
    --store-dir run/scores --data-version v1
# 2b. ... or from the matrix
translonscorer score-matrix --events-dir run/events \
    --partitions matrix/p0 --partitions matrix/p1 \
    --store-dir run/scores --data-version v1

# 3. compose the per-translon report (+ default consequentiality)
translonscorer report --store-dir run/scores --events-dir run/events \
    --out run/report.parquet --data-version v1

# 4. (optional) re-gate under a stricter policy without re-scoring
translonscorer consequential --report run/report.parquet \
    --out run/strict.parquet --min-tier-confidence 1.0
```

Key options: `--offset-method {global,file,metagene}` and `--global-offset` /
`--offsets-file` control P-site offsets (BAM mode); `--site {A,P}` picks the
queried site; `--multimap unique` keeps only unique mappers; `--min-tier-confidence`
and `--min-expression-percentile` set the consequentiality floors.

### Output

- `events/` — `events`, `feature_event`, `event_overlap` Parquet trees (the
  deduplicated genomic events and translon→event membership).
- `scores/` — `fact_event_score` store, append-only, partitioned by
  `data_version`/`tier`. Re-scoring with new data/methods adds partitions; it
  never mutates existing ones (reproducible, A/B-able).
- `report.parquet` — one row per translon: `{init,elongation,term}_{call,metric,
  n_reads,supported_frac}`, `junction_*`, `total_reads`, `tier_confidence`,
  `expression_percentile`, `consequentiality_score`, `consequential`.

---

## Package layout

The package is organised as a **functional core** (pure transforms) wrapped by
an **imperative shell** (I/O + orchestration). See
[`docs/architecture.md`](docs/architecture.md) for the design rationale.

```
TranslonScorer/
├── cli.py              Click CLI — every command lives here
├── workflows.py        orchestration shell: pipeline / extract / score / report
├── model.py            dataclasses: Region, OffsetParams, ScoreThresholds,
│                       FrameSupportParams, ConsequentialityPolicy, ScoreRecord
│
│   ── functional core (pure, no I/O) ──
├── events.py           annotation → deduplicated genomic events + membership
├── offsets.py          P/A-site offset calibration (global | file | metagene)
├── scoring/            the event scorer
│   ├── aspects.py        init / elongation / termination / junction scorers
│   ├── evidence.py       evidence → eligibility/call decisions (thresholds)
│   └── run.py            vectorised scorer; score_events()
├── report.py           compose per-translon report from event scores
├── consequential.py    dynamic consequentiality policy (no hard gates)
├── qc.py               periodicity / frame-dominance QC (pure)
├── clustering.py       profile clustering for confusing/contended loci
├── frame_support.py    frame-posterior estimation (linear/latent/HMM)
│
│   ── coverage providers + profile building (the stateful edge) ──
├── coverage/
│   ├── base.py           capability Protocols (CoverageProvider, SupportsSites…)
│   ├── bam.py            BamSetProvider (1–20 BAMs; offsets; transcriptome→genome)
│   ├── bigwig.py         BigwigSetProvider (coverage-only)
│   ├── matrix.py         MatrixProvider (sparse annotation-scale matrix)
│   ├── profile.py        pure P/A-site profile from reads + offset table
│   ├── transcriptome.py  transcript→genome projection (exon walk)
│   └── profiles.py, locus_profiles.py, locus_features.py, transcript_coords.py,
│       mapped_index.py, index_from_bam.py, junctions.py, junction_model.py
│                         profile/index/junction construction from sources
│
│   ── the sparse-matrix engine ──
├── matrix_rollup.py    tabulate coverage/junctions over the matrix (region_coverage)
├── matrix_scoring.py   score ORFs/profiles off the matrix
├── matrix_normalisation.py  per-locus matrix normalisation
├── matrix_qc.py        per-sample length/periodicity QC straight off the matrix
│
│   ── I/O adapters ──
├── io/
│   ├── annotation.py     GTF → CDS / exon blocks (build_cds_blocks, build_exon_blocks)
│   ├── annotation_bundle.py  build/load an annotation bundle directory
│   ├── bam.py            read-id parsing, chrom normalisation, CIGAR, junctions
│   ├── matrix.py         sparse-matrix partition manifests / readers
│   ├── store.py          event + fact_event_score Parquet store read/write
│   └── inspect.py        Parquet inspection helper (the `inspect` command)
│
│   ── frame analysis ──
├── frame/
│   ├── bleed.py, deblur.py, hmm.py, latent.py   frame-correction models
│   ├── frame_crosstalk.py, frame_method_compare.py
│   ├── frame_disambiguation.py   ambiguous-read frame assignment
│   └── read_assignment.py        unique/fractional/EM/frame-aware assignment
│
│   ── ORF-composite path + score-table schema (research utilities) ──
├── orf/
│   ├── score_schema.py, score_gates.py, panel_manifest.py, profile_compare.py
│   ├── orfs_import.py (BED12→transcripts), map_orfs.py, assemble.py
│   └── rdg_flux_export.py   RDG-Flux per-position frame-posterior export
│
│   ── legacy single-sample path (kept for `process-bam`) ──
├── config.py           Config dataclass + validate_config
├── legacy_workflow.py  process-bam / zarr orchestration
├── core/               coordinates.py, orffinder.py, scoring.py (classic ORF scoring)
│
│   ── lower-level format readers & misc ──
├── file_handlers/      bam.py, bed.py, bigwig.py, sparse_parquet.py, zarr.py
├── utils/              io.py, logging.py
└── visualization/      plots.py, report.py (HTML), riboseq_profile_style.py
```

### A note on the `matrix_*.py` modules at the package root

`matrix_rollup.py`, `matrix_scoring.py`, `matrix_normalisation.py` and
`matrix_qc.py` sit at the package root (not in a subpackage). This is the
**matrix engine** — it's deliberately kept as flat top-level peers, alongside
the other top-level core modules (`clustering.py`, `frame_support.py`, etc.),
which is how the codebase was consolidated when the old `pipeline/` package was
removed. They are cohesive (one subject each) and import cleanly from `io/`,
`coverage/` and `scoring/`. Grouping them into a `matrix/` subpackage would be a
reasonable future tidy, but it is purely cosmetic and is intentionally not done
here.

---

## Documentation

Design/reference docs live in [`docs/`](docs/):

- [`architecture.md`](docs/architecture.md) — functional-core/imperative-shell design + module map.
- [`event_scoring_model.md`](docs/event_scoring_model.md) — the event model, dedup measurements, scoring spec.
- [`matrix_query_engine.md`](docs/matrix_query_engine.md) / [`matrix_normalisation_strategy.md`](docs/matrix_normalisation_strategy.md) — matrix querying + normalisation.
- [`RELEASE.md`](docs/RELEASE.md) — release machinery; [`REFACTOR_TASKS.md`](docs/REFACTOR_TASKS.md) — the migration record.

> **Housekeeping note:** `docs/` also contains ~35 `read_assignment_*.md` files
> (most are `read_assignment_ribometric_10pct.*.md` per-locus review notes from a
> May 2025 investigation). These are **research lab-notes, not user
> documentation**, and are effectively archival. Consider moving them to an
> `docs/notes/` or `archive/` subfolder (or pruning) so the docs index reflects
> only current reference material.

---

## Development

```sh
make gate        # the regression gate (golden + BAM-provider tests)
make test        # full test suite
make lint        # ruff + black --check
make format      # black + ruff --fix
make typecheck   # mypy (advisory)
```

## License

MIT.

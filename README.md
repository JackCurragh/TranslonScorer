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
pip install "TranslonScorer @ git+https://github.com/jackcurragh/translonscorer"
# optional extras: [bigwig] [zarr] [viz] [fastbam] [dev] or [full]
pip install "TranslonScorer[full] @ git+https://github.com/jackcurragh/translonscorer"
```

Core needs `pysam`, `polars`, `numpy`. `scipy` is optional (hierarchical
clustering of confusing loci). `pyBigWig` (the `[bigwig]` extra) is needed for
bigWig coverage and bigBed feature input.

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

The **feature source** (what to score) and the **coverage source** (the reads)
are independent — mix any of them: e.g. score GTF CDSs against your BAMs, or your
own BED12 ORFs against the matrix.

**You do not need a pre-built annotation database.** `extract-events` builds the
features to score from whatever you have — a GTF/GFF (score annotated CDSs), a
BED12/bigBed (your own ORF set), a FASTA (de-novo ORF finding), or the
annotation sqlite:

```sh
translonscorer extract-events --gtf anno.gtf --feature-type CDS --out-dir run/events  # annotated CDSs
translonscorer extract-events --bed12 orfs.bed                   --out-dir run/events  # your ORFs
translonscorer extract-events --fasta contigs.fa --start-codons ATG,CTG --out-dir run/events  # de-novo
translonscorer extract-events --sqlite anno.translons.sqlite     --out-dir run/events  # annotation DB
```

The end-to-end flow is four steps:

```
extract-events        feature source (GTF/BED/FASTA/sqlite) ─► genomic event store
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

Score annotated CDSs (from a GTF) against a handful of genome-aligned Ribo-seq
BAMs — no annotation database needed:

```sh
# one-shot: feature source = --gtf, coverage source = --bam
translonscorer pipeline \
  --gtf annotation.gtf --feature-type CDS \
  --bam sampleA.bam --bam sampleB.bam \
  --out-dir results/ \
  --offset-method global --global-offset 12

# results/report.parquet now has one row per translon with per-aspect calls
```

Swap the feature source freely — `--bed12 my_orfs.bed` to score your own ORF
set, or `--fasta contigs.fa` to find ORFs de-novo.

Transcriptome-aligned BAMs (reads positioned along the mRNA) are projected to
genome coordinates — UTRs included, so 5′UTR uORFs land correctly:

```sh
translonscorer pipeline \
  --gtf annotation.gtf \
  --bam tx_aligned.bam --transcriptome --annotation annotation.gtf \
  --out-dir results/
```

### B. The matrix (annotation scale)

The sparse matrix is **one logical matrix sharded into partition directories by
read sequence**. Reads for any locus are spread across *all* partitions, so the
matrix must **always be used in full** — there is no single-partition or subset
usage. You point at the matrix **root** and every partition under it is scanned
together.

Matrix scoring needs a **P-site index** built once per matrix. It supplies the
per-(sample, length) P-site offsets: a single flat offset across read lengths
smears the P-site across frames and collapses in-frame fraction to the ~0.33
random floor, so the index is required, not optional.

```sh
# once per matrix — calibrates offsets and bakes the index
translonscorer build-psite-index \
  --matrix-dir matrix/global_partitioned \
  --cds-gtf annotation.gtf --calibration-gtf annotation.gtf \
  --out-dir matrix/psite_index

translonscorer pipeline \
  --sqlite annotation.translons.sqlite \
  --matrix-dir matrix/global_partitioned \
  --psite-index matrix/psite_index \
  --chrom chr12 \
  --out-dir results_chr12/ \
  --data-version matrix_v1
```

`--chrom` restricts which **events** are extracted (a tractable scope); it does
not restrict the matrix — all partitions are still read. Omit it for the whole
genome. See [`runs/run_matrix_scoring.sh`](runs/run_matrix_scoring.sh) for a
reproducible matrix run you can copy.

BAM/bigWig mode needs no index — it calibrates its own offsets per file (see
`--offset-method`). That is the only difference between the two modes: the
events, the scorer, the report and the store are identical.

### C. Step by step (full control)

```sh
# 1. events (once per feature set; reused across all scoring runs).
#    Source can be --gtf / --bed12 / --bigbed / --fasta / --sqlite.
translonscorer extract-events --gtf annotation.gtf --feature-type CDS \
    --out-dir run/events --chrom chr12

# 2a. score from BAMs ...
translonscorer score-bams  --events-dir run/events --bam a.bam --bam b.bam \
    --store-dir run/scores --data-version v1
# 2b. ... or from the whole matrix (all partitions under the root).
#     --psite-index is required; build it once with build-psite-index.
translonscorer score-matrix --events-dir run/events \
    --matrix-dir matrix/global_partitioned \
    --psite-index matrix/psite_index \
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

Reference/design docs live in [`docs/`](docs/):

- [`architecture.md`](docs/architecture.md) — functional-core/imperative-shell design + module map.
- [`event_scoring_model.md`](docs/event_scoring_model.md) — the event model, dedup measurements, scoring spec.
- [`matrix_query_engine.md`](docs/matrix_query_engine.md) / [`matrix_normalisation_strategy.md`](docs/matrix_normalisation_strategy.md) — matrix querying + normalisation.
- [`rdg_flux_export.md`](docs/rdg_flux_export.md) — RDG-Flux export format.
- [`RELEASE.md`](docs/RELEASE.md) — release machinery.

(Research lab-notes and internal tracking are kept in a git-ignored `notes/`
directory, not in `docs/`.)

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

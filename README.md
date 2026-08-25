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
| 1–N bigWigs (coverage only) | `BigwigSetProvider` | `score-bigwig` |
| annotation-scale sparse matrix | `MatrixProvider` | `score-matrix` |

All three feed one scoring core and produce the same `fact_event_score` store
and per-translon report, and all three run under `pipeline` as a single call.

**bigWig is the simplest source and the most lossy.** No offsets to calibrate
and no index to build, but per-base depth carries no read-level splice
information — so **junction events cannot be scored at all** and come back
`INSUFFICIENT` with a null call (unassessable, *not* unsupported). The run warns
with the count. Score junctions with `--bam` or matrix mode.

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
score-bams / score-bigwig / score-matrix   events + coverage ─► score store
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

### B. Coverage tracks (bigWig)

If all you have is coverage, pass bigWigs. Ribo-seq is stranded, so prefer the
forward/reverse pair — an unstranded track mixes +/- signal at every position:

```sh
# one-shot: feature source = --gtf, coverage source = stranded bigWig pair
translonscorer pipeline \
  --gtf annotation.gtf --feature-type CDS --chrom 12 \
  --forward-bigwig fwd.bw --reverse-bigwig rev.bw \
  --sample mysample \
  --out-dir results/
```

Repeat `--forward-bigwig`/`--reverse-bigwig` for more samples (Nth pairs with
Nth). `--bigwig` takes unstranded tracks instead. Junction events are not
scorable from bigWig — see the caveat above.

### C. The matrix (annotation scale)

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

### D. Step by step (full control)

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

The spine is one line of dependency:

```
sources → events → provider.coverage() → scorer → report → consequentiality → store
```

Everything else hangs off it or sits beside it.

```
TranslonScorer/
├── cli.py              Click CLI (top level = the workflow; research/ legacy groups)
├── workflows.py        orchestration shell: pipeline / extract / score / report
├── model.py            dataclasses: Region, OffsetParams, ScoreThresholds,
│                       FrameSupportParams, ConsequentialityPolicy
│
│   ── the spine (pure, no I/O) ──
├── events.py           features → deduplicated genomic events + membership
├── io/feature_sources.py  GTF / BED12 / bigBed / FASTA / sqlite → features
├── offsets.py          P/A-site offset calibration (global | file | metagene)
├── scoring/            the ONE event scorer
│   ├── aspects.py        init / elongation / termination / junction scorers
│   ├── evidence.py       evidence → eligibility/call decisions (thresholds)
│   └── run.py            score_events() — batched elongation + per-event loop
├── report.py           compose per-translon report from event scores
├── consequential.py    dynamic consequentiality policy (no hard gates)
│
│   ── coverage providers: the ONE point where sources differ ──
├── coverage/
│   ├── bam.py            BamSetProvider (1–20 BAMs; offsets; transcriptome→genome)
│   ├── bigwig.py         BigwigSetProvider (coverage-only, lossy by design)
│   ├── profile.py        pure P/A-site profile from reads + offset table
│   ├── transcriptome.py  transcript→genome projection (exon walk)
│   ├── base.py           MAPPABILITY_LEDGER_SCHEMA + the duck-typing convention
│   └── profiles.py, locus_profiles.py, locus_features.py, transcript_coords.py,
│       mapped_index.py, index_from_bam.py, junctions.py
│                         a separate cohort profile BUILDER (the `profiles` and
│                         `index-from-bam` commands) — not the scorer's providers
│
│   ── the sparse-matrix engine ──
├── matrix/
│   ├── provider.py       MatrixProvider (requires a P-site index)
│   ├── psite_index.py    build-psite-index + query_genomic_coverage
│   ├── rollup.py         partition scanning; FrameRollup; junction tabulation
│   ├── qc.py             per-sample length/periodicity QC off the matrix
│   ├── normalisation.py  per-locus matrix normalisation
│   └── clustering.py     profile clustering for confusing/contended loci
│
│   ── I/O adapters ──
├── io/
│   ├── annotation.py     GTF → CDS / exon blocks
│   ├── annotation_bundle.py  build/load an annotation bundle directory
│   ├── bam.py            read-id parsing, chrom normalisation, CIGAR, junctions
│   ├── matrix.py         sparse-matrix partition manifests / readers
│   ├── store.py          event + fact_event_score Parquet store read/write
│   └── inspect.py        Parquet inspection helper (the `inspect` command)
│
│   ── in-flight: a read/position index (nothing reads it yet) ──
├── alignments/         builder.py, schema.py, provenance.py  (one row per alignment)
├── counts/             builder.py, schema.py      (5′-end + junction counts per sample)
│
│   ── research: frame assignment (the `research` command group) ──
├── frame_support.py    frame-posterior estimation (linear/latent/HMM)
├── frame/
│   ├── bleed.py, deblur.py, hmm.py, latent.py   frame-correction models
│   ├── validation.py             posteriors vs annotated CDS frame
│   ├── rdg_flux_export.py        RDG-Flux per-position frame-posterior export
│   └── read_assignment.py        unique/fractional/EM/frame-aware assignment
│
│   ── ORF table handling (the `legacy` command group) ──
├── orf/
│   ├── orfs_import.py (BED12→transcripts), map_orfs.py
│   └── score_schema.py, panel_manifest.py
├── assemble.py         overlap resolution over scored candidates
│
│   ── single-sample path (kept for `process-bam` / `profiles`) ──
├── config.py           Config dataclass + validate_config
├── legacy_workflow.py  process-bam / zarr orchestration
├── qc.py               periodicity / frame-dominance QC (pure)
├── core/               coordinates.py, orffinder.py
│                       (orffinder is NOT legacy — it is the shared ORF-from-
│                        sequence rule engine behind feature_sources.from_fasta)
│
│   ── lower-level format readers & misc ──
├── file_handlers/      bam.py, bed.py, bigwig.py, sparse_parquet.py, zarr.py
└── utils/              io.py, logging.py
```

### Reading the CLI

`translonscorer --help` shows the workflow plus build/inspect infrastructure.
Two groups hold everything that is not on the scoring path:

- **`research`** — `compare-read-assignment`, `validate-panel`,
  `export-rdg-flux`. Nothing under `frame/` is imported by the scorer; these
  support the frame-assignment analyses.
- **`legacy`** — `orfs-import`, `assemble`, `map-orfs`. ORF table handling;
  `orf/` is not imported by the spine.

Both still run (`translonscorer research export-rdg-flux --help`). They are
grouped, not deprecated — the split exists so the top level reads as the
product.

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

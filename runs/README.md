# Reproducible matrix event-scoring runs

End-to-end demonstration of the new TranslonScorer event-scoring pipeline driven
entirely by the CLI — no bespoke scripts, just the four subcommands chained in
[`run_matrix_scoring.sh`](run_matrix_scoring.sh).

## What it does

```
extract-events  annotation sqlite        ->  genomic event store (deduplicated)
score-matrix    events + sparse matrix   ->  append-only fact_event_score store
report          scores + events          ->  per-translon report + consequentiality
consequential   report (strict re-gate)  ->  full-chain-supported subset
```

## Run it

```bash
# defaults: chr12, first 16 matrix partitions (~30 s)
./run_matrix_scoring.sh

# a different chromosome, full cohort depth (all partitions)
CHROM=chr7 N_PARTITIONS=0 ./run_matrix_scoring.sh
```

Configurable via environment (see the CONFIG block in the script): `SQLITE`,
`MATRIX_DIR`, `CHROM`, `DATA_VERSION`, `N_PARTITIONS`, `STRICT_TIER_CONFIDENCE`,
`OUTDIR`.

`N_PARTITIONS` is a speed/depth knob: the matrix is partitioned by read sequence,
so a subset is a lower-depth read sample of the whole cohort. Use `0` for full
depth.

## Outputs (under `matrix_<chrom>/`, git-ignored)

| path | contents |
|---|---|
| `events/` | event store: `events/`, `feature_event/`, `event_overlap/` (one Parquet per chrom) |
| `scores/` | append-only `fact_event_score` store, partitioned `data_version=…/tier=…` |
| `report.parquet` | one row per translon: per-aspect calls/metrics, `consequentiality_score`, `consequential` |
| `consequential_strict.parquet` | report re-gated at `min-tier-confidence` (full translation chain supported) |

## Reading the report

`report.parquet` columns include, per chain aspect, `{init,elongation,term}_call`
(SUPPORTED/UNSUPPORTED), `_metric`, `_n_reads`, `_supported_frac`; plus
`junction_n`/`junction_supported_frac`, `total_reads`, `tier_confidence`,
`expression_percentile`, `consequentiality_score`, and the boolean `consequential`.

Sort by `consequentiality_score` for the most-supported, most-translated
translons. The default `report` policy uses zero floors (nothing hard-gated);
`consequential` with `--min-tier-confidence` / `--min-expression-percentile`
re-gates an existing report without rescoring.

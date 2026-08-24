# Matrix query engine — optimality observations & a generalisable API

Status: design note (2026-06-12). Backed by benchmarks on the 30-sample
`global_partitioned` matrix (PRJNA604580). Companion code:
`TranslonScorer/pipeline/matrix_qc.py`; demos `scripts/demo_frame_dominance.py`,
`scripts/bench_*.py`.

Related design note: `docs/matrix_normalisation_strategy.md` defines the
normalisation contract for locus-level sample clustering, cluster aggregation,
and handoff to event-level scoring.

---

## 1. The matrix and what a query costs

The sparse-Parquet matrix stores each **unique read once** (companion BAM:
`read_id → genomic locus`) and **per-sample multiplicities sparsely**
(`read_id → sample_id → count`), partitioned by read-sequence prefix
(`read_bucket`). It deduplicates aggressively: 65.4M total reads → 3.8M unique
(~17× sharing, 0.61 duplicate rate) across this 30-sample cohort.

A periodicity/QC query over the matrix has **two cost phases**:

| phase | work | scales with | original impl |
|---|---|---|---|
| **A — annotate** | scan unique-read BAMs, assign CDS frame per read | `N_unique` (sample- & query-independent) | pysam Python loop over ~54M alignments |
| **B — tabulate** | join counts, fan out to per-sample frame tallies | `nnz` (reads × samples present) | Python `iter_rows` dict accumulation over ~154M rows |

### Measured cost model (256 partitions, this matrix)

Two data points (1 sample = 722 s, 30 samples = 3120 s) solve to:

```
cost ≈ 640 s  (Phase A floor, fixed)  +  83 s/sample  (Phase B, Python loop)
```

So the original `periodicity_qc` is **not** flat in sample count — the Python
accumulation loop costs ~83 s per sample. That single loop, not the BAM scan,
dominates at cohort scale.

### Comparison to per-sample RiboMetric

- RiboMetric on one transcriptome BAM (SRR11005875, 8 threads): **503 s/sample**.
- Old matrix path, 30 samples: 3139 s → **4.8× faster** than 30× RiboMetric,
  but asymptotic ceiling only ~6× (503/83) because the per-sample marginal is
  so heavy.
- Equivalence is sound: matrix vs RiboMetric per-read-length agree at
  **dominance r = 0.92, periodicity r = 0.93**, same 29 nt frame-0 peak.

---

## 2. The optimisation: split static index from dynamic tabulation

Phase A is a pure function of `(unique reads, annotation, offset)` — it does not
depend on which samples you select, nor on aggregate/cluster/individual
granularity. **Pay it once; cache it.**

```
build_read_loci_index(partitions, cds_df)  ->  read_loci.parquet
    read_id | length | frame            # (+ extensible: chrom, pos, tran_id, biotype…)

tabulate_frames(partitions, read_loci, group_level=…)   # pure Polars, streaming
    counts.scan_parquet()
      .join(read_loci, on=read_id)       # attribute attach
      .join(sample_map, on=sample_id)
      .group_by(GROUPKEY, length, frame) # GROUPKEY = ∅ | sample | cluster
      .agg(count.sum())
      .collect(engine="streaming")
```

`GROUPKEY` is the only thing that changes between **aggregate**, **cluster**,
and **per-sample** — one materialised `counts ⋈ read_loci` scan answers all
three. Swap `frame` for `pos`/`tran_pos` and the identical pipeline yields the
positional **profiles** used for clustering. This is the
"aggregate/cluster/individual in one pass" primitive.

### Measured impact (6-partition equivalence + micro-benchmark)

- Fast path is **bit-identical** to the old loop (max count diff = 0 / 1524 cells).
- Phase B tabulation: ~20 s → **0.2 s** (~100×, Python loop → Polars group-by).
- Repeat query (index cached): 54 s → **0.2 s** (~270×).

Result: first run amortises Phase A once (~640 s floor, isolated and
parallelisable); **every subsequent query is seconds**, and the per-sample slope
collapses from 83 s to native — the matrix becomes ~flat in sample count, as the
architecture implies it should be.

### Full-scale numbers (256 partitions, 30 samples, index cached)

| query | old | new (warm) | speedup |
|---|---|---|---|
| `periodicity_qc` (30 samples) | 3120 s | **7.7 s** | **~405×** |
| tabulate 1 sample | — | 6.9 s | |
| tabulate aggregate | — | 7.7 s | |
| tabulate 30 samples | (3120) | 9.7 s | |
| index build (once, 35.6M reads) | — | 514 s | |

New cost model: `514 s (one-time build) + ~7 s + 0.1 s/sample`. The per-sample
marginal fell from **83 s → ~0.1 s** (1 sample 6.9 s vs 30 samples 9.7 s) — flat
in sample count, confirming the architecture. vs per-sample RiboMetric
(503 s × 30 = 15,090 s) the warm query is **~1550×** faster; the index is reused
across every clustering / feature-scoring pass.

### Phase-A wins

1. **Parallelise partitions** — DONE. `build_read_loci_index(n_workers=…)` scans
   the independent partitions across processes. **514 s → 50.6 s (10.2×)** at
   8 workers on the full 256-partition / 35.6M-read build, **bit-identical** to
   the serial baseline. Per-partition work is the exact serial
   `_scan_partition_loci`, so the result is identical regardless of worker count.
   Combined with the ~8 s cached tabulate, a full *cold* periodicity run is now
   ~60 s vs the original 3120 s (~52×); warm queries stay ~8 s. Uses the **spawn**
   start
   method: `fork` deadlocks because polars' Rust threadpool is already live in
   the parent. Caveat: spawn re-imports the entry module per worker, so a script
   that calls this must guard execution under `if __name__ == "__main__":`.
2. **oxbow** BAM→Arrow read instead of the pysam Python loop — pending, and now
   *gated on the unique-mapper change below*: oxbow's row order differs from
   pysam's, which only matters because of order-dependent multimapper handling.
   Once we keep unique mappers only, every retained read has exactly one
   alignment, order is irrelevant, and oxbow becomes a safe drop-in.
3. **Vectorised interval frame-join** (Polars `join_where`) instead of the
   per-partition Python sweep, so frame can be (re)assigned at query time for
   any offset without an index rebuild.
4. Store `pos`/`length`/`strand` (fully static) in the index and derive `frame`
   at query time — changing offset must not force a rebuild.

---

## 2-impl. Status: implemented in `pipeline/matrix_rollup.py`

§2a (multimapping) and §2b (offset-independent map-reduce rollup) below are
**implemented and validated**:

- `build_alignment_index(…, ref_offset, n_workers)` — per in-CDS alignment
  `(read_id, length, strand, phase0, gene_code, n_align)`, all candidate loci
  retained (no last-wins), parallel build **bit-identical to serial** (3.76M
  rows / 16 partitions verified).
- `tabulate_rollup(…, multimap_mode="unique"|"gba")` — vectorised reduction to
  `(sample, length, strand, phase0) → count`.
- `calibrate_offsets()` / `rollup_to_periodicity()` — per-(sample, length)
  offset calibrated **centrally on the reduced rollup** (frame is analytic in
  offset), then applied; output matches `periodicity_qc` schema.
- `periodicity_qc_scalable(…)` — end-to-end wrapper.

Validation: **unique-mode vs stock RiboMetric dominance r = 0.924**; GBA
conserves counts (unique 11.2M → GBA 14.9M = unique + redistributed multimapper
counts) and rescues periodic multimappers (SRR11005875 29 nt: 100k→134k reads at
identical 0.93 dominance); central calibration picks per-sample offsets (29 nt →
12, not the default 15). `build_cds_blocks` now carries `gene_id` and lives in
the library (the old scalar-CDS `getexons_and_cds` is no longer on this path).

Both former pending items are now **done**:

- **oxbow reader** (`reader="oxbow"`): vectorised BAM→Arrow scan replacing the
  per-record pysam loop. Coordinate mapping verified (`reference_start = pos-1`,
  `reference_end = end` — oxbow's `end` already accounts for CIGAR N/D, so no
  CIGAR parsing). **Set-identical to pysam** in unique mode (order no longer
  matters), ~1.5× faster on the scan; compounds with the 10× parallelism.
- **Fold-in** (`build_matrix_rollup`): the count-join is done *inside* each
  worker, which reduces its own partition's nnz to small `(sample,length,strand,
  phase0)` + `(sample,gene)` density partials. The parent sums the partials; for
  GBA the (16%) multimapper rows are spilled to per-partition parquet shards and
  weighted in a single streaming pass against the reduced global density. The
  full nnz join is **never materialised in one place** — this is the 6k path.
  Verified identical to the two-step `tabulate_rollup` for both unique and GBA.
  `build_alignment_index`/`tabulate_rollup` are kept for read-level / region
  queries; `build_matrix_rollup` is the QC/cohort path.

**Full-scale folded benchmark** (256 partitions, 30 samples, oxbow, 8 workers,
end-to-end scan + count-join + reduce in one pass):

| mode | rollup build | rollup rows | total count | cohort 29 nt dom |
|---|---|---|---|---|
| unique | **52.8 s** | 3 379 | 228 M | 0.907 |
| gba | **60.9 s** | 3 408 | 305 M (+77 M rescued) | 0.909 |

vs the original `periodicity_qc` (3120 s) this is **~59×**, in a single pass with
no separate index materialisation, and calibrate+score on the ~3.4k-row rollup
is instant. GBA rescues ~35% more periodic reads (SRR11005875 29 nt 2.1M→2.84M)
at identical 0.92 dominance. The 52.8 s is the full nnz processed in sharded
parallel; per-query cost afterwards is O(rollup), independent of cohort size.

Still pending: run GBA pass-2 as a second *parallel* worker pass (each shard is
self-contained per read, so weighting can parallelise once global density is
broadcast) — currently a single streaming pass, adequate to mid cohorts.

## 2a. Multimapping strategy (replaces "last-wins")

The current scan does `read_pos[rid] = …` — **last alignment in BAM order wins**
for any read with multiple loci (8.5% of read_ids here, up to 10 loci each, plus
~18% secondary/supplementary records). This is arbitrary and biologically
uninformed, and it is the sole reason the scan is row-order-dependent.

Supported strategies to implement:

1. **`unique` (default).** Keep only reads with exactly one alignment
   (`NH == 1` / single primary record); drop multimappers entirely. Deterministic,
   conservative, defensible — and order-independent, so it also unblocks oxbow.
2. **`guilty_by_association` (rescue / EM, opt-in).** Distribute a multimapper's
   count fractionally across its candidate loci, weighted by the local
   **unique-read coverage** prior at each locus (cf. RSEM / salmon / ribotricer).
   For frame/periodicity this means a multimapper contributes a weighted mixture
   to each locus's `(length, frame)` tally rather than 1.0 to one arbitrary locus.
   Requires retaining all candidate loci per read + a unique-coverage pass first.

"last-wins" should be removed, not kept as an option.

---

## 2b. Scaling to large cohorts (6k+ samples)

The 30-sample numbers above hide a Phase-B cliff: `tabulate_frames` scans `nnz`
((read × sample) non-zeros), which grows ~linearly with sample count. At ~30
samples nnz ≈ 154M and the streaming group-by is ~8 s; at **6,000 samples** nnz
is ~200× larger (billions of rows) — minutes per query and memory-bound. The
"~0.1 s/sample" flatness held only because the *cached-index tabulate* was still
small; it does not survive 200× more nnz.

The output, though, stays tiny, and the frame distribution is an **associative
reduction** (sum over partitions), so the fix is to **fold tabulation into the
build as a map-reduce**.

**Offset must not be baked in.** Offset is per-(sample, read length) and a
sample's reads are spread across every partition, so no single worker can
calibrate it. But frame is *linear* in offset:

```
frame = (phase0 + offset) % 3        phase0 = (cds_phase + 5'pos) % 3   (sign flips on − strand)
```

`phase0` (frame at offset 0) depends only on the read's genomic position and the
annotation — **sample- and offset-independent**, computable per partition with
no cross-partition info. So the rollup is keyed on `phase0`, not on a
pre-offset frame:

- Each worker scans its partition, computes `(read_id, length, strand, phase0)`
  (offset 0), joins *that partition's* counts, and emits
  `(sample, length, strand, phase0) → partial_count`.
- The parent **sums** the partials → `frame_rollup.parquet`
  (≈ 6k × 15 × 2 × 3 ≈ 540k rows).
- **Offset calibration runs centrally on the reduced rollup**, after partitions
  are summed — so it sees every partition for each sample. Pick per-(sample,
  length) offset that puts the dominant frame at 0, then apply analytically
  `frame = (phase0 ± offset) % 3` and collapse strand.

This is strictly better than the current code (which bakes a uniform `offset=15`
per partition and nudges post-hoc): offsets become first-class and
per-sample/per-length, and the index is **offset-independent** — changing
offsets never forces a rebuild.

The full nnz scan then happens **once at build time, distributed and never
materialised whole** (bounded memory regardless of nnz). Afterwards, periodicity
QC for *any* sample subset (1 or 6,000) is a filter + analytic offset shift over
the ~540k-row rollup — **milliseconds, independent of nnz**. Per-query cost
becomes a function of the *output* size (groups × attributes), not the matrix
size.

*Caveat:* `(phase0 + offset) % 3` is exact for reads solidly inside a CDS
interval; for reads within ~`offset` nt of an exon boundary or start codon the
shifted A-site can cross into a different-phase interval or out of CDS, so the
analytic frame is slightly off for that small boundary-proximal fraction —
negligible for statistical periodicity (offset search range is only ±a few nt).
`read_loci` retains the genomic position if an exact per-read recompute is ever
needed.

This promotes the engine to three artefacts:

| artefact | grows with | built | queried |
|---|---|---|---|
| `read_loci` (read → attrs) | unique reads | scan, once | region-scoped read-level queries |
| `frame_rollup` (group × attr → metric) | groups (samples/clusters) | associative reduce during scan | QC / per-sample / cluster — O(output) |
| counts (sparse) | nnz | matrix build | only via pushed-down region queries |

Read-level queries that genuinely need per-read resolution (positional profiles
for clustering specific loci) stay tractable at 6k because they are
**region-scoped** — predicate pushdown bounds the nnz touched. The unbounded
all-samples × all-reads scan only ever runs at build time.

---

## 3. Generalising: a track-matrix query engine

The pattern above is not RiboSeq-specific. Any **sparse multi-sample track
matrix** — expression, coverage, methylation, ATAC, RiboSeq frame — is:

```
features stored once (id → coordinates/attributes)   ×   counts[feature, sample]
```

and every analytical query is the same three-operator pipeline:

```
SELECT   features        (all | region | id-set)         → predicate pushdown
ANNOTATE feature → attrs (frame, position, biotype, …)   → cached static index
GROUP    by level        (aggregate | cluster | sample | study)
REDUCE   metric          (sum | mean | frame_dist | periodicity | entropy | …)
```

### Engine design principles (for hyper-optimised querying)

1. **Static/dynamic split.** Feature-attribute indices are built once and
   versioned next to the matrix; queries are joins. Never recompute annotation
   per query.
2. **LazyFrame all the way to `collect`.** Predicate + projection pushdown,
   streaming, multithreaded; zero Python row loops. Region/sample/`read_bucket`
   filters push to the scan.
3. **Group key is a parameter, not a code path.** aggregate/cluster/sample/study
   are one operator with different keys; compute several at once via `collect_all`
   sharing the subplan (one scan, many granularities).
4. **Composable reductions.** Each metric is a Polars expression over the grouped
   frame, registered in a table. New metric ⇒ new expression, not a new scan.
5. **Materialised common aggregates.** Study- and cluster-level rollups cached so
   hot queries skip even the join.
6. **Cost guards.** Streaming + region pushdown bound worst-case so a naive
   caller can't scan the whole matrix unintentionally.

### Public API / MCP surface (track data → LLMs)

A thin verb set over the engine, track-type agnostic:

- `query(features, cohort|samples, level, metric)` — the general primitive.
- `expression(genes, samples)` / `profile(region, level)` — sugar over `query`.
- `compare(groupA, groupB, metric)` — differential at any granularity.
- `qc(samples)` — per-sample assessment (periodicity, length dist, dominance).

Only the **attribute registry** (how a feature maps to coordinates/attributes)
and the **reduction registry** differ per track type; the SELECT/ANNOTATE/GROUP/
REDUCE engine and its caching are shared. This makes a single backend able to
expose many heterogeneous track datasets to an LLM/MCP client with predictable,
cache-warm latency.

---

## 4. Open items

- Move `build_cds_blocks` (currently in `scripts/demo_frame_dominance.py`) into
  the library as the canonical CDS-block builder — `getexons_and_cds` returns
  *scalar* CDS bounds and silently breaks frame assignment (see
  `project-periodicity-qc-deconfliction` memory).
- Apply Phase-A wins (oxbow + parallel partitions) and re-benchmark the build.
- Generalise `tabulate_frames` to a positional `tabulate_profiles` (group by
  `pos`) to feed `score_clustered` / profile clustering directly.
- Persist `read_loci.parquet` as a first-class matrix sidecar with a version/
  offset/annotation provenance stamp.

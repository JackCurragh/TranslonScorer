# TranslonScorer architecture

Design target for the tool: a **functional core** (pure transforms over Polars
DataFrames + frozen dataclasses) with **state confined to coverage providers**
(the only place that owns file handles, computed offsets, and cached indices).
Everything else is `data → function → data`.

## Principles

- **Pure functions** for all transformation logic (events, offsets, profiles,
  scoring, consequentiality, reporting). Same inputs → same outputs; trivially
  testable; this is what made the golden/`verify_vectorised` gates work.
- **Frozen dataclasses** for config and immutable records (`OffsetParams`,
  `ScoreThresholds`, `ConsequentialityPolicy`, `Region`, score rows). No hidden
  mutable defaults.
- **Stateful objects only where warranted** — exactly one kind: `CoverageProvider`
  subclasses. They legitimately hold open BAM/bigwig handles, per-BAM computed
  offsets, and cached indices, and expose a pure interface so the rest stays
  functional. I/O adapters open/close within a call (no persistent state).

## Layered module map

```
TranslonScorer/
  model.py          dataclasses: Region, Event, Feature, OffsetParams,
                    ScoreThresholds, ConsequentialityPolicy, Tier, ScoreRecord
  io/               adapters — open→read→close, return DataFrames (pure-ish)
    bam.py            alignments via oxbow/pysam; CIGAR blocks; read_id parse
    bigwig.py         bigwig coverage
    matrix.py         sparse matrix: manifest, counts, unique-read BAM
    annotation.py     GTF/BED → cds_blocks, gene_spans
    store.py          parquet read/write: events, scores, ledger, profiles
  events.py         extract_events, frame intervals, contention, feature_event   [pure]
  offsets.py        metagene | file | global offset calc + plausibility           [pure]
  coverage/
    base.py           capability PROTOCOLS (not one fat ABC) — see below          [interface]
    profile.py        reads → P/A-site coverage, fold/merge, size factors         [pure]
    bam.py            BamSetProvider     (state: handles, per-BAM offsets)        [stateful]
    bigwig.py         BigwigSetProvider                                           [stateful]
    matrix.py         MatrixProvider                                              [stateful]
  scoring/
    aspects.py        init/elong/term/junction scorers + identifiability          [pure]
    evidence.py       evidence builders, eligibility/call (thresholds)            [pure]
    run.py            score_events / vectorised batch orchestration               [pure]
  clustering.py     translation-mode clustering (cluster tier)                    [pure]
  consequential.py  tier × confidence × expression × context → consequentiality   [pure + policy]
  report.py         per-feature composition, annotation-set stratification        [pure]
  qc.py             periodicity / frame-dominance (RiboMetric-equivalent)         [pure]
  workflows.py      end-to-end compositions (score-bams, score-matrix, …)
  cli.py            thin click wrappers over workflows
```

Data flow (all arrows are pure functions except the provider boundary):
```
annotation ──events.py──▶ events + feature_event + contention
                                   │
features/BAMs/bigwigs ─CoverageProvider─▶ P/A-site coverage + junctions + size_factors + mappability
                                   │
                          scoring/run.py ──▶ score records (fact_event_score)
                                   │
              report.py ⨝ feature_event ──▶ per-feature evidence
                                   │
              consequential.py ──▶ ranked / stratified sets
```

## Where each wishlist item lives

| functionality | home | kind |
|---|---|---|
| 1–20 BAMs / bigwigs / matrix inputs | `coverage/{bam,bigwig,matrix}.py` (one ABC) | stateful providers |
| offset calc: metagene \| file \| global | `offsets.py` (pure) → provider holds chosen offsets | pure + state |
| **per-BAM, per-read-length offsets before aggregation** | `BamSetProvider`: compute offsets per BAM → per-sample P/A profile → *then* fold/merge | provider |
| plausible-offset constraint (25 nt ≠ 18 offset) | `offsets.plausible_offset_range` | pure |
| P-site (init) / A-site (elong, term) | `coverage(site=…)` param; scorers request site; `A = P+3` | pure |
| genome-native + transcriptome→genome (isoform multimap) | `io/bam.py` projection + `BamSetProvider` | adapter + provider |
| multimapper unique / GBA / copy-number | read selection in `coverage/*`; ledger output | provider |
| per-event mappability ledger (2nd pass) | `coverage/*` emits → `io/store.py`; relaxation in `consequential.py` | provider + policy |
| event dedup (genomic events) | `events.py` (extract_events) | pure |
| per-aspect scoring + identifiability | `scoring/aspects.py` | pure |
| aggregate default; cluster on confusing loci | **trigger** (is locus confusing?) in `scoring/aspects.py` next to identifiability; **response** (cluster+rescore) in `clustering.py`; provider `by_sample` | pure |
| size-factor expression normalisation | `coverage/profile.size_factors`; used in `consequential.py` | pure |
| consequentiality (tier/confidence/expr/context) | `consequential.py` + `ConsequentialityPolicy` | pure + policy |
| per-feature reports + annotation sets | `report.py` | pure |
| periodicity/frame QC | `qc.py` | pure |

## Migration: where the demoed code goes

> **STATUS — COMPLETE (T15, 2026-06).** `pipeline/` has been fully removed; every
> module now lives in the new tree. Final homes: matrix engine + QC at top level
> (`matrix_rollup.py`, `matrix_scoring.py`, `matrix_normalisation.py`,
> `matrix_qc.py`); frame subsystem in `frame/`; profiles/coords/indexing in
> `coverage/`; the (deprecated) ORF-composite path + score-table schema in `orf/`;
> `config.py`/`legacy_workflow.py` at top level; annotation/inspect adapters in
> `io/`. No re-export shims remain. The table below is kept as historical record
> of the intended mapping.

| current (demo) | →  target | notes |
|---|---|---|
| `pipeline/event_extract.py` | `events.py` + `io/annotation.py` (`build_cds_blocks`, `build_gene_spans`) | already mostly pure |
| `pipeline/event_score.py` | `scoring/aspects.py` + `scoring/evidence.py` + `scoring/run.py`; `ScoreThresholds`→`model.py`; `persist_scores`→`io/store.py`; prefix-sum helpers→`coverage/profile.py` | the v1 scorers |
| `pipeline/matrix_rollup.py` | `region_coverage`/`tabulate_*`→`coverage/matrix.py` + `coverage/profile.py`; oxbow reader→`io/bam.py`; junction logic→`coverage`; `profile_matrix`→`report.py`/`clustering.py` | split by concern |
| `pipeline/matrix_qc.py` | `qc.py`; `_build_frame_intervals`/`_deconflict_intervals`→`events.py`; manifest/count helpers→`io/matrix.py` | frame intervals are annotation-derived |
| `pipeline/coverage_providers.py` (skeleton) | `coverage/base.py` + split into `coverage/{bam,bigwig,matrix}.py`; `offsets`→`offsets.py` | the new home |
| `pipeline/profiles.py` (`_compute_asite_profiles`) | `coverage/profile.py` | profile gen |
| `pipeline/profile_clustering.py`, `score_clustered` | `clustering.py` | cluster tier |
| `pipeline/feature_metrics.py`, `orf_composite.py` | **deprecated** (superseded by events+scoring) | old composite path |
| `file_handlers/{bam,bigwig,sparse_parquet}.py` | `io/{bam,bigwig,matrix}.py` | consolidate |
| `scripts/*` (run/score/refine/expression/finalize) | `workflows.py` functions + `cli.py` subcommands | scripts become first-class |
| `scripts/verify_vectorised.py`, GAPDH golden | `tests/` (pytest regression gates) | CI-enforced |

## Capability protocols (resolved: NOT one fat ABC)

The three providers have genuinely different capabilities — a single ABC forces
`junction_support` to raise on bigwig and `site` to be silently ignored. Instead,
`coverage/base.py` defines **runtime-checkable Protocols** so capability is in the
type system, not in `NotImplementedError`:

```python
class CoverageProvider(Protocol):          # minimal — all sources
    def coverage(self, regions, *, by_sample=False) -> pl.DataFrame: ...
    def size_factors(self) -> dict[str, float]: ...

class SupportsSites(Protocol):             # BAM, matrix (P/A distinction)
    def coverage(self, regions, *, site: str = "A", by_sample=False) -> pl.DataFrame: ...

class SupportsJunctions(Protocol):         # BAM, matrix (needs CIGAR)
    def junction_support(self, junctions, *, by_sample=False) -> pl.DataFrame: ...

class SupportsMappability(Protocol):       # BAM, matrix
    def mappability_ledger(self, events) -> pl.DataFrame: ...
```

Scoring code requiring junctions is typed `SupportsJunctions` → a bigwig provider
is a **static type error if passed**, not a runtime crash. BAM/matrix satisfy all
four; bigwig satisfies only `CoverageProvider`. Init scorer requires
`SupportsSites`; elongation/term need only `coverage` (A-site default).

## Migration order (strangler-fig — the order is the whole game)

The new tree and `pipeline/` **coexist** until each piece is proven; flip
module-by-module, golden-gate green at every step, leaf-deps first:

0. **Make the gate one command.** Convert `verify_vectorised.py` + the GAPDH/chr12
   golden into `tests/test_golden.py` (pytest). `make gate` must be GREEN on
   current `pipeline/` code before anything moves.
1. **Leaf deps:** `model.py` (dataclasses, zero deps) + `io/` adapters. Gate green.
2. **Pure transforms, bit-identical:** `events.py`, `offsets.py`, `scoring/*`,
   `consequential.py`, `report.py`, `qc.py`, `clustering.py` — moved but producing
   byte-identical output (GAPDH golden Δ=0). `pipeline/` re-imports from the new
   homes (shim) so nothing breaks. Providers untouched.
3. **Coverage providers last** (the only new *behaviour*): `coverage/base.py`
   protocols, `profile.py`, then `MatrixProvider` must reproduce the GAPDH golden
   *through the provider path* before `BamSetProvider`/`BigwigSetProvider` (new
   behaviour → new tests, not the old golden).
4. **Orchestration + cleanup:** `workflows.py`, `cli.py` (deprecate
   `orf-composite`/`score-orfs`/`feature-metrics`); delete `pipeline/` shims once
   green.

Rule: **never advance a step with a red gate.** Phase 1–2 are pure moves (Δ=0
required); only Phase 3 introduces behaviour and earns new tests.

## Two decisions noted

1. **Offset calc reliability is a test, not an assumption.** `offsets.py` supports
   metagene/file/global; a `tests/` benchmark must compare them on a truth set
   (which method's per-length offsets maximise downstream CDS in-frame /
   periodicity, or best match known reading frame) and *select the default*. We do
   not pick a canonical method by assertion.
2. **Offsets are per-BAM, per-read-length, pre-aggregation.** `BamSetProvider`
   never shares offsets across BAMs: each BAM → its own metagene offsets → its own
   P/A profile → only then folded/merged into the aggregate (and kept per-sample
   for tiers). Bigwigs carry no read length, so they bypass this (coverage only).
3. **`A = P + 3` is a validated default, not a law.** The reliability test
   (decision 1) validates the **A-site** offset on the truth set alongside P-site;
   we don't assume the P→A gap is a constant 3 nt across read lengths — if the
   data disagrees for some lengths, the provider carries an explicit A-offset
   table rather than `P+3`.

## Open method sub-choices (still yours)

- Metagene offset-pick rule within the plausible window (argmax of 5′ pile-up vs
  changepoint vs riboWaltz) — to be decided by the reliability test above.
- Consequentiality policy weights (how tier/confidence/expression/context combine).

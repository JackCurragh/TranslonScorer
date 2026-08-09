# Event-level scoring model

Score genomic **events** once, compose them into per-**translon** (feature)
reports. No annotation decisions — just evidence. Mirrors the matrix engine's
"score unique things once, fan out" philosophy, applied to scoring.

## Coverage & method spec (LOCKED 2026-06-13)

Decisions that apply app-wide (matrix *and* user-BAM/bigwig paths), via a shared
`CoverageProvider` interface and one scoring core:

- **Inputs / audience.** General users: 1–20 BAMs or bigwigs scoring a feature
  set; annotation scale: the sparse matrix. Same scoring strategy; only the
  coverage provider differs. BAM is first-class; **bigwig is a lossy mode**
  (coverage only → no junction spanning, no multimapper resolution, no P/A
  distinction).
- **Coordinate space: genome-native.** Transcriptome BAMs are accepted by
  projecting reads to the genome and resolving **isoform multimappers** (same-gene
  isoform hits collapse to one genomic locus; only truly distinct genomic loci
  remain multimappers).
- **Ribosome site (app-wide):** **initiation = P-site**, **elongation = A-site**,
  **termination = A-site**. A-site = P-site + 3 nt (one codon) → same frame, so
  periodicity/identifiability are unchanged; only the initiation step shifts one
  codon to the P-site (sharper start signal). Generate P-site and A-site profiles
  (offsets 3 nt apart). **This corrects the current code, which uses A-site for
  init — re-validate init against the GAPDH golden.**
- **Offset determination:** metagene (canonical) + offsets-file + global. Metagene
  = aggregate **5′ ends of unique reads around annotated start codons**, per read
  length → P-site offset per length. **Constrain to plausible offsets**: offset ∈
  [offset_min, min(offset_max, ⌊read_len × max_frac⌋)] (e.g. 8…min(20, 0.667·L)) —
  a 25 nt read cannot have an 18 nt offset. Restrict to a sensible read-length
  range. (Drops the frame-nudge for offset determination.)
- **Aggregation:** score on the **aggregate** by default; on *confusing* loci
  (ambiguous/heterogeneous) **cluster and score per-cluster aggregate**
  (`score_clustered`).
- **Multimapper:** **unique by default**, but keep a **per-event mappability
  ledger** (`{unique_reads, multimapper_reads, n_loci_per_multiread}`, second-pass
  OK) so we can relax without re-reading: **copy-number-aware** for paralogs
  (distribute to expected N), **guilty-by-association** when no clear paralogy.
- **Consequentiality (redesign, not hard gates):** dynamic, context-aware —
  measurability *tier* (coverage/support/**signal concentration**, à la the
  notebook's `support_tier`/`top_position_fraction`) × per-aspect **confidence** ×
  **size-factor-normalised expression** × context attributes (biotype, novelty,
  small-RNA overlap as *attributes*, never length/biotype hard-deletes — a
  start-stop ORF or translated small RNA stays in; spike-concentrated pileups are
  down-ranked by concentration).

## Why events (measured on the real translon DB, 8.85M translons)

Genome-wide via `extract_events` / `run_extract` (277 s, whole genome):

| event type | raw instances | unique events | dedup |
|---|---|---|---|
| elongation | 54,465,691 blocks | 687,274 | 79× |
| initiation | 8,852,481 | 387,756 | 23× |
| termination | 8,852,481 | 176,397 | 50× |
| junction | ~45,623,000 | 267,028 | **171×** |
| **total** | **~117.8M** | **1,518,455** | **~78×** |

~117.8M raw structural units → **1.52M unique events**. Scoring per-event instead
of per-translon is ~78× less work *and* makes shared evidence consistent: two
isoforms sharing a CDS exon get the **same** elongation score, never two
contradictory ones. Phase-aware elongation (687k) is within 1% of the
phase-agnostic count (681k) — almost no segment is translated in two frames, so
the dedup is clean. Outputs: `event_store/{events,feature_event,event_overlap}/`
(26 MB / 311 MB / 7.8 MB).

## Event taxonomy & identity

| type | genomic key | shared by | aspect scored |
|---|---|---|---|
| `init` | (chrom, start_pos, strand) | translons starting there | initiation peak / in-frame establishment |
| `term` | (chrom, stop_pos, strand) | translons ending there | stop dwell / step-down / readthrough |
| `elongation` | (chrom, seg_start, seg_end, strand, **phase**) | isoforms using that CDS segment in that frame | periodicity / in-frame fraction / uniformity |
| `junction` | (chrom, donor, acceptor, strand) | isoforms using that junction | spanning-read support / frame continuity |

**Deterministic id** (stable across annotation rebuilds → stable joins,
incremental updates, score reuse):

```
event_id = xxh3_64(f"{type}|{chrom}|{strand}|{start}|{end}|{phase}")
```

Never autoincrement — a rebuilt annotation must reproduce the same id for the
same genomic event so old scores still join.

## Overlap / contention (the convolution problem)

A uORF overlapping a CDS has clean signal where only it claims the region and
**convoluted** signal where another frame's feature also claims it. Key moves:

1. **Contention is structural — computed once, not per sample.** Which positions
   are claimed by ≥2 reading frames is a property of the annotation (or of all
   sequence ORFs). Compute a contention map once at extraction; cache it. The
   expensive "all ORFs" pass is here, *once*, never per score.
2. **Known overlaps deconvolve by frame, no ORF enumeration needed.** In a
   contended window the signal splits into f0/f1/f2; the uORF's contribution is
   simply its frame component — which the matrix engine already produces. So
   attributing signal between two *known* translons is cheap; full deconvolution
   (`frame_disambiguation` / `deblur`) is reserved for the escalation path.

Each elongation event therefore carries `clean_bp` / `contended_bp` and links to
the overlapping events, so "does this need a second pass?" is a lookup.

**Measured (chr22, `extract_events`):** of 16,340 elongation events, **48.9% are
frame-contended** (overlap another event in a different phase), covering **34% of
elongation bp**. So in a comprehensive translon catalog you cannot naively
score per-segment for roughly half the events — contention handling is load-
bearing, not a nice-to-have. (Junctions dedup even harder than segments:
957k raw → 5,996 unique on chr22, **160×**.)

## Scoring tiers (default cheap, escalate on suspicion)

| tier | what | when |
|---|---|---|
| **0 aggregate** | score events from cohort-summed coverage | always; resolves clean, well-covered events |
| **1 per-sample / cluster** | per-sample rollup, dispersion, clustering | signal concentration high (one sample/study driving it), borderline aggregate |
| **2 deconvolved** | frame-resolved attribution / `deblur` on contended events | event is contended and identifiability is low |

Escalation **triggers are cheap** (the rollup gives per-sample dispersion almost
free); only the gated work (clustering, deconvolution) is expensive, and it runs
only on the flagged minority. Tiers are additive rows, never recomputation.

## Storage: structure (relational) ⟂ scores (columnar)

Static, annotation-versioned **structure** → reference DB (extends
`translons.sqlite`). Large, sparse, recompute-heavy **scores** → append-only
partitioned Parquet star schema. Joined on `event_id`.

### Reference DB (DDL sketch)

```sql
CREATE TABLE events (
  event_id            INTEGER PRIMARY KEY,   -- xxh3_64 of canonical key
  type                TEXT NOT NULL,         -- init|term|elongation|junction
  chrom               TEXT NOT NULL,
  strand              INTEGER NOT NULL,
  start               INTEGER NOT NULL,      -- 0-based genomic
  "end"               INTEGER NOT NULL,
  phase               INTEGER,               -- 0/1/2 for elongation; NULL else
  attrs               TEXT,                  -- JSON: start_codon_class, donor/acceptor…
  clean_bp            INTEGER,               -- contention summary
  contended_bp        INTEGER,
  annotation_version  TEXT NOT NULL
);
CREATE INDEX events_loc ON events(chrom, start, "end");

-- features = translons (existing table); thin link to a stable feature_id
CREATE TABLE feature_event (
  feature_id  TEXT    NOT NULL,   -- translon_id
  event_id    INTEGER NOT NULL,
  role        TEXT    NOT NULL,   -- init|elongation|junction|term
  rank        INTEGER,            -- translation order
  phase       INTEGER,            -- frame the feature uses this event in
  PRIMARY KEY (feature_id, event_id, role)
);
CREATE INDEX fe_event ON feature_event(event_id);

-- contention detail (only for elongation events that overlap another frame)
CREATE TABLE event_overlap (
  event_id        INTEGER NOT NULL,   -- this elongation event
  other_event_id  INTEGER NOT NULL,   -- overlapping event in a different phase
  overlap_start   INTEGER NOT NULL,
  overlap_end     INTEGER NOT NULL,
  PRIMARY KEY (event_id, other_event_id)
);
```

### Score store (Parquet, append-only star schema)

```
fact_event_score/
  data_version=<matrix gen>/ group_kind=<aggregate|sample|cluster|study>/ part-*.parquet

columns:
  event_id            uint64
  group_id            str        # 'aggregate' | sample_name | cluster_id | study_id
  aspect              str        # init | elongation | term | splice
  tier                str        # aggregate | per_sample | deconvolved
  status              str        # NA | INSUFFICIENT | SUPPORTED | UNSUPPORTED | AMBIGUOUS
  score               f32        # within-eligible [0,1]
  identifiability     f32        # frame-attributable fraction in contended regions
  n_reads             f32
  metrics             struct/json
  annotation_version  str
  method_version      str
```

Grain = `(event_id, group_id, aspect, tier, *_version)`. Append-only: re-scoring
(new method/data) writes new partitions; the reference DB and prior scores are
untouched — reproducible, A/B-able, and the second pass is just more rows.

Cardinality stays bounded: Tier 0 is ~1 row per event per aspect; the fan-out to
per-sample/deconvolved happens only for escalated events.

## Build phases

1. **`extract_events`** (annotation only, once per annotation version):
   `translons` + `translon_blocks` → `events`, `feature_event`, contention.
   Reuses block ranks for phase; deterministic ids.
2. **`score_events`** (matrix engine → `fact_event_score`): Tier 0 aggregate by
   default via region-scoped `tabulate_profiles` / `feature_metrics`; escalate
   per triggers.
3. **`report_translon`** (`feature_event ⨝ fact_event_score`): per-translon,
   per-aspect, per-group evidence. Elongation = length-weighted over its segment
   events. No collapse to one number, **no annotation decision**.

## Score semantics (v0) — measurement vs decision

Scoring is split so the **continuous evidence** is primary and the
**SUPPORTED/UNSUPPORTED call is a thin, re-derivable function** of that evidence
plus an explicit, versioned `ScoreThresholds`. Every event-aspect is one
long-form row: `event_id, aspect, group, tier, n_reads, metric, metric_name,
eligibility, call, evidence(json), thresholds_version`. `eligibility ∈ {NA,
INSUFFICIENT, ELIGIBLE}` is never conflated with `call ∈ {SUPPORTED,
UNSUPPORTED, AMBIGUOUS}` — "couldn't test" ≠ "failed". Change a threshold and
re-derive calls; never re-measure.

| aspect | headline metric (v1) | range | meaning | continuous evidence | call logic |
|---|---|---|---|---|---|
| `init` | `init_rise` = log2((body+α)/(leader+α)) | ℝ | log2 fold-change of body vs leader (robust median codon levels, α=1), median over flanks 9/18/30/60 | `rise_by_flank`, `stability`, `flank_peakiness` (review flags, not gates) | SUPPORTED if rise≥1 (2×); UNSUPPORTED if ≤0; AMBIGUOUS if borderline |
| `elongation` | `elong_in_frame` | 0…1 | in-frame fraction where **unconfounded**; paired with `breadth` | `clean_in_frame`, `overall_in_frame` (= Chothani PIF), `cif`, `breadth` (covered_nt/span), `identifiability`, `competitor_share`, `noise_share` | AMBIGUOUS if contended & identifiability<0.5; SUPPORTED if in-frame≥0.5 **and** breadth≥0.2; else UNSUPPORTED |
| `term` | `term_drop` = log2((body+α)/(UTR+α)) | ℝ | log2 fold-change of body before vs UTR after the stop | same step machinery as init, plus `dropoff` (Chothani, bounded) | SUPPORTED if drop≥1 (2×); UNSUPPORTED if ≤0; AMBIGUOUS if borderline |
| `junction` | `junc_spanning` | 0…∞ | reads whose CIGAR intron matches (donor,acceptor) | `spanning` | SUPPORTED if spanning≥20 |

`metric` is the natural per-aspect quantity — **not** normalised across aspects
(no composite). Golden = `outputs/reference_scores_gapdh.parquet` (22 rows);
`verify_vectorised.py` gates scalar≡vectorised; `run_score_chromosome.py` gates
run≡golden — all currently Δ=0.

**v1 (2026-06-12):** init/term moved ratio→log2 fold-change (de-saturated:
canonical GAPDH start 5.36, stop 7.68; internal start −10.36); elongation gained
a `breadth` guard to kill single-spike false positives; peakiness/stability
demoted to evidence. Still to iterate: junction normalisation (spanning/crossing
+ overhang-weight) & frame continuity; null-model anchoring; term stop-peak +
readthrough; per-sample tier escalation.

**v1.2 (2026-08-09) — translation signature scores as evidence:** the three
components of the Chothani et al. (*Mol Cell* 2022) signature now ride in the
evidence blob rather than as new headline metrics, because the record schema
carries exactly one metric per row and the incumbents already own it.

- **PIF was already there under another name.** `overall_in_frame` is exactly
  `sum(psites[seq(1,n,3)]) / sum(psites)`. No field was added. The headline
  `elong_in_frame` is *not* PIF — it is the same ratio restricted to
  uncontended positions, so the two diverge precisely where events overlap in
  different frames.
- **`cif`** (elongation) reproduces two quirks of the published R on purpose:
  the threshold is the literal `33.33`, not `100/3`, so an exactly-uniform codon
  counts as frame-0 dominant; and the denominator is *every* codon, so
  zero-coverage codons lower the score instead of being excluded.
- **`dropoff`** (termination) is a bounded [0,1] ratio over a fixed 33-nt window
  using frame-0 positions only — deliberately *not* the same statistic as
  `term_drop`, which is an unbounded multi-flank log2 fold-change over all
  positions. Neither supersedes the other. It assumes the translon span
  includes its stop codon; if the upstream annotation excludes it, every value
  is 3 nt out of register.

Validated against the R *formulas* (fuzz-tested against a transliteration in
`tests/reference_signature.py`), **not** yet against published values on real
data.

**v1.1 (2026-07-14) — splice-aware init/term flanks, junction wiring fixed:**
`init_rise`/`term_drop` previously read the leader/UTR flank as raw flanking
*genomic* bases (`aspects._flank_positions`), which is wrong whenever a splice
site falls inside the flank window — the true leader/UTR sits on the far side
of the intron, and the flat window instead reads mostly zero-coverage
intronic sequence, spuriously inflating the step metric toward a false
SUPPORTED. `aspects._project_flank` now walks in transcript order, jumping
across a nearby intron from an optional `splice_context` (built by
`events.build_splice_context(context_gtf)`, threaded through
`score_events(_vectorised)` → `_score_events_over_provider` →
`score_matrix_workflow`/`score_bams_workflow`/`pipeline_workflow` →
`--context-gtf` on the corresponding CLI commands). Evidence gained a
`flank_spliced` bool marking whether a jump was applied. This follows the
same genome-wide, strand-scoped, not-gene-scoped intron lookup already used
by `events._context_junction_events`; where isoforms disagree on the
immediate upstream/downstream exon it follows exactly one path — isoform-of-
origin remains deferred (§6.3 of reannotation_engine_design.md), this only
fixes the single-path case, which is the common one.

Separately: `_score_events_over_provider` was never calling
`provider.junction_support()` — every junction event scored against an empty
support dict regardless of provider capability, so junction events always
came out INSUFFICIENT/n_reads=0 even though both `MatrixProvider` and
`BamSetProvider` fully implement `junction_support()`. This was a pure
wiring gap (already flagged as "Junction scoring dead… Undiagnosed" in
project_overview.md §10) — now fixed via `_junction_support_for_chrom`. That
fix also exposed a second footgun in the same area: a `runtime_checkable`
Protocol (`SupportsJunctions`) only checks that a method exists, not that it
works — `BigwigSetProvider.junction_support()` exists and always raises
`NotImplementedError`, so `isinstance(provider, SupportsJunctions)` was
`True` for it too. `_junction_support_for_chrom` now also catches
`NotImplementedError` explicitly rather than relying on `isinstance` alone.

**v1.2 (2026-07-14) — bigwig coverage restored, `--offset-method metagene`
fixed on score-bams:** the bigwig path described in the "LOCKED" spec above
was previously scaffolded but never reachable: `BigwigSetProvider` existed
but `workflows.py` never imported it, and its `stranded` flag was a stub
that always sum-merged forward+reverse into one `count` regardless. Both are
now fixed — `BigwigSetProvider.coverage()` emits a real `strand` column when
`stranded=True` (forward/reverse dict entries required in that mode), and
`score_bigwigs_workflow` (mirroring `score_bams_workflow`, same
`_score_events_over_provider` core, same `--context-gtf` support) is wired
to a new `score-bigwig` CLI command (`--bigwig` for unstranded, or
`--forward-bigwig`/`--reverse-bigwig` pairs for stranded — strongly
preferred, since Ribo-seq is stranded). Still lossy as documented: no P/A
site, no junction spanning, no mappability.

Separately, `--offset-method metagene` on `score-bams`/`pipeline` was dead in
practice — `BamSetProvider`'s metagene branch needs genomic start-codon
coordinates it was never given, so it silently fell back to the fixed
global offset every time (only a log warning, despite the CLI presenting
`metagene` as a normal working choice and this doc calling it "canonical").
`score_bams_workflow` now derives start-codon coordinates from the
extracted `init` events (free — already read for scoring) and passes them
through. Caveat: the metagene histogram fetches directly from the BAM by
genomic coordinate, so this only works for genome-aligned BAMs
(`transcriptome=False`) — transcriptome-aligned BAMs still fall back to
global, unchanged from before.

Mappability-based eligibility is still unimplemented (see
project_translonscorer_bigwig_mappability_gap memory / audit notes) — next
candidate if bigwig/BAM scoring correctness needs a repetitive-region caveat.

**v1.3 (2026-07-14) — mappability-track diagnostic annotation
(`map_track_mean`/`map_track_low`):** a user-supplied precomputed mappability
track (Umap/GEM-mappability style bigwig, values ~0..1) can now be attached
via `--mappability-bigwig` on `score-matrix`/`score-bams`/`score-bigwig`/
`pipeline`. Deliberately named `map_track_*`, not `mappability_*`, to avoid
conflating it with the pre-existing (and still unimplemented) BAM-NH-tag
`SupportsMappability`/`mappability_ledger()` concept in `coverage/base.py` —
those are different ideas that happened to share a name.

`workflows._map_track_for_chrom` reads the track over each event's padded
window (reusing `_merge_event_regions`, same pad as coverage) via a second,
independent `BigwigSetProvider` instance and computes a mean; `map_track_low`
is `mean < thr.mappability_low` (new `ScoreThresholds` field, default 0.5).
Threaded through `score_events(_vectorised)` exactly like `junction_support` —
merged into the raw evidence dict right before `event_record(...)`, uniformly
across all event types. **Never affects eligibility/call** — same "review
flag, not a gate" precedent as `flank_peakiness`/`stability`.

Promoted to first-class `_RECORD_SCHEMA` columns (not left in the `evidence`
JSON blob) specifically because `report.py`'s `_compose_per_translon` only
ever selects a fixed column subset from the score store — anything left in
`evidence` never reaches the report a user actually reads. `report.py`
aggregates per aspect using an **unweighted** mean (not the read-weighted
mean used for `metric`) — read-weighting would wash out the signal for
exactly the zero/low-read events this annotation exists to explain — and
`map_track_low` as `.any()` across the aspect's events (a single flagged
event is worth surfacing). Fully backward compatible: `mappability_bigwig=
None` everywhere by default, new columns stay null, golden gate unaffected.

**v1.4 (2026-07-14) — one coherent matrix-scoring path: `MatrixProvider`
reads the P-site index directly.** Closes the "two matrix pipelines" gap
named earlier this session: `score_matrix_workflow` (`MatrixProvider`, flat
`ref_offset`) already shared `_score_events_over_provider`/
`score_events_vectorised`/`scoring/aspects.py` with score-bams/score-bigwig,
but had a known accuracy bug (one flat P-site offset regardless of read
length → `elong_in_frame` ~0.33 floor). `score_matrix_rollup_workflow`
(the `--gtf`/FrameRollup path) fixed the offset accuracy but did it by
hand-rolling a second, transcript-coordinate scoring pipeline that bypasses
`_score_events_over_provider` entirely and only ever scores elongation —
duplicated orchestration, permanently forfeiting init/term/junction/
mappability on the "accurate" path.

Investigation found `build-psite-index`'s Phase 2 already writes exactly the
missing substrate: a genomic, per-(sample, length)-offset-corrected,
chrom-sharded P-site index (`_SHARD_SCHEMA`: `p_site, strand, length,
sample_id, count`) — it was just never read as ordinary genomic coverage,
only re-projected into transcript space to feed FrameRollup. New
`psite_index.query_genomic_coverage` reads it directly (`pos, count[,
strand][, sample_name]` — the `CoverageProvider` contract); `MatrixProvider`
gained a second coverage strategy (`psite_index_dir` constructor param) that
uses it instead of the flat `ref_offset` path when set —
`.junction_support()`/`.mappability_ledger()` are untouched (already
offset-independent, read raw partition BAMs directly). No new provider
class, no new scorer.

CLI: `score-matrix --psite-index` no longer requires `--gtf` — without it,
this new mode activates (recommended); with it, unchanged FrameRollup
behaviour (backwards compatible, not removed). Also added to `pipeline`
(previously absent entirely). `--gtf`/FrameRollup is now documented as
elongation-only and legacy — a candidate for deletion
(`score_matrix_rollup_workflow`/`score_frame_rollup`/
`score_elongation_from_rollup`) once the new path is validated on real
cohort-scale data (a manual/HPC step, not done as part of this pass).
`build_frame_rollup`/`calibrate_offsets` stay regardless (Phase 1 of
`build-psite-index` calls them directly); `build_coverage_index`/
`profile_from_index` (CoverageIndex, "FR6") stay too — a separate
downstream clustering/differential-translation feature, not part of event
scoring, not redundant with this.

## Open choices

- Junction event extraction (same dedup pattern; not yet measured).
- Elongation combination: report per-segment *and* a length-weighted aggregate.
- Whether `feature_id` reuses `translon_id` or a content hash (`feature_key`
  already exists on `translons` — check if it pre-dedupes isoform-identical ORFs).
- Score store engine: plain partitioned Parquet vs DuckDB-over-Parquet for the
  report joins at scale.

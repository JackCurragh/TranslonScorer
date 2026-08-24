# TranslonScorer — Implementation Roadmap

Companion to `reannotation_engine_design.md`. That document says *what* the
system should be; this one maps the path from the current code to it, with an
**optimal-implementation assessment for every task** — the realistic best way to
build it, not just a checkbox.

> **Naming note.** The L0–L4 labels below are a conceptual layering, not module
> names. Two of these layers exist as packages and were renamed to say what they
> hold: **L0b → `TranslonScorer/alignments/`** (one row per alignment) and
> **L1 → `TranslonScorer/counts/`** (5′-end and junction counts per sample).
> `VersionKey.l0a_version` is now `source_data_version`. The remaining layers
> (L0a upstream counts, L2 QC gate, L3 aggregate, L4 scoring) have no package of
> their own, so the labels are retained here as design vocabulary only.

---

## Current state (ground truth)

The code is further along than a two-track reading suggests, but in a shape that
conflicts with the design in two specific ways.

| Layer | Exists today | Gap to target |
|---|---|---|
| **L0a** read matrix | upstream: count Parquets `(read_id, sample_id, count)` per partition + manifests | not version-stamped |
| **L0b** alignment | 257 genome BAMs | **no persisted loci table; no multimapping metadata captured** |
| **L1** raw 5′ cache | `build_psite_index` persists per-chrom `(p_site, strand, length, sample_id)→count` shards | **offset baked in at build (not raw); CDS-restricted; no aln metadata; full-BAM rescan per build** |
| **L2** QC gate | `calibrate_cohort` → offsets + rich per-sample QC; `usable_sample_lengths` | no explicit include/exclude + reason + offset_confidence as first-class fields |
| **L3** aggregate profile | `query_frame_rollup` / `query_coverage_index` aggregate at query | reads offset-baked index; "passing" filter not built in |
| **L4** scoring | `score_events_vectorised` (init/term/elong/contention); `score_frame_rollup` (per-sample offset, elong only) | **two tracks; no formal event keys; elongation frame not splice-path-aware** |

Plus a structural debt: **three overlapping scan/index implementations** —
`psite_index.build_psite_index`, `matrix_rollup.build_frame_rollup` /
`build_coverage_index`, `matrix_qc.build_read_loci_index`. Each independently
full-scans the BAMs.

### The one cross-cutting efficiency principle

Every current builder full-scans the 257 BAMs. **The single highest-leverage
move is to scan each BAM exactly once into a persisted L0b loci table, then never
touch BAMs again** — L1, QC, and scoring all read from L0b + counts. This alone
removes the dominant I/O cost and is the precondition for caching. Treat "BAMs
are touched once, ever" as an invariant the whole roadmap serves.

A note on physical layout: the **scan** is naturally by the 257 sequence-prefix
partitions, but **queries** are by genomic locus. So building L0b/L1 is a
one-time **transpose** (scan-by-prefix → write-by-chrom). Pay that shuffle once;
everything downstream gets locus-local reads for free. **Do not underestimate
this step:** at cohort scale (billions of alignment rows across 257 partitions)
it is an *out-of-core external shuffle* with real memory-pressure and spill
behaviour. The logic is simple; making it robust at scale is the engineering.

> **Effort estimates throughout are logic-complexity, not scale-hardening.**
> Every "small / N days" tag reflects how hard the *code* is, ignoring the
> scaling tail (6k+ samples, billions of reads). At this scale, simple features
> routinely grow days of memory/throughput tuning. Use the estimates for
> *sequencing*, not scheduling.

---

## Phase 0 — Decisions that gate sequencing

### 0.1 Canonical scan/index path
**Optimal:** adopt `psite_index` as the spine — it already produces per-chrom,
per-sample, per-length sparse shards with a working query layer, which is 70% of
L1+L3. Refactor *it* rather than `build_frame_rollup` (whose per-feature framing
bakes in transcript-coupling) or `matrix_qc.build_read_loci_index`. Mark the
other two for deletion in 7.3 once parity is proven by golden-output tests.
**Effort:** decision + a parity test harness (½ day).

### 0.2 L1 scope — genome-wide, minus a contamination mask
**Settled: genome-wide.** CDS-restriction blocks novel-ORF scoring (the whole
point of the tool); it is not a real option. The genuine decision is **what to
exclude**, and it is a **blacklist mask, not a coordinate restriction.**

Ribo-Seq is heavily rRNA/tRNA-contaminated — these are frequently the *majority*
of reads. Dropping the CDS filter naively does **not** give a small uplift; it
risks indexing enormous rRNA/tRNA pileups that dwarf the real signal and inflate
L1. So: index genome-wide but subtract a **mask** of rRNA, tRNA, and known
artifact/blacklist regions (ENCODE-style). "Don't restrict to CDS" is right;
"index every rRNA stack" is its own mistake.

Do **not** trust any a-priori size multiplier — the footprint depends entirely
on how contamination-heavy the cohort is and must be **measured** on one
chromosome before committing cohort-wide. Drop the `_assign_frames_sweep` in-CDS
filter.

**DECIDED (2026-06-29): no mask for now — allow rRNA/tRNA to align.** Build L1
genome-wide and unmasked initially; accept the size risk. Keep masking as a clean
insertion point (a region-exclusion filter at the L0b write in 1.2) so it can be
switched on later without redesign. Still worth a 1-chrom size probe to know how
big unmasked L1 actually is.

### 0.3 Multimapper accumulation semantics (dictates the L1 schema)
**DECIDED (2026-06-29): (a) unique-only L1.** L1 holds unique mappers only; all
multimapper handling is deferred to an L0b pass (6.3). L1 sums are always valid;
no `weight` column needed in L1. L0b still retains *all* alignments + metadata so
the multimapper work (filter / guilty-by-association) can be added later without
a rescan.

(Original analysis retained for context:) A multimapper contributes
its *full* count to **every** locus it maps to, so naively summing L1 over
multimappers **overcounts**. A `unique`/`multi` flag only cleanly supports the
*exclude* case. The real fork:
- **(a) Unique-only L1.** L1 holds unique mappers; *all* multimapper handling
  goes through L0b separately (filter/EM at query). Simplest, keeps L1 sums
  always valid, but multimapper-inclusive queries always pay an L0b pass.
- **(b) Fractional-weight L1.** L1 stores a weight (1/NH, or local-density-
  weighted) so sums stay valid with multimappers included. Richer default
  queries, but bakes a weighting choice into L1 and complicates re-weighting.

This is genuinely a requirements call (how much do multimappers matter to the
default path vs only to hard loci?) and it sets the L1 schema in 2.2.

### 0.4 Versioning key
**Optimal:** a single content fingerprint
`{genome_id, aligner_cfg_hash, junction_set_id, l0a_version, multimap_policy}`
stamped into a `_meta.json` at every artifact root and echoed in directory names
(`l1/genome=GRCh38/jset=v44/l0a=2026-06/…`). Offset version is a *sub*-key on L3
only (so recalibration invalidates L3, not L1 — see 2.1). **Effort:** small;
define once, thread everywhere.

---

## Phase 1 — L0b: multimapping-aware loci (critical path)

### 1.1 L0b schema
**Optimal:** one row per (unique_read × alignment), columns
`read_id, chrom, pos5, end, strand, length, cigar, mapq, nh, is_secondary,
aln_score, mismatches, junctions_crossed (list[donor,acceptor]), weight`.
Read `nh` from the **aligner's `NH` tag**, *not* by counting BAM records:
aligners cap reported secondary alignments (STAR `outSAMmultNmax` etc.), so a
record count *underestimates* true multiplicity. Prefix-partitioning still helps
— all *reported* alignments of a read are co-located, so enumerating a read's
loci needs no cross-partition join — but the multiplicity *value* must come from
the tag. `weight` is left null at build (a policy input, not a baked decision).
Store `cigar`/`junctions_crossed` only when non-trivial (spliced/clipped) to keep
the table lean.
**Effort:** schema + extractor (below).

### 1.2 One-time extractor — merge-then-linear-scan
**Optimal: offload the transpose to `samtools`, don't hand-roll it.** The
scan-by-prefix → write-by-chrom transpose is an out-of-core shuffle; rather than
implement that in Polars, let the mature external-sort tool do it:

```
257 prefix BAMs → samtools merge/sort → per-chrom coordinate-sorted BAM(s)
              → single streaming pysam pass per chrom → L0b Parquet shards
```

This makes L0b construction an **embarrassingly-parallel, per-chrom linear pass**
(no shuffle in our code), and the coordinate sort gives **position-sorted output
for free** → exactly the `pos5` row-group statistics / predicate pushdown 1.1
wants. The per-chrom sorted BAMs are a reusable artifact. Split by chrom — do
**not** keep one monolithic 10s-of-TB BAM.

The per-record transform (extract `pos5`, length, strand, `NH` tag, junctions
from CIGAR) is the same work regardless of input order; coordinate-sorted input
just lets each chrom shard stream out without holding the cohort in memory.

> **Precondition: globally-unique `read_id` — VERIFIED (2026-06-29).** Merging is
> only safe if read IDs don't collide across partitions. Checked on the
> `data/global_partitioned` cohort: IDs are encoded `partition_ordinal × 10⁹ +
> local`, giving disjoint billion-wide blocks; all sampled pairwise intersections
> = 0, total == unique. **No re-keying pass needed.** (If a future cohort breaks
> this, add one offset-by-partition re-key applied to *both* BAMs and counts
> before merge.)

Apply the **contamination mask (0.2) at this stage** — masking definite
rRNA/tRNA before writing L0b is the biggest size win, since those reads can be
the majority. Keep all *alignments* (multimapping) for non-masked regions.
**Effort:** medium; the `samtools` merge/sort is the heavy machine-time job
(*logic only* — see scale caveat).

### 1.3 Version stamping
**Optimal:** write `_meta.json` (0.4 key) into the L0b root; refuse to build L1
from an L0b whose key disagrees with the requested config. **Effort:** small.

---

## Phase 2 — L1: raw 5′ cache, offset decoupled

### 2.1 Strip offset-baking
**Optimal:** delete the `p_site = pos5 ± offset` step from `_build_worker`
([psite_index.py:436](../TranslonScorer/psite_index.py)); store raw `pos5`.
Offsets move to query/L3 as a join + arithmetic (`score_frame_rollup` already
does exactly this analytically — reuse that code, don't reinvent). Net: a
deletion plus a moved join. The payoff is large — offset recalibration becomes a
cheap L3 recompute instead of a full cohort rescan.
**Effort:** small, high value.

### 2.2 Build L1 from L0b + counts, factored
**Schema is set by decision 0.3, not here.** Materialize L1 per-chrom as
`(pos5, strand, length, sample_id, …) → count`, where the `…` and the meaning of
`count` depend on the 0.3 fork:
- **Unique-only L1:** L1 carries unique mappers; sums are always valid;
  multimappers are handled exclusively via an L0b pass (6.3). No multimapper rows
  in L1.
- **Fractional-weight L1:** L1 stores a `weight`-adjusted count so multimapper-
  inclusive sums stay valid; full per-alignment metadata (MAPQ, score, CIGAR)
  still stays in **L0b** for re-weighting / EM (6.3).

In **both** cases, full metadata lives in L0b, not L1 — L1 stays the small, fast
hot-path substrate. Whatever 0.3 chooses, do **not** store raw full-count
multimapper rows in L1 with only a flag: that makes "include multi" silently
overcount.

Build by joining the L0b loci (chrom-partitioned) to the L0a count Parquets on
`read_id`, then group. This replaces the BAM rescan with a Parquet join — the
"scan once" payoff. Keep it lazy/streaming in Polars to bound memory.
**Effort:** medium (2–3 days), *logic only* — see the scale caveat above.

### 2.3 Apply 0.2 scope
**Optimal:** falls out of 1.2 (no CDS filter) — L1 inherits genome-wide coverage
automatically. Just confirm the query layer no longer assumes CDS-only.
**Effort:** small.

---

## Phase 3 — L2: explicit QC artifacts

### 3.1 First-class decision fields
**Optimal:** `calibrate_cohort` already computes the hard parts
(periodicity_score, length/ligation metrics, recommended offsets/lengths). Add a
thin **policy layer** on top that emits, per `(sample, length)`:
`authenticity_score, periodicity_score, offset, offset_confidence,
include_exclude, reason`. Keep the policy *declarative* (thresholds in one
config object) so the conservative bias is tunable, not hardcoded.
`offset_confidence` = margin between best and second-best frame fraction (cheap,
already computable from the frame counts). **Effort:** small–medium.

### 3.2 Name the circularity
**Optimal:** add a `reason` enum value `CDS_TRAINED_QC_EXCLUSION` and document in
the artifact's schema doc that unusual-but-real biology may be filtered. Make the
QC pass also emit the *excluded* rows (not silently drop them) so downstream can
override. **Effort:** small.

### 3.3 Persist keyed to L0b version
**Optimal:** QC artifact root carries the 0.4 key + an `offset_version` it
defines. L3 consumes this. **Effort:** trivial once 0.4 exists.

---

## Phase 4 — L3: aggregate profile (offsets live here)

### 4.1 Apply offsets at aggregation
**Optimal:** repoint `query_frame_rollup` / `query_coverage_index` at raw L1;
join the `(sample_id, length) → offset` table and compute `p_site = pos5 ±
offset` (or `frame = (tx_pos + offset) % 3`) in the query expression. The pivot
and frame logic already exist — only the offset source changes (from
"pre-baked" to "joined"). **Effort:** small (mirrors 2.1).

### 4.2 Build the passing-data default aggregate
**Optimal:** join the L2 include/exclude table as a semi-join filter before
aggregation, then sum to the default genomic P-site profile. Cache this profile
keyed by `(0.4 key, offset_version, qc_policy_version)` so routine re-scores skip
it entirely. **Effort:** small–medium.

---

## Phase 5 — L4: one scorer + formal event keys

### 5.1 Collapse to one scorer
**Optimal:** keep `score_events_vectorised` (prefix-sum + contention — the
*complete* scorer) and **delete `score_frame_rollup`**. Its only unique value was
per-sample offsets, which now live upstream in L1/L3. Feed the vectorised scorer
the L3 aggregate profile per chrom/strand (the `_score_events_over_provider`
plumbing already does per-chrom, per-strand fetch — reuse it, swap the provider's
source from BAM-rescan to L3). This *is* the resolution of the original two-track
problem. **Effort:** medium; mostly deletion + re-pointing.

### 5.2 Formal event keys
**Optimal:** compute a stable key at **extraction** time (cheap, once), not at
score time. `splice_path_id = hash(ordered (donor,acceptor) chain)`; event key =
hash of the tuple in §6.1 of the design. Dedup is then a `group_by(key)`. Store
the key as a column in the events Parquet so dedup and join are trivial
downstream. **Effort:** medium.

### 5.3 Splice-path-aware frame
**Optimal:** precompute, per event, the genomic→transcript coordinate map
(cumulative exon-length prefix offsets). The existing prefix-sum elongation
kernel already operates in *a* coordinate — feed it transcript-projected
positions for spliced events; single-exon events are unchanged. Reuse the kernel;
only the coordinate projection is new — a cumulative exon-length prefix map. The
projection itself is a cumsum, but it is **not** trivial: reads must be assigned
to the correct exon, and reads straddling an exon boundary need explicit handling.
Single-exon events (the majority) skip all of this. **Effort:** medium.

### 5.4 Junction as first-class event
**Optimal:** already supported (`tabulate_junctions`, `score_junction_event`) —
verify it flows through the unified scorer and is keyed by `junction_event_key`.
The junction-spanning counts come from L1's junction stream (built in 1.2/2.2 from
`junctions_crossed`), so no new scan. **Effort:** small (wiring + tests).

---

## Phase 6 — Dynamic aggregation + guardrails

### 6.1 Arbitrary subset/length-band query API
**Optimal:** one function over L1: `(sample_ids?, length_band?, nh_class?,
region) → profile`. Because L1 is per-chrom sparse with sample/length columns,
this is predicate pushdown + group_by — no new machinery, just a clean public
signature. **Effort:** small.

### 6.2 Provenance reporting
**Optimal:** every dynamic-subset result carries a sidecar struct: selection
rule, n_samples, n_studies (join `sample_id → study` from the sample map),
read-length composition, and a cross-study replication flag (signal present in
≥k independent studies). Compute inline from the same query — cheap. Make it
**non-optional** in the return type so it can't be forgotten. **Effort:** small–medium.

### 6.3 Multimap filter + guilty-by-association
**Optimal, staged:** (a) immediate — `nh_class` filter in 6.1 (free). (b) later —
fractional/EM read assignment using L0b metadata: distribute a multimapper's
count across its loci by local unique-read density (guilty-by-association), an
EM-style iteration over the L0b loci table. Self-contained, operates on the
persisted L0b — no rescan, no impact on the hot path. **Effort:** (a) trivial;
(b) a contained research task, post-MVP.

---

## Phase 7 — Orchestration, caching, cleanup

### 7.1 Load-or-build wrappers
**Optimal:** a small `materialize(artifact, key)` helper per layer: hash the 0.4
key (+ layer-specific sub-keys), check the artifact root, build only on miss.
Uniform across L0b/L1/L3/QC. **Effort:** small.

### 7.2 Per-release workflow
**Optimal:** a dispatcher that diffs the requested config against existing
artifact `_meta.json`s: ORF-only change → reuse L0b/L1/L2/L3, re-extract events +
re-score; junction change → rebuild from L0b down. The diff is just key
comparison. **Effort:** small–medium once 7.1 exists.

### 7.3 Delete deprecated paths
**Optimal:** once the parity harness (0.1) is green, remove
`build_frame_rollup`, `build_coverage_index`, `matrix_qc.build_read_loci_index`
and their dead helpers. Reduces three scanners to one. **Effort:** small;
satisfying.

---

## Critical path & sequencing

```
0.1 0.2 0.3 0.4              (decisions — do together, ~1 day)
   └─► 1.1 1.2 1.3          (L0b: the one BAM scan — unblocks everything)
          └─► 2.1 2.2 2.3   (L1 raw + factored, multimap class)
                 ├─► 3.x    (L2 policy layer — parallelisable with Phase 2)
                 └─► 4.x    (L3 offsets-at-query + passing aggregate)
                        └─► 5.1 5.2 5.3 5.4  (one scorer + keys)
                               └─► 6.x 7.x   (dynamic agg, caching, cleanup)
```

**Load-bearing, expensive-to-retrofit, do not defer:** 1.1/1.2 (multimapping-aware
L0b) and 5.2 (event keys). Everything else layers on cleanly or is additive.

**Highest value-to-effort:** 2.1 (offset decoupling — small diff, turns cohort
rescans into cheap L3 recomputes) and 1.2 (scan-once L0b — removes the dominant
I/O cost and enables all caching).

**Realistic MVP cut** (correct, single-track, cached, multimap-capable):
0.x → 1.x → 2.x → 3.1 → 4.x → 5.1 → 5.2 → 5.4 → 7.1. Defer 5.3 splice-frame to
loci that need it, 6.2/6.3b provenance+EM, and the covariate work.

---

## Two parallel tracks: cache build (machine) ‖ scoring dev (human)

The dependency graph above is serial, but the *wall-clock* plan is not: the L0b
build is long machine time, and the scorer can be developed against a stub in
parallel because the `CoverageProvider` seam already exists
(`DictCoverageProvider`, `coverage/matrix.py:167`).

**Step A — freeze three contracts (cheap, unblocks both tracks):**
1. L0b Parquet schema (1.1);
2. the **L3-profile / `CoverageProvider` interface** the scorer reads;
3. the 0.4 version key.
Also settle 0.3 (multimap policy) + the 0.2 mask — but note these gate **L1, not
L0b**, so they need not block the build kickoff.

**Step B — launch the L0b build immediately (machine, runs unattended).** Needs
only 1.1 + 0.4. Retains everything, so independent of 0.3 and of the scorer.
Route: `samtools merge/sort → per-chrom → L0b Parquet` (1.2). This is the long
pole — start it first.

**Step C — develop scoring against fixtures, in parallel (human).** Unified
scorer (5.1), event keys (5.2), splice-frame (5.3), junction events (5.4), QC
policy layer (3.x) — all built/tested through `DictCoverageProvider` over
synthetic L3 profiles. Zero dependency on the running build.

**Step D — when L0b lands:** build L1 (now needs 0.3 + mask) as a second machine
job while the scorer work finishes.

**Step E — converge:** real QC → L3 → run the finished scorer on real L3; run
the golden-parity + sample-count scaling tests (below).

---

## How we know it's efficient enough

**Headline gate (architecture pass/fail), make it a test:** an **ORF-only
re-assessment opens zero BAMs and finishes in minutes** (assert no
`pysam.AlignmentFile` open on a cached re-score; target < ~15–30 min genome-wide).
If a re-score still costs hours, the design failed — no micro-benchmark redeems it.

**The diagnostic tell — slope, not absolute time.** The old `periodicity_qc` was
`640s + 83s/sample` (a per-sample Python loop) → ~140 h at 6k samples. The target
is a **flat/sublinear slope** in sample count: vectorised Polars join/groupby, and
**no per-sample or per-read Python loop in any hot path.** "Is there a Python loop
over samples/reads here?" is the single best efficiency smell-test.

**Supporting targets:**
- **L0b build:** I/O-bound — within ~2× of `samtools view > /dev/null` over the
  same BAMs; linear in workers; **resumable** (per-chrom/per-partition checkpoints).
- **L1 build:** Parquet-join bound, faster than L0b, no BAM access.
- **Per-locus dynamic query:** sub-second to low-seconds from chrom-partitioned,
  `pos5`-sorted L1 (predicate pushdown).
- **Storage:** L1 footprint *measured* after masking (the 0.2 probe), within budget.

**Gate on correctness:** every number above is meaningless without the
golden-parity test (fast path bit-identical to the scalar reference, task 0.1).

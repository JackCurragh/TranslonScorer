# TranslonScorer — Design Descriptor

**Status:** design specification (principles + requirements). This document
describes *what the system should do and why*, independent of the current code.
It is written to be read cold by someone with no prior exposure to the project.

TranslonScorer is **not** an ORF caller over BAMs. It is a **re-interpretation
engine over a durable evidence substrate**: evidence is collected once, and each
transcriptome release is a cheap re-assessment of that evidence — not a fresh
analysis.

---

## 1. Purpose

There is a very large and growing body of published Ribo-Seq data. Ribo-Seq
captures **Ribosome Protected Fragments (RPFs)** — the short stretches of RNA
physically shielded by a translating ribosome — and is therefore direct
evidence of *which regions of the transcriptome are being translated*.

The goal of this system is to **use that evidence in aggregate to annotate
translated regions**, including rare events that no single experiment could
support, and to do so in a way that can be **cheaply re-assessed against a
changing transcriptome** (e.g. for a recurring GENCODE-style release of
translated regions).

Two facts shape the entire design:

1. **Translation annotation is transcriptome-dependent.** The same physical
   evidence is *interpreted* differently under different transcript models. A
   read is only evidence for an ORF if it falls within that ORF's exons in that
   transcriptome. Therefore a re-release of the transcriptome requires
   re-assessing the evidence — but it must not require re-collecting it.

2. **Aligning directly to a transcriptome hides uncertainty.** If you align
   reads to a transcriptome, the aligner silently decides which isoform each
   read belongs to and which positions are exonic. Those decisions are
   transcriptome-dependent, irreversible, and invisible downstream. Aligning to
   the **genome** instead keeps each read's position a fixed, transcriptome-
   independent fact, and turns "which transcript does this support?" into an
   **explicit, inspectable** projection step where ambiguity can be represented
   rather than buried. This is the difference between *accounted-for* and
   *hidden* uncertainty — a theme the rest of this document holds itself to.

The system is therefore a **transcriptome-relative re-interpretation engine
over a fixed body of genomic translation evidence.**

---

## 2. Upstream data preparation (context, assumed to exist)

The scoring engine consumes the products of an upstream pipeline. That pipeline
is summarised here because its outputs and its fallibilities are part of the
design contract.

1. **Sample discovery.** Candidate Ribo-Seq runs are detected from published
   data, largely by sample-name heuristics. *This is fallible:* non-Ribo-Seq
   samples (RNA-seq, etc.) can contaminate the set. The scoring engine must not
   assume every input sample is genuine Ribo-Seq (see §7, QC gate).

2. **RPF cleaning.** UMIs and adapters are trimmed to recover "clean RPFs."
   *This is fallible:* mis-trimming distorts read length and 5′ position.

3. **Deduplication + counting.** Within each sample, identical clean RPF
   sequences are collapsed to a single unique sequence with an occurrence
   **count**. This is the key space reduction: thousands of samples share many
   identical fragments.

4. **Per-study → global count matrix.** Per-study matrices (unique fragment ×
   sample → count) are aggregated into a global matrix.

5. **Alignment.** The unique fragment sequences are aligned to the **genome**.
   The unique-read store is **partitioned by nucleotide prefix** — the current
   upstream pipeline uses **257 prefix partitions**, giving 257 FASTAs and, after
   alignment, 257 BAMs. *The partition count is an implementation detail of the
   current read store, not an architectural invariant:* the scorer consumes
   partitioned unique-read alignments and must not depend on the number 257.

The two durable products handed to the scoring engine are therefore a **count
matrix** (`unique_read → (sample → count)`) and **genomic alignments** of those
unique reads.

> **Everything above is versioned.** The unique reads and counts are
> observations *after a particular preprocessing pipeline*. Trimming, UMI
> handling, contamination detection, sample inclusion, and deduplication logic
> can all improve over time. See §5 (L0a) — these products are immutable *within
> a named preprocessing/data-release version*, not eternal.

---

## 3. Core principle: the cache boundary is the *junction set*

The single most important design decision is **where the boundary sits between
work that is done once and work that is repeated per transcriptome release.**

Naively one might place that boundary at "the genome alignment" and treat
everything genomic as transcriptome-invariant. **This is wrong**, for a reason
specific to Ribo-Seq:

- Aligners are **splice-aware**: the set of annotated splice **junctions** fed
  to the aligner index determines which reads map and where.
- RPFs are **short** (~26–34 nt). A read straddling a junction may have only a
  handful of bases on the short arm, leaving it nearly unalignable unless the
  junction is already known to the aligner.
- Therefore adding newly annotated junctions can change *which reads align at
  all* and *where they land* — the BAM itself changes.

So among the inputs that come from a transcriptome **release**, the junction set
is the one that moves the alignment. But it is not the *only* thing the evidence
substrate depends on. The substrate is also coupled to the **genome
assembly/version**, **aligner parameters**, **multimapping policy**,
**trimming/UMI preprocessing**, **read-length filters**, and the **canonical CDS
set** used for QC/offset learning. Those are held fixed across routine releases,
but they are real dependencies. The precise, non-overclaimed statement is:

> **For a fixed genome assembly, preprocessing/data-release version, aligner
> configuration, multimapping policy, and unique-read set, the junction set is
> the only transcriptome-derived input that changes the alignment evidence.**

Within that frame, the invariance that drives the cache design is:

> **Invariant to ORF / translon annotation changes.** New candidate translated
> regions move no reads. Reuse everything; only re-extract events and re-score.
>
> **NOT invariant to junction-set changes.** A changed junction set changes the
> aligner index, hence the alignment, hence everything downstream of it. The
> unique-read FASTAs must be re-aligned (only the unique reads — never raw
> data), and the profile rebuilt where junctions changed.

This also resolves an otherwise puzzling question — *what does an updated
transcriptome actually buy us?* The answer is: **the junctions.** They improve
alignment sensitivity (especially for short, junction-spanning RPFs) and they
are required to interpret frame across introns (§6.2).

---

## 4. Coordinate and identity principles

The pipeline of abstractions is:

```
observed reads
  → genomic evidence            (alignment, multimapping-aware)
  → QC / offset calibration     (per sample × length)
  → aggregate P-site profile
  → transcriptome-relative event extraction
  → event-level scoring         (genomic, per de-duplicated event, once)
  → feature-level reporting
```

Three principles cut across the layers:

- **Genome-native evidence.** All evidence is held in genomic coordinates.
  Transcript coordinates are derived on demand, never the storage frame.
- **Multimapping is first-class.** Genomic placement of a short RPF is itself
  uncertain; that uncertainty is represented, not resolved away (§5, L0b).
- **Score each event once, under a formal key.** Identity of "the same event"
  is defined explicitly (§6.1), not left to coordinate coincidence.

---

## 5. Layered data model

```
 L0a  Read matrix          unique_read → sequence ; unique_read → (sample → count)
        └─ immutable WITHIN a named preprocessing/data-release version

 L0b  Alignment            unique_read → {all compatible alignments, with metadata}
        └─ rebuilt PER RELEASE iff the junction set changes; over unique reads only

 L1   Dynamic-aggregation cache (genomic, per-sample) — RAW 5′ evidence
        positional stream:  (genomic 5′ pos, length, sample) → count [+ aln metadata]
        junction stream:    (donor, acceptor, length, sample) → count
        └─ pre-offset; preserves every dimension one might aggregate over

 L2   QC gate              periodicity / length → per-(sample,length) artefacts:
        authenticity, periodicity, offset, offset confidence, include/exclude, reason

 L3   Default aggregate profile
        offset-applied P-site profile, summed over the passing (sample, length)

 L4   Scoring             genomic, per de-duplicated event (formal key), once
        └─ dynamic re-aggregation from L1 for complex / ambiguous loci
```

### L0a — Read matrix (immutable within a version)

`unique_read → sequence` and `unique_read → (sample → count)`. These are
observations — the sequences seen and how often, per sample — and do **not**
depend on the genome or transcriptome. They are **immutable within a named
preprocessing/data-release version**. If trimming, UMI handling, contamination
detection, sample inclusion, or deduplication logic changes, that is a **new
L0a version**, not a mutation of the old one. Everything downstream records the
L0a version it was built from.

### L0b — Alignment (per release iff junctions change), multimapping-aware

The unique-read sequences aligned to the genome. **Multimapping is first-class
here, and this is one of the most important design commitments in the
document.** Short RPFs map to multiple places — paralogues, repeats,
pseudogenes, homologous CDS, rRNA/contamination, immunoglobulin regions — and
collapsing each read to a single chosen alignment would simply replace *hidden
transcriptome uncertainty* with *hidden genomic-placement uncertainty*, which
this design exists to avoid.

L0b therefore stores, for each unique read, **all compatible alignments with the
metadata needed to apply any placement policy downstream**:

```
unique_read_id
chrom, start, end, strand
CIGAR
MAPQ
NH / number_of_alignments
primary_or_secondary
mismatch_count / alignment_score
splice_junctions_crossed        (donor/acceptor coords, if any)
alignment_weight                 (policy input, not a baked-in decision)
```

L0b makes **no placement decision**. It records the placement *evidence*. L1
and the scorer decide whether to use unique mappers only, fractional counts,
MAPQ-weighted counts, or locus-specific policies. Because alignment is
splice-aware, L0b is a function of the junction set and is the only expensive
artefact that may be rebuilt per release — and only over unique reads.

### L1 — Dynamic-aggregation cache (genomic, per-sample, RAW 5′)

The central substrate. It is **raw, pre-offset, 5′-end alignment evidence** —
*not* P-site, and never "5′ or P-site". P-site projection depends on L2 offset
calibration and belongs strictly downstream (L3). Keeping L1 pre-offset means
that if offset calibration changes, **L1 does not need rebuilding** — only L3
does. This removes any ambiguity about what invalidates the cache.

L1 deliberately **preserves the dimensions a query might aggregate over** —
genomic 5′ position, read length, and **sample** — rather than pre-collapsing
them, and carries the relevant **alignment metadata** from L0b so that placement
policy can be applied at query time. It has **two parallel streams**:

- **Positional stream:** `(genomic 5′ pos, length, sample) → count`.
- **Junction stream:** `(donor, acceptor, length, sample) → count` — counts of
  reads spanning each splice junction.

It is held **factored** wherever possible — the version-pinned per-sample counts
(L0a) joined to the per-release unique-read loci (L0b) — so that per-sample
resolution at thousands of samples remains affordable. Materialising the full
dense tensor is not required.

**Why per-sample (not metadata groups).** Translation/expression programs do
not map cleanly onto sample metadata (study, tissue, condition). At a hard
locus the signal may live in a *data-defined* subset of samples that no
metadata label captures. Dynamic aggregation must therefore be able to select
samples **on the data**, which requires the cache to retain individual-sample
resolution.

### L2 — QC gate → per-(sample, length) artefacts

The QC gate decides which `(sample, length)` pairs carry genuine RPF signal and
calibrates their P-site offsets. Because the whole design is about making
uncertainty explicit, **"usable" must not be a black-box verdict.** L2 emits
*separate, inspectable artefacts* per `(sample, length)`:

```
authenticity_score
periodicity_score
offset
offset_confidence
include_exclude_decision
reason_for_exclusion
```

The authenticity signal uses the signatures of real Ribo-Seq — a characteristic
**read-length** distribution and strong **triplet periodicity** over
known-coding regions. The offset is the length-dependent shift from the 5′ end
to the P-site that best aligns reads into one frame over canonical CDS.

> **Named conservative bias / circularity.** Authenticity and offsets are learnt
> on annotated **canonical CDS**. This makes the gate *circular* in a specific
> way: unusual but real biology — tissue-specific translation, short ORFs,
> non-AUG initiation, stress-state reprogramming — can be scored as "bad
> Ribo-Seq" and filtered out. This is acceptable for a **conservative release**,
> but it is a deliberate, named bias, not a neutral filter. The separate
> artefacts above exist precisely so a downstream consumer can see *why* a
> `(sample, length)` was excluded and override the policy if a use case needs to.

L2 runs on L1 and must complete in reasonable time, since it gates everything
downstream.

### L3 — Default aggregate profile

For the common case, the passing `(sample, length)` data is offset-shifted to
**P-site** positions and summed into a **single genomic profile**. This is the
substrate for the bulk of event scoring; most loci want "the aggregate of all
read lengths that pass QC." L3 is the *only* layer where offsets are applied, so
it is the only layer invalidated by an offset-calibration change.

### L4 — Scoring (see §6)

---

## 6. Scoring model

### 6.1 Score genomically, per de-duplicated event, once — with a formal key

Events are scored in **genomic coordinates**, and each distinct genomic event is
scored **exactly once** — never once per transcript that contains it. This
avoids redundant work and avoids inflating or fragmenting evidence across
overlapping isoforms.

But "the same event" is not always obvious from coordinates alone. Two ORFs can
cover the same genomic bases in different **phases**; a spliced ORF and a
retained-intron ORF can share a genomic span but follow different **splice
paths**. De-duplication is therefore defined over **formal event keys**, not raw
spans:

```
initiation_event_key  = (chrom, pos, strand, start_codon, context)
termination_event_key = (chrom, pos, strand, stop_codon)
junction_event_key    = (donor, acceptor, strand, phase_in, phase_out)
elongation_event_key  = (chrom, span(s), strand, phase, splice_path_id)
```

An elongation event that crosses a junction is **path-aware**: its key includes
the splice path, so two ORFs sharing genomic bases but differing in phase or
splice path are correctly treated as **different events**. Without these keys,
"de-duplicated event" could be read too loosely and silently merge distinct
biology.

An **event** is thus a keyed, scorable genomic feature derived from candidate
translated regions. Event extraction is the transcriptome-dependent part of
scoring — the candidate ORF/translon set determines which events exist — but the
*evidence* each event is scored against is the transcriptome-independent genomic
profile. A separate composition step fans per-event scores back to the
per-translon (feature) level for reporting. **Scoring makes no annotation
decisions; it only quantifies evidence.**

### 6.2 Frame is genomic arithmetic — except across introns

For elongation, "in-frame" is determined by a read's position relative to the
event's start and phase. Within a single exon this is pure genomic arithmetic.
**The exception is an intron inside the scored window:** genomic distance then
diverges from transcript distance, and frame must be carried across the junction
using the splice path in the event key. This is the scoring-time role of
junctions, and it is why "nearby splicing affects scores" is in scope.

The same divergence applies to the init/term leader/UTR flanks, not just
elongation frame continuity: a start or stop codon near a splice site has its
leader/UTR flank on the far side of an intron. This is implemented (as of
2026-07-14, see event_scoring_model.md v1.1) as a single-path flank
projection (`scoring.aspects._project_flank`) — genuinely fixes the common
case, but does not resolve the isoform-of-origin ambiguity when candidate
transcripts disagree on the immediate flanking exon (§6.3 below still
applies: exactly one path is followed, not enumerated).

### 6.3 Isoform-of-origin is deferred

Determining *which* isoform produced an observed signal is hard — especially
where no junction-spanning reads are present to disambiguate. The design
deliberately **defers exact isoform attribution.** The committed scope is to
score translation support **genomically** and to **account for nearby splicing**
that would affect a score (frame continuity across introns, and reads arriving
via a junction). Exact isoform attribution is a later target, not part of the
core engine.

### 6.4 Junctions: first-class event *and* covariate

Junction-spanning reads are direct evidence that a junction is used.

- **Primary (required): a first-class event type.** A junction is scored on its
  own spanning-read support, exactly as initiation/termination/elongation events
  are scored on their evidence. Keyed by `junction_event_key`.
- **Secondary (desirable, deferred): a covariate.** It is also valuable to see
  junction evidence *modulate* nearby elongation/frame scores — e.g. a strongly
  supported junction inside an elongation window informs how positions on either
  side are stitched into one frame and how confident the elongation call is.
  This is layered on top of the event-type treatment, not a replacement, and is
  out of scope for the MVP (§8).

### 6.5 Dynamic aggregation for hard loci — with guardrails

Because L1 retains per-sample, per-length resolution, complex or ambiguous loci
can be **re-aggregated on demand** to sharpen signal — e.g. restricting to a
data-defined subset of samples, or to a specific read-length band — rather than
relying solely on the default L3 aggregate. This is how rare/contested events
get a more careful read-out than the bulk path provides.

This power is also a hazard: searching arbitrary sample subsets until a rare ORF
*looks* periodic can manufacture signal out of batch effects, read-length
artefacts, study-specific trimming errors, or contamination. It is statistically
p-hacking. The design therefore imposes a guardrail:

> **Dynamic aggregation may be used to diagnose or rescue ambiguous loci, but
> any data-defined subset must be reported with its selection rule, number of
> samples, number of independent studies, read-length composition, and whether
> the signal replicates across independent studies.**

This matters especially for a GENCODE-style release: a human reviewer must be
able to tell whether support comes from one weird sample or ten independent
studies.

---

## 7. Per-release re-assessment workflow

The cost of re-assessing translation against a new transcriptome release depends
on **whether the junction set changed** (under the fixed genome / preprocessing
/ aligner / multimapping frame of §3):

**ORF-only release (no junction change).**
Reuse L0a, L0b, L1, L2, L3 unchanged. Re-extract events from the new ORF set and
re-score (L4). Cheap — the target case for routine releases.

**Junction-changing release.**
1. Re-align the unique-read FASTAs against the genome with the new junction set
   → new L0b. (Unique reads only; raw data and L0a counts untouched.)
2. Rebuild L1 where junctions changed (positional + junction streams).
3. Re-run the L2 QC gate (offsets/authenticity can shift as alignment shifts).
4. Re-aggregate L3 over the passing data.
5. Re-extract events and re-score (L4).

Work is bounded by the unique-read count, not raw data volume, and is
concentrated where junctions actually changed.

**Preprocessing / genome / aligner change.** Treated as a new L0a (and therefore
L0b) **version**, not a release re-assessment — a full rebuild, by definition.

---

## 8. Minimum viable implementation

The sections above describe the elegant final system. To prevent the project
stalling on the "correct" full L1/L4 design, the committed MVP is:

1. **L0a** unique-read matrix exists (version-pinned).
2. **L0b** genome alignment with **all alignments retained** and their metadata
   (multimapping-aware from day one — retrofitting it later is far harder).
3. **L1** raw 5′ positional cache, per sample and length, with alignment
   metadata; pre-offset.
4. **L2** CDS-based periodicity + offset calibration, emitting the separate
   per-(sample, length) artefacts of §5 (L2).
5. **L3** aggregate P-site profile over passing data.
6. **L4** score initiation / elongation / termination events under formal keys.
7. **Junctions** scored as first-class events.
8. **Deferred:** junction-as-covariate (§6.4); exact isoform attribution
   (§6.3); incremental (rather than wholesale) L1 rebuild after junction change.

The two MVP commitments most easily skipped but most expensive to retrofit are
**multimapping-aware L0b** and **formal event keys** — both should be in from
the start even if their *policies* (placement weighting, dedup tie-breaks) start
simple.

---

## 9. Design invariants and non-goals

**Invariants**
- One evidence substrate (the profile), one scorer. There are **not** two
  scoring tracks; all scores — and any future scores — are computed on the
  resulting profiles regardless of how a locus is reached.
- Evidence is collected once per L0a version and re-interpreted cheaply; the
  only per-release recomputation is triggered by junction-set changes.
- Scoring is genomic and per-keyed-event-once; transcript membership is a
  composition/reporting concern, not a scoring concern.
- Uncertainty is **represented, not hidden** — at every layer: genomic placement
  (multimapping, L0b), authenticity/offset (L2 artefacts), and data-defined
  aggregation (selection-rule reporting, §6.5).

**Non-goals (for the core engine)**
- Exact isoform-of-origin attribution (deferred — §6.3).
- Trusting sample metadata for aggregation (rejected — §5, L1).
- Treating every input sample as genuine Ribo-Seq (rejected — the QC gate, §7,
  exists precisely because discovery is fallible).
- Collapsing multimappers to a single placement (rejected — §5, L0b).

---

## 10. Glossary

- **RPF (Ribosome Protected Fragment):** the short RNA fragment shielded by a
  translating ribosome; the unit of Ribo-Seq evidence.
- **Unique read:** a distinct clean-RPF sequence, stored once with a per-sample
  occurrence count.
- **Count matrix:** `unique_read → (sample → count)`.
- **Partition:** one prefix-defined slice of the unique reads; the current store
  uses 257, but that number is not architectural.
- **Multimapper:** a unique read with more than one compatible genomic
  alignment; represented with all alignments + metadata, never pre-resolved.
- **Profile:** per-position (and per-junction) read counts that quantify
  translation evidence; the single substrate all scores consume. Held raw and
  pre-offset at `(5′ pos, length, sample)` in L1; offset-applied and aggregated
  in L3.
- **P-site offset:** the length-dependent distance from a read's 5′ end to the
  ribosome P-site; calibrated per `(sample, length)` in L2.
- **Periodicity:** the triplet (every-3-nt) read pattern produced by ribosomes
  stepping one codon at a time; the primary signal of genuine translation.
- **Frame / phase:** position modulo 3 relative to an ORF's reading frame.
- **Event:** a keyed, scorable genomic feature — initiation, termination,
  elongation span, or junction — derived from candidate translated regions
  (§6.1).
- **Splice path:** the ordered exon chain an elongation event follows; part of
  the elongation event key, so path-distinct ORFs over shared bases stay
  distinct.
- **Translon / feature:** a candidate translated region; events compose back up
  to translons for reporting.
- **Junction set:** the annotated splice junctions supplied to the aligner; the
  sole *transcriptome-derived* alignment coupling under a fixed
  genome/preprocessing/aligner/multimapping frame (§3).
- **L0a version:** a named preprocessing/data-release version; L0a (and hence
  everything downstream) is immutable within it.

---

## 11. Open questions

- **L1 storage strategy.** How factored vs materialised the cache is, and whether
  L1 is partitioned by chromosome and/or by the upstream prefix partitions, is a
  build-vs-query trade-off left to implementation.
- **Multimapping placement policy.** L0b retains all alignments; the default
  policy L1/L4 apply (unique-only vs fractional vs MAPQ-weighted vs
  locus-specific) is not yet fixed.
- **Junction-change detection granularity.** Whether L1 rebuild after a
  junction-set change can be made incremental (only reads near changed
  junctions) or is done wholesale per release.
- **Covariate formulation.** The precise way junction support modulates nearby
  elongation/frame scores (§6.4) is desirable but not yet specified.
- **Event-key canonicalisation.** Exact normalisation rules for each key
  (e.g. how much initiation context, how splice paths are identified and
  compared) need pinning down before dedup is implemented.

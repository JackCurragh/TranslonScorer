# TranslonScorer Score-First Frame Plan

Branch: `local/read-assignment-prototype`
Status: local-only working branch; do not push until the design gates below are agreed.

Progress log: see `docs/gate_progress.md`.

## 1. Design reset

TranslonScorer should stay score-first. The central product is not per-read assignment by itself; it is uncertainty-aware translon evidence:

- ORF/translon scores: SRU, HRF, coverage, NZC, composite score.
- Frame posteriors: p0, p1, p2, entropy, frame-weighted coverage.
- Assignment and identifiability labels only where transcript ambiguity affects the scores.
- Position-level frame posterior exports for downstream path models such as RDG-Flux.

Read assignment is a means to avoid contaminated transcript profiles. It becomes a headline method only if the data show that frame compatibility resolves enough ambiguous signal to change scores or transcript-level conclusions.

## 2. Minimal combinations to test

The combinations should isolate one added assumption at a time.

| Label | Input evidence | Frame model | Isoform model | Output being tested | Why it exists |
| --- | --- | --- | --- | --- | --- |
| A | BigWig/profile | none | none | existing raw scores | Baseline: current TranslonScorer behavior. |
| B | BigWig/profile | linear+HMM | none | frame-weighted scores and entropy | Tests whether frame support improves scoring without read-level complexity. |
| C | BAM/read profiles | linear+HMM | fixed convention | same as B, with read-length offsets retained | Tests read-level input and memory without introducing joint EM. |
| D | BAM/read equivalence classes | none or fixed frame | frame-blind EM | transcript-weighted profiles and scores | RSEM-style control for isoform ambiguity. |
| E | BAM/read equivalence classes | linear+HMM or local frame posterior | frame-aware EM | isoform-weighted, frame-weighted scores | Tests the only true joint-inference claim. |
| F | Any transcript-level output | same as selected model | same as selected model | aliased/data-limited/identifiable labels | Required guardrail for transcript-specific claims. |

This is the working matrix. We should not expand it until each row has a clear failure condition.

## 3. Immediate gates

### Gate 0: Define score outputs

Before more modelling, define the output schema for scores:

- `orf_id`, `tran_id`, `gene_id`, coordinates, ORF type.
- Raw score terms: `rise_up`, `step_down`, `hrf`, `avg`, `nzc`, `score`.
- Frame terms: `frame_posterior_mean`, `frame_entropy_mean`, `frame_weighted_count`, `frame_method`.
- Ambiguity terms: `assignment_weight`, `assignment_entropy`, `identifiability_class`.
- Provenance terms: input type, offset source, annotation version, method version.

Pass condition: one schema supports A-F without special-case output files.

### Gate 1: Does frame support improve scoring?

Run A vs B on a small frozen panel:

- annotated CDS sanity loci,
- known uORF loci,
- frameshift or overlapping ORF examples if available,
- negative/background UTR examples.

Metric: improvement in rank separation or calibration of known translated vs background candidates, not only `delta_p0`.

Pass condition: frame-weighted scores improve or clarify calls on the panel. If not, frame support remains a browser/QC track rather than a scoring default.

Frame-assignment comparison inside this gate:

- `linear` learns global or read-length-specific leakage from CDS interiors, such as `[0.70, 0.10, 0.20]`, builds the cyclic confusion matrix, and inverts it per codon. Local observed `[80, 10, 10]` can therefore become adjusted support close to `[100, 0, 0]`.
- `linear+hmm` is the same linear correction followed by HMM smoothing along each transcript. It is a smoothing layer, not a separate leakage model.
- `latent` is the frame-only EM model. It starts from the learned leakage matrix and now fits count-space per-codon latent frame support with that calibration held fixed by default. Adaptive leakage updates are intentionally regularized/opt-in because re-learning leakage from all codons can let sparse or non-canonical regions corrupt the CDS calibration.

The executable comparison is `translonscorer compare-frame-methods --methods linear,latent`. Add `linear+hmm` when the question is whether smoothing improves calibration beyond direct correction.

This comparison now writes CDS-interior validation summaries. The truth label is the transcript frame implied by annotated CDS start, `cds_start % 3`, not always frame 0.

Important real-data distinction: `p0/p1/p2` are conditional frame probabilities, not a calibrated translated/background posterior. Frame-support output therefore carries `support_evidence`, `support_gated_p0..2`, and `secondary_frame_mass`. This lets low-support high-confidence calls and possible overlapping translation remain visible instead of being hidden by HMM smoothing.

### Gate 2: Does read-level processing reproduce BigWig/profile results?

Run B vs C on the same samples where BigWig and BAM are available.

Pass condition: C agrees with B within a defined tolerance, while preserving read length, offset, and mapping metadata. This validates the read-level substrate before EM.

### Gate 3: Is frame-informative ambiguity common enough?

Cheap empirical check first:

- Build read equivalence classes for multi-mappers or isoform-shared loci.
- For each candidate assignment, compute transcript position and candidate frame.
- Count reads/loci where candidate assignments are frame-discordant.
- Stratify by gene, ORF type, coding structure difference, and read support.

Pass condition: enough ambiguous signal is frame-discordant that row E could change scores or transcript calls. If the rate is low, joint EM becomes optional and the paper centers on frame-aware scoring plus identifiability.

### Gate 4: Does frame-aware EM beat frame-blind EM where it should?

Run D vs E first on synthetic frame-informative mixtures, then curated loci.

Pass condition: E improves assignment or score calibration on frame-informative cases and does not regress frame-uninformative cases. If not, report the negative result and keep frame-blind weighting plus identifiability.

## 4. Development sequence

1. Stabilize score schema.
2. Wire linear frame correction, optional HMM smoothing, and latent frame EM into score outputs cleanly.
3. Add comparison runner for A vs B on a small frozen panel.
4. Make BAM/profile generation reproduce BigWig/profile outputs with retained read-length metadata.
5. Implement the cheap frame-disambiguation analysis.
6. Decide whether frame-aware EM is a core deliverable or an optional branch.
7. If positive, implement frame-blind EM first, then add the frame compatibility term.
8. Add identifiability labels before emitting transcript-specific claims as final outputs.

## 5. First local deliverables

These are the first concrete branch tasks:

- `docs/score_first_frame_plan.md`: this design plan.
- Score schema implemented in `TranslonScorer.pipeline.score_schema`.
- Small panel manifest validation implemented as `translonscorer validate-panel`; `score-compare-frame` can merge it before scoring.
- Paired raw vs frame-weighted score runner implemented as `translonscorer score-compare-frame`.
- Frame-only correction comparison implemented as `translonscorer compare-frame-methods`.
- RDG-Flux v1 substrate export implemented as `translonscorer export-rdg-flux`.
- Cheap frame-disambiguation runner implemented as `translonscorer frame-disambiguation`.
- Read-assignment prototype implemented as `translonscorer compare-read-assignment`; see `docs/read_assignment_plan.md`.

## 6. Non-goals for this branch

- No full multi-sample joint inference.
- No differential translation statistics.
- No ORF discovery redesign.
- No commitment that joint EM is the paper centerpiece until Gate 3 passes.

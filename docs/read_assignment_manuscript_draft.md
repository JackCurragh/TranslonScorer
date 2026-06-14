# TranslonScorer Read-Origin Assignment Draft

## Working Title

Uncertainty-preserving read-origin assignment for frame-aware Ribo-seq substrate generation.

## Thesis

Ribo-seq read assignment is not just transcript quantification with shorter reads. The origin of a ribosome footprint can be ambiguous at two levels: which genomic locus generated the read, and which transcript isoform within that locus carried the translating ribosome. Standard RNA-seq quantifiers already solve an important part of this problem by assigning multi-compatible fragments probabilistically, but they do not use the most distinctive evidence in Ribo-seq: reading-frame structure.

TranslonScorer should therefore produce posterior evidence, not hard calls:

```text
P(locus, transcript, position, frame | read/profile evidence)
```

The v1 read-assignment layer does not need to be a complete ORF caller. It needs to generate disambiguated, uncertainty-preserving signal that downstream scoring and graph/path models can consume.

## Why This Is Needed

RNA-seq established the core assignment problem: a read can be compatible with multiple genes or isoforms, and transcript abundance must be inferred rather than counted naively. RSEM made this explicit with EM-based transcript quantification, including multi-mapping reads and isoform ambiguity ([Li and Dewey 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3163565/)). Salmon and kallisto later pushed the same idea into fast equivalence-class or mapping-based quantification, with bias models and pseudo/quasi-mapping speedups ([Patro et al. 2017](https://www.nature.com/articles/nmeth.4197), [Bray et al. 2016](https://www.nature.com/articles/nbt.3519)).

Ribo-seq adds a second source of information. A footprint compatible with two transcript origins may land in different transcript-coordinate frames or different translated structures. That frame compatibility should influence assignment when local frame evidence is strong. Existing Ribo-seq tools focus mainly on translated ORF detection, periodicity, or ORF-level quantification rather than a general posterior substrate layer for transcript/locus origin. Relevant reference points include Rp-Bp's Bayesian use of ribosome periodicity ([Malone et al. 2017](https://academic.oup.com/nar/article/45/6/2960/2953491)), RiboCode's de novo translatome annotation from Ribo-seq ([Xiao et al. 2018](https://academic.oup.com/nar/article/46/10/e61/4925760)), and ORFquant's ORF-level translation quantification across transcript mixtures ([Calviello et al. 2020](https://www.nature.com/articles/s41594-020-0450-4)).

The gap TranslonScorer occupies is narrower and useful: before ORF/path inference, assign read-origin and frame evidence probabilistically, keep uncertainty visible, and report intrinsic aliasing rather than forcing a best origin.

## Model Sketch

For a read `r`, define candidate origins:

```text
c = (locus_id, transcript_id, transcript_pos)
```

The prototype estimates:

```text
P(c | r, data) proportional to
  theta_target
  * alignment_likelihood(r | c)
  * frame_likelihood(transcript_pos | local frame posterior)
```

where `theta_target` is the abundance of the selected target level:

- `tran_id`: transcript isoform assignment
- `locus_id`: genomic-origin assignment
- `locus_id,tran_id`: joint locus/transcript assignment

The frame term is deliberately conservative:

```text
frame_likelihood = P(signal frame = transcript_pos mod 3 | local frame evidence)
```

It uses transcript-coordinate signal frame because TranslonScorer frame-support tables store `p0/p1/p2` by transcript coordinate residue class. It does not use CDS-relative frame. When `support_evidence` is low, frame likelihood is blended back to neutral so weak frame calls do not create false assignment confidence.

For target-level assignment, candidate-origin rows are collapsed to `(read, assignment_target)` before posterior inference. The target likelihood uses the best internal candidate likelihood, and the inferred target posterior is then distributed back over candidate-origin rows. This avoids an annotation-multiplicity artifact where a locus with many transcript models would otherwise become more likely simply because it has more candidate rows.

## Method Ladder

The comparison ladder isolates assumptions one at a time.

| Method | Interpretation | Expected role |
| --- | --- | --- |
| `unique` | only assign reads with one candidate origin | conservative lower bound |
| `fractional` | split by alignment/candidate likelihood within read | frame-blind heuristic |
| `frame_fractional` | split by alignment likelihood times frame compatibility | tests frame evidence without abundance EM |
| `em` | abundance EM over the chosen target | RSEM-style control |
| `frame_em` | abundance EM with frame compatibility | TranslonScorer read-origin model |

The key falsification conditions are:

- `frame_em` must match `em` when frame evidence is uninformative.
- `frame_em` must beat `em` in frame-discordant mixtures.
- `em` must beat fractional assignment when unique reads anchor origin abundance.
- all methods must remain uncertain when origins are intrinsically aliased.

## Evidence Incorporated So Far

### Synthetic Gate

The deterministic simulation in `notebooks/read_assignment_prototype_simulation.py` covers four regimes.

| Regime | Result | Interpretation |
| --- | --- | --- |
| frame-discordant isoforms | `frame_em` reaches near-complete posterior mass on the true origin | frame compatibility works where it should |
| abundance-resolvable isoforms | `em` and `frame_em` outperform fractional assignment | unique/high-confidence reads anchor abundance |
| genomic multi-locus reads with unique anchors | `em` and `frame_em` outperform fractional assignment at `locus_id` target level | same model handles locus ambiguity |
| aliased/no-evidence origins | methods stay near chance and entropy remains high | identifiability guardrail works |

Current synthetic summary:

- `frame_em` on frame-discordant isoforms: soft true posterior approximately 1.000.
- `em` on abundance-resolvable isoforms: soft true posterior approximately 0.911.
- `em` on genomic multi-locus unique-anchor simulation: soft true posterior approximately 0.948.
- aliased no-evidence case: posterior remains 0.500 with normalized entropy 1.000.

### Real Annotation-Side Evidence

The TransCODE pilot iRibo candidate-transcript inference output is not read-level evidence, but it identifies a real ambiguity pressure point. From 132,419 iRibo candidate records:

- 126,695 have transcript context assigned, 95.68%.
- 13,918 have ambiguous inferred feature type, 10.51%.
- 18,568 unique transcripts are represented.
- compatible transcript count is exactly 1 for 56.59% of candidates.
- compatible transcript count is greater than 1 for 39.09% of candidates.
- per-sample ambiguous feature-type rates range from 2.62% to 12.62%.

This supports the story that transcript context is frequently not a single deterministic label, even before read-level multi-mapping is considered.

### Real Candidate-Origin Probe

The first real read-substrate probe uses a transcriptome-aligned RiboMetric sample BAM:

```text
/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/subsampled_1_percent.bam
```

The BAM reference names encode transcript and gene IDs, so the probe can construct candidate origins without a full genomic projection. This is not final biological validation, but it is a real candidate-origin table with repeated read names across compatible transcripts and genes.

Current result:

- 35,804 read keys and 68,417 weighted reads.
- mean candidate transcript targets per read: 7.799.
- mean candidate locus targets per read: 1.863.
- transcript-ambiguous read keys: 88.50%.
- locus-ambiguous read keys: 20.78%.
- transcript-level EM reduces normalized assignment entropy from 0.846 to 0.384.
- locus-level EM reduces normalized assignment entropy from 0.263 to 0.084.
- `frame_em` exactly matches `em` without a frame-support join, confirming neutral frame behaviour.
- transcript-level EM classes by weighted reads: 15.38% unique, 27.73% resolved, 24.15% ambiguous, and 32.75% data-limited.
- locus-level EM classes by weighted reads: 73.69% unique, 15.66% resolved, 6.13% ambiguous, and 4.52% data-limited.

The real probe changes the emphasis of the story. Transcript ambiguity is not a rare corner case in this substrate; it is the dominant state. Locus ambiguity is smaller but still substantial. That makes uncertainty reporting a core deliverable even before we know how often frame evidence resolves the ambiguity.

Evidence artifact:

- `docs/read_assignment_story_evidence.md`
- `notebooks/read_assignment_story_evidence.py`
- `notebooks/read_assignment_story_evidence.png`

### Matched Frame-Support Probe

The first matched frame-aware run joins the real candidate-origin table to two frame-support tracks:

- `unique_transcript`: independent support from reads with exactly one compatible transcript target.
- `fractional_all`: exploratory support from all candidate origins with fractional weights, useful as an upper-bound sensitivity check but not as independent evidence.

The independent track is sparse:

- 3,933 supported frame rows across 1,623 transcripts.
- total unique-support count: 10,521.0.
- only 9.8% of weighted candidate reads have matched independent frame support.

With this independent support, `frame_em` changes transcript-level assignment only modestly relative to `em`:

- weighted mean target-posterior TV distance: 0.0018.
- weighted reads with TV distance at least 0.10: 0.65%.
- after a minimal quality gate (`total_count >= 10`, `support_evidence >= 0.2`), weighted mean TV distance is 0.0016 and the weighted TV>=0.10 fraction remains 0.65%.

The exploratory `fractional_all` support moves more mass:

- transcript-level weighted mean TV distance: 0.0318.
- weighted reads with TV distance at least 0.10: 6.17%.
- after the same gate, weighted mean TV distance is 0.0281 and the weighted TV>=0.10 fraction is 5.76%.

Because `fractional_all` reuses ambiguous reads to build the frame support, it should be treated as a sensitivity bound rather than proof. The independent high-shift cases are also concentrated: 57 high-shift read keys, 443.0 weighted reads, mean TV 0.205, and mean 7 candidate transcripts all cluster around `LINC01783` / neighbouring pseudogene-like transcript models. Targeted case review shows independent frame support for only one transcript in the cluster, `ENST00000415386.2`, with 29 unique-support counts across 3 codons.

This does not falsify frame-aware assignment. It says the current real substrate is support-limited. The honest claim is that the model is ready, frame terms behave neutrally when support is absent, weak frame evidence can now be gated, and a production matched sample is needed before claiming broad real isoform deconvolution.

Evidence artifacts:

- `docs/read_assignment_matched_frame_probe.md`
- `docs/read_assignment_linc01783_case_review.md`
- `notebooks/read_assignment_matched_frame_probe.py`
- `notebooks/read_assignment_locus_case_review.py`

### Scaled 10% BAM Probe

The larger local RiboMetric `subsampled_10_percent.bam` probe asks whether the sparse-support conclusion survives a roughly 10x input scale-up.

Current result:

- 358,510 read keys and 828,006 weighted reads.
- transcript-ambiguous read keys: 88.53%.
- locus-ambiguous read keys: 21.25%.
- independent frame support: 33,269 rows across 5,540 transcripts.
- weighted candidate reads with independent frame support: 12.8%.

At this scale, independent frame support moves more posterior mass:

- ungated independent `frame_em` versus `em`: weighted mean TV distance 0.0257; 3.00% of weighted reads move by at least 0.10.
- gated independent `frame_em` versus `em`: weighted mean TV distance 0.0231; 2.93% of weighted reads move by at least 0.10.

The important caveat is where that movement goes. The largest high-shift cluster remains the `ENSG00000228549` / `LINC01783` lncRNA ambiguity block, accounting for 97.3% of weighted independent high-shift reads in the exported high-shift locus table. A coding-associated scan finds 111 high-shift read keys involving protein-coding candidates, 1,989 weighted reads total. Of these, 63 read keys and 1,254 weighted reads have gated local frame support; the rest are abundance-propagated EM shifts. The top coding-associated signal (`CIC`) is broad multi-locus ambiguity: 4 reads, 671 weighted reads, 46 candidate gene labels, and 28.5 compatible transcript rows per read, with only 1 weighted read carrying gated local frame support.

A seed coding review panel now ranks direct local-frame-supported candidates separately from abundance-propagated shifts. This gives a concrete curation path: start with lower-locus, direct-support candidates such as `POGK`, `H1-2`, `OAZ2`, `H1-4`, `MT-CO1`, and `NDUFS5`, while routing broad multi-locus cases to mappability review. Initial reviews of `POGK` and `OAZ2` show both are still cross-gene cases rather than clean within-gene isoform validation examples.

The scaled result sharpens the manuscript story. The frame term is not inert; it produces measurable assignment shifts with more data. But the first real shifts are dominated by mapping/annotation ambiguity, not yet by clean isoform-origin resolution. That makes identifiability and locus triage a central contribution rather than a footnote.

Evidence artifacts:

- `docs/read_assignment_ribometric_10pct.matched_frame_probe.md`
- `docs/read_assignment_ribometric_10pct.ENSG00000228549.high_shift_examples.locus_review.md`
- `docs/read_assignment_ribometric_10pct.coding_shift_scan.md`
- `docs/read_assignment_ribometric_10pct.CIC.coding_shift_scan.locus_review.md`
- `docs/read_assignment_ribometric_10pct.PTCH1.coding_shift_scan.locus_review.md`
- `docs/read_assignment_ribometric_10pct.POGK.coding_shift_scan.locus_review.md`
- `docs/read_assignment_ribometric_10pct.OAZ2.coding_shift_scan.locus_review.md`
- `docs/read_assignment_ribometric_10pct.coding_locus_panel_seed.md`
- `notebooks/read_assignment_scaled_bam_probe.py`
- `notebooks/read_assignment_coding_shift_scan.py`
- `notebooks/read_assignment_high_shift_locus_review.py`

## What Would Flesh Out The Story

The story needs one more layer of real data: read-level candidate-origin tables. The current model proves the inference logic, and the TransCODE pilot results show annotation-side ambiguity, but the central empirical question is how often real reads carry useful disambiguating evidence.
The first candidate-origin table has now been obtained from a transcriptome BAM, a first matched frame-support join has been run, and a scaled 10% local BAM probe confirms measurable but biologically concentrated frame-induced movement. The remaining missing layer is curated coding/dual-coding locus validation and a true production matched BAM/profile pair, so `frame_em` can be evaluated beyond ambiguity triage.

The evidence to obtain is:

| Evidence | Question answered | Minimum viable source |
| --- | --- | --- |
| read candidate-origin table | How many reads are ambiguous at transcript or locus level? | obtained for one transcriptome BAM; still needed for matched production samples |
| frame-discordant candidate rate | How often can frame evidence help assignment? | first sparse matched probe complete; repeat on production matched sample |
| posterior entropy distribution | How many loci are identifiable versus aliased? | assignment outputs grouped by locus/transcript |
| real abundance shifts | Does EM materially alter profiles compared with fractional assignment? | `compare-read-assignment` on real candidates |
| structural aliasing map | Which transcript pairs are indistinguishable independent of read depth? | transcript models, ORF/frame structure, candidate-origin equivalence classes |

## Incorporation Plan

### 1. Candidate Table Acquisition

Use the existing mapped-index/frame-disambiguation machinery to produce a table with:

```text
read_key
locus_id
gene_id
tran_id
tran_start_bam
count or read_weight
candidate_likelihood or mapq
```

This table is the substrate for both transcript and locus assignment. If the source is genomic BAM projection, `read_key` should come from `qname` when possible; if the source is a collapsed/Zarr index, it should come from the stable read row id.

### 2. First Real Assignment Run

Run:

```bash
python3 -m TranslonScorer.cli compare-read-assignment \
  --candidates <real_candidates.parquet> \
  --frame-support <frame_support.parquet> \
  --out-prefix <prefix> \
  --methods unique,fractional,frame_fractional,em,frame_em \
  --abundance-key tran_id \
  --max-iter 300 \
  --write-assignments
```

Then repeat with:

```bash
--abundance-key locus_id
--abundance-key locus_id,tran_id
```

The first pass does not need truth labels. It should inspect:

- ambiguous read fraction
- mean candidate origins per read
- mean max posterior
- normalized assignment entropy
- high-confidence fraction
- EM convergence
- transcript/locus abundance shifts relative to fractional assignment

### 3. Identifiability Labels

Add per-target labels:

| Label | Operational definition |
| --- | --- |
| unique | candidate set has one origin |
| abundance_resolved | EM posterior high, entropy low, supported by unique/high-confidence reads |
| frame_resolved | frame-aware posterior high and frame-blind posterior ambiguous |
| aliased | candidate origins have indistinguishable likelihoods and frame support |
| data_limited | posterior remains broad because depth/support is insufficient |
| unsupported_frame | frame posterior is weak; frame term was neutralized |

These labels should land in the score schema and later in the sparse RDG-Flux joint posterior export.

## Draft Claim Wording

TranslonScorer extends transcript-level assignment ideas from RNA-seq to Ribo-seq by treating reading-frame compatibility as an additional likelihood term, while preserving posterior uncertainty. In simulations, the model reduces to ordinary abundance EM when frame evidence is uninformative, improves over frame-blind assignment when candidate origins are frame-discordant, and preserves ambiguity for intrinsically aliased origins. Real TransCODE pilot annotation outputs show that transcript context and candidate feature class are frequently not unique, motivating read-level posterior assignment before downstream scoring. The remaining empirical gate is to quantify how often real Ribo-seq read ambiguity is frame-informative rather than intrinsically or data-limited aliased.

## Current Limitation

This draft is deliberately honest about the boundary. The first matched real probe shows that frame-aware assignment can move posterior mass, but the independent support is sparse and the high-shift cases are concentrated in one small ambiguity cluster. We have not yet shown that real reads are often frame-discordant enough for `frame_em` to matter genome-wide. If production matched samples show the same pattern, the paper should not claim broad isoform deconvolution. It should claim a calibrated substrate generator with an explicit negative result: many transcript/locus assignments are better reported as ambiguous than over-resolved.

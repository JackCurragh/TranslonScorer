# Read Assignment Story Evidence

![Story evidence](../notebooks/read_assignment_story_evidence.png)

## What We Can Say Now

The branch now has four kinds of evidence.

First, deterministic read-origin simulations show that the method ladder behaves as intended: frame-aware EM wins when candidate origins are frame-discordant, abundance EM wins when unique reads anchor origin abundance, frame-aware EM collapses back to abundance EM when frame evidence is uninformative, and aliased origins remain unresolved.

Second, the TransCODE pilot candidate-transcript inference output gives a real annotation-side reason to care about assignment uncertainty. It is not read-level evidence, but it shows that transcript context and feature class are not always uniquely determined from candidate ORF coordinates.

Third, the real candidate-origin probe converts a transcriptome-aligned RiboMetric BAM into read-to-origin candidates. This adds a read-level substrate sanity check: transcript ambiguity is common in a real alignment file, and EM materially reduces assignment entropy even without frame evidence.

Fourth, the first matched frame-support probe joins that candidate table to independent and exploratory frame-support tracks. Independent support is sparse and shifts only a small fraction of weighted reads, while exploratory all-candidate support shifts more but is not independent validation.

Fifth, the scaled 10% local BAM probe shows that frame-aware assignment has a measurable real-data effect at larger depth, but the dominant movement is still annotation/mappability ambiguity rather than clean transcript isoform resolution.

## Local Evidence Obtained

- TransCODE pilot iRibo candidate records: 132,419.
- Transcript context assigned: 126,695 (95.68%).
- Ambiguous inferred feature types: 13,918 (10.51%).
- Unique transcripts represented: 18,568.
- Per-sample ambiguous feature-type rate ranges from 2.62% to 12.62% with mean 7.91%.
- Compatible transcript count is 1 for 56.59% of candidates and greater than 1 for 39.09%.
- Real transcriptome-BAM candidate probe: 35,804 read keys, 88.50% transcript-ambiguous, and 20.78% locus-ambiguous.
- In that probe, transcript-level EM reduces normalized assignment entropy from 0.846 to 0.384; locus-level EM reduces it from 0.263 to 0.084.
- Transcript-level EM weighted classes: 15.38% unique, 27.73% resolved, 24.15% ambiguous, 32.75% data-limited.
- Locus-level EM weighted classes: 73.69% unique, 15.66% resolved, 6.13% ambiguous, 4.52% data-limited.
- Matched independent frame support covers 3,933 frame rows across 1,623 transcripts and 9.8% of weighted candidate reads.
- Independent frame-aware EM shifts transcript assignment modestly: weighted mean TV distance 0.0018 versus frame-blind EM, with 0.65% of weighted reads moving by at least 0.10.
- A minimal frame-support quality gate changes that to weighted mean TV 0.0016, still with 0.65% of weighted reads moving by at least 0.10.
- High independent-support shifts concentrate around `LINC01783`, not genome-wide.
- Scaled 10% transcriptome-BAM probe: 358,510 read keys, 828,006 weighted reads, 88.53% transcript-ambiguous, and 21.25% locus-ambiguous.
- In the scaled probe, independent frame-aware EM shifts transcript assignment by weighted mean TV 0.0257; 3.00% of weighted reads move by at least 0.10.
- With the minimal frame-support gate, the scaled independent shift is weighted mean TV 0.0231; 2.93% of weighted reads move by at least 0.10.
- The scaled high-shift signal is dominated by `ENSG00000228549` / `LINC01783` lncRNA ambiguity. Coding-associated high-shift reads exist: 111 read keys and 1,989 weighted reads, of which 63 read keys and 1,254 weighted reads have gated local frame support.
- The top coding-associated case (`CIC`) spans 46 candidate gene labels and is broad multi-locus ambiguity; `PTCH1` is lower-locus but has no matched local frame-support rows, showing an abundance-propagated EM shift.
- A seed coding review panel now separates direct local-frame-supported candidates from abundance-propagated shifts. Initial `POGK` and `OAZ2` reviews show direct local support but cross-gene ambiguity, so they remain triage examples rather than validation wins.

## Synthetic Assignment Gate

| scenario | best method | soft true posterior | normalized entropy |
|---|---:|---:|---:|
| abundance_resolvable_isoforms | em | 0.911 | 0.321 |
| aliased_no_evidence | fractional | 0.500 | 1.000 |
| frame_discordant_isoforms | frame_em | 1.000 | 0.000 |
| genomic_multilocus_unique_anchors | em | 0.948 | 0.237 |

## What This Identifies

- The transcript isoform question should be framed as an uncertainty-preserving substrate-generation problem, not as a hard transcript caller.
- The genomic multi-locus question is the same mathematical problem with a different abundance target; the prototype already supports `locus_id` and composite `locus_id,tran_id` targets.
- Frame evidence is only a useful disambiguator when candidate origins differ in transcript-coordinate frame and local frame support is strong. Otherwise it must be neutral.
- The real missing evidence is not more model variants. It is candidate-origin tables from real read alignments, with enough metadata to separate transcript ambiguity, locus ambiguity, and intrinsic aliasing.

## What Remains To Obtain

| evidence | why it matters | immediate source | incorporation point |
|---|---|---|---|
| Real read candidate-origin table | Tests how often ambiguity exists at read level | obtained for one transcriptome BAM; repeat for matched production samples | `compare-read-assignment` real-data gate |
| Real frame support joined to candidates | Tests whether frame evidence actually resolves candidate origins | first sparse and scaled local probes complete; repeat on curated coding loci and production matched samples | frame-aware likelihood in assignment EM |
| Locus/transcript equivalence-class counts | Makes genome-scale EM streamable | mapped-index grouped by read key | out-of-core assignment runner |
| Identifiability labels | Prevents overconfident transcript claims | posterior entropy, max posterior, candidate structural equivalence | score schema and RDG sparse export |
| Orthogonal transcript evidence | Separates Ribo-seq ambiguity from true isoform usage | matched RNA-seq/long-read/proteomics where available | priors and validation, not hard truth |
| Curated coding locus review queue | Focuses manual validation on plausible frame-supported biology | `docs/read_assignment_ribometric_10pct.coding_locus_panel_seed.md` | Panel B / Stage 3 candidate triage |

## High-Ambiguity Examples

| feature_id | sample | locus | strand | feature_type | compatible_transcripts | compatible_feature_types |
|---|---|---|---:|---:|---:|---|
| orf_001104442 | Ribo_Fib_pooled | chr3:114389026-114389077 | - | ambiguous | 350 | noncoding_ORF;uORF |
| orf_001124582 | Ribo_Pancreas_pooled | chr3:114389026-114389077 | - | ambiguous | 350 | noncoding_ORF;uORF |
| orf_001156015 | SRR11005880_to_84 | chr3:114389026-114389077 | - | ambiguous | 350 | noncoding_ORF;uORF |
| orf_001174655 | SRR11005885_to_89 | chr3:114389026-114389077 | - | ambiguous | 350 | noncoding_ORF;uORF |
| orf_001195762 | SRR11005890_to_94 | chr3:114389026-114389077 | - | ambiguous | 350 | noncoding_ORF;uORF |
| orf_001237722 | SRR11005900_to_04 | chr3:114389026-114389077 | - | ambiguous | 350 | noncoding_ORF;uORF |
| orf_001104441 | Ribo_Fib_pooled | chr3:114380907-114380940 | - | ambiguous | 342 | noncoding_ORF;uORF |
| orf_001124581 | Ribo_Pancreas_pooled | chr3:114380907-114380940 | - | ambiguous | 342 | noncoding_ORF;uORF |

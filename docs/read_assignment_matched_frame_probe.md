# Matched Frame-Aware Assignment Probe

![Matched frame probe](../notebooks/read_assignment_matched_frame_probe.png)

## Frame Support Tracks

- `none`: neutral frame control.
- `unique_transcript`: primary independent frame support, built only from reads with one compatible transcript target.
- `fractional_all`: exploratory support, built from all candidate origins with fractional read weights. This can be circular and should be treated as a sensitivity check.

| frame_support | rows | transcripts | total_count | mean_support | weighted_reads_with_support |
|---|---:|---:|---:|---:|---:|
| unique_transcript | 3933 | 1623 | 10521.0 | 0.059 | 0.098 |
| fractional_all | 268625 | 71156 | 68417.0 | 0.004 | 0.247 |

## Main Readout

| support | target | method | entropy | high_confidence | final_delta |
|---|---|---:|---:|---:|---:|
| none | tran_id | em | 0.384 | 0.460 | 1.71e-04 |
| none | tran_id | frame_em | 0.384 | 0.460 | 1.71e-04 |
| none | locus_id | em | 0.084 | 0.903 | 1.83e-05 |
| none | locus_id | frame_em | 0.084 | 0.903 | 1.83e-05 |
| unique_transcript | tran_id | em | 0.384 | 0.460 | 1.71e-04 |
| unique_transcript | tran_id | frame_em | 0.387 | 0.453 | 1.72e-04 |
| unique_transcript | locus_id | em | 0.084 | 0.903 | 1.83e-05 |
| unique_transcript | locus_id | frame_em | 0.088 | 0.896 | 1.84e-05 |
| unique_gate10_e02 | tran_id | em | 0.384 | 0.460 | 1.71e-04 |
| unique_gate10_e02 | tran_id | frame_em | 0.386 | 0.453 | 1.71e-04 |
| unique_gate10_e02 | locus_id | em | 0.084 | 0.903 | 1.83e-05 |
| unique_gate10_e02 | locus_id | frame_em | 0.088 | 0.896 | 1.83e-05 |
| fractional_all | tran_id | em | 0.384 | 0.460 | 1.71e-04 |
| fractional_all | tran_id | frame_em | 0.375 | 0.453 | 1.54e-04 |
| fractional_all | locus_id | em | 0.084 | 0.903 | 1.83e-05 |
| fractional_all | locus_id | frame_em | 0.071 | 0.913 | 1.75e-05 |
| fractional_gate10_e02 | tran_id | em | 0.384 | 0.460 | 1.71e-04 |
| fractional_gate10_e02 | tran_id | frame_em | 0.376 | 0.447 | 1.45e-04 |
| fractional_gate10_e02 | locus_id | em | 0.084 | 0.903 | 1.83e-05 |
| fractional_gate10_e02 | locus_id | frame_em | 0.071 | 0.913 | 1.87e-05 |

## Frame-Induced Posterior Shift

| support | target | comparison | weighted_mean_tv | weighted_reads_tv_ge_0.10 |
|---|---|---:|---:|---:|
| none | tran_id | em->frame_em | 0.0000 | 0.0000 |
| none | locus_id | em->frame_em | 0.0000 | 0.0000 |
| unique_transcript | tran_id | em->frame_em | 0.0018 | 0.0065 |
| unique_transcript | locus_id | em->frame_em | 0.0018 | 0.0065 |
| unique_gate10_e02 | tran_id | em->frame_em | 0.0016 | 0.0065 |
| unique_gate10_e02 | locus_id | em->frame_em | 0.0016 | 0.0065 |
| fractional_all | tran_id | em->frame_em | 0.0318 | 0.0617 |
| fractional_all | locus_id | em->frame_em | 0.0103 | 0.0232 |
| fractional_gate10_e02 | tran_id | em->frame_em | 0.0281 | 0.0576 |
| fractional_gate10_e02 | locus_id | em->frame_em | 0.0101 | 0.0232 |

## Interpretation

The independent `unique_transcript` frame support changes transcript-level EM modestly: weighted mean target-posterior TV distance is 0.0018, with 0.65% of weighted reads moving by at least 0.10.
Applying a minimal quality gate (`total_count >= 10` and `support_evidence >= 0.2`) changes that to weighted mean TV 0.0016, with 0.65% of weighted reads moving by at least 0.10.
The exploratory `fractional_all` support produces a larger shift: weighted mean TV distance is 0.0318. Because it reuses ambiguous reads to build frame evidence, this should be treated as an upper-bound/sensitivity result rather than proof.
With the same gate, the exploratory support shifts weighted mean TV 0.0281, with 5.76% of weighted reads moving by at least 0.10.

The next useful step is not another global model variant. It is locus-level inspection of reads with high `em->frame_em` TV distance under the independent support track, checking whether the frame-driven moves are biologically plausible or merely sparse-support artifacts.

Those independent high-shift reads are concentrated rather than genome-wide: the top shifted gene label is `LINC01783`, with 57 exported reads, 443.0 weighted reads, mean TV 0.205, and mean 7.0 candidate transcripts. That pattern argues for targeted locus review before treating the global frame-aware shift as biological.

Targeted review confirms the caution. The LINC01783 cluster has 57 high-shift read keys compatible with 8 transcript targets, but independent frame support exists for only one transcript (`ENST00000415386.2`), with 29 unique-support counts across 3 codons and max support evidence 0.804. The frame term is technically doing what it was designed to do, but this probe does not yet justify a broad claim that frame evidence resolves real isoform origin genome-wide.

Follow-up artifact: `docs/read_assignment_linc01783_case_review.md`.

## Top Independent-Support Transcript Shifts

| read_key | tv_distance | baseline_top | frame_top | largest_shift | frame_top_gene | candidate_transcripts | candidate_loci |
|---|---:|---|---|---|---|---:|---:|
| read959365_x2 | 0.526 | ENST00000415386.2 | ENST00000415386.2 | ENST00000415386.2 -> ENST00000662795.1 | LINC01783 | 7 | 2 |
| read2158018_x3 | 0.526 | ENST00000415386.2 | ENST00000415386.2 | ENST00000415386.2 -> ENST00000660907.1 | LINC01783 | 7 | 2 |
| read2912207_x1 | 0.526 | ENST00000415386.2 | ENST00000415386.2 | ENST00000415386.2 -> ENST00000453554.1 | LINC01783 | 7 | 2 |
| read815391_x3 | 0.526 | ENST00000415386.2 | ENST00000415386.2 | ENST00000415386.2 -> ENST00000662795.1 | LINC01783 | 7 | 2 |
| read671063_x3 | 0.526 | ENST00000415386.2 | ENST00000415386.2 | ENST00000415386.2 -> ENST00000666785.1 | LINC01783 | 7 | 2 |
| read1130542_x1 | 0.198 | ENST00000415386.2 | ENST00000415386.2 | ENST00000415386.2 -> ENST00000660907.1 | LINC01783 | 7 | 2 |
| read1006763_x2 | 0.198 | ENST00000415386.2 | ENST00000415386.2 | ENST00000415386.2 -> ENST00000453554.1 | LINC01783 | 7 | 2 |
| read1285681_x2 | 0.192 | ENST00000415386.2 | ENST00000415386.2 | ENST00000415386.2 -> ENST00000453554.1 | LINC01783 | 7 | 2 |

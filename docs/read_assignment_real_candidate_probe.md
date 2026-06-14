# Real Candidate-Origin Probe

Source BAM: `/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/subsampled_1_percent.bam`

![Real candidate probe](read_assignment_real_candidate_probe.png)

## Candidate-Origin Summary

- Reads: 35,804
- Weighted reads: 68,417
- Mean candidate transcript origins per read: 7.799
- Mean candidate locus origins per read: 1.863
- Transcript-ambiguous reads: 88.50%
- Locus-ambiguous reads: 20.78%
- Weighted transcript-ambiguous fraction: 84.62%
- Weighted locus-ambiguous fraction: 26.31%

## Assignment Summary

| target | method | assigned_fraction | mean_max_posterior | entropy | high_confidence | targets | iterations | final_delta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| tran_id | unique | 0.154 | 0.154 | 0.000 | 0.154 | 71156 | 0 | 0.00e+00 |
| tran_id | fractional | 1.000 | 0.331 | 0.846 | 0.154 | 71156 | 0 | 0.00e+00 |
| tran_id | em | 1.000 | 0.652 | 0.384 | 0.460 | 71156 | 300 | 1.71e-04 |
| tran_id | frame_em | 1.000 | 0.652 | 0.384 | 0.460 | 71156 | 300 | 1.71e-04 |
| locus_id | unique | 0.737 | 0.737 | 0.000 | 0.737 | 14928 | 0 | 0.00e+00 |
| locus_id | fractional | 1.000 | 0.820 | 0.263 | 0.737 | 14928 | 0 | 0.00e+00 |
| locus_id | em | 1.000 | 0.947 | 0.084 | 0.903 | 14928 | 300 | 1.83e-05 |
| locus_id | frame_em | 1.000 | 0.947 | 0.084 | 0.903 | 14928 | 300 | 1.83e-05 |

## Identifiability Summary

See `read_assignment_real_candidate_probe.identifiability_summary.csv` for per-method class fractions. The most important contrast is transcript-level versus locus-level ambiguity: transcript assignment remains much harder than locus assignment in this transcriptome-aligned substrate.

## Interpretation

This is a substrate sanity check, not a final biological validation. The BAM is already transcriptome-aligned, so locus identity is inferred from encoded gene IDs in transcript names. It nevertheless gives a real candidate-origin table with repeated read names across compatible transcripts and genes, which is the immediate input shape required by the assignment model.

Because no frame-support table is joined here, `frame_em` should behave like `em`. Any difference would indicate a bug in neutral frame handling.

# Scaled Matched Frame-Aware Assignment Probe: read_assignment_ribometric_10pct

Source BAM: `/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/subsampled_10_percent.bam`

![Scaled matched frame probe](../notebooks/read_assignment_ribometric_10pct.matched_frame_probe.png)

## Candidate-Origin Summary

- Read keys: 358,510.
- Weighted reads: 828,006.
- Mean candidate transcript targets per read: 7.739.
- Mean candidate locus targets per read: 1.872.
- Transcript-ambiguous read keys: 88.53%.
- Locus-ambiguous read keys: 21.25%.

## Frame Support Tracks

| frame_support | rows | transcripts | total_count | mean_support | weighted_reads_with_support |
|---|---:|---:|---:|---:|---:|
| unique_transcript | 33269 | 5540 | 127010.0 | 0.054 | 0.128 |
| fractional_all | 2182334 | 135587 | 828006.0 | 0.004 | 0.395 |

## Assignment Readout

| support | target | method | entropy | high_confidence | final_delta |
|---|---|---:|---:|---:|---:|
| none | tran_id | em | 0.324 | 0.541 | 1.92e-04 |
| none | tran_id | frame_em | 0.324 | 0.541 | 1.92e-04 |
| none | locus_id | em | 0.144 | 0.829 | 3.82e-05 |
| none | locus_id | frame_em | 0.144 | 0.829 | 3.82e-05 |
| unique_transcript | tran_id | em | 0.324 | 0.541 | 1.92e-04 |
| unique_transcript | tran_id | frame_em | 0.349 | 0.516 | 1.88e-04 |
| unique_transcript | locus_id | em | 0.144 | 0.829 | 3.82e-05 |
| unique_transcript | locus_id | frame_em | 0.151 | 0.828 | 3.49e-05 |
| unique_gate10_e02 | tran_id | em | 0.324 | 0.541 | 1.92e-04 |
| unique_gate10_e02 | tran_id | frame_em | 0.349 | 0.516 | 1.88e-04 |
| unique_gate10_e02 | locus_id | em | 0.144 | 0.829 | 3.82e-05 |
| unique_gate10_e02 | locus_id | frame_em | 0.157 | 0.828 | 3.52e-05 |
| fractional_all | tran_id | em | 0.324 | 0.541 | 1.92e-04 |
| fractional_all | tran_id | frame_em | 0.281 | 0.584 | 2.21e-04 |
| fractional_all | locus_id | em | 0.144 | 0.829 | 3.82e-05 |
| fractional_all | locus_id | frame_em | 0.126 | 0.840 | 3.69e-05 |
| fractional_gate10_e02 | tran_id | em | 0.324 | 0.541 | 1.92e-04 |
| fractional_gate10_e02 | tran_id | frame_em | 0.280 | 0.584 | 2.19e-04 |
| fractional_gate10_e02 | locus_id | em | 0.144 | 0.829 | 3.82e-05 |
| fractional_gate10_e02 | locus_id | frame_em | 0.126 | 0.840 | 3.24e-05 |

## Frame-Induced Posterior Shift

| support | target | comparison | weighted_mean_tv | weighted_reads_tv_ge_0.10 |
|---|---|---:|---:|---:|
| none | tran_id | em->frame_em | 0.0000 | 0.0000 |
| none | locus_id | em->frame_em | 0.0000 | 0.0000 |
| unique_transcript | tran_id | em->frame_em | 0.0257 | 0.0300 |
| unique_transcript | locus_id | em->frame_em | 0.0257 | 0.0299 |
| unique_gate10_e02 | tran_id | em->frame_em | 0.0231 | 0.0293 |
| unique_gate10_e02 | locus_id | em->frame_em | 0.0231 | 0.0293 |
| fractional_all | tran_id | em->frame_em | 0.0478 | 0.1190 |
| fractional_all | locus_id | em->frame_em | 0.0394 | 0.0693 |
| fractional_gate10_e02 | tran_id | em->frame_em | 0.0470 | 0.1174 |
| fractional_gate10_e02 | locus_id | em->frame_em | 0.0396 | 0.0697 |

## Interpretation

Independent frame support shifts transcript-level EM by weighted mean TV 0.0257, with 3.00% of weighted reads moving by at least 0.10.
With the minimal quality gate (`total_count >= 10`, `support_evidence >= 0.2`), the independent shift is weighted mean TV 0.0231, with 2.93% of weighted reads moving by at least 0.10.
The exploratory all-candidate support shifts more: weighted mean TV 0.0478, or 0.0470 after the same gate.
The high independent-support shifts remain concentrated: `ENSG00000228549` accounts for 97.3% of weighted high-shift reads, with 99 exported reads and mean TV 0.899.

## Top Independent-Support High-Shift Loci

| gene_label | reads | weighted_reads | mean_tv | max_tv | mean_transcripts | mean_loci |
|---|---:|---:|---:|---:|---:|---:|
| ENSG00000228549 | 99 | 13588.0 | 0.899 | 0.925 | 7.02 | 2.02 |
| ESRRAP2 | 1 | 374.0 | 0.912 | 0.912 | 2.00 | 2.00 |

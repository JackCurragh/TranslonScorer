# High-Shift Locus Review: PTCH1

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.PTCH1.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 3.
- Weighted high-shift reads: 3.0.
- Mean transcript-level TV distance: 0.686.
- Most frequent compatible transcript: `ENST00000331920.11` with 3 read keys and 3.0 weighted candidate rows.
- Strongest independent frame-support transcript: `` with total unique support 0.0 and max support evidence 0.000.
- Largest individual read shift: `read1754897_x1` TV 0.686, `ENST00000375290.6` to `ENST00000331920.11`.
- Candidate gene labels represented: 1.
- Candidate biotypes represented: `nonsense_mediated_decay;processed_transcript;protein_coding;retained_intron`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated and relatively locus-local, but none of the reviewed candidate transcripts has matched independent local frame support in this probe. The shift is therefore abundance-propagated through the EM rather than directly supported by local frame evidence at the candidate positions.

Operationally, the reviewed reads average 8.3 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000331920.11 | PTCH1 | protein_coding | 3 | 3.0 | 2930-4222 | 0/2/1 |
| ENST00000375290.6 | PTCH1 | nonsense_mediated_decay | 3 | 3.0 | 1794-3086 | 1/0/2 |
| ENST00000429896.6 | PTCH1 | protein_coding | 3 | 3.0 | 2013-3305 | 1/0/2 |
| ENST00000430669.6 | PTCH1 | protein_coding | 3 | 3.0 | 2413-3705 | 2/1/0 |
| ENST00000437951.6 | PTCH1 | protein_coding | 3 | 3.0 | 2173-3465 | 2/1/0 |
| ENST00000690194.1 | PTCH1 | nonsense_mediated_decay | 3 | 3.0 | 1848-3140 | 1/0/2 |
| ENST00000692981.1 | PTCH1 | protein_coding | 3 | 3.0 | 1713-3005 | 1/0/2 |
| ENST00000693534.1 | PTCH1 | processed_transcript | 2 | 2.0 | 516-648 | 2/0/0 |
| ENST00000375271.4 | PTCH1 | protein_coding | 1 | 1.0 | 1020-1020 | 1/0/0 |
| ENST00000549678.1 | PTCH1 | retained_intron | 1 | 1.0 | 214-214 | 0/1/0 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|

# High-Shift Locus Review: RCE1

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.RCE1.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 9.
- Weighted high-shift reads: 9.0.
- Mean transcript-level TV distance: 0.272.
- Most frequent compatible transcript: `ENST00000309657.8` with 9 read keys and 9.0 weighted candidate rows.
- Strongest independent frame-support transcript: `` with total unique support 0.0 and max support evidence 0.000.
- Largest individual read shift: `read3618810_x1` TV 0.272, `ENST00000524849.5` to `ENST00000533277.5`.
- Candidate gene labels represented: 1.
- Candidate biotypes represented: `nonsense_mediated_decay;protein_coding;retained_intron`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated and relatively locus-local, but none of the reviewed candidate transcripts has matched independent local frame support in this probe. The shift is therefore abundance-propagated through the EM rather than directly supported by local frame evidence at the candidate positions.

Operationally, the reviewed reads average 5.6 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000309657.8 | RCE1 | protein_coding | 9 | 9.0 | 95-950 | 1/4/4 |
| ENST00000524849.5 | RCE1 | nonsense_mediated_decay | 9 | 9.0 | 85-936 | 2/3/4 |
| ENST00000533277.5 | RCE1 | retained_intron | 9 | 9.0 | 701-2233 | 4/3/2 |
| ENST00000524506.5 | RCE1 | protein_coding | 8 | 8.0 | 99-891 | 4/1/3 |
| ENST00000525356.1 | RCE1 | protein_coding | 8 | 8.0 | 51-782 | 2/4/2 |
| ENST00000532775.5 | RCE1 | retained_intron | 3 | 3.0 | 126-399 | 2/1/0 |
| ENST00000530610.1 | RCE1 | retained_intron | 2 | 2.0 | 26-146 | 0/0/2 |
| ENST00000534822.1 | RCE1 | retained_intron | 2 | 2.0 | 1-121 | 0/2/0 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|

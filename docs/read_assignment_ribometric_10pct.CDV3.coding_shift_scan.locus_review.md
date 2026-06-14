# High-Shift Locus Review: CDV3

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.CDV3.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 2.
- Weighted high-shift reads: 250.0.
- Mean transcript-level TV distance: 0.391.
- Most frequent compatible transcript: `ENST00000264993.8` with 2 read keys and 250.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000216268.6` with total unique support 480.0 and max support evidence 0.949.
- Largest individual read shift: `read11008_x249` TV 0.674, `ENST00000216268.6` to `ENST00000264993.8`.
- Candidate gene labels represented: 15.
- Candidate biotypes represented: `lncRNA;processed_transcript;protein_coding;retained_intron`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 13.5 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000264993.8 | CDV3 | protein_coding | 2 | 250.0 | 495-495 | 2/0/0 |
| ENST00000431519.6 | CDV3 | protein_coding | 2 | 250.0 | 205-205 | 0/2/0 |
| ENST00000688838.1 | CDV3 | protein_coding | 2 | 250.0 | 495-495 | 2/0/0 |
| ENST00000216268.6 | ZBED4 | protein_coding | 1 | 249.0 | 20-20 | 0/0/1 |
| ENST00000315073.10 | TRIM41 | protein_coding | 1 | 249.0 | 2191-2191 | 0/1/0 |
| ENST00000322945.11 | MAZ | protein_coding | 1 | 249.0 | 26-26 | 0/0/1 |
| ENST00000408973.3 | LCNL1 | protein_coding | 1 | 249.0 | 900-900 | 1/0/0 |
| ENST00000482657.1 | LCNL1 | retained_intron | 1 | 249.0 | 2770-2770 | 0/1/0 |
| ENST00000508930.1 | TRIM41 | retained_intron | 1 | 249.0 | 2382-2382 | 1/0/0 |
| ENST00000510072.1 | TRIM41 | processed_transcript | 1 | 249.0 | 288-288 | 1/0/0 |
| ENST00000515223.1 | TRIM41 | retained_intron | 1 | 249.0 | 4258-4258 | 0/1/0 |
| ENST00000530595.1 | ENSG00000254602 | lncRNA | 1 | 249.0 | 111-111 | 1/0/0 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000216268.6 | 7 | 480.0 | 0.949 | 0.175 | 0.143 |
| ENST00000649141.2 | 17 | 47.0 | 0.428 | 0.052 | 0.333 |
| ENST00000262367.10 | 4 | 5.0 | 0.082 | 0.040 | 0.143 |
| ENST00000530595.1 | 1 | 3.0 | 0.145 | 0.145 | 0.067 |
| ENST00000322945.11 | 1 | 1.0 | 0.026 | 0.026 | 0.143 |
| ENST00000431519.6 | 1 | 1.0 | 0.026 | 0.026 | 0.143 |

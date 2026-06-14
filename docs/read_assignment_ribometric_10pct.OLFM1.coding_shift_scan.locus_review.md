# High-Shift Locus Review: OLFM1

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.OLFM1.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 4.
- Weighted high-shift reads: 144.0.
- Mean transcript-level TV distance: 0.625.
- Most frequent compatible transcript: `ENST00000371793.8` with 4 read keys and 144.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000400454.6` with total unique support 2786.0 and max support evidence 0.968.
- Largest individual read shift: `read1280659_x1` TV 0.648, `ENST00000400454.6` to `ENST00000392991.8`.
- Candidate gene labels represented: 10.
- Candidate biotypes represented: `nonsense_mediated_decay;processed_pseudogene;processed_transcript;protein_coding`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 6.8 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000371793.8 | OLFM1 | protein_coding | 4 | 144.0 | 158-167 | 0/1/3 |
| ENST00000392991.8 | OLFM1 | protein_coding | 4 | 144.0 | 152-161 | 0/1/3 |
| ENST00000400454.6 | DSCAM | protein_coding | 3 | 143.0 | 364-365 | 0/1/2 |
| ENST00000331628.8 | SYNDIG1L | protein_coding | 2 | 39.0 | 104-105 | 1/0/1 |
| ENST00000251809.4 | SPAG1 | protein_coding | 1 | 1.0 | 1356-1356 | 1/0/0 |
| ENST00000326043.5 | MAF | protein_coding | 1 | 1.0 | 1501-1501 | 0/1/0 |
| ENST00000388798.7 | SPAG1 | protein_coding | 1 | 1.0 | 1302-1302 | 1/0/0 |
| ENST00000393350.1 | MAF | protein_coding | 1 | 1.0 | 1478-1478 | 0/0/1 |
| ENST00000413684.1 | KMT5AP2 | processed_pseudogene | 1 | 1.0 | 29-29 | 0/0/1 |
| ENST00000420124.4 | KMT2B | protein_coding | 1 | 1.0 | 17-17 | 0/0/1 |
| ENST00000431232.7 | PGAP6 | protein_coding | 1 | 1.0 | 93-93 | 1/0/0 |
| ENST00000523302.1 | SPAG1 | processed_transcript | 1 | 1.0 | 177-177 | 1/0/0 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000400454.6 | 3 | 2786.0 | 0.968 | 0.478 | 0.440 |
| ENST00000682841.1 | 4 | 602.0 | 0.988 | 0.532 | 0.143 |
| ENST00000431232.7 | 2 | 2.0 | 0.026 | 0.026 | 0.143 |
| ENST00000331628.8 | 1 | 1.0 | 0.026 | 0.026 | 0.143 |

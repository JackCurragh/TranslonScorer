# High-Shift Locus Review: DSCAM

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.DSCAM.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 4.
- Weighted high-shift reads: 74.0.
- Mean transcript-level TV distance: 0.185.
- Most frequent compatible transcript: `ENST00000400454.6` with 4 read keys and 74.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000400454.6` with total unique support 2786.0 and max support evidence 0.968.
- Largest individual read shift: `read74731_x43` TV 0.229, `ENST00000400454.6` to `ENST00000262367.10`.
- Candidate gene labels represented: 4.
- Candidate biotypes represented: `processed_transcript;protein_coding`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 2.8 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000400454.6 | DSCAM | protein_coding | 4 | 74.0 | 364-365 | 0/3/1 |
| ENST00000262367.10 | CREBBP | protein_coding | 2 | 69.0 | 295-296 | 0/1/1 |
| ENST00000399451.6 | ANKRD28 | protein_coding | 2 | 5.0 | 176-176 | 0/0/2 |
| ENST00000251809.4 | SPAG1 | protein_coding | 1 | 43.0 | 1356-1356 | 1/0/0 |
| ENST00000388798.7 | SPAG1 | protein_coding | 1 | 43.0 | 1302-1302 | 1/0/0 |
| ENST00000523302.1 | SPAG1 | processed_transcript | 1 | 43.0 | 177-177 | 1/0/0 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000400454.6 | 3 | 2786.0 | 0.968 | 0.478 | 0.440 |
| ENST00000399451.6 | 6 | 10.0 | 0.145 | 0.051 | 0.333 |
| ENST00000262367.10 | 4 | 5.0 | 0.082 | 0.040 | 0.143 |

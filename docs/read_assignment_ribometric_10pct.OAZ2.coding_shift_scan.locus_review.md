# High-Shift Locus Review: OAZ2

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.OAZ2.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 1.
- Weighted high-shift reads: 1.0.
- Mean transcript-level TV distance: 0.692.
- Most frequent compatible transcript: `ENST00000326005.10` with 1 read keys and 1.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000682841.1` with total unique support 602.0 and max support evidence 0.988.
- Largest individual read shift: `read424386_x1` TV 0.692, `ENST00000682841.1` to `ENST00000326005.10`.
- Candidate gene labels represented: 2.
- Candidate biotypes represented: `protein_coding`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 2.0 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000326005.10 | OAZ2 | protein_coding | 1 | 1.0 | 18-18 | 1/0/0 |
| ENST00000682841.1 | IL17D | protein_coding | 1 | 1.0 | 35-35 | 0/0/1 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000682841.1 | 4 | 602.0 | 0.988 | 0.532 | 0.143 |
| ENST00000326005.10 | 2 | 2.0 | 0.026 | 0.026 | 0.143 |

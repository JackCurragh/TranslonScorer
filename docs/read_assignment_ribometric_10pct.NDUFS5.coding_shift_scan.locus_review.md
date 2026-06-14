# High-Shift Locus Review: NDUFS5

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.NDUFS5.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 1.
- Weighted high-shift reads: 2.0.
- Mean transcript-level TV distance: 0.208.
- Most frequent compatible transcript: `ENST00000372967.3` with 1 read keys and 2.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000372967.3` with total unique support 40.0 and max support evidence 0.576.
- Largest individual read shift: `read2135854_x2` TV 0.208, `ENST00000372967.3` to `ENST00000372969.8`.
- Candidate gene labels represented: 1.
- Candidate biotypes represented: `protein_coding`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 2.0 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000372967.3 | NDUFS5 | protein_coding | 1 | 2.0 | 32-32 | 0/0/1 |
| ENST00000372969.8 | NDUFS5 | protein_coding | 1 | 2.0 | 0-0 | 1/0/0 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000372967.3 | 6 | 40.0 | 0.576 | 0.222 | 0.474 |
| ENST00000372969.8 | 4 | 33.0 | 0.816 | 0.265 | 0.263 |

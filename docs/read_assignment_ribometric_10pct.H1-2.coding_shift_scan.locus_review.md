# High-Shift Locus Review: H1-2

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.H1-2.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 2.
- Weighted high-shift reads: 8.0.
- Mean transcript-level TV distance: 0.347.
- Most frequent compatible transcript: `ENST00000304218.6` with 2 read keys and 8.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000343677.4` with total unique support 336.0 and max support evidence 0.961.
- Largest individual read shift: `read1367346_x2` TV 0.349, `ENST00000304218.6` to `ENST00000343677.4`.
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
| ENST00000304218.6 | H1-4 | protein_coding | 2 | 8.0 | 63-64 | 1/1/0 |
| ENST00000343677.4 | H1-2 | protein_coding | 2 | 8.0 | 43-44 | 0/1/1 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000343677.4 | 88 | 336.0 | 0.961 | 0.118 | 0.474 |
| ENST00000304218.6 | 84 | 294.0 | 0.868 | 0.095 | 0.481 |

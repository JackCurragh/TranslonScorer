# High-Shift Locus Review: POGK

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.POGK.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 1.
- Weighted high-shift reads: 18.0.
- Mean transcript-level TV distance: 0.365.
- Most frequent compatible transcript: `ENST00000253928.14` with 1 read keys and 18.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000449930.5` with total unique support 611.0 and max support evidence 0.992.
- Largest individual read shift: `read512906_x18` TV 0.365, `ENST00000449930.5` to `ENST00000253928.14`.
- Candidate gene labels represented: 2.
- Candidate biotypes represented: `processed_transcript;protein_coding`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 6.0 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000253928.14 | FLYWCH1 | protein_coding | 1 | 18.0 | 103-103 | 0/1/0 |
| ENST00000416288.6 | FLYWCH1 | protein_coding | 1 | 18.0 | 35-35 | 0/0/1 |
| ENST00000449930.5 | POGK | protein_coding | 1 | 18.0 | 166-166 | 0/1/0 |
| ENST00000570425.5 | FLYWCH1 | protein_coding | 1 | 18.0 | 27-27 | 1/0/0 |
| ENST00000571140.5 | FLYWCH1 | processed_transcript | 1 | 18.0 | 53-53 | 0/0/1 |
| ENST00000573525.1 | FLYWCH1 | protein_coding | 1 | 18.0 | 10-10 | 0/1/0 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000449930.5 | 4 | 611.0 | 0.992 | 0.502 | 0.333 |

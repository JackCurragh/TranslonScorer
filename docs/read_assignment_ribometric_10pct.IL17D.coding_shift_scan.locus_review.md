# High-Shift Locus Review: IL17D

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.IL17D.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 4.
- Weighted high-shift reads: 51.0.
- Mean transcript-level TV distance: 0.597.
- Most frequent compatible transcript: `ENST00000253928.14` with 4 read keys and 51.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000449930.5` with total unique support 611.0 and max support evidence 0.992.
- Largest individual read shift: `read3354147_x1` TV 0.706, `ENST00000682841.1` to `ENST00000253928.14`.
- Candidate gene labels represented: 6.
- Candidate biotypes represented: `nonsense_mediated_decay;processed_transcript;protein_coding`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 8.5 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000253928.14 | FLYWCH1 | protein_coding | 4 | 51.0 | 104-105 | 3/0/1 |
| ENST00000416288.6 | FLYWCH1 | protein_coding | 4 | 51.0 | 36-37 | 1/3/0 |
| ENST00000570425.5 | FLYWCH1 | protein_coding | 4 | 51.0 | 28-29 | 0/1/3 |
| ENST00000571140.5 | FLYWCH1 | processed_transcript | 4 | 51.0 | 54-55 | 1/3/0 |
| ENST00000573525.1 | FLYWCH1 | protein_coding | 4 | 51.0 | 11-12 | 3/0/1 |
| ENST00000682841.1 | IL17D | protein_coding | 4 | 51.0 | 32-33 | 3/0/1 |
| ENST00000317991.10 | GRAMD1A | protein_coding | 2 | 10.0 | 20-20 | 0/0/2 |
| ENST00000600231.5 | GRAMD1A | nonsense_mediated_decay | 2 | 10.0 | 20-20 | 0/0/2 |
| ENST00000680623.1 | GRAMD1A | protein_coding | 2 | 10.0 | 20-20 | 0/0/2 |
| ENST00000281419.8 | ASAP2 | protein_coding | 1 | 40.0 | 178-178 | 0/1/0 |
| ENST00000315273.4 | ASAP2 | protein_coding | 1 | 40.0 | 225-225 | 1/0/0 |
| ENST00000449930.5 | POGK | protein_coding | 1 | 40.0 | 167-167 | 0/0/1 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000449930.5 | 4 | 611.0 | 0.992 | 0.502 | 0.333 |
| ENST00000682841.1 | 4 | 602.0 | 0.988 | 0.532 | 0.143 |
| ENST00000354300.5 | 8 | 9.0 | 0.082 | 0.033 | 0.143 |

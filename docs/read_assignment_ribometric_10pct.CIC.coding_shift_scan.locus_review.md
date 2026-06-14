# High-Shift Locus Review: CIC

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.CIC.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 4.
- Weighted high-shift reads: 671.0.
- Mean transcript-level TV distance: 0.147.
- Most frequent compatible transcript: `ENST00000160740.7` with 4 read keys and 671.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000160740.7` with total unique support 1280.0 and max support evidence 0.997.
- Largest individual read shift: `read59144_x255` TV 0.166, `ENST00000160740.7` to `ENST00000418437.1`.
- Candidate gene labels represented: 46.
- Candidate biotypes represented: `lncRNA;nonsense_mediated_decay;processed_transcript;protein_coding;retained_intron;transcribed_processed_pseudogene`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 28.5 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000160740.7 | CIC | protein_coding | 4 | 671.0 | 416-418 | 1/2/1 |
| ENST00000534404.1 | LRP4 | protein_coding | 4 | 671.0 | 31-33 | 2/1/1 |
| ENST00000258098.6 | RAB11FIP5 | protein_coding | 3 | 670.0 | 24-25 | 1/2/0 |
| ENST00000286548.9 | GNAQ | protein_coding | 3 | 670.0 | 218-219 | 2/0/1 |
| ENST00000314888.10 | TLN1 | protein_coding | 3 | 670.0 | 83-84 | 2/0/1 |
| ENST00000418437.1 | ESRRAP2 | transcribed_processed_pseudogene | 3 | 670.0 | 130-131 | 0/1/2 |
| ENST00000573070.5 | PLSCR3 | nonsense_mediated_decay | 3 | 670.0 | 259-260 | 0/1/2 |
| ENST00000573331.5 | TMEM256-PLSCR3 | nonsense_mediated_decay | 3 | 670.0 | 545-546 | 2/0/1 |
| ENST00000574401.5 | PLSCR3 | protein_coding | 3 | 670.0 | 224-225 | 2/0/1 |
| ENST00000215730.12 | SNAP29 | protein_coding | 2 | 415.0 | 17-17 | 0/0/2 |
| ENST00000262160.11 | SMAD2 | protein_coding | 2 | 415.0 | 127-127 | 0/2/0 |
| ENST00000263257.6 | NOVA2 | protein_coding | 2 | 415.0 | 1617-1617 | 2/0/0 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000160740.7 | 1 | 1280.0 | 0.997 | 0.997 | 0.000 |
| ENST00000314888.10 | 164 | 199.0 | 0.145 | 0.035 | 0.455 |
| ENST00000371711.4 | 25 | 28.0 | 0.082 | 0.028 | 0.455 |
| ENST00000505490.3 | 12 | 17.0 | 0.207 | 0.051 | 0.143 |
| ENST00000367607.8 | 6 | 6.0 | 0.026 | 0.026 | 0.143 |
| ENST00000338343.10 | 5 | 5.0 | 0.026 | 0.026 | 0.143 |
| ENST00000266077.5 | 4 | 4.0 | 0.026 | 0.026 | 0.143 |
| ENST00000286548.9 | 4 | 4.0 | 0.026 | 0.026 | 0.143 |
| ENST00000215730.12 | 3 | 3.0 | 0.026 | 0.026 | 0.143 |
| ENST00000418437.1 | 1 | 2.0 | 0.082 | 0.082 | 0.091 |
| ENST00000546324.1 | 2 | 2.0 | 0.026 | 0.026 | 0.143 |
| ENST00000470383.1 | 1 | 1.0 | 0.026 | 0.026 | 0.143 |

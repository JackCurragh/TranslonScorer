# High-Shift Locus Review: TMEM271

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.TMEM271.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 23.
- Weighted high-shift reads: 235.0.
- Mean transcript-level TV distance: 0.268.
- Most frequent compatible transcript: `ENST00000610212.3` with 23 read keys and 283.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000373210.4` with total unique support 3227.0 and max support evidence 0.801.
- Largest individual read shift: `read1736764_x1` TV 0.551, `ENST00000400454.6` to `ENST00000610212.3`.
- Candidate gene labels represented: 195.
- Candidate biotypes represented: `lncRNA;nonsense_mediated_decay;processed_pseudogene;processed_transcript;protein_coding;retained_intron;transcribed_unprocessed_pseudogene`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 63.7 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000610212.3 | TMEM271 | protein_coding | 23 | 283.0 | 27-178 | 6/19/0 |
| ENST00000262752.5 | RPS6KA6 | protein_coding | 20 | 186.0 | 185-198 | 16/0/4 |
| ENST00000399451.6 | ANKRD28 | protein_coding | 18 | 228.0 | 165-177 | 17/0/1 |
| ENST00000286201.3 | FZD7 | protein_coding | 17 | 198.0 | 24-35 | 16/0/1 |
| ENST00000361132.9 | RASGEF1C | protein_coding | 17 | 198.0 | 75-86 | 16/0/1 |
| ENST00000396432.7 | ARHGAP21 | protein_coding | 17 | 198.0 | 188-199 | 0/1/16 |
| ENST00000525908.6 | C11orf80 | protein_coding | 17 | 198.0 | 84-92 | 16/0/1 |
| ENST00000527368.5 | C11orf80 | processed_transcript | 17 | 198.0 | 46-54 | 1/16/0 |
| ENST00000527634.5 | C11orf80 | protein_coding | 17 | 198.0 | 49-57 | 1/16/0 |
| ENST00000531400.6 | C11orf80 | nonsense_mediated_decay | 17 | 198.0 | 0-8 | 16/0/1 |
| ENST00000540737.7 | C11orf80 | protein_coding | 17 | 198.0 | 51-59 | 16/0/1 |
| ENST00000615330.4 | RASGEF1C | protein_coding | 17 | 198.0 | 9-20 | 16/0/1 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000373210.4 | 4 | 3227.0 | 0.801 | 0.294 | 0.143 |
| ENST00000400454.6 | 3 | 2786.0 | 0.968 | 0.478 | 0.440 |
| ENST00000682841.1 | 4 | 602.0 | 0.988 | 0.532 | 0.143 |
| ENST00000341744.8 | 38 | 266.0 | 0.881 | 0.170 | 0.455 |
| ENST00000396432.7 | 5 | 91.0 | 0.889 | 0.279 | 0.143 |
| ENST00000262644.9 | 32 | 58.0 | 0.238 | 0.050 | 0.455 |
| ENST00000395699.5 | 37 | 53.0 | 0.145 | 0.041 | 0.474 |
| ENST00000244745.4 | 44 | 51.0 | 0.082 | 0.033 | 0.455 |
| ENST00000332707.10 | 25 | 51.0 | 0.542 | 0.077 | 0.333 |
| ENST00000221419.10 | 23 | 41.0 | 0.421 | 0.057 | 0.455 |
| ENST00000374580.10 | 23 | 28.0 | 0.082 | 0.032 | 0.455 |
| ENST00000525234.1 | 4 | 27.0 | 0.481 | 0.227 | 0.455 |

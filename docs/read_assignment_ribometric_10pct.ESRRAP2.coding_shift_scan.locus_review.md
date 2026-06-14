# High-Shift Locus Review: ESRRAP2

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.ESRRAP2.coding_shift_scan.png)

## Summary

- High-shift read keys reviewed: 15.
- Weighted high-shift reads: 388.0.
- Mean transcript-level TV distance: 0.493.
- Most frequent compatible transcript: `ENST00000418437.1` with 15 read keys and 388.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000160740.7` with total unique support 1280.0 and max support evidence 0.997.
- Largest individual read shift: `read36209_x374` TV 0.913, `ENST00000160740.7` to `ENST00000418437.1`.
- Candidate gene labels represented: 8.
- Candidate biotypes represented: `nonsense_mediated_decay;processed_pseudogene;protein_coding;retained_intron;transcribed_processed_pseudogene`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are coding-associated, but the candidate set is broad and multi-locus rather than a clean within-gene isoform choice. The frame term changes posterior mass for these reads, yet this is still better treated as ambiguity triage than as transcript-origin validation until genomic mappability and curated locus context are checked.

Operationally, the reviewed reads average 6.7 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000418437.1 | ESRRAP2 | transcribed_processed_pseudogene | 15 | 388.0 | 128-1493 | 8/4/3 |
| ENST00000000442.11 | ESRRA | protein_coding | 13 | 13.0 | 240-1438 | 4/2/7 |
| ENST00000405666.5 | ESRRA | protein_coding | 13 | 13.0 | 249-1447 | 4/2/7 |
| ENST00000406310.6 | ESRRA | protein_coding | 13 | 13.0 | 267-1462 | 4/2/7 |
| ENST00000677967.1 | ESRRA | protein_coding | 13 | 13.0 | 240-1435 | 4/2/7 |
| ENST00000400596.2 | ESRRAP1 | processed_pseudogene | 11 | 11.0 | 45-1198 | 5/4/2 |
| ENST00000468670.2 | ESRRA | protein_coding | 5 | 5.0 | 428-652 | 0/4/1 |
| ENST00000539594.5 | ESRRA | protein_coding | 5 | 5.0 | 230-573 | 2/1/2 |
| ENST00000545035.1 | ESRRA | protein_coding | 4 | 4.0 | 83-555 | 1/1/2 |
| ENST00000467987.1 | ESRRA | retained_intron | 3 | 3.0 | 253-443 | 0/2/1 |
| ENST00000160740.7 | CIC | protein_coding | 1 | 374.0 | 415-415 | 0/1/0 |
| ENST00000536050.5 | WDR59 | protein_coding | 1 | 1.0 | 80-80 | 0/0/1 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000160740.7 | 1 | 1280.0 | 0.997 | 0.997 | 0.000 |
| ENST00000418437.1 | 1 | 2.0 | 0.082 | 0.082 | 0.091 |

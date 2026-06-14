# High-Shift Locus Review: ENSG00000228549

Scaled probe label: `read_assignment_ribometric_10pct`

![High-shift locus review](../notebooks/read_assignment_ribometric_10pct.ENSG00000228549.high_shift_examples.png)

## Summary

- High-shift read keys reviewed: 99.
- Weighted high-shift reads: 13588.0.
- Mean transcript-level TV distance: 0.899.
- Most frequent compatible transcript: `ENST00000415386.2` with 99 read keys and 13588.0 weighted candidate rows.
- Strongest independent frame-support transcript: `ENST00000415386.2` with total unique support 327.0 and max support evidence 0.617.
- Largest individual read shift: `read72852_x141` TV 0.925, `ENST00000415386.2` to `ENST00000438002.1`.
- Candidate gene labels represented: 3.
- Candidate biotypes represented: `lncRNA`.

## Interpretation

This locus review asks whether the scaled frame-aware shift is a broad deconvolution signal or a concentrated ambiguity case.

The high-shift reads are compatible with a tight noncoding or pseudogene-like transcript block. The frame term redistributes mass among neighbouring compatible transcript models, but the candidate set is biologically weak as a validation target because the labels are not a clean translated isoform case.

Operationally, the reviewed reads average 7.0 compatible transcript rows per read in this case subset.

For the project gate, this supports keeping frame-aware assignment as an uncertainty-preserving model, while requiring production matched samples and curated coding/dual-coding loci before claiming real isoform-origin resolution.

## Candidate Transcript Summary

| tran_id | gene_name | biotype | read_keys | weighted_candidate_rows | pos_range | frame_rows_0/1/2 |
|---|---|---|---:|---:|---|---|
| ENST00000415386.2 | LINC01783 | lncRNA | 99 | 13588.0 | 336-344 | 50/13/36 |
| ENST00000438002.1 | ENSG00000228549 | lncRNA | 99 | 13588.0 | 1675-1683 | 36/50/13 |
| ENST00000453554.1 | ENSG00000228549 | lncRNA | 99 | 13588.0 | 336-344 | 50/13/36 |
| ENST00000660907.1 | ENSG00000228549 | lncRNA | 99 | 13588.0 | 339-347 | 50/13/36 |
| ENST00000662795.1 | ENSG00000228549 | lncRNA | 99 | 13588.0 | 394-402 | 36/50/13 |
| ENST00000666785.1 | ENSG00000228549 | lncRNA | 99 | 13588.0 | 61-69 | 36/50/13 |
| ENST00000668460.1 | ENSG00000228549 | lncRNA | 99 | 13588.0 | 394-402 | 36/50/13 |
| ENST00000599640.6 | ENSG00000227733 | lncRNA | 2 | 4.0 | 1965-1970 | 1/0/1 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence | max_secondary_mass |
|---|---:|---:|---:|---:|---:|
| ENST00000415386.2 | 5 | 327.0 | 0.617 | 0.209 | 0.326 |

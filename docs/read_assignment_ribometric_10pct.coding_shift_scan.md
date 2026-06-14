# Coding High-Shift Scan: read_assignment_ribometric_10pct

## Summary

- Total read keys scanned: 358,510.
- High-shift read keys at TV >= 0.10: 669.
- Coding-associated high-shift read keys: 111.
- Coding-associated weighted reads: 1989.0.
- Coding-associated read keys with gated local frame support: 63.
- Coding-associated weighted reads with gated local frame support: 1254.0.
- Top coding-associated frame-top gene label: `CIC` with 671.0 weighted reads and mean TV 0.147.

## Interpretation

This scan excludes the dominant lncRNA-only ambiguity story by asking whether any high frame-induced shifts involve protein-coding candidates. These are not truth-labelled validations, but they are the next candidates for curated inspection because they are closer to the biological claim of translated isoform disambiguation.

## Gene Summary

| frame_top_gene | frame_top_biotype | reads | weighted_reads | local_supported_weight | mean_tv | max_tv | mean_transcripts | mean_loci |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CIC | protein_coding | 4 | 671.0 | 1.0 | 0.147 | 0.166 | 28.50 | 18.50 |
| ESRRAP2 | transcribed_processed_pseudogene | 15 | 388.0 | 374.0 | 0.493 | 0.913 | 6.73 | 2.93 |
| CDV3 | protein_coding | 2 | 250.0 | 249.0 | 0.391 | 0.674 | 13.50 | 8.00 |
| TMEM271 | protein_coding | 23 | 235.0 | 235.0 | 0.268 | 0.551 | 63.74 | 28.57 |
| OLFM1 | protein_coding | 4 | 144.0 | 144.0 | 0.625 | 0.648 | 6.75 | 4.00 |
| DSCAM | protein_coding | 4 | 74.0 | 74.0 | 0.185 | 0.229 | 2.75 | 2.25 |
| IL17D | protein_coding | 4 | 51.0 | 51.0 | 0.597 | 0.706 | 8.50 | 3.25 |
| PGP | protein_coding | 1 | 34.0 | 34.0 | 0.145 | 0.145 | 25.00 | 11.00 |
| TMEM256-PLSCR3 | nonsense_mediated_decay | 2 | 28.0 | 26.0 | 0.424 | 0.739 | 9.00 | 5.00 |
| POGK | protein_coding | 1 | 18.0 | 18.0 | 0.365 | 0.365 | 6.00 | 2.00 |
| TRIM28 | protein_coding | 1 | 16.0 | 0.0 | 0.147 | 0.147 | 31.00 | 16.00 |
| AGO4 | protein_coding | 4 | 9.0 | 9.0 | 0.248 | 0.304 | 65.25 | 24.25 |
| RCE1 | retained_intron | 9 | 9.0 | 0.0 | 0.272 | 0.272 | 5.56 | 1.00 |
| MT-CO1 | protein_coding | 2 | 9.0 | 9.0 | 0.140 | 0.140 | 2.00 | 2.00 |
| H1-2 | protein_coding | 2 | 8.0 | 8.0 | 0.347 | 0.349 | 2.00 | 2.00 |
| LINC01783 | lncRNA | 3 | 7.0 | 0.0 | 0.773 | 0.776 | 42.33 | 10.33 |
| FLYWCH2 | protein_coding | 7 | 7.0 | 0.0 | 0.113 | 0.113 | 8.00 | 2.00 |
| SIX1 | protein_coding | 1 | 5.0 | 5.0 | 0.309 | 0.309 | 5.00 | 3.00 |
| TLN1 | protein_coding | 3 | 3.0 | 1.0 | 0.534 | 0.740 | 7.33 | 6.00 |
| PTCH1 | protein_coding | 3 | 3.0 | 0.0 | 0.686 | 0.686 | 8.33 | 1.00 |

## Top Reads

| read_key | weight | tv | local_support | baseline_top | frame_top | frame_top_gene | frame_top_biotype | candidate_biotypes |
|---|---:|---:|---:|---|---|---|---|---|
| read36209_x374 | 374.0 | 0.913 | 0.997 | ENST00000160740.7 | ENST00000418437.1 | ESRRAP2 | transcribed_processed_pseudogene | transcribed_processed_pseudogene;protein_coding |
| read3669034_x1 | 1.0 | 0.776 | 0.000 | ENST00000415386.2 | ENST00000415386.2 | LINC01783 | lncRNA | lncRNA;retained_intron;protein_coding |
| read484734_x4 | 4.0 | 0.774 | 0.000 | ENST00000415386.2 | ENST00000415386.2 | LINC01783 | lncRNA | protein_coding;retained_intron;nonsense_mediated_decay;processed_transcript;lncRNA;processed_pseudogene |
| read758335_x2 | 2.0 | 0.768 | 0.000 | ENST00000415386.2 | ENST00000415386.2 | LINC01783 | lncRNA | lncRNA;protein_coding;retained_intron;processed_transcript;processed_pseudogene |
| read57823_x1 | 1.0 | 0.748 | 0.997 | ENST00000160740.7 | ENST00000376759.8 | RBM3 | protein_coding | protein_coding;unprocessed_pseudogene;nonsense_mediated_decay;processed_transcript;retained_intron;lncRNA;TEC |
| read255017_x1 | 1.0 | 0.740 | 0.997 | ENST00000160740.7 | ENST00000314888.10 | TLN1 | protein_coding | lncRNA;protein_coding;retained_intron |
| read357901_x26 | 26.0 | 0.739 | 0.997 | ENST00000160740.7 | ENST00000573331.5 | TMEM256-PLSCR3 | nonsense_mediated_decay | nonsense_mediated_decay;protein_coding |
| read3354147_x1 | 1.0 | 0.706 | 0.988 | ENST00000682841.1 | ENST00000682841.1 | IL17D | protein_coding | protein_coding;processed_transcript;nonsense_mediated_decay |
| read424386_x1 | 1.0 | 0.692 | 0.988 | ENST00000682841.1 | ENST00000326005.10 | OAZ2 | protein_coding | protein_coding |
| read823778_x9 | 9.0 | 0.689 | 0.988 | ENST00000682841.1 | ENST00000682841.1 | IL17D | protein_coding | protein_coding;processed_transcript;nonsense_mediated_decay |
| read1754897_x1 | 1.0 | 0.686 | 0.000 | ENST00000375290.6 | ENST00000331920.11 | PTCH1 | protein_coding | protein_coding;nonsense_mediated_decay;processed_transcript |
| read3653463_x1 | 1.0 | 0.686 | 0.000 | ENST00000375290.6 | ENST00000331920.11 | PTCH1 | protein_coding | protein_coding;nonsense_mediated_decay;processed_transcript |
| read2754214_x1 | 1.0 | 0.686 | 0.000 | ENST00000375290.6 | ENST00000331920.11 | PTCH1 | protein_coding | protein_coding;nonsense_mediated_decay;retained_intron |
| read11008_x249 | 249.0 | 0.674 | 0.949 | ENST00000216268.6 | ENST00000264993.8 | CDV3 | protein_coding | protein_coding;retained_intron;processed_transcript;lncRNA |
| read1280659_x1 | 1.0 | 0.648 | 0.968 | ENST00000400454.6 | ENST00000392991.8 | OLFM1 | protein_coding | protein_coding;processed_transcript;nonsense_mediated_decay |
| read115159_x38 | 38.0 | 0.645 | 0.968 | ENST00000400454.6 | ENST00000392991.8 | OLFM1 | protein_coding | protein_coding |
| read156842_x104 | 104.0 | 0.645 | 0.968 | ENST00000400454.6 | ENST00000371793.8 | OLFM1 | protein_coding | protein_coding |
| read1102481_x1 | 1.0 | 0.644 | 0.000 | ENST00000573331.5 | ENST00000418437.1 | ESRRAP2 | transcribed_processed_pseudogene | transcribed_processed_pseudogene;protein_coding;nonsense_mediated_decay |
| read3903415_x1 | 1.0 | 0.643 | 0.968 | ENST00000400454.6 | ENST00000676846.1 | NOMO3 | nonsense_mediated_decay | protein_coding;processed_transcript;nonsense_mediated_decay;retained_intron |
| read589355_x1 | 1.0 | 0.601 | 0.988 | ENST00000682841.1 | ENST00000682841.1 | IL17D | protein_coding | protein_coding;processed_transcript |

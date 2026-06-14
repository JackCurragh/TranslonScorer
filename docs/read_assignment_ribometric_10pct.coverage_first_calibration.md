# Coverage-First Calibration Scan: read_assignment_ribometric_10pct

## Takeaway

Yes: there are better-covered features to calibrate on than the previous high-posterior-shift examples.

The earlier scan asked where frame-aware assignment moved reads the most. That over-selected sparse or odd annotation cases. This scan asks a better first question: where do we have many ambiguous reads in protein-coding transcript candidates, and do those reads land in codons with independent frame/profile support?

The answer is that same-locus transcript ambiguity is available at useful coverage in genes such as the top entries below. Those are the right first calibration cases. Multi-locus/paralog-heavy cases are still useful, but as controls for mappability and locus-origin ambiguity rather than as isoform validation wins.

Important caveat: this still uses the cached RiboMetric 10 percent BAM candidate table, which failed the global P-site/frame audit. So these cases are good for assignment-method calibration and profile inspection, not yet for biological claims about true translated frame until the same loci are rerun with validated P-site coordinates.

## Gates Used

- Local frame support gate: `total_count >= 10` and `support_evidence >= 0.2`.
- Frame-informative read gate: gated frame-likelihood spread >= `0.15` across that read's candidates.
- Main calibration class: more than one transcript candidate, exactly one locus candidate, and at least one protein-coding candidate.
- Multi-locus candidate rows are kept separately as controls.

## Global Counts

- Read keys in candidate table: 358,510.
- Transcript-ambiguous read keys: 317,390.
- Same-locus protein-coding ambiguous read keys: 238,949.
- Same-locus protein-coding weighted reads: 393880.0.
- Same-locus weighted reads with gated local frame support: 7.0.
- Same-locus weighted reads with frame-informative gated support: 6.0.
- Multi-locus protein-coding read keys: 73,443.
- Multi-locus protein-coding weighted reads: 193396.0.

## Panel Seed

| role | gene | read_keys | weighted_reads | candidate_rows_weight | transcripts | loci | gated_fraction | informative_fraction | max_support | score |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| multi_locus_mappability_control | ACTB | 1,212 |  | 111431.0 | 15 | 1 | 0.000 | 0.000 | 0.000 | 11.621 |
| multi_locus_mappability_control | ACTG1 | 963 |  | 93538.0 | 12 | 1 | 0.000 | 0.000 | 0.000 | 11.446 |
| multi_locus_mappability_control | ZNF354B | 171 |  | 38513.0 | 2 | 1 | 0.000 | 0.000 | 0.000 | 10.559 |
| multi_locus_mappability_control | ANXA2 | 329 |  | 20148.0 | 30 | 1 | 0.000 | 0.000 | 0.000 | 9.911 |
| multi_locus_mappability_control | EEF1A1 | 508 |  | 18522.0 | 14 | 1 | 0.000 | 0.000 | 0.000 | 9.827 |
| multi_locus_mappability_control | UBC | 230 |  | 13811.0 | 11 | 1 | 0.000 | 0.000 | 0.000 | 9.533 |
| multi_locus_mappability_control | CSNK1E | 178 |  | 12924.0 | 10 | 1 | 0.043 | 0.038 | 0.000 | 9.467 |
| multi_locus_mappability_control | POTEE | 640 |  | 12371.0 | 3 | 1 | 0.000 | 0.000 | 0.000 | 9.423 |
| same_locus_isoform_calibration | TGFBI | 1,980 | 16106.0 | 88882.0 | 8 | 1 | 0.000 | 0.000 | 0.949 | 11.914 |
| same_locus_isoform_calibration | CDK6 | 36 | 18082.0 | 54235.0 | 3 | 1 | 0.000 | 0.000 | 0.000 | 9.803 |
| same_locus_isoform_calibration | FN1 | 2,632 | 8879.0 | 99294.0 | 16 | 1 | 0.000 | 0.000 | 0.000 | 9.092 |
| same_locus_isoform_calibration | TMC5 | 54 | 6706.0 | 26824.0 | 4 | 1 | 0.000 | 0.000 | 0.000 | 8.811 |
| same_locus_isoform_calibration | S100A6 | 334 | 5419.0 | 25478.0 | 5 | 1 | 0.000 | 0.000 | 0.000 | 8.598 |
| same_locus_isoform_calibration | ANPEP | 1,385 | 4608.0 | 22981.0 | 7 | 1 | 0.000 | 0.000 | 0.000 | 8.436 |
| same_locus_isoform_calibration | TIMP1 | 510 | 3864.0 | 12569.0 | 4 | 1 | 0.000 | 0.000 | 0.000 | 8.260 |
| same_locus_isoform_calibration | HRAS | 485 | 2998.0 | 21855.0 | 8 | 1 | 0.000 | 0.000 | 0.000 | 8.006 |
| same_locus_isoform_calibration | ANGPTL4 | 680 | 2895.0 | 19560.0 | 10 | 1 | 0.000 | 0.000 | 0.000 | 7.971 |
| same_locus_isoform_calibration | LGALS1 | 295 | 2023.0 | 9029.0 | 7 | 1 | 0.000 | 0.000 | 0.000 | 7.959 |
| same_locus_isoform_calibration | LAMB3 | 1,047 | 2536.0 | 8382.0 | 4 | 1 | 0.000 | 0.000 | 0.000 | 7.839 |
| same_locus_isoform_calibration | CTC1 | 89 | 2267.0 | 6813.0 | 6 | 1 | 0.000 | 0.000 | 0.000 | 7.727 |
| same_locus_isoform_calibration | KRT7 | 424 | 1677.0 | 6080.0 | 6 | 1 | 0.000 | 0.000 | 0.000 | 7.599 |
| same_locus_isoform_calibration | LAMC2 | 862 | 1935.0 | 4029.0 | 3 | 1 | 0.000 | 0.000 | 0.000 | 7.568 |
| same_locus_isoform_calibration | ITGA5 | 841 | 1889.0 | 5974.0 | 6 | 1 | 0.000 | 0.000 | 0.000 | 7.544 |

## Top Same-Locus Isoform Calibration Genes

| gene | read_keys | weighted_reads | candidate_rows_weight | mean_tx | max_tx | gated_fraction | informative_fraction | max_support | max_support_count | score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TGFBI | 1,980 | 16106.0 | 88882.0 | 5.56 | 8 | 0.000 | 0.000 | 0.949 | 165.0 | 11.914 |
| CDK6 | 36 | 18082.0 | 54235.0 | 2.69 | 3 | 0.000 | 0.000 | 0.000 | 0.0 | 9.803 |
| FN1 | 2,632 | 8879.0 | 99294.0 | 11.31 | 16 | 0.000 | 0.000 | 0.000 | 0.0 | 9.092 |
| TMC5 | 54 | 6706.0 | 26824.0 | 4.00 | 4 | 0.000 | 0.000 | 0.000 | 0.0 | 8.811 |
| S100A6 | 334 | 5419.0 | 25478.0 | 4.37 | 5 | 0.000 | 0.000 | 0.000 | 0.0 | 8.598 |
| ANPEP | 1,385 | 4608.0 | 22981.0 | 4.90 | 7 | 0.000 | 0.000 | 0.000 | 0.0 | 8.436 |
| TIMP1 | 510 | 3864.0 | 12569.0 | 3.27 | 4 | 0.000 | 0.000 | 0.000 | 0.0 | 8.260 |
| HRAS | 485 | 2998.0 | 21855.0 | 7.23 | 8 | 0.000 | 0.000 | 0.000 | 0.0 | 8.006 |
| ANGPTL4 | 680 | 2895.0 | 19560.0 | 6.59 | 10 | 0.000 | 0.000 | 0.000 | 0.0 | 7.971 |
| LGALS1 | 295 | 2023.0 | 9029.0 | 4.66 | 7 | 0.000 | 0.000 | 0.000 | 3.0 | 7.959 |
| LAMB3 | 1,047 | 2536.0 | 8382.0 | 3.34 | 4 | 0.000 | 0.000 | 0.000 | 0.0 | 7.839 |
| CTC1 | 89 | 2267.0 | 6813.0 | 3.13 | 6 | 0.000 | 0.000 | 0.000 | 0.0 | 7.727 |

## Top Candidate Transcripts Inside Same-Locus Cases

| gene | transcript | biotype | read_keys | candidate_rows_weight | gated_weight_fraction | unique_frame_count | max_unique_support | local_max_support | pos_range | score |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| TGFBI | ENST00000442011.7 | protein_coding | 1,980 | 16106.0 | 0.000 | 440.0 | 0.984 | 0.949 | 3-2522 | 11.889 |
| TGFBI | ENST00000506699.5 | retained_intron | 1,960 | 15991.0 | 0.000 | 0.0 | 0.000 | 0.000 | 0-2971 | 9.680 |
| TGFBI | ENST00000507018.5 | nonsense_mediated_decay | 1,791 | 14118.0 | 0.000 | 0.0 | 0.000 | 0.000 | 0-2432 | 9.555 |
| TGFBI | ENST00000514554.5 | protein_coding | 999 | 7856.0 | 0.000 | 0.0 | 0.000 | 0.000 | 0-1501 | 8.969 |
| TGFBI | ENST00000515433.1 | retained_intron | 922 | 7271.0 | 0.000 | 0.0 | 0.000 | 0.000 | 582-4297 | 8.892 |
| TGFBI | ENST00000504185.5 | processed_transcript | 450 | 4541.0 | 0.000 | 0.0 | 0.000 | 0.000 | 3-531 | 8.421 |
| TGFBI | ENST00000604555.5 | protein_coding | 535 | 4103.0 | 0.000 | 0.0 | 0.000 | 0.000 | 0-606 | 8.320 |
| TGFBI | ENST00000509485.5 | nonsense_mediated_decay | 531 | 3984.0 | 0.000 | 0.0 | 0.000 | 0.000 | 1-701 | 8.290 |
| TGFBI | ENST00000514242.5 | retained_intron | 409 | 3451.0 | 0.000 | 0.0 | 0.000 | 0.000 | 179-542 | 8.147 |
| TGFBI | ENST00000508767.5 | protein_coding | 405 | 3004.0 | 0.000 | 0.0 | 0.000 | 0.000 | 4-553 | 8.008 |
| TGFBI | ENST00000508076.5 | protein_coding | 326 | 2680.0 | 0.000 | 0.0 | 0.000 | 0.000 | 69-847 | 7.894 |
| TGFBI | ENST00000513497.1 | retained_intron | 221 | 2068.0 | 0.000 | 0.0 | 0.000 | 0.000 | 320-562 | 7.635 |
| TGFBI | ENST00000503087.1 | protein_coding | 192 | 1593.0 | 0.000 | 0.0 | 0.000 | 0.000 | 2-362 | 7.374 |
| TGFBI | ENST00000509749.1 | retained_intron | 159 | 1381.0 | 0.000 | 0.0 | 0.000 | 0.000 | 1-182 | 7.231 |
| TGFBI | ENST00000504411.1 | retained_intron | 124 | 735.0 | 0.000 | 0.0 | 0.000 | 0.000 | 306-858 | 6.601 |
| CDK6 | ENST00000424848.3 | protein_coding | 36 | 18082.0 | 0.000 | 1.0 | 0.026 | 0.000 | 94-1420 | 9.968 |
| CDK6 | ENST00000265734.8 | protein_coding | 36 | 18082.0 | 0.000 | 0.0 | 0.000 | 0.000 | 43-1369 | 9.803 |
| CDK6 | ENST00000491250.1 | processed_transcript | 18 | 18063.0 | 0.000 | 0.0 | 0.000 | 0.000 | 268-495 | 9.802 |
| CDK6 | ENST00000467166.1 | retained_intron | 5 | 5.0 | 0.000 | 0.0 | 0.000 | 0.000 | 69-329 | 1.792 |
| CDK6 | ENST00000473078.1 | retained_intron | 2 | 3.0 | 0.000 | 0.0 | 0.000 | 0.000 | 24-96 | 1.386 |

## Multi-Locus Controls

These are deliberately not mixed into same-locus isoform validation. They are useful for checking whether the model reports genomic-origin ambiguity instead of pretending it resolved an isoform question.

| gene | biotype | read_keys | candidate_rows_weight | transcripts | loci | gated_fraction | informative_fraction | max_support |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| ACTB | protein_coding | 1,212 | 111431.0 | 15 | 1 | 0.000 | 0.000 | 0.000 |
| ACTG1 | protein_coding | 963 | 93538.0 | 12 | 1 | 0.000 | 0.000 | 0.000 |
| ACTG1 | nonsense_mediated_decay | 963 | 54802.0 | 7 | 1 | 0.000 | 0.000 | 0.000 |
| ACTG1 | retained_intron | 963 | 46766.0 | 6 | 1 | 0.000 | 0.000 | 0.000 |
| RNU1-1 | snRNA | 180 | 39590.0 | 1 | 1 | 0.000 | 0.000 | 0.000 |
| RNU1-2 | snRNA | 180 | 39590.0 | 1 | 1 | 0.000 | 0.000 | 0.000 |
| RNU1-3 | snRNA | 180 | 39590.0 | 1 | 1 | 0.000 | 0.000 | 0.000 |
| RNVU1-29 | snRNA | 180 | 39590.0 | 1 | 1 | 0.000 | 0.000 | 0.000 |
| RNU1-27P | snRNA | 180 | 39590.0 | 1 | 1 | 0.000 | 0.000 | 0.000 |
| RNVU1-18 | snRNA | 180 | 39590.0 | 1 | 1 | 0.000 | 0.000 | 0.000 |
| RNU1-4 | snRNA | 180 | 39590.0 | 1 | 1 | 0.000 | 0.000 | 0.000 |
| RNVU1-7 | snRNA | 180 | 39590.0 | 1 | 1 | 0.000 | 0.000 | 0.000 |

## Decision Change

The validation order should change from `largest model movement -> inspect whether the case is meaningful` to `coverage and biological class gate -> inspect model movement`. That directly addresses the problem that sparse positions are often uninformative once ambiguity is resolved.

Immediate next step: generate before/after profile plots for the panel-seed same-locus genes, using the standard phase-coloured Ribo-seq profile style, and use those plots to decide which features deserve a deeper read-level rerun with validated P-site offsets.

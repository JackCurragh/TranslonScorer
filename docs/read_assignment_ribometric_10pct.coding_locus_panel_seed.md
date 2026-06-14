# Coding Locus Panel Seed: read_assignment_ribometric_10pct

## Purpose

This is not a frozen validation panel. It is a review queue generated from the scaled coding-shift scan, intended to identify which coding-associated frame-induced assignment shifts deserve manual curation next.

## Summary

- Candidate gene labels: 35.
- Lower-locus-ambiguity candidates: 10.
- Broad multi-locus ambiguity candidates: 16.
- Candidates with direct local frame support: 27.

## Top Review Candidates

| gene_label | biotype | support | ambiguity | weighted_reads | local_supported_weight | max_tv | mean_transcripts | mean_loci | priority |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| ESRRAP2 | transcribed_processed_pseudogene | direct_local_frame_support | moderate_locus_ambiguity | 388.0 | 374.0 | 0.913 | 6.73 | 2.93 | 2.235 |
| OLFM1 | protein_coding | direct_local_frame_support | moderate_locus_ambiguity | 144.0 | 144.0 | 0.648 | 6.75 | 4.00 | 1.236 |
| CDV3 | protein_coding | direct_local_frame_support | broad_multilocus_ambiguity | 250.0 | 249.0 | 0.674 | 13.50 | 8.00 | 1.161 |
| IL17D | protein_coding | direct_local_frame_support | moderate_locus_ambiguity | 51.0 | 51.0 | 0.706 | 8.50 | 3.25 | 1.141 |
| TMEM256-PLSCR3 | nonsense_mediated_decay | direct_local_frame_support | moderate_locus_ambiguity | 28.0 | 26.0 | 0.739 | 9.00 | 5.00 | 0.843 |
| TMEM271 | protein_coding | direct_local_frame_support | broad_multilocus_ambiguity | 235.0 | 235.0 | 0.551 | 63.74 | 28.57 | 0.686 |
| POGK | protein_coding | direct_local_frame_support | lower_locus_ambiguity | 18.0 | 18.0 | 0.365 | 6.00 | 2.00 | 0.512 |
| DSCAM | protein_coding | direct_local_frame_support | moderate_locus_ambiguity | 74.0 | 74.0 | 0.229 | 2.75 | 2.25 | 0.453 |
| H1-2 | protein_coding | direct_local_frame_support | lower_locus_ambiguity | 8.0 | 8.0 | 0.349 | 2.00 | 2.00 | 0.365 |
| SIX1 | protein_coding | direct_local_frame_support | moderate_locus_ambiguity | 5.0 | 5.0 | 0.309 | 5.00 | 3.00 | 0.232 |
| OAZ2 | protein_coding | direct_local_frame_support | lower_locus_ambiguity | 1.0 | 1.0 | 0.692 | 2.00 | 2.00 | 0.229 |
| TLN1 | protein_coding | direct_local_frame_support | broad_multilocus_ambiguity | 3.0 | 1.0 | 0.740 | 7.33 | 6.00 | 0.174 |
| AGO4 | protein_coding | direct_local_frame_support | broad_multilocus_ambiguity | 9.0 | 9.0 | 0.304 | 65.25 | 24.25 | 0.165 |
| H1-4 | protein_coding | direct_local_frame_support | lower_locus_ambiguity | 1.0 | 1.0 | 0.474 | 2.00 | 2.00 | 0.157 |
| MT-CO1 | protein_coding | direct_local_frame_support | lower_locus_ambiguity | 9.0 | 9.0 | 0.140 | 2.00 | 2.00 | 0.154 |
| PGP | protein_coding | direct_local_frame_support | broad_multilocus_ambiguity | 34.0 | 34.0 | 0.145 | 25.00 | 11.00 | 0.148 |
| ZBED4 | protein_coding | direct_local_frame_support | moderate_locus_ambiguity | 1.0 | 1.0 | 0.577 | 17.00 | 5.00 | 0.143 |
| PTCH1 | protein_coding | abundance_propagated_shift | lower_locus_ambiguity | 3.0 | 0.0 | 0.686 | 8.33 | 1.00 | 0.140 |
| GNAS | protein_coding | direct_local_frame_support | broad_multilocus_ambiguity | 2.0 | 2.0 | 0.451 | 20.00 | 12.00 | 0.139 |
| NDUFS5 | protein_coding | direct_local_frame_support | lower_locus_ambiguity | 2.0 | 2.0 | 0.208 | 2.00 | 1.00 | 0.135 |

## Lower-Locus-Ambiguity Subset

| gene_label | biotype | support | weighted_reads | local_supported_weight | max_tv | mean_transcripts | mean_loci | reason |
|---|---|---|---:|---:|---:|---:|---:|---|
| POGK | protein_coding | direct_local_frame_support | 18.0 | 18.0 | 0.365 | 6.00 | 2.00 | weighted_reads=18.0; max_tv=0.365; mean_candidate_loci=2.0 |
| H1-2 | protein_coding | direct_local_frame_support | 8.0 | 8.0 | 0.349 | 2.00 | 2.00 | weighted_reads=8.0; max_tv=0.349; mean_candidate_loci=2.0 |
| OAZ2 | protein_coding | direct_local_frame_support | 1.0 | 1.0 | 0.692 | 2.00 | 2.00 | weighted_reads=1.0; max_tv=0.692; mean_candidate_loci=2.0 |
| H1-4 | protein_coding | direct_local_frame_support | 1.0 | 1.0 | 0.474 | 2.00 | 2.00 | weighted_reads=1.0; max_tv=0.474; mean_candidate_loci=2.0 |
| MT-CO1 | protein_coding | direct_local_frame_support | 9.0 | 9.0 | 0.140 | 2.00 | 2.00 | weighted_reads=9.0; max_tv=0.14; mean_candidate_loci=2.0 |
| PTCH1 | protein_coding | abundance_propagated_shift | 3.0 | 0.0 | 0.686 | 8.33 | 1.00 | weighted_reads=3.0; max_tv=0.686; mean_candidate_loci=1.0 |
| NDUFS5 | protein_coding | direct_local_frame_support | 2.0 | 2.0 | 0.208 | 2.00 | 1.00 | weighted_reads=2.0; max_tv=0.208; mean_candidate_loci=1.0 |
| RCE1 | retained_intron | abundance_propagated_shift | 9.0 | 0.0 | 0.272 | 5.56 | 1.00 | weighted_reads=9.0; max_tv=0.272; mean_candidate_loci=1.0 |
| FLYWCH2 | protein_coding | abundance_propagated_shift | 7.0 | 0.0 | 0.113 | 8.00 | 2.00 | weighted_reads=7.0; max_tv=0.113; mean_candidate_loci=2.0 |
| PLSCR3 | retained_intron | abundance_propagated_shift | 1.0 | 0.0 | 0.232 | 5.00 | 2.00 | weighted_reads=1.0; max_tv=0.232; mean_candidate_loci=2.0 |

## Interpretation

The best immediate manual-review targets are lower-locus-ambiguity coding candidates, because broad multi-locus cases mostly test genomic mappability rather than isoform deconvolution. High-priority broad cases are still useful, but they should be routed to mappability/paralog review before being used as evidence for transcript-origin assignment.

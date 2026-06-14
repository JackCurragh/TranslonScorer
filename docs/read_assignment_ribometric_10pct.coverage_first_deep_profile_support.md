# Deep-Profile Frame Support For Coverage-First Panel: read_assignment_ribometric_10pct

## Takeaway

There are better-covered features, but the coverage has to come from the high-depth transcript profiles rather than the sparse unique-read frame table produced by the current read-level probe.

This run takes the coverage-first same-locus panel, pulls those transcripts from the full forward/reverse RiboCrypt profile parquets, builds frame posteriors from those aggregate profiles, then asks how much of the ambiguous read table is covered by those stronger frame posteriors.

This is still an interim calibration layer: aggregate profiles can tell us whether a transcript has strong frame structure, but they cannot by themselves prove which ambiguous read originated from which transcript. The final validation still needs a read-level rerun with validated P-site offsets and multi-mappers preserved.

Design rule learned from this pass: frame evidence is only used when every candidate target for a read has gated frame support. Otherwise the model can accidentally prefer a no-support transcript over a supported transcript simply because the supported transcript has locally incompatible frame evidence.

## Inputs

- Profile parquets: `/Users/jackt/projects/all-RiboSeq/ribocrypt_fwd/ribocrypt_fwd_transcript_profiles.parquet`, `/Users/jackt/projects/all-RiboSeq/ribocrypt_rev/ribocrypt_rev_transcript_profiles.parquet`.
- CDS annotation: `/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/gencode.v25.annotation_RiboMetric.tsv`.
- Same-locus genes selected: TGFBI, CDK6, FN1, TMC5, S100A6, ANPEP, TIMP1, HRAS, ANGPTL4, LGALS1.
- Transcripts pulled: 97.
- Candidate/read overlap gate: profile-frame `total_count >= 100`, `support_evidence >= 0.2`, complete frame support for every candidate target of the read, and candidate frame-likelihood spread >= `0.15`.

## Candidate-Read Coverage By Deep Profile Frame Support

| gene | read_keys | weighted_reads | supported_weight | supported_fraction | informative_weight | informative_fraction | max_support | max_support_count | mean_spread |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TGFBI | 1,986 | 16112.0 | 14510.0 | 0.901 | 559.0 | 0.035 | 1.000 | 423694.0 | 0.029 |
| TIMP1 | 511 | 3865.0 | 3836.0 | 0.992 | 195.0 | 0.050 | 1.000 | 488087.0 | 0.029 |
| FN1 | 2,640 | 8887.0 | 3696.0 | 0.416 | 3696.0 | 0.416 | 1.000 | 217371.0 | 0.439 |
| TMC5 | 63 | 9035.0 | 2325.0 | 0.257 | 2.0 | 0.000 | 1.000 | 225478.0 | 0.015 |
| LGALS1 | 295 | 2023.0 | 1993.0 | 0.985 | 64.0 | 0.032 | 1.000 | 1156232.0 | 0.033 |
| ANGPTL4 | 680 | 2895.0 | 1883.0 | 0.650 | 289.0 | 0.100 | 0.992 | 30579.0 | 0.043 |
| ANPEP | 1,392 | 4629.0 | 1307.0 | 0.282 | 1294.0 | 0.280 | 1.000 | 49206.0 | 0.169 |
| HRAS | 497 | 3091.0 | 157.0 | 0.051 | 157.0 | 0.051 | 0.887 | 32011.0 | 0.034 |
| S100A6 | 337 | 5422.0 | 145.0 | 0.027 | 144.0 | 0.027 | 1.000 | 1842220.0 | 0.035 |
| CDK6 | 38 | 18085.0 | 0.0 | 0.000 | 0.0 | 0.000 | 0.000 | 24812.0 | 0.000 |

## Target Transcript Profile Coverage

| gene | transcript | biotype | profile_count | profile_rows | candidate_rows_weight | has_cds | profile_pos_range |
|---|---|---|---:|---:|---:|---|---|
| ANGPTL4 | ENST00000301455.7 | protein_coding | 1247307.0 | 1,802 | 2889.0 | False | 7-1864 |
| ANGPTL4 | ENST00000593998.5 | nonsense_mediated_decay | 1245276.0 | 1,741 | 2888.0 | True | 0-1790 |
| ANGPTL4 | ENST00000595079.5 | nonsense_mediated_decay | 1202768.0 | 1,677 | 2722.0 | True | 6-1716 |
| ANGPTL4 | ENST00000393962.6 | protein_coding | 1126884.0 | 1,660 | 2592.0 | True | 0-1704 |
| ANGPTL4 | ENST00000594348.1 | retained_intron | 860560.0 | 1,004 | 2021.0 | False | 0-1020 |
| ANGPTL4 | ENST00000598255.5 | retained_intron | 742097.0 | 892 | 1725.0 | False | 0-908 |
| ANGPTL4 | ENST00000597137.5 | retained_intron | 485276.0 | 549 | 945.0 | False | 2-581 |
| ANGPTL4 | ENST00000599192.5 | protein_coding | 331354.0 | 521 | 1045.0 | True | 2-570 |
| ANGPTL4 | ENST00000594875.1 | protein_coding | 320929.0 | 486 | 847.0 | True | 0-523 |
| ANGPTL4 | ENST00000601886.1 | protein_coding | 313360.0 | 497 | 960.0 | True | 0-551 |
| ANGPTL4 | ENST00000601770.1 | protein_coding | 294684.0 | 436 | 926.0 | True | 84-812 |
| ANPEP | ENST00000679248.1 | protein_coding | 3808814.0 | 3,806 | 4608.0 | False | 8-3842 |
| ANPEP | ENST00000560137.2 | protein_coding | 3808461.0 | 3,783 | 4608.0 | False | 8-3812 |
| ANPEP | ENST00000300060.7 | protein_coding | 3806015.0 | 3,633 | 4608.0 | False | 8-3661 |
| ANPEP | ENST00000559874.2 | protein_coding | 3799700.0 | 3,613 | 4600.0 | False | 8-3677 |
| ANPEP | ENST00000560030.1 | nonsense_mediated_decay | 1169472.0 | 873 | 1164.0 | True | 0-895 |
| ANPEP | ENST00000558177.5 | processed_transcript | 723387.0 | 745 | 615.0 | False | 0-915 |
| ANPEP | ENST00000558740.1 | retained_intron | 617494.0 | 546 | 664.0 | False | 0-572 |
| ANPEP | ENST00000560028.1 | nonsense_mediated_decay | 590860.0 | 522 | 1016.0 | True | 0-581 |
| ANPEP | ENST00000559887.1 | retained_intron | 432496.0 | 415 | 721.0 | False | 0-533 |
| ANPEP | ENST00000559761.5 | retained_intron | 431198.0 | 368 | 377.0 | False | 0-467 |
| CDK6 | ENST00000424848.3 | protein_coding | 2722547.0 | 11,265 | 18082.0 | False | 8-11662 |
| CDK6 | ENST00000265734.8 | protein_coding | 2713337.0 | 11,214 | 18082.0 | True | 8-11611 |
| CDK6 | ENST00000467166.1 | retained_intron | 588665.0 | 723 | 5.0 | False | 15-775 |
| CDK6 | ENST00000491250.1 | processed_transcript | 502289.0 | 511 | 18063.0 | False | 0-557 |
| CDK6 | ENST00000473078.1 | retained_intron | 295011.0 | 235 | 3.0 | False | 37-385 |
| FN1 | ENST00000354785.11 | protein_coding | 56823160.0 | 8,172 | 8869.0 | False | 0-8389 |
| FN1 | ENST00000323926.10 | protein_coding | 56642514.0 | 8,187 | 8863.0 | True | 5-8707 |
| FN1 | ENST00000359671.5 | protein_coding | 56387948.0 | 8,007 | 8847.0 | True | 1-8523 |
| FN1 | ENST00000336916.8 | protein_coding | 56207082.0 | 7,914 | 8841.0 | True | 5-8434 |

## Linear+HMM Frame Support Summary

| gene | transcript | biotype | frame_count | rows | max_support | mean_support | mean_pmax | mean_entropy |
|---|---|---|---:|---:|---:|---:|---:|---:|
| ANGPTL4 | ENST00000301455.7 | protein_coding | 1247307.0 | 617 | 0.980 | 0.229 | 0.949 | 0.171 |
| ANGPTL4 | ENST00000593998.5 | nonsense_mediated_decay | 1245276.0 | 596 | 0.990 | 0.235 | 0.955 | 0.147 |
| ANGPTL4 | ENST00000595079.5 | nonsense_mediated_decay | 1202768.0 | 570 | 0.992 | 0.236 | 0.953 | 0.156 |
| ANGPTL4 | ENST00000393962.6 | protein_coding | 1126884.0 | 568 | 0.990 | 0.223 | 0.949 | 0.168 |
| ANGPTL4 | ENST00000594348.1 | retained_intron | 860560.0 | 340 | 0.990 | 0.281 | 0.976 | 0.088 |
| ANGPTL4 | ENST00000598255.5 | retained_intron | 742097.0 | 302 | 0.990 | 0.272 | 0.972 | 0.099 |
| ANGPTL4 | ENST00000597137.5 | retained_intron | 485276.0 | 193 | 0.802 | 0.286 | 0.972 | 0.086 |
| ANGPTL4 | ENST00000599192.5 | protein_coding | 331354.0 | 189 | 0.980 | 0.191 | 0.963 | 0.131 |
| ANGPTL4 | ENST00000594875.1 | protein_coding | 320929.0 | 172 | 0.777 | 0.225 | 0.966 | 0.118 |
| ANGPTL4 | ENST00000601886.1 | protein_coding | 313360.0 | 179 | 0.990 | 0.188 | 0.955 | 0.161 |
| ANGPTL4 | ENST00000601770.1 | protein_coding | 294684.0 | 158 | 0.992 | 0.195 | 0.947 | 0.168 |
| ANPEP | ENST00000679248.1 | protein_coding | 3808814.0 | 1,279 | 1.000 | 0.317 | 0.973 | 0.101 |
| ANPEP | ENST00000560137.2 | protein_coding | 3808461.0 | 1,269 | 1.000 | 0.320 | 0.973 | 0.100 |
| ANPEP | ENST00000300060.7 | protein_coding | 3806015.0 | 1,219 | 1.000 | 0.332 | 0.974 | 0.096 |
| ANPEP | ENST00000559874.2 | protein_coding | 3799700.0 | 1,221 | 1.000 | 0.331 | 0.975 | 0.094 |
| ANPEP | ENST00000560030.1 | nonsense_mediated_decay | 1169472.0 | 295 | 1.000 | 0.408 | 0.997 | 0.014 |
| ANPEP | ENST00000558177.5 | processed_transcript | 723387.0 | 280 | 0.926 | 0.283 | 0.978 | 0.078 |
| ANPEP | ENST00000558740.1 | retained_intron | 617494.0 | 188 | 0.864 | 0.362 | 0.990 | 0.034 |
| ANPEP | ENST00000560028.1 | nonsense_mediated_decay | 590860.0 | 186 | 0.983 | 0.334 | 0.977 | 0.074 |
| ANPEP | ENST00000559887.1 | retained_intron | 432496.0 | 159 | 0.983 | 0.272 | 0.959 | 0.122 |
| ANPEP | ENST00000559761.5 | retained_intron | 431198.0 | 133 | 0.833 | 0.357 | 0.990 | 0.044 |
| CDK6 | ENST00000424848.3 | protein_coding | 2722547.0 | 3,869 | 1.000 | 0.057 | 0.870 | 0.430 |
| CDK6 | ENST00000265734.8 | protein_coding | 2713337.0 | 3,852 | 1.000 | 0.058 | 0.871 | 0.428 |
| CDK6 | ENST00000467166.1 | retained_intron | 588665.0 | 249 | 0.992 | 0.203 | 0.903 | 0.314 |
| CDK6 | ENST00000491250.1 | processed_transcript | 502289.0 | 182 | 1.000 | 0.100 | 0.905 | 0.315 |
| CDK6 | ENST00000473078.1 | retained_intron | 295011.0 | 92 | 0.806 | 0.341 | 0.976 | 0.065 |
| FN1 | ENST00000354785.11 | protein_coding | 56823160.0 | 2,734 | 1.000 | 0.798 | 0.987 | 0.047 |
| FN1 | ENST00000323926.10 | protein_coding | 56642514.0 | 2,786 | 1.000 | 0.778 | 0.987 | 0.049 |
| FN1 | ENST00000359671.5 | protein_coding | 56387948.0 | 2,719 | 1.000 | 0.789 | 0.984 | 0.054 |
| FN1 | ENST00000336916.8 | protein_coding | 56207082.0 | 2,695 | 1.000 | 0.788 | 0.986 | 0.050 |

## Files

- `target_transcripts`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.target_transcripts.csv`
- `target_profiles`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.target_profiles.parquet`
- `profile_summary`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.profile_summary.csv`
- `support_summary`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.frame_support_summary.csv`
- `candidate_overlap_reads`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.candidate_overlap_reads.parquet`
- `candidate_overlap_by_gene`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.candidate_overlap_by_gene.csv`
- `linear_frame_support`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.linear_frame_support.parquet`
- `linear_hmm_frame_support`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.linear_hmm_frame_support.parquet`
- `latent_frame_support`: `/Users/jackt/projects/all-RiboSeq/translonscorer/notebooks/read_assignment_ribometric_10pct.coverage_first.deep_profile.latent_frame_support.parquet`

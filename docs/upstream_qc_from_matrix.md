# Investigation: What can the matrix tell us about upstream processing quality?

**Status**: Open — actively investigating  
**Branch**: `local/read-assignment-prototype`

## Question

Can we learn what worked and what didn't in upstream steps (trimming, alignment, library preparation) directly from the multi-sample sparse Parquet matrix — without re-running per-sample pipelines?

## Motivation

The 6k-sample cohort contains a wide spread of Ribo-seq quality. Some samples with weak periodicity scores may be salvageable if the root cause is an upstream processing failure (e.g. poor adapter trimming, wrong library protocol) rather than biology. Identifying those cases from the matrix avoids re-running RiboMetric on 6k BAMs.

## What we can compute from the matrix (no BAM)

All of these are now implemented in `calibrate_cohort()` (Phase 1 of `build-psite-index`):

| Metric | Upstream step it reflects | Where computed |
|--------|--------------------------|----------------|
| `rpf_28_32_prop` | Trimming — fraction of reads in the expected RPF window | `fast_reads_qc` |
| `peak_length`, `mean_length` | Trimming — are reads the right length? | `fast_reads_qc` |
| `rld_CV_metric`, `rld_IQR_metric` | Trimming quality — tight vs broad length peak | `fast_reads_qc` |
| `rld_bimodality` | Two distinct length populations (e.g. mixed protocols) | `fast_reads_qc` |
| `ligation_bias_KL_5p` / `ligation_bias_score_5p` | Ligation protocol — 5' end nucleotide bias | `fast_reads_qc` |
| `ligation_bias_KL_3p` / `ligation_bias_score_3p` | Adapter trimming — 3' adapter remnants surviving | `fast_reads_qc` |
| `periodicity_score` | Overall Ribo-seq signal quality | matrix rollup |
| `prop_cds` | Alignment / rRNA depletion — fraction of reads landing in CDS | rollup + reads QC |
| `library_type` | Composite: `elongation` vs `low_quality` | derived |
| `recommended_lengths` | Which read lengths carry clean periodicity | derived |

## What we cannot compute without extra data

| Metric | Why not feasible |
|--------|-----------------|
| mRNA region distribution (5'UTR/CDS/3'UTR proportions) | Needs dense positional data across all features |
| Metagene profiles | Same — needs position-level coverage across many genes |
| CDS coverage uniformity (Gini, Theil, autocorrelation) | Same |
| FLOSS heterogeneity | Needs reference transcript set |
| Start codon enrichment ratio | Needs initiation-enriched metagene at codon resolution |

## Hypotheses to test

1. **Trimming failure** → `rpf_28_32_prop` low + `ligation_bias_KL_3p` high (adapter remnants at 3' end shift read lengths and bias 3' dinucleotides).
2. **Ligation bias only** → `ligation_bias_score_5p` low but `periodicity_score` and `rpf_28_32_prop` are still good. These samples are usable.
3. **Not Ribo-seq** → `periodicity_score` < 0.1 + `prop_cds` < 0.3. Likely RNA-seq or degraded.
4. **Mixed protocols in one study** → `rld_bimodality` high within a study. Some samples at 28-32nt, others at 35-40nt.
5. **rRNA contamination / poor depletion** → `prop_cds` low despite good length distribution and periodicity at the reads that do align.

## Next steps

- [ ] Inspect `qc_per_sample.parquet` after the new Phase 1 run completes (includes all metrics above)
- [ ] Cluster samples by QC profile to identify failure modes
- [ ] Cross-reference study metadata to see if failure modes correlate with protocol/lab/year
- [ ] Define thresholds for each failure mode → downstream filtering flags
- [ ] Consider whether any failure modes are correctable (e.g. re-trimming or length filtering on the matrix)

## Related docs

- [Matrix query engine](matrix_query_engine.md)
- [Event scoring model](event_scoring_model.md)

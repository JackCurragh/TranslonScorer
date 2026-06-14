# LINC01783 Frame-Shift Case Review

![LINC01783 case review](../notebooks/read_assignment_linc01783_case.png)

## Summary

- High-shift read keys reviewed: 57.
- Compatible transcript targets among those reads: 8.
- Most frequent candidate transcript: `ENST00000415386.2` with 57 read keys.
- Strongest independent frame-support transcript by total unique support: `ENST00000415386.2` with total unique support 29.0 and max support evidence 0.804.

## Interpretation

The independent frame-aware shift is concentrated in a small lncRNA/pseudogene-like ambiguity cluster rather than distributed across many coding loci. The top baseline and frame-aware transcript often remain the same, while posterior mass is redistributed among several compatible transcript models. That makes this a cautionary case: frame evidence is detecting a local frame-pattern difference, but the biological interpretation is weak without stronger independent support and annotation review.

For the manuscript story, this supports the conservative framing: frame-aware assignment can move reads, but the trustworthy independent signal in this probe is sparse and concentrated. The main deliverable remains uncertainty and identifiability reporting until matched production frame support shows broader, biologically interpretable shifts.

## Candidate Transcript Summary

| tran_id | gene_name | read_keys | weighted_candidate_rows | pos_range |
|---|---|---:|---:|---|
| ENST00000415386.2 | LINC01783 | 57 | 443.0 | 298-342 |
| ENST00000438002.1 | ENSG00000228549 | 57 | 443.0 | 1637-1681 |
| ENST00000453554.1 | ENSG00000228549 | 57 | 443.0 | 298-342 |
| ENST00000660907.1 | ENSG00000228549 | 57 | 443.0 | 301-345 |
| ENST00000662795.1 | ENSG00000228549 | 57 | 443.0 | 356-400 |
| ENST00000666785.1 | ENSG00000228549 | 57 | 443.0 | 23-67 |
| ENST00000668460.1 | ENSG00000228549 | 57 | 443.0 | 356-400 |
| ENST00000599640.6 | ENSG00000227733 | 1 | 3.0 | 1965-1965 |

## Independent Frame Support

| tran_id | frame_rows | total_unique_support | max_support_evidence | mean_support_evidence |
|---|---:|---:|---:|---:|
| ENST00000415386.2 | 3 | 29.0 | 0.804 | 0.325 |

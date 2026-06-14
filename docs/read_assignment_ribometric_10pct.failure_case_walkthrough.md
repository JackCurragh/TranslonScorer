# Read Assignment Failure Case Walkthrough: read_assignment_ribometric_10pct

Notebook: `../notebooks/read_assignment_failure_case_walkthrough.ipynb`
Diagnostic plot: `../notebooks/read_assignment_ribometric_10pct.failure_case.diagnostics.png`

## Takeaway

The current evidence supports ordinary EM more clearly than frame-aware EM. Frame-aware EM is useful only when frame evidence is present at the covered positions and separates the compatible transcripts in the right direction.

The failure cases are not mysterious: the ambiguity graphs can contain tens to hundreds of compatible transcripts, often with very little unique read anchoring. In that setting, the frame term can move mass toward a transcript with a better local frame likelihood even if abundance evidence favored another transcript.

## Case Summary

| scenario | observed effect | EM gain | frame gain | frame comparable | frame favors assumed | frame favors competitor | explanation |
|---|---:|---:|---:|---:|---:|---:|---|
| ANPEP:frame_top | frame increased assumed score | 0.212 | 0.082 | 0.282 | 0.066 | 0.207 | ordinary EM gives a moderate gain; the frame gain is an EM-level redistribution, not a clean read-level vote for the assumed transcript; the ambiguity graph is broad and almost unanchored by unique reads. |
| FN1:frame_top | frame increased assumed score | 0.276 | 0.052 | 0.417 | 0.071 | 0.318 | ordinary EM gives a moderate gain; the frame gain is an EM-level redistribution, not a clean read-level vote for the assumed transcript; the ambiguity graph is broad and almost unanchored by unique reads. |
| HRAS:frame_top | frame mostly unchanged | 0.871 | -0.000 | 0.045 | 0.004 | 0.042 | ordinary EM is doing useful work relative to fractional splitting; very little covered signal passed the frame-evidence gate; the ambiguity graph is broad and almost unanchored by unique reads. |
| HRAS:em_top | frame mostly unchanged | 0.871 | -0.000 | 0.045 | 0.004 | 0.042 | ordinary EM is doing useful work relative to fractional splitting; very little covered signal passed the frame-evidence gate; the ambiguity graph is broad and almost unanchored by unique reads. |
| LGALS1:frame_top | frame mostly unchanged | 0.755 | -0.000 | 0.983 | 0.260 | 0.380 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative. |
| LGALS1:em_top | frame mostly unchanged | 0.755 | -0.000 | 0.983 | 0.260 | 0.380 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative. |
| S100A6:frame_top | frame mostly unchanged | 0.626 | -0.006 | 0.027 | 0.010 | 0.017 | ordinary EM is doing useful work relative to fractional splitting; very little covered signal passed the frame-evidence gate; the ambiguity graph is broad and almost unanchored by unique reads. |
| TGFBI:em_top | frame mostly unchanged | 0.780 | -0.007 | 0.903 | 0.000 | 0.514 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative. |
| TGFBI:frame_top | frame mostly unchanged | 0.780 | -0.007 | 0.903 | 0.000 | 0.514 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative. |
| TIMP1:em_top | frame mostly unchanged | 0.685 | -0.012 | 0.992 | 0.177 | 0.586 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative. |
| TIMP1:frame_top | frame mostly unchanged | 0.685 | -0.012 | 0.992 | 0.177 | 0.586 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative. |
| ANGPTL4:em_top | frame reduced assumed score | 0.381 | -0.042 | 0.650 | 0.141 | 0.393 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative. |
| ANGPTL4:frame_top | frame reduced assumed score | 0.381 | -0.042 | 0.650 | 0.141 | 0.393 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative. |
| ANPEP:em_top | frame reduced assumed score | 0.349 | -0.214 | 0.282 | 0.038 | 0.219 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative; the ambiguity graph is broad and almost unanchored by unique reads. |
| FN1:em_top | frame reduced assumed score | 0.463 | -0.410 | 0.415 | 0.005 | 0.389 | ordinary EM is doing useful work relative to fractional splitting; among comparable reads, the strongest frame likelihood often belongs to an alternative; the ambiguity graph is broad and almost unanchored by unique reads. |
| S100A6:em_top | frame reduced assumed score | 0.747 | -0.623 | 0.035 | 0.000 | 0.024 | ordinary EM is doing useful work relative to fractional splitting; very little covered signal passed the frame-evidence gate; the ambiguity graph is broad and almost unanchored by unique reads. |

## Selected Allocation Plots

![read_assignment_ribometric_10pct.failure_case.ANPEP_frame_top.allocation](../notebooks/read_assignment_ribometric_10pct.failure_case.ANPEP_frame_top.allocation.png)

![read_assignment_ribometric_10pct.failure_case.FN1_frame_top.allocation](../notebooks/read_assignment_ribometric_10pct.failure_case.FN1_frame_top.allocation.png)

![read_assignment_ribometric_10pct.failure_case.S100A6_em_top.allocation](../notebooks/read_assignment_ribometric_10pct.failure_case.S100A6_em_top.allocation.png)

![read_assignment_ribometric_10pct.failure_case.FN1_em_top.allocation](../notebooks/read_assignment_ribometric_10pct.failure_case.FN1_em_top.allocation.png)

![read_assignment_ribometric_10pct.failure_case.HRAS_frame_top.allocation](../notebooks/read_assignment_ribometric_10pct.failure_case.HRAS_frame_top.allocation.png)


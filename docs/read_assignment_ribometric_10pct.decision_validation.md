# Read Assignment Decision Validation: read_assignment_ribometric_10pct

Notebook: `../notebooks/read_assignment_decision_validation.ipynb`
Long summary: `../notebooks/read_assignment_ribometric_10pct.decision_validation.summary.csv`
Wide decision table: `../notebooks/read_assignment_ribometric_10pct.decision_validation.wide.csv`
Stress-test metadata: `../notebooks/read_assignment_ribometric_10pct.decision_validation.real_stress_metadata.csv`

## Takeaway

The current evidence supports ordinary EM more clearly than frame-aware EM. EM is useful because it keeps ambiguous reads while using unique/compatibility structure to avoid raw-compatible overcounting. Frame-aware EM is only useful in regimes where candidate transcripts put the same read into different supported frames.

The real-locus stress tests are not biological truth. They show which loci could benefit from frame evidence if a particular transcript were truly active.

![Toy truth validation](../notebooks/read_assignment_ribometric_10pct.decision_validation.toy_truth.png)

![Real compatibility stress](../notebooks/read_assignment_ribometric_10pct.decision_validation.real_stress.png)

## Largest Real-Compatibility Frame Gains

| scenario | assumed transcript | EM | frame EM | frame gain |
|---|---:|---:|---:|---:|
| ANPEP:frame_top | ENST00000560137.2 | 0.419 | 0.501 | 0.082 |
| FN1:frame_top | ENST00000359671.5 | 0.367 | 0.418 | 0.052 |
| HRAS:frame_top | ENST00000397596.6 | 1.000 | 1.000 | -0.000 |
| HRAS:em_top | ENST00000397596.6 | 1.000 | 1.000 | -0.000 |
| LGALS1:em_top | ENST00000215909.10 | 1.000 | 1.000 | -0.000 |
| LGALS1:frame_top | ENST00000215909.10 | 1.000 | 1.000 | -0.000 |
| S100A6:frame_top | ENST00000368719.9 | 0.842 | 0.836 | -0.006 |
| TGFBI:em_top | ENST00000442011.7 | 1.000 | 0.993 | -0.007 |

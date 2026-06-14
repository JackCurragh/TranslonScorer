# RDG-Gated Frame Assignment Experiment: read_assignment_ribometric_10pct

Notebook: `../notebooks/read_assignment_rdg_gated_assignment_experiment.ipynb`
Plot: `../notebooks/read_assignment_ribometric_10pct.rdg_gated_assignment.validation.png`

## Question

Can the same logic RDG-Flux uses for path identifiability constrain TranslonScorer read assignment, so frame evidence is used only when it actually separates candidate origins?

## Gate Tested

For each read, candidate transcript origins are treated like a small RDG path set. The local gate allows frame evidence to influence EM only when three conditions are met:

1. all candidate transcript origins have comparable frame support at the candidate position;
2. frame likelihood differs by at least 0.25 across candidate origins;
3. exactly one transcript origin has the best frame compatibility.

The stricter consensus gate adds one more requirement: at least 25% of read-weighted evidence must pass the local gate, and the uniquely frame-separable fraction must exceed the alias-group fraction. If this locus-level gate fails, all reads fall back to ordinary EM.

## Readout

![RDG-gated frame assignment takeaway](../notebooks/read_assignment_ribometric_10pct.rdg_gated_assignment.takeaway.png)

Takeaway table:

| pattern | frame gain | local gate gain | consensus gate gain | unique frame target | alias-group contrast | consensus pass rate |
|---|---:|---:|---:|---:|---:|---:|
| Clean separable, correct frame | 0.667 | 0.667 | 0.667 | 1.000 | 0.000 | 1.000 |
| Mixed candidates, correct frame | -0.004 | -0.000 | -0.000 | 0.000 | 1.000 | 0.000 |
| Same-frame alias | 0.000 | -0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Dual-frame overlap | 0.059 | -0.000 | -0.000 | 0.000 | 0.417 | 0.000 |
| Misleading frame evidence | -0.567 | -0.040 | 0.000 | 0.008 | 0.992 | 0.000 |

Detailed scenario plot:

![RDG-gated frame assignment](../notebooks/read_assignment_ribometric_10pct.rdg_gated_assignment.validation.png)

Detailed table:

| validation | frame model | category | EM | frame EM | local gate | consensus gate | frame gain | local gain | consensus gain | scenario gate | reads locally gated | alias-group contrast |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| alias_rich_same_gene | neutral | same_frame_aliased | 0.333 | 0.333 | 0.333 | 0.333 | 0.000 | 0.000 | 0.000 | False | 0.000 | 0.000 |
| alias_rich_same_gene | single_true_frame | same_frame_aliased | 0.333 | 0.333 | 0.333 | 0.333 | 0.000 | -0.000 | 0.000 | False | 0.000 | 0.000 |
| broad_mixed_all_origin | dual_frame_overlap | mixed_realistic | 1.000 | 1.000 | 1.000 | 1.000 | 0.000 | -0.000 | -0.000 | False | 0.000 | 0.455 |
| broad_mixed_all_origin | misleading_competitor_frame | mixed_realistic | 1.000 | 0.000 | 0.877 | 1.000 | -1.000 | -0.123 | 0.000 | False | 0.014 | 0.986 |
| broad_mixed_all_origin | neutral | mixed_realistic | 1.000 | 1.000 | 1.000 | 1.000 | -0.000 | 0.000 | 0.000 | False | 0.000 | 0.000 |
| broad_mixed_all_origin | single_true_frame | mixed_realistic | 1.000 | 1.000 | 1.000 | 1.000 | 0.000 | -0.000 | -0.000 | False | 0.000 | 1.000 |
| clean_rich_all_origin | dual_frame_overlap | clean_frame_informative | 0.333 | 0.452 | 0.333 | 0.333 | 0.119 | 0.000 | 0.000 | False | 0.000 | 0.471 |
| clean_rich_all_origin | misleading_competitor_frame | clean_frame_informative | 0.333 | 0.000 | 0.333 | 0.333 | -0.333 | 0.000 | 0.000 | False | 0.000 | 1.000 |
| clean_rich_all_origin | neutral | clean_frame_informative | 0.333 | 0.333 | 0.333 | 0.333 | -0.000 | 0.000 | 0.000 | False | 0.000 | 0.000 |
| clean_rich_all_origin | single_true_frame | clean_frame_informative | 0.333 | 1.000 | 1.000 | 1.000 | 0.667 | 0.667 | 0.667 | True | 1.000 | 0.000 |
| clean_rich_same_gene | dual_frame_overlap | clean_frame_informative | 0.333 | 0.452 | 0.333 | 0.333 | 0.119 | 0.000 | 0.000 | False | 0.000 | 0.471 |
| clean_rich_same_gene | misleading_competitor_frame | clean_frame_informative | 0.333 | 0.000 | 0.333 | 0.333 | -0.333 | 0.000 | 0.000 | False | 0.000 | 1.000 |
| clean_rich_same_gene | neutral | clean_frame_informative | 0.333 | 0.333 | 0.333 | 0.333 | -0.000 | 0.000 | 0.000 | False | 0.000 | 0.000 |
| clean_rich_same_gene | single_true_frame | clean_frame_informative | 0.333 | 1.000 | 1.000 | 1.000 | 0.667 | 0.667 | 0.667 | True | 1.000 | 0.000 |
| mixed_rich_same_gene | dual_frame_overlap | mixed_realistic | 0.602 | 0.602 | 0.602 | 0.602 | 0.000 | 0.000 | 0.000 | False | 0.000 | 0.271 |
| mixed_rich_same_gene | misleading_competitor_frame | mixed_realistic | 0.602 | 0.002 | 0.565 | 0.602 | -0.600 | -0.037 | 0.000 | False | 0.017 | 0.980 |
| mixed_rich_same_gene | neutral | mixed_realistic | 0.602 | 0.602 | 0.602 | 0.602 | -0.000 | 0.000 | 0.000 | False | 0.000 | 0.000 |
| mixed_rich_same_gene | single_true_frame | mixed_realistic | 0.602 | 0.594 | 0.602 | 0.602 | -0.008 | 0.000 | 0.000 | False | 0.000 | 1.000 |

## What Worked

- Clean frame-informative cases kept the frame benefit. RDG-gated frame EM reached 1.000 mean true-transcript posterior, matching unrestricted frame EM at 1.000.
- Mixed realistic cases stopped using frame as a false discriminator. Unrestricted frame EM averaged -0.004 gain over EM, while the RDG-gated version averaged -0.000.
- Same-frame alias cases stayed unresolved. The gate fraction was 0.000, which is the intended behaviour because frame cannot distinguish those transcripts.
- Dual-frame overlap was not collapsed into a single frame. Reads where frame evidence supported more than one candidate group were treated as ambiguous assignment evidence.

## What Still Fails

- The local gate alone was not enough. In mixed misleading cases, a small fraction of uniquely frame-separable reads could still pull the abundance EM and move many ambiguous reads. The consensus gate blocked that by requiring enough uniquely informative evidence at the locus level.
- The consensus gate still cannot prove that a high-quality frame posterior is biologically correct. It only decides when the candidate structure makes frame useful for assignment. Independent frame calibration and support-quality filters remain necessary.
- A global RDG rank test was not sufficient for this read-assignment problem. Some candidate graphs were globally identifiable even when a particular read's frame evidence only separated an alias group. The useful unit is the read-local candidate set plus a locus-level sufficiency check.

## Decision

This is a better formulation than unrestricted frame EM and better than a local gate alone. The production version should expose this as a consensus-gated frame-assignment mode and report the gate diagnostics alongside the read assignments. The reported fields should include how many reads were frame-separable, how many were frame-aliased, whether the locus-level gate passed, and how much posterior mass came from ordinary EM fallback.

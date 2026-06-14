# Locus-Shaped Truth Validation: read_assignment_ribometric_10pct

Notebook: `../notebooks/read_assignment_locus_shaped_truth_validation.ipynb`
Plot: `../notebooks/read_assignment_ribometric_10pct.locus_shaped_truth.validation.png`

## Formulation Fix

The previous real-locus stress test asked whether a method could recover an assumed transcript on real data. That is not a clean validation question because the assumed transcript may not be biologically true and the real frame-support track may already contain assignment artefacts.

This validation keeps real read-to-transcript ambiguity graphs, but controls the truth and the frame evidence. The question becomes: does adding an independent frame likelihood improve over ordinary EM only when the candidate origins are frame-separable?

## Readout

![Locus-shaped truth validation](../notebooks/read_assignment_ribometric_10pct.locus_shaped_truth.validation.png)

Summary table:

| validation | frame model | graph | reads | transcripts | fractional | EM | frame EM | EM gain | frame gain |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| alias_rich_same_gene | neutral | TGFBI:same_gene | 70 | 3 | 0.333 | 0.333 | 0.333 | 0.000 | -0.000 |
| alias_rich_same_gene | single_true_frame | TGFBI:same_gene | 70 | 3 | 0.333 | 0.333 | 0.333 | -0.000 | 0.000 |
| broad_mixed_all_origin | dual_frame_overlap | HRAS:all_origin | 192 | 261 | 0.140 | 1.000 | 1.000 | 0.860 | 0.000 |
| broad_mixed_all_origin | misleading_competitor_frame | HRAS:all_origin | 192 | 261 | 0.140 | 1.000 | 0.000 | 0.860 | -1.000 |
| broad_mixed_all_origin | neutral | HRAS:all_origin | 192 | 261 | 0.140 | 1.000 | 1.000 | 0.860 | -0.000 |
| broad_mixed_all_origin | single_true_frame | HRAS:all_origin | 192 | 261 | 0.140 | 1.000 | 1.000 | 0.860 | 0.000 |
| clean_rich_all_origin | dual_frame_overlap | ANGPTL4:all_origin | 157 | 7 | 0.208 | 0.333 | 0.452 | 0.125 | 0.119 |
| clean_rich_all_origin | misleading_competitor_frame | ANGPTL4:all_origin | 157 | 7 | 0.208 | 0.333 | 0.000 | 0.125 | -0.333 |
| clean_rich_all_origin | neutral | ANGPTL4:all_origin | 157 | 7 | 0.208 | 0.333 | 0.333 | 0.125 | -0.000 |
| clean_rich_all_origin | single_true_frame | ANGPTL4:all_origin | 157 | 7 | 0.208 | 0.333 | 1.000 | 0.125 | 0.667 |
| clean_rich_same_gene | dual_frame_overlap | ANGPTL4:same_gene | 157 | 7 | 0.208 | 0.333 | 0.452 | 0.125 | 0.119 |
| clean_rich_same_gene | misleading_competitor_frame | ANGPTL4:same_gene | 157 | 7 | 0.208 | 0.333 | 0.000 | 0.125 | -0.333 |
| clean_rich_same_gene | neutral | ANGPTL4:same_gene | 157 | 7 | 0.208 | 0.333 | 0.333 | 0.125 | -0.000 |
| clean_rich_same_gene | single_true_frame | ANGPTL4:same_gene | 157 | 7 | 0.208 | 0.333 | 1.000 | 0.125 | 0.667 |
| mixed_rich_same_gene | dual_frame_overlap | FN1:same_gene | 400 | 27 | 0.091 | 0.602 | 0.602 | 0.511 | 0.000 |
| mixed_rich_same_gene | misleading_competitor_frame | FN1:same_gene | 400 | 27 | 0.091 | 0.602 | 0.002 | 0.511 | -0.600 |
| mixed_rich_same_gene | neutral | FN1:same_gene | 400 | 27 | 0.091 | 0.602 | 0.602 | 0.511 | -0.000 |
| mixed_rich_same_gene | single_true_frame | FN1:same_gene | 400 | 27 | 0.091 | 0.602 | 0.594 | 0.511 | -0.008 |

## Interpretation

- In clean frame-informative graphs with a correct single-frame model, frame EM adds an average 0.667 posterior over ordinary EM.
- In neutral-frame controls, the largest absolute frame gain is 0.000; this is the expected negative control.
- With a deliberately misleading frame model, frame EM can lose as much as -1.000 posterior. This is the failure mode we must guard against in real data.
- Dual-frame overlap is not the same as bad periodicity. If two frames are genuinely supported, the model should avoid forcing all mass into one frame simply because one frame is dominant.
- The practical gate should therefore be conditional: use frame-aware assignment where frame support is independent, sufficiently deep, and actually separates candidate origins; otherwise report EM assignment plus uncertainty.

Selected real-shaped validation cases:

| validation | case | read category | selected because |
|---|---|---|---:|
| clean_rich_same_gene | ANGPTL4:em_top:same_gene | clean_frame_informative | 398 reads |
| mixed_rich_same_gene | FN1:em_top:same_gene | mixed_realistic | 2633 reads |
| alias_rich_same_gene | TGFBI:em_top:same_gene | same_frame_aliased | 188 reads |
| clean_rich_all_origin | ANGPTL4:em_top:all_origin | clean_frame_informative | 398 reads |
| broad_mixed_all_origin | HRAS:em_top:all_origin | mixed_realistic | 507 reads |

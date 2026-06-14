# Score-First Gate Progress

Branch: `local/read-assignment-prototype`

## Current decision state

The implementation now supports Gates 0-3 as executable local workflows. The frame-assignment comparison is executable for linear correction versus latent EM. Gate 4 now has a local prototype for read-origin assignment plus first real candidate-origin and sparse matched frame-support probes. It remains a decision-gated method rather than a committed core deliverable until production matched samples show enough independent frame evidence to matter.

## Gate 0: Unified score schema

Status: implemented.

Decision: score outputs are the stable product surface. Frame and assignment models add columns to the same ORF/translon score table rather than producing separate method-specific result formats.

Implemented:

- `TranslonScorer.pipeline.score_schema.ensure_score_schema`
- `TranslonScorer.pipeline.score_schema.add_frame_score_columns`
- `TranslonScorer.pipeline.score_schema.compare_score_tables`

Schema groups:

- raw score terms: `rise_up`, `step_down`, `hrf`, `avg`, `nzc`, `score`
- frame terms: `frame_posterior_mean`, `frame_entropy_mean`, `frame_weighted_count`, `frame_support_codons`, `frame_method`
- ambiguity terms: `assignment_weight`, `assignment_entropy`, `assignment_model`, `identifiability_class`
- provenance terms: `score_mode`, `input_type`, `offset_source`, `annotation_version`, `method_version`

Important correction made during Gate 0:

- CDS intervals from `getexons_and_cds()` are genomic, while frame support operates on transcript-space profiles. Added `cds_to_transcript_space()` and routed profile/frame-support paths through it.
- HMM smoothing now recomputes entropy after posterior smoothing.
- Existing ORF classification order was corrected so CDS-extending ORFs ending at the CDS stop are labelled `extORF` rather than `uoORF`.

## Gate 1: Raw scores vs frame-weighted scores

Status: executable, awaiting real panel run.

Implemented:

- CLI: `translonscorer score-compare-frame`
- CLI: `translonscorer score-compare-existing`
- CLI: `translonscorer validate-panel`
- Python: `TranslonScorer.pipeline.score_gates.compare_raw_frame_scoring`
- Python: `TranslonScorer.pipeline.panel_manifest.validate_panel_manifest`

Output:

- `<prefix>.raw_scores.parquet`
- `<prefix>.frame_scores.parquet`
- `<prefix>.comparison.parquet`
- `<prefix>.summary.csv`

Panel manifest:

- Recommended columns: `panel_id`, `category`, `label`, `tran_id`, `start`, `stop`, `type`, `gene_id`, `evidence`
- Required key: either `orf_id`, or `tran_id,start,stop`
- `validate-panel` reports duplicate keys and a SHA256 hash for freezing.

Decision rule:

- If frame-weighted scoring improves rank separation/calibration on the frozen panel, frame support becomes part of the default scoring evidence.
- If not, frame support remains a QC/browser track and not a scoring default.

## Frame-assignment correction comparison

Status: executable, awaiting real profile/CDS run.

Implemented:

- asymmetric cyclic confusion learning in `TranslonScorer.frame.bleed.learn_confusion`
- adjusted true-frame count output via `TranslonScorer.frame.bleed.apply_confusion_counts`
- observed and adjusted columns in every frame-support table: `observed_f0..2`, `adjusted_f0..2`, `total_count`
- support diagnostics in every frame-support table: `p_max`, `secondary_frame_mass`, `frame_periodicity_score`, `depth_score`, `support_evidence`, and `support_gated_p0..2`
- annotated-CDS validation in `TranslonScorer.pipeline.frame_method_compare.validate_frame_support_on_cds`
- CLI: `translonscorer compare-frame-methods`
- Python: `TranslonScorer.pipeline.frame_method_compare.compare_frame_methods`

Decision logic:

- `linear` is the direct correction: learn global/per-length leakage from CDS interiors and invert it per codon.
- `linear+hmm` is not a different correction model. It is linear correction followed by HMM smoothing along the transcript to reduce isolated codon noise.
- `latent` is frame-assignment EM: initialise from the learned leakage matrix, then fit count-space per-codon latent frame support with that leakage calibration held fixed by default. This makes EM a probabilistic refinement of the linear correction instead of an under-regularized re-learning of the leakage matrix from every codon.
- `p0/p1/p2` are conditional frame probabilities. They should not be read as a calibrated probability that a weakly supported position is translated.
- `secondary_frame_mass` is kept as an overlap/mixed-frame diagnostic. A drop in secondary mass after HMM smoothing is treated as a possible over-cleaning warning, not automatic improvement.

Output:

- `<prefix>.linear.frame_support.parquet`
- `<prefix>.latent.frame_support.parquet`
- `<prefix>.linear_vs_latent.frame_compare.parquet`
- `<prefix>.frame_compare.summary.csv`
- `<prefix>.frame_validation.summary.csv`
- optional `<prefix>.<method>.frame_validation.rows.parquet`

Decision rule:

- If latent EM materially improves annotated-CDS calibration or held-out frame recovery over linear correction, keep it as a supported correction mode.
- If latent EM mostly agrees with linear correction but adds instability, use linear or linear+hmm as the default and report latent as exploratory.
- Validation compares against the transcript frame implied by annotated CDS start (`cds_start % 3`), so transcripts whose CDS starts in frame 1 or 2 are handled correctly.

Report:

- Executable notebook: `notebooks/frame_assignment_comparison_report.ipynb`
- Contents: method definitions, synthetic benchmark, depth sweep, pairwise posterior deltas, runtime timing, transcript profile visual panels, and a guarded real-data subset demo for `ribocrypt_fwd_full`.
- Real-data plots now separate raw conditional frame posterior from support-gated evidence and show entropy plus secondary-frame mass so overlapping or weak translation is visible.
- Targeted real-data triage: `notebooks/frame_assignment_real_targeted_investigation.py` writes `notebooks/frame_assignment_real_targeted_cases.csv`, `.md`, and `.png`. After switching latent EM to count-space fixed-calibration fitting, the first global-top coding transcript pass produced 375 windows where HMM smoothing suppressed secondary-frame mass retained by both linear and latent, 119 windows with supported secondary-frame signal, 116 mixed/low-support windows, and 9 consistent single-frame windows. HMM should remain an auxiliary denoising view rather than the primary biological posterior.

## RDG-Flux v1 substrate export

Status: implemented and smoke-tested.

Implemented:

- CLI: `translonscorer export-rdg-flux`
- Python: `TranslonScorer.pipeline.rdg_flux_export.export_rdg_flux_v1`
- Python: `TranslonScorer.pipeline.rdg_flux_export.rdg_flux_position_table`

Contract:

- Input coordinates are the same transcript profile convention: `tran_id,pos`, 0-based transcript positions.
- Output uses RDG-Flux names: `sample_id, transcript_id, pos, count`.
- Frame evidence is exported as posterior mass: `p_frame0, p_frame1, p_frame2`.
- If `support_evidence` is present in frame support, `p_frame0..2` are scaled by it and the residual mass goes to `p_background`; otherwise export falls back to the old supported-codon behavior.
- `p_translated`, `frame_entropy`, `effective_depth`, and `local_periodicity_score` are included.
- Sidecar metadata JSON records coordinate system, offset model, annotation source, FASTA, TranslonScorer version, model stage, normalization, sample ids, and frame method.
- When `--frame-support` is supplied as Parquet and single-file output is requested, export uses a Polars lazy scan/sink path rather than materialising the full RDG table in Python.

Output:

- single-file mode: `<out>.parquet` and `<out>.metadata.json`
- partitioned mode: `<out_dir>/sample_id=<sample>/part.parquet` and `<out_dir>/metadata.json`

Decision:

- RDG-Flux should consume this per-position posterior table as its v1 emission layer.
- Stage 3 joint transcript-frame posterior export remains future work and should be added as a second sparse posterior layer rather than overloading this v1 table.

## Gate 2: BAM/read-level profiles reproduce BigWig/profile evidence

Status: executable, awaiting matched BAM/BigWig run.

Implemented:

- CLI: `translonscorer compare-profiles`
- Python: `TranslonScorer.pipeline.profile_compare.compare_profiles`

Output:

- `<prefix>.profile_compare.summary.csv`
- optional `<prefix>.profile_compare.deltas.parquet`

Decision rule:

- BAM-derived transcript profiles must agree with BigWig-derived profiles within a defined tolerance before read-level EM work begins.
- This gate validates read-level infrastructure and offset handling without claiming isoform deconvolution.

## Gate 3: Frame-disambiguation rate

Status: executable, awaiting real BAM/candidate run.

Implemented:

- CLI: `translonscorer frame-disambiguation`
- Python: `TranslonScorer.pipeline.frame_disambiguation.frame_disambiguation_from_bam`
- Python: `TranslonScorer.pipeline.frame_disambiguation.frame_disambiguation_from_candidates`

Output:

- `<prefix>.cds_transcript_space.parquet` for BAM mode
- `<prefix>.candidate_frames.parquet`
- `<prefix>.read_classes.parquet`
- `<prefix>.summary.csv`

Decision rule:

- If a meaningful fraction of ambiguous read assignments are frame-discordant, frame-aware EM becomes a core development branch.
- If not, Stage 3 is downgraded to frame-blind weighting plus identifiability reporting.

## Gate 4: Frame-aware EM

Status: prototype implemented, first real candidate-origin and sparse matched frame-support probes complete. RDG-gated frame assignment is now implemented in the reusable pipeline API; awaiting production matched sample run.

Implemented:

- CLI: `translonscorer compare-read-assignment`
- Python: `TranslonScorer.pipeline.read_assignment.assign_reads`
- Python: `TranslonScorer.pipeline.read_assignment.compare_read_assignment_methods`
- Python: `TranslonScorer.pipeline.read_assignment.frame_assignment_gate_diagnostics`
- Python: `TranslonScorer.pipeline.read_assignment.assignment_identifiability`
- frame-support quality gates: `min_frame_support_count` and `min_frame_support_evidence`, exposed in the CLI as `--min-frame-support-count` and `--min-frame-support-evidence`
- RDG-style frame assignment gates: `frame_gate_min_likelihood_range` and `frame_gate_min_read_fraction`, exposed in the CLI as `--frame-gate-min-likelihood-range` and `--frame-gate-min-read-fraction`
- tests: `tests/test_read_assignment.py`
- simulation: `notebooks/read_assignment_prototype_simulation.py`
- real candidate probe: `notebooks/read_assignment_real_candidate_probe.py`

Supported assignment methods:

- `unique`: assigns only reads with one candidate origin.
- `fractional`: normalizes alignment/candidate likelihood within each read.
- `frame_fractional`: fractional assignment with local frame-posterior compatibility.
- `em`: RSEM-style abundance EM over `abundance_key`.
- `frame_em`: abundance EM with frame compatibility.
- `rdg_local_frame_em`: EM where frame compatibility is used only for reads whose candidate set has one unique best frame-compatible assignment target.
- `rdg_gated_frame_em`: consensus-gated EM where local frame-compatible reads are used only when enough read-weighted evidence passes the local gate and that evidence exceeds alias-group frame contrast.

Assignment targets:

- `abundance_key=tran_id` for transcript isoform assignment.
- `abundance_key=locus_id` for genomic locus-origin assignment.
- `abundance_key=locus_id,tran_id` for a joint locus/transcript target.

Important correction:

- Frame likelihood uses the candidate P-site's transcript-coordinate signal frame (`pos % 3`). It does not use CDS-relative frame, because frame-support tables store posterior mass by transcript-coordinate residue class.
- Low `support_evidence` makes the frame-compatibility term neutral rather than overconfident, so weak frame calls do not drive read assignment.
- Candidate-origin rows are collapsed to `(read_key, assignment_target)` before target-level posterior inference. This prevents annotation multiplicity from inflating locus-level assignment when many transcript models share one locus. The target posterior is then distributed back over candidate-origin rows.
- RDG-style gating treats each read's candidate origins like a small path set. Frame evidence can move assignment only when it separates one assignment target, not merely a group of frame-compatible aliases. The consensus gate then prevents a tiny number of locally frame-separable reads from pulling the abundance EM for the whole locus.

Decision: no read-assignment EM implementation should be treated as core until Gate 3 is run on real candidate tables. The prototype exists so the model assumptions can be tested quickly and falsified early.

Prototype simulation result:

- `frame_em` beats frame-blind EM in a frame-discordant isoform scenario.
- `em` beats fractional assignment when unique reads anchor isoform or locus abundance.
- `frame_em` matches `em` when frame evidence is uninformative.
- deliberately aliased reads remain at chance, which is the desired identifiability behaviour.
- `rdg_gated_frame_em` keeps the clean frame-disambiguation gain, matches ordinary EM for same-frame aliases and mixed alias groups, and blocks sparse misleading frame support in focused tests.

Real candidate-origin probe:

- Source: `/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/subsampled_1_percent.bam`
- 35,804 read keys and 68,417 weighted reads.
- 88.50% of read keys are transcript-ambiguous; 20.78% are locus-ambiguous.
- Transcript-level EM reduces normalized assignment entropy from 0.846 to 0.384.
- Locus-level EM reduces normalized assignment entropy from 0.263 to 0.084.
- Without frame support, `frame_em` exactly matches `em`, confirming neutral frame handling.
- Weighted transcript-level EM classes: 15.38% unique, 27.73% resolved, 24.15% ambiguous, 32.75% data-limited.
- Weighted locus-level EM classes: 73.69% unique, 15.66% resolved, 6.13% ambiguous, 4.52% data-limited.

Matched frame-support probe:

- `unique_transcript` independent support: 3,933 frame rows across 1,623 transcripts, total count 10,521.0, matched support for 9.8% of weighted candidate reads.
- `fractional_all` exploratory support: 268,625 frame rows across 71,156 transcripts, total count 68,417.0, matched support for 24.7% of weighted candidate reads.
- With independent support, `frame_em` changes transcript assignment modestly relative to `em`: weighted mean target-posterior TV distance 0.0018; 0.65% of weighted reads move by at least 0.10.
- Applying a minimal quality gate (`total_count >= 10`, `support_evidence >= 0.2`) changes the independent-support shift only slightly: weighted mean TV distance 0.0016; 0.65% of weighted reads still move by at least 0.10.
- With exploratory all-fractional support, transcript-level weighted mean TV distance rises to 0.0318; 6.17% of weighted reads move by at least 0.10. This is sensitivity analysis, not independent validation.
- With the same gate, exploratory all-fractional support shifts weighted mean TV distance 0.0281; 5.76% of weighted reads move by at least 0.10.
- Independent high-shift cases concentrate around `LINC01783`: 57 read keys, 443.0 weighted reads, mean TV 0.205, and roughly 7 compatible transcripts per read.
- LINC01783 case review found independent frame support for only one transcript (`ENST00000415386.2`), with 29 unique-support counts across 3 codons.

Decision from the sparse matched probe:

- The model behavior is coherent: `frame_em` matches `em` when frame support is absent and shifts posterior mass only where matched support exists.
- The current real substrate is too sparse and too concentrated to justify a broad claim that frame resolves isoform origin genome-wide.
- The next Gate 4 run should use a production matched sample and calibrate the frame-support quality thresholds before transcript-specific shifts are interpreted biologically.

Scaled 10% BAM probe:

- Source: `/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/subsampled_10_percent.bam`
- 358,510 read keys and 828,006 weighted reads.
- 88.53% of read keys are transcript-ambiguous; 21.25% are locus-ambiguous.
- Independent frame support: 33,269 frame rows across 5,540 transcripts; matched support for 12.8% of weighted candidate reads.
- With independent support, `frame_em` shifts transcript assignment relative to `em`: weighted mean TV distance 0.0257; 3.00% of weighted reads move by at least 0.10.
- With the minimal quality gate (`total_count >= 10`, `support_evidence >= 0.2`), weighted mean TV distance is 0.0231; 2.93% of weighted reads move by at least 0.10.
- The high-shift independent mass remains concentrated: `ENSG00000228549` accounts for 97.3% of weighted high-shift reads in the exported high-shift locus summary.
- Targeted review classifies that cluster as `ENSG00000228549` / `LINC01783` lncRNA ambiguity, not clean translated isoform validation.
- Coding-associated scan: 111 high-shift read keys, 1,989 weighted reads. Of these, 63 read keys and 1,254 weighted reads have gated local frame support; the remaining coding-associated shifts are abundance-propagated through EM.
- The top coding-associated frame-top label is `CIC`, but the reviewed reads span 46 candidate gene labels, average 28.5 compatible transcript rows per read, and only 1 weighted read has gated local frame support.
- The lower-locus `PTCH1` review shows the opposite failure mode: one gene, high TV shift, but no matched local frame-support rows, so the shift is abundance-propagated rather than directly frame-supported.
- A seed coding review panel is written at `docs/read_assignment_ribometric_10pct.coding_locus_panel_seed.md`; it prioritizes direct local-frame-supported candidates over abundance-propagated shifts.
- Targeted reviews of `POGK` and `OAZ2` show direct local frame support but cross-gene ambiguity (`POGK`/`FLYWCH1`, `OAZ2`/`IL17D`), so they remain triage candidates rather than validation wins.
- Additional targeted reviews are now written for `ESRRAP2`, `OLFM1`, `CDV3`, `TMEM271`, `DSCAM`, `IL17D`, `H1-2`, `NDUFS5`, and `RCE1`.
- The executed notebook `notebooks/read_assignment_failure_case_walkthrough.ipynb` groups the failures into noncoding annotation-block ambiguity, broad multi-locus ambiguity, pseudogene/parent-gene ambiguity, same-gene ambiguity without local frame support, paralog ambiguity, underpowered same-gene candidates, and direct-frame but cross-gene reads.
- The input-assumption audit `docs/read_assignment_ribometric_10pct.input_assumption_audit.md` fails this BAM as a biological frame-validation substrate: best offset from alignment start is 13, but best annotated-CDS frame fraction is only 0.368, below the 0.55 gate.
- The executed notebook `notebooks/read_assignment_profile_disambiguation_walkthrough.ipynb` now shows transcript profiles for selected ambiguous reads before and after `raw_compatible`, `unique_only`, `fractional`, `em`, and gated `frame_em` disambiguation.
- Coverage routes and local high-depth paths are documented in `docs/read_assignment_coverage_paths.md`. The key distinction is that deep P-site profiles already exist locally, but a deep read-level BAM with valid P-site candidate coordinates is still missing.

Decision from the scaled probe:

- Frame-aware assignment has measurable real-data effect at larger local depth.
- The observed effect is dominated by annotation/mappability ambiguity, and the input BAM does not meet the frame-signal assumption. This strengthens the identifiability-reporting story and weakens any broad claim of immediate genome-wide isoform deconvolution.
- The next validation gate should manually curate the direct local-frame-supported coding seed panel and then rerun on a true production matched BAM/profile pair.

RDG-gated assignment experiment:

- Notebook: `notebooks/read_assignment_rdg_gated_assignment_experiment.ipynb`
- Report: `docs/read_assignment_ribometric_10pct.rdg_gated_assignment_experiment.md`
- Clean separable, correct-frame cases: unrestricted `frame_em` and `rdg_gated_frame_em` both add mean true-transcript posterior gain of 0.667 over ordinary EM.
- Mixed candidates with correct frame: unrestricted `frame_em` gives mean gain -0.004; `rdg_gated_frame_em` falls back to EM with gain 0.000.
- Same-frame aliases: all methods remain unresolved, as expected.
- Dual-frame overlap: unrestricted `frame_em` gives apparent mean gain 0.059, but the gate treats the signal as alias-group frame contrast rather than transcript-origin evidence and falls back to EM.
- Misleading frame evidence: unrestricted `frame_em` loses mean true-transcript posterior -0.567, the local-only RDG gate still loses -0.040, and the consensus RDG gate blocks the failure with gain 0.000.

Decision from the RDG-gated experiment:

- `frame_em` should remain an exploratory method.
- `rdg_gated_frame_em` is the safer production default for frame-aware read assignment because it only uses frame when candidate structure makes frame informative for origin assignment.
- Gate diagnostics must be exported with comparisons: the important values are the local frame-separable read fraction, alias-group frame contrast fraction, and whether the consensus gate passed.

## Verification so far

Focused unit tests:

```bash
python3 -m pytest tests/test_read_assignment.py tests/test_score_first_gates.py -q
```

Result: 23 passed.

Latest result after adding `rdg_gated_frame_em`: 31 passed.

Full-suite status:

```bash
python3 -m pytest -q
```

Result: collection still fails on pre-existing environment/fixture issues: missing `data/annotationsubset.gtf`, missing `orfipy_core`, and missing `oxbow`.

Actionable subset:

```bash
python3 -m pytest tests/test_coordinates_unit.py tests/test_orffinder_unit.py tests/test_score_first_gates.py -q
```

Result: 14 passed.

Compile check:

```bash
python3 -m compileall -q TranslonScorer tests notebooks/read_assignment_prototype_simulation.py notebooks/read_assignment_story_evidence.py notebooks/read_assignment_real_candidate_probe.py notebooks/read_assignment_matched_frame_probe.py notebooks/read_assignment_locus_case_review.py
```

Result: passed.

CLI smoke check:

```bash
python3 -m TranslonScorer.cli --help
```

Result: gate commands are registered.

RDG-Flux export smoke:

```bash
python3 -m TranslonScorer.cli export-rdg-flux --profiles <profiles.parquet> --cds <cds.parquet> --sample-id ribocrypt_fwd_full --out <rdg.parquet> --frame-method linear --no-frame-by-length
```

Result: passed; wrote a 12-column RDG-Flux position table and sidecar metadata.

Precomputed frame-support lazy smoke:

```bash
python3 -m TranslonScorer.cli export-rdg-flux --profiles <profiles.parquet> --frame-support <frame_support.parquet> --sample-id ribocrypt_fwd_full --out <rdg.parquet>
```

Result: passed; metadata records `lazy_precomputed_export: true`.

Immediate RDG test inputs checked:

- `/Users/jackt/projects/all-RiboSeq/ribocrypt_fwd/ribocrypt_fwd_transcript_profiles.parquet` exists.
- Profile schema is `tran_id,pos,count`.
- Profile size is 157,586,237 rows with total count approximately `9.1204e10`.
- GENCODE v45 GTF and hg38 FASTA paths exist.

Decision: do not run the full 157M-row export interactively until the curated locus subset or precomputed frame-support parquet is available; the command is documented in `docs/rdg_flux_export.md`.

Frame method comparison smoke:

```bash
python3 -m TranslonScorer.cli compare-frame-methods --profiles <profiles.parquet> --cds <cds.parquet> --out-prefix <prefix> --methods linear,latent --no-frame-by-length
```

Result: passed; wrote linear support, latent support, comparison, and summary outputs.

Read assignment prototype smoke:

```bash
python3 notebooks/read_assignment_prototype_simulation.py
python3 notebooks/read_assignment_real_candidate_probe.py
python3 notebooks/read_assignment_matched_frame_probe.py
python3 notebooks/read_assignment_locus_case_review.py
python3 -m TranslonScorer.cli compare-read-assignment --candidates <candidates.parquet> --frame-support <frame_support.parquet> --out-prefix <prefix> --methods unique,fractional,em,frame_em --abundance-key tran_id
```

Result: passed; simulation summary, visual reports, real candidate-origin probe, matched frame-support probe, LINC01783 case review, CLI summary, convergence, abundance, and identifiability outputs were written.

Gated read-assignment CLI smoke:

```bash
python3 -m TranslonScorer.cli compare-read-assignment --candidates notebooks/read_assignment_real_candidate_probe.candidates.parquet --frame-support notebooks/read_assignment_matched_frame.unique_transcript_frame_support.parquet --out-prefix /tmp/read_assignment_cli_gated_smoke --methods em,frame_em --abundance-key tran_id --max-iter 20 --min-frame-support-count 10 --min-frame-support-evidence 0.2
```

Result: passed; summary, convergence, abundance, and identifiability outputs were written. The smoke intentionally used `--max-iter 20`, so it checks CLI plumbing rather than final EM convergence.

Full legacy suite:

```bash
python3 -m pytest -q
```

Result: blocked during collection by pre-existing local test environment gaps:

- `data/annotationsubset.gtf` is missing.
- `orfipy_core` is not installed.
- `oxbow` is not installed in the `python3` environment used for the run.

# Coverage Paths For Read-Assignment Validation

## Decision

We need higher coverage, but more specifically we need higher coverage that preserves the right coordinate semantics.

The current RiboMetric transcriptome BAM is not enough for biological validation of frame-aware read assignment because it fails the P-site/frame-signal audit. It is still useful for exercising candidate-origin construction, EM behaviour, and ambiguity reporting.

## Local Data Paths

### Deep P-site profiles

These are the best local resources for frame-profile inspection and RDG-Flux-style substrate work:

- `/Users/jackt/projects/all-RiboSeq/ribocrypt_fwd/ribocrypt_fwd_transcript_profiles.parquet`
  - 157,586,237 rows.
  - Total count about `9.12e10`.
  - 249,066 transcripts.
- `/Users/jackt/projects/all-RiboSeq/ribocrypt_rev/ribocrypt_rev_transcript_profiles.parquet`
  - 152,619,699 rows.
  - Total count about `9.56e10`.
  - 249,021 transcripts.
- `/Users/jackt/projects/all-RiboSeq/bench/hg38/ribocrypt_forward.bw`
- `/Users/jackt/projects/all-RiboSeq/bench/hg38/ribocrypt_reverse.bw`

Use these for:

- high-coverage frame inference;
- profile visualisation;
- Stage 1 validation;
- RDG-Flux emission export;
- selecting high-depth loci to revisit at read level.

Do not use these alone for:

- read-origin validation;
- transcript-isoform read assignment;
- genomic-locus read assignment.

They are already aggregated profiles, so read-level ambiguity is gone.

### Curated high-coverage benchmark subsets

These are smaller, local, profile-level subsets that are useful for fast iteration:

- `/Users/jackt/projects/all-RiboSeq/bench/hg38/full/ribocrypt_fwd/ribocrypt_fwd_transcript_profiles.parquet`
  - 990,776 rows.
  - Total count about `1.27e8`.
  - 238 transcripts.
- `/Users/jackt/projects/all-RiboSeq/bench/hg38/full/ribocrypt_rev/ribocrypt_rev_transcript_profiles.parquet`
- `/Users/jackt/projects/all-RiboSeq/bench/hg38/mid/ribocrypt_fwd/ribocrypt_fwd_transcript_profiles.parquet`
- `/Users/jackt/projects/all-RiboSeq/bench/hg38/small/ribocrypt_fwd/ribocrypt_fwd_transcript_profiles.parquet`

Use these for:

- fast plotting;
- method comparison;
- selecting candidate loci before running expensive read-level processing.

### Current read-level BAMs

The only non-empty local read-level BAMs suitable for the current prototype are:

- `/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/subsampled_10_percent.bam`
- `/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/subsampled_5_percent.bam`
- `/Users/jackt/projects/all-RiboSeq/RiboMetric/sample_data/subsampled_1_percent.bam`

The 10 percent BAM failed the P-site/frame-signal audit:

- best offset from alignment start: 13;
- best annotated-CDS frame fraction: 0.368;
- offset-0 frame fraction: 0.365;
- current gate: 0.55.

Use these for:

- candidate-origin table construction;
- assignment algorithm plumbing;
- ambiguity/failure-mode plotting.

Do not use these for:

- positive biological claims about frame-aware read assignment.

The Nextflow `work/` BAMs currently found under `/Users/jackt/projects/all-RiboSeq/work` are zero-byte outputs in this checkout, so they are not usable as a coverage route.

## Practical Routes To More Coverage

### Route 1: Profile-first locus selection, then read-level rerun

This is the best next step.

1. Use the deep profile parquet or BigWigs to identify high-depth curated loci.
2. Require strong annotated-CDS periodicity in the profile data.
3. Build a short list of loci where alternative transcripts differ in translated structure.
4. Reprocess the underlying reads for those loci with multi-mappers preserved and length-specific P-site offsets recorded.
5. Run read assignment only on those gene groups first.

This avoids spending compute on low-depth or intrinsically uninformative loci.

### Route 2: Production read-level run with P-site offsets

This is the necessary route for the real Stage 3 validation.

Required output:

- read-level candidate origins;
- multi-mapping candidates retained;
- transcript coordinate for each candidate;
- length-specific P-site offset applied or stored;
- `NH`/alignment multiplicity or equivalent mapping likelihood;
- annotation and transcriptome version recorded.

Acceptance gate before assignment:

- annotated-CDS frame concentration should pass the audit, provisionally `max_frame_fraction >= 0.55`;
- frame concentration should be length-specific, not only global;
- enough reads should remain in the curated loci after filtering.

### Route 3: Synthetic high-depth mixtures

This should run in parallel with the real-data route.

Use synthetic read-origin tables where truth is known, varying:

- depth;
- isoform similarity;
- frame-discordant versus frame-preserving isoforms;
- multi-locus ambiguity;
- P-site offset error;
- overlapping translation.

This is where we can prove that EM, frame-aware EM, and linear/fractional strategies behave correctly before relying on messy real data.

### Route 4: Aggregate profiles for frame support only

For an interim hybrid test, use the deep P-site profiles to estimate frame support, but use read-level candidate tables only for loci where reads are available.

This is not a final validation because the frame evidence and read evidence may come from different processing paths. It is useful as a sensitivity check.

## Plotting Standard Going Forward

Ribo-seq profile plots should use:

- vertical bars, not connected lines or scatter;
- bars coloured by transcript-coordinate phase (`pos % 3`);
- y-axis starting at zero;
- shared y-scale across methods for the same transcript;
- transcript position on x-axis;
- compact track-like axes with minimal grid;
- explicit caveat when coordinates are raw read starts rather than P-sites.

The reusable style helper is:

- `/Users/jackt/projects/all-RiboSeq/translonscorer/TranslonScorer/visualization/riboseq_profile_style.py`

The profile-disambiguation walkthrough now uses this style.

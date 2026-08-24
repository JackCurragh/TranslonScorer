# Matrix normalisation strategy for locus clustering and event scoring

Status: design note (2026-06-12). This document defines the normalisation
contract for matrix-based translon scoring. It is intended to be implementable
in parallel with event-level scoring work.

## Goal

The matrix engine gives positional Ribo-seq profiles for thousands of bulk
samples at each candidate locus. The goal is to cluster samples by the *shape*
of their local ribosome profile, then score each cluster aggregate. This should
separate alternative translon configurations used in different cellular
programs without letting expression level or sequencing depth dominate the
clustering.

Normalisation therefore has two distinct jobs:

1. Build a shape-normalised matrix for clustering samples at a locus.
2. Preserve count-scale signal for cluster aggregation and scoring.

The same transformed matrix must not be used for both jobs.

## Core principle

Maintain two matrices for every locus or genomic event window:

| matrix | purpose | count scale |
|---|---|---|
| `X_score` | aggregation and event scoring | raw or sample-depth-normalised counts |
| `X_cluster` | sample clustering only | within-locus shape-normalised values |

`X_cluster` is allowed to remove magnitude because clustering should ask
"where are reads distributed?". `X_score` must retain interpretable abundance
because scoring should ask "how much evidence supports this event?".

## Inputs

For a candidate locus or event window, start from a dense or sparse profile:

```text
X_raw[sample, position] = A-site count
```

Required metadata:

- sample identifier
- total usable RPF depth, or another matrix-level size factor
- raw locus count total per sample
- genomic/transcript position vector
- optional sample/batch metadata for diagnostics

The locus profile can be transcript-space or genomic-space, but the coordinate
system must be explicit and stable. Event-level scoring should prefer genomic
events as the canonical storage layer, then compose feature/translon reports
from those events.

## Coverage gate

Coverage filtering is performed before shape normalisation and uses raw or
depth-normalised counts, never row-normalised values.

Recommended fields:

```text
locus_total_raw[sample] = sum_positions X_raw[sample, position]
locus_total_depth_norm[sample] = sum_positions X_depth_norm[sample, position]
```

Recommended initial policy:

- exclude samples with `locus_total_raw < min_raw_locus_counts`
- default threshold should be conservative and configurable, e.g. 20-100 A-site
  counts depending on window size and score type
- retain excluded samples in the label table with `cluster_id = -1` and an
  exclusion reason

This prevents one- or two-read profiles from being stretched into apparently
strong shapes by row normalisation.

## Sample-depth correction

Apply sample-level depth correction before deriving either clustering or scoring
matrices when samples have different usable RPF depth:

```text
X_depth_norm[s, p] = X_raw[s, p] / size_factor[s]
```

The size factor should be based on usable Ribo-seq signal, not total file reads
when possible. Candidate choices, in order of preference:

1. matrix-wide CDS-assigned A-site counts passing QC
2. high-confidence coding-region RPF counts
3. total usable A-site counts after read-length and mapping filters
4. library/read depth only as a fallback

Depth normalisation must be global per sample, not estimated from the candidate
locus itself. Locus totals are biological signal and should not define the
sample size factor.

## Shape normalisation for clustering

The default clustering transform should remove expression magnitude while
preserving positional structure:

```text
X0 = X_depth_norm[covered_samples, :]
X1 = log1p(X0 / profile_scale)
X_cluster = row_shape_normalise(X1)
```

Where `profile_scale` is optional and can be `1.0` initially. The first
implementation should support at least these row transforms:

### `l1_log`

```text
u[s, p] = X_depth_norm[s, p] / sum_p X_depth_norm[s, p]
X_cluster[s, p] = log(u[s, p] + pseudocount)
```

Then row-center or L2-normalise before cosine clustering.

### `shifted_clr`

Single-cell depth-normalisation work motivates a shifted centered log-ratio
transform for count compositions:

```text
u[s, p] = X_depth_norm[s, p] / sum_p X_depth_norm[s, p]
z[s, p] = log(u[s, p] + c) - mean_p log(u[s, p] + c)
X_cluster = z
```

This is the preferred default to test because it explicitly removes row-depth
geometry after the log step. The pseudocount `c` must be configurable. Use a
stable default, record it in outputs, and benchmark sensitivity.

### `l2_log`

```text
X1 = log1p(X_depth_norm)
X_cluster[s, :] = X1[s, :] / ||X1[s, :]||_2
```

This is simple, works well with cosine distance, and is useful as a baseline.

### Methods to avoid as defaults

- raw Euclidean distance: clusters mostly by locus count
- min-max scaling: over-amplifies low-count spike profiles
- z-score alone: zeros and sparse spikes can become artificial structure

These can remain diagnostic options, but should not be the default.

## Distance and clustering

Default clustering should use cosine or correlation-like geometry on
`X_cluster`, not Euclidean distance on raw counts.

Initial supported methods:

- cosine k-means for fixed `k`
- agglomerative clustering with cosine distance
- optional HDBSCAN/DBSCAN only after coverage and transform behaviour is stable

The number of clusters should be conservative. For implementation, expose
`n_clusters` but record enough diagnostics to support later model selection.

## Cluster aggregation for scoring

Cluster aggregates are built from `X_score`, not `X_cluster`.

Recommended default:

```text
cluster_profile[c, p] = sum_{s in cluster c} X_depth_norm[s, p]
```

Alternative supported summaries:

- `sum_depth_norm`: default; preserves evidence scale while correcting depth
- `sum_raw`: useful if upstream counts are already comparable
- `mean_depth_norm`: diagnostic only unless scoring is proven scale-invariant

The aggregation output must include:

- cluster id
- sample count in cluster
- total raw counts
- total depth-normalised counts
- excluded sample count
- normalisation method and parameters

## Relationship to event-level scoring

The event-level scoring plan is compatible with this normalisation strategy.
The intended division of responsibility is:

1. The normalisation layer prepares comparable sample profiles and cluster
   aggregates.
2. The scoring layer stores and scores canonical genomic events once.
3. Feature/translon reports compose event scores without rescoring shared
   evidence.

The event database should treat a genomic event as the canonical evidence unit,
for example:

```text
event_id
event_type
genome interval / strand
position vector or window definition
cluster_id / aggregate_id
normalisation_id
score fields
```

Feature-level translons should then reference one or more event IDs. This avoids
scoring the same genomic evidence repeatedly for overlapping translons or
shared exons, and it prevents contradictory scores for the same underlying
elongation or initiation evidence.

No annotation decision is made at this layer. The output is evidence:

```text
genomic event scores -> composed translon-level reports
```

Classification or annotation policy can be added later.

## Required outputs

Every clustering/scoring run should write enough metadata to make scores
reproducible:

```text
normalisation_id
input matrix id / manifest id
coordinate system
window definition
sample size-factor source
coverage gate threshold
shape transform
pseudocount
row centering / L2 flag
distance metric
cluster method
cluster parameters
aggregation method
```

Label table:

```text
sample_id
cluster_id
included
exclusion_reason
locus_total_raw
locus_total_depth_norm
size_factor
```

Cluster summary:

```text
cluster_id
n_samples
raw_count_total
depth_norm_count_total
mean_locus_total_raw
median_locus_total_raw
within_cluster_cosine_similarity
```

## Diagnostics

Normalisation is acceptable only if it removes depth-driven geometry without
destroying event signal.

Minimum diagnostics per locus/event window:

- correlation between raw locus totals and first two clustering PCs
- correlation between raw locus totals and cluster labels
- per-cluster raw and depth-normalised count distributions
- within-cluster and between-cluster cosine distances
- cluster stability after downsampling high-depth samples
- cluster stability after modest window padding/trimming
- fraction of samples excluded for low coverage

Failure modes to flag:

- clusters ordered almost entirely by `locus_total_raw`
- clusters dominated by a single study/batch after depth correction
- many low-count samples assigned to confident clusters
- cluster aggregate scores driven by one or two high-depth samples

## Implementation sketch

The normalisation API can be independent of the event scorer:

```python
normalise_locus_profiles(
    matrix,
    sample_ids,
    pos_vector,
    size_factors,
    *,
    min_raw_locus_counts=50,
    score_count_scale="depth_norm",
    cluster_transform="shifted_clr",
    pseudocount=1e-6,
    row_l2=False,
)
```

Return:

```text
X_score
X_cluster
included_mask
normalisation_metadata
sample_qc_table
```

Then:

```python
cluster_locus_profiles(X_cluster, included_sample_ids, ...)
aggregate_cluster_profiles(X_score, labels, method="sum_depth_norm")
```

The existing `profile_clustering.normalise_profiles` can be extended rather
than replaced, but it should distinguish coverage filtering, sample-depth
correction, shape transformation, and scoring aggregation explicitly.

## Open parameters to benchmark

- default coverage threshold
- size-factor source
- shifted-CLR pseudocount
- `shifted_clr` versus `l2_log`
- sum versus mean cluster aggregation for each score type
- fixed `k` versus adaptive cluster selection
- whether to compute transforms on the full event window or only informative
  positions

These are tuning choices. The architectural requirement is fixed: cluster on a
shape-normalised matrix and score on raw/depth-normalised cluster aggregates.

# Annotation decisions

Decisions taken, with the reason for each. Scoring measures; this document is
where annotation choices are recorded. Where a decision was made on pancreas data
it is marked, because scope is part of the decision.

## Unit of annotation

**Decision:** the annotation unit is the locus, keyed by
`chromosome + terminal-stop interval + strand`. Each distinct start is kept as an
alternative record under that locus. The annotation carries both `locus_id` and
`start_variant_id`.

**Reason:** tools frequently disagree about starts while agreeing about the stop.
Different start selection is not grounds for creating separate biological loci.

**Correction on record:** a valid continuous ORF cannot have a different frame at
the same stop. Alternative starts are retained and audited for frame consistency.

## Scoring unit

**Decision:** all scores are per event or per block. Translons are not rescored
as whole units.

**Reason:** events are shared across features, so scoring them once is what makes
annotation-scale work tractable. Rescoring per translon discards that and
recomputes the same signal many times.

**Consequence:** belief in a translon is an assessment over the set of blocks it
claims, derived from per-block evidence.

## Missing evidence

**Decision:** `NA` is missing evidence. Never zero, never a negative score.

**Decision:** missing mappability is an unknown, not evidence of low mappability.
A call with no mappability information does not enter the high-confidence tier by
default. **Enforced in code**: `_map_track_for_chrom` distinguishes a confirmed
0.0 from a true gap in the mappability bigwig and returns
`map_track_mean=None`/`map_track_low=None` when nothing is known, rather than
defaulting an unknown position into the low-mappability band.

**Reason:** a data gap and a confirmed-bad measurement are different claims, and
collapsing them silently converts absence into evidence against.

## Profile evidence

**Decision:** profile support is evidence, not truth.

**Reason:** profile-derived evidence is part of the scoring workflow and is not
independent validation.

**Decision:** observability — mappability and profile depth — modifies confidence.
It does not define biological class.

## Tool combination

**Decision:** no unweighted majority vote.

**Reason:** the tools have different purposes, thresholds, and available metrics,
and their call counts differ by an order of magnitude (RibORF2 produces far more
calls than the others). A vote weights that imbalance rather than the evidence.

**Decision:** native tool scores are compared only within a tool, by within-tool
percentile. Raw scores are not compared between tools.

## Spliced features (pancreas-scoped)

**Decision:** policy development starts with non-spliced loci.

**Reason:** this removes transcript/genomic projection complexity while basic tool
and profile behaviour is being established.

**Decision:** fully spliced and mixed-splicing calls are reported, but treated as
a separate interpretation problem.

**Reason:** a genomic profile and a transcript-projected ORF are not equivalent
views at a splice junction.

## Data sources (pancreas-scoped)

**Decision:** prefer pooled pancreas BED12 annotations over earlier local or
JBrowse-only representations where splicing matters.

**Reason:** the pooled source preserves the block structure needed to interpret
spliced calls.

**Decision:** use the Hoffman/Umap k24 mappability track.

**Reason:** the earlier sparse track had missing or partial values that could be
read as zero mappability.

## Hard exclusions

A call is rejected from the release annotation — and retained in a
rejected/candidate table — if:

1. It is not frame-consistent from start to terminal stop, or its blocks cannot be
   projected coherently onto a transcript.
2. It has no usable transcript/exon mapping and no independent genomic profile
   evidence strong enough to support a discovery candidate.
3. The pooled profile has zero total signal, zero usable events, or no observable
   signal in any expected initiation, elongation, or termination region.
4. Mappability is missing over the feature and there is no independent support.
5. A required score is unavailable or not applicable.

These are evidence and representation failures. None of them is a claim that the
underlying biology is impossible.

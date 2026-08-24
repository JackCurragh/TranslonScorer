# Locus-level attribution: interpretability and scale — exploration notes

Not a plan, not yet implemented. Captures a design discussion so it survives
past one conversation, in the same spirit `significance_testing_plan.md` was
written before anything in it existed. Nothing here should be treated as
decided.

## The problem

A score currently reads as if it belongs entirely to its named candidate
event. `competitor_share`/attribution machinery computes with awareness of
neighbors, but the neighbor information itself doesn't survive past the
computation — it collapses into a confidence number or an outcome label
(`clean_drop`, `multi_way_overlap`, etc.) and disappears. A caller looking at
one event's score today has no way to ask "what else nearby could plausibly
have contributed to this," short of re-deriving it themselves.

That's the actual requirement: not just a faster way to detect overlap, but
a retained, queryable record — for any scored event, the set of other
candidate events within its window, their spans, their frames, their
distance, and their own battery status — surfaced in `report.py`'s output,
not folded away inside an aggregate.

## Performance: don't recompute per-locus

Whatever the interpretability layer looks like, it must be computed once per
some natural grouping unit and joined in at scoring time — never
recomputed on demand per locus. `events.py`'s existing `_contention`
function (which produces the `event_overlap` table) already does exactly
this shape of thing for exact same-frame overlaps: sort candidates by start,
one linear sweep, O(n log n). The gap flagged in
`docs/significance_testing_results.md` is that this table is written to
disk during `extract_events_workflow` and never read back in by
`score_bams_workflow`/`score_matrix_workflow` — so the fix for "make this
fast" and the fix for "make this happen at all in the standard CLI path"
are the same fix. Extending `_contention` to add a proximity slop (not just
exact overlap) is a small, well-scoped addition to something that already
works, independent of the larger question below.

## The open question: what's the right grouping unit?

Raised directly: is per-chromosome the right boundary, or per-gene — all
transcripts of one gene together (a splice graph), with an RDG layered on
top representing decision points across all isoforms, as "the whole
translatable context" of a gene?

**Per-chromosome** (what `_contention` does today): computationally clean,
a single sweep, no boundary-effects at gene edges. But the output is a flat
table you query into, not itself a "context" object — "what's nearby this
locus" means slicing a chromosome-wide table, not walking a self-contained
per-gene structure.

**Per-gene, splice graph + RDG**: closer to what "the whole translatable
context" actually means intuitively — one graph object per gene that
already encodes exon structure, per-path frame, and decision points, so
containment/independence/overlap become properties readable off that graph
rather than externally computed facts bolted on afterward. This is
appealing precisely because it's what RDG's stated design goal already is
(see RDG's README: "Ribosome Decision Graph... to represent RNA transcripts
that encode multiple proteoforms").

**But gene-level doesn't fully dissolve the boundary problem — it just
moves it up one level.** Genes themselves overlap: antisense pairs, nested
genes, readthrough transcripts spanning what's annotated as two genes. The
GAPDH validation locus itself sits among several neighboring genes in the
same GTF window. A "per-gene" object still needs a way to reach across into
a neighboring gene's object when a candidate event's window crosses that
boundary — the same problem restated at a coarser grain, not solved.

**A tension with existing, considered design work, worth reconciling
explicitly rather than picking gene-level because it feels more natural:**
the L0b/L1 substrate architecture (see project memory
`project_l0b_builder`/`project_l1_builder`/`project_reannotation_engine_design`)
already made a deliberate decision that the cache/scoring boundary is the
**junction set**, not the transcript or gene — specifically to avoid
transcript-centric assumptions and to score genomic events once regardless
of how many transcripts they're annotated under. A gene-level
splice-graph+RDG object would reintroduce exactly the transcript/gene-centric
boundary that decision moved away from. Before building toward "gene =
splice graph + RDG," this needs to be reconciled with why that boundary was
rejected before — either the junction-set decision was scoped narrower than
this new question (plausible — it was about caching scored *events*, not
about an interpretability/context layer) or the two genuinely conflict and
one has to give. Not resolved here.

## What RDG actually offers today, concretely (from direct investigation)

- `RDG/index.py` has a real interval tree (`PointIntervalNode`/
  `PointIntervalIndex`, O(log n + k) queries) but it's **point-containment
  only** — no range/overlap query, no proximity/slop. Not directly reusable
  for "what overlaps or is near this candidate."
- `RDGDocument.enumerate_paths()`/`EnumeratedPath` (schema.py) exposes
  per-candidate interval + frame cleanly, and `CandidateORF` maps 1:1 onto
  TranslonScorer's own event fields — the ingestion side is easy.
- But RDG's flat candidate builder (`rdg_document_from_candidates`) does
  **not** merge overlapping candidates into one shared graph — every
  candidate gets its own independent sibling path. The builder that
  actually merges/shares structure (`insert_translon`/`build_orf` in
  `RDG.py`) only works from raw FASTA-sequence scanning, not from an
  arbitrary candidate list. So "gene = splice graph + RDG, built from our
  own candidate calls" is **not something RDG does today** — it would be
  real new construction work, not a matter of calling existing library
  functions in the right order.

## Interpretability requirement, restated (the part to hold onto regardless
of how the grouping-unit question resolves)

Whatever computes "what's nearby," the output contract should be: retain a
per-event neighbor list (ids, spans, frames, distances, each neighbor's own
battery outcome), and surface it in the report — not just consume it into a
single confidence scalar or outcome label. That's true whether the
underlying computation ends up being a per-chromosome sweep, a per-gene
graph, or something else entirely; it's the actual interpretability
requirement independent of implementation.

## Near-term, smaller, independently-scoped step

Regardless of how the larger question above eventually resolves: extending
`_contention`/`event_overlap` with a proximity slop and wiring it into
`score_bams_workflow`/`score_matrix_workflow` is well-scoped, chromosome-
level, and already close to done (the machinery exists, it's just
disconnected). This can proceed without waiting on the gene/graph question.

## Explicitly open, not decided here

- Chromosome-level flat table vs. gene-level splice-graph+RDG object vs.
  something else, as the grouping unit for attribution/interpretability.
- Reconciling any gene-level design with the L0b/L1 junction-set cache
  boundary decision.
- Whether RDG's graph-merging construction path needs real new engineering
  to support arbitrary multi-candidate input, or whether a lighter, custom,
  non-RDG structure is the more honest choice given how much of RDG's
  "decision graph" framing (branch-point exclusivity) was already found to
  be the wrong mental model for population-scale Ribo-seq data.

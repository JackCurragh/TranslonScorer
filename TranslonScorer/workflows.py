"""Orchestration layer — the imperative shell that wires the functional core.

These functions compose the pure pieces (events extraction, coverage providers,
the event scorer, the consequentiality policy, the score store) into end-to-end
runs. They own *all* the I/O and sequencing; the modules they call stay pure.

Public API
----------
extract_events_workflow  — annotation sqlite → genomic event Parquet store
score_matrix_workflow    — events + sparse matrix partitions → score store
score_bams_workflow      — events + genome BAMs → score store
report_workflow          — score store + events → per-translon report + policy
consequential_workflow   — per-translon report → consequentiality labels
pipeline_workflow        — one-shot extract → score → report (matrix or BAMs)

Scoring contract (shared by score_matrix/score_bams, mirrors the golden gate):
a single per-position A-site coverage table per chromosome is fed to
score_events_vectorised, which scores init/term scalar and elongation via prefix
sums. Coverage is queried per chromosome so genomic positions never collide
across chromosomes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import polars as pl

from TranslonScorer.consequential import apply_policy
from TranslonScorer.events import SpliceContext, build_splice_context, run_extract
from TranslonScorer.io.store import (
    persist_scores,
    read_events,
    read_feature_event,
    read_scores,
)
from TranslonScorer.model import (
    ConsequentialityPolicy,
    OffsetParams,
    Region,
    ScoreThresholds,
)
from TranslonScorer.report import compose_report
from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS, score_events_vectorised

# ---------------------------------------------------------------------------
# extract-events
# ---------------------------------------------------------------------------


def extract_events_workflow(
    out_dir: str,
    *,
    sqlite_path: Optional[str] = None,
    gtf_path: Optional[str] = None,
    feature_type: str = "CDS",
    bed12_path: Optional[str] = None,
    bigbed_path: Optional[str] = None,
    fasta_path: Optional[str] = None,
    start_codons: Optional[List[str]] = None,
    stop_codons: Optional[List[str]] = None,
    min_len: int = 0,
    max_len: int = 1_000_000,
    annotation_version: str = "",
    chroms: Optional[List[str]] = None,
    context_gtf: Optional[str] = None,
    context_flank: int = 200,
) -> dict:
    """Extract deduplicated genomic events from any feature source.

    Exactly one source must be given: ``sqlite_path`` (annotation DB),
    ``gtf_path`` (GTF/GFF, scoring ``feature_type``), ``bed12_path``,
    ``bigbed_path``, or ``fasta_path`` (de-novo ORF finding). Writes events/,
    feature_event/, event_overlap/ Parquet trees under ``out_dir``.
    """
    from TranslonScorer.events import write_events
    from TranslonScorer.io import feature_sources as fs

    sources = [sqlite_path, gtf_path, bed12_path, bigbed_path, fasta_path]
    if sum(s is not None for s in sources) != 1:
        raise ValueError(
            "provide exactly one feature source: sqlite_path | gtf_path | "
            "bed12_path | bigbed_path | fasta_path"
        )

    # sqlite keeps its streaming per-chrom reader (the 8.85M-translon scale).
    if sqlite_path is not None:
        return run_extract(
            sqlite_path, out_dir, annotation_version=annotation_version, chroms=chroms
        )

    if gtf_path is not None:
        blocks, translons = fs.from_gtf(gtf_path, feature_type=feature_type)
    elif bed12_path is not None:
        blocks, translons = fs.from_bed12(bed12_path)
    elif bigbed_path is not None:
        blocks, translons = fs.from_bigbed(bigbed_path)
    else:
        blocks, translons = fs.from_fasta(
            fasta_path,
            start_codons=start_codons,
            stop_codons=stop_codons,
            min_len=min_len,
            max_len=max_len,
        )
    return write_events(
        blocks,
        translons,
        out_dir,
        annotation_version=annotation_version,
        chroms=chroms,
        context_gtf=context_gtf,
        context_flank=context_flank,
    )


# ---------------------------------------------------------------------------
# shared scoring helper
# ---------------------------------------------------------------------------


# Flank pad (nt) around each event when fetching coverage. The init/term
# scorers read at most ±60 nt (flanks 9/18/30/60); 128 covers that plus the
# P/A-site offset with margin, so scores are identical to a whole-chromosome
# fetch while pulling only reads near events.
_EVENT_FLANK_PAD = 128


def _merge_event_regions(
    chrom: str, starts: List[int], ends: List[int], pad: int = _EVENT_FLANK_PAD
) -> List[Region]:
    """Padded, merged genomic intervals covering every event's scoring window.

    Querying these instead of one chromosome-spanning region pulls only reads
    near events — orders of magnitude less on a sparse event set over a deep
    matrix — without changing any score (scorers read within ±60 nt of an event).
    """
    ivs = sorted((max(0, s - pad), e + pad) for s, e in zip(starts, ends))
    merged: List[Region] = []
    for s, e in ivs:
        if merged and s <= merged[-1].end:
            merged[-1] = Region(chrom, merged[-1].start, max(merged[-1].end, e))
        else:
            merged.append(Region(chrom, s, e))
    return merged


def _junction_support_for_chrom(provider, ev_chrom: pl.DataFrame) -> Dict[int, dict]:
    """Spanning-read support for every junction event on one chromosome.

    Was previously never called: score_matrix_workflow/score_bams_workflow
    always scored junction events against an empty support dict (every
    junction event fell out as INSUFFICIENT/n_reads=0), even though both
    MatrixProvider.junction_support() and BamSetProvider.junction_support()
    are fully implemented — this was a pure wiring gap, not a missing
    capability.

    A runtime-checkable Protocol only checks that a `junction_support` METHOD
    exists, not that it works — BigwigSetProvider defines one that always
    raises NotImplementedError (bigwigs have no CIGAR), so `isinstance(...,
    SupportsJunctions)` is True for it too. NotImplementedError is therefore
    caught explicitly below rather than relied on the isinstance check alone.

    Returns {junction_event_id: {"span_conf": n, "span_short": n, "unspliced": n}}.
    """
    from TranslonScorer.coverage.base import SupportsJunctions

    junc_ev = ev_chrom.filter(pl.col("type") == "junction")
    if junc_ev.is_empty() or not isinstance(provider, SupportsJunctions):
        return {}

    junctions = list(
        zip(
            junc_ev["chrom"].to_list(),
            junc_ev["start"].to_list(),
            junc_ev["end"].to_list(),
            junc_ev["strand"].to_list(),
            junc_ev["event_id"].to_list(),
        )
    )
    try:
        supp_df = provider.junction_support(junctions)
    except NotImplementedError:
        return {}
    out: Dict[int, dict] = {}
    if supp_df is None or supp_df.is_empty():
        return out
    for row in supp_df.iter_rows(named=True):
        d = out.setdefault(int(row["junction_id"]), {})
        d[row["kind"]] = d.get(row["kind"], 0.0) + float(row["count"])
    return out


def _map_track_for_chrom(
    map_provider, ev_chrom: pl.DataFrame, thr: ScoreThresholds
) -> Dict[int, dict]:
    """Region-context mappability annotation for every event on one chromosome.

    Not the BAM-NH-tag SupportsMappability/mappability_ledger concept
    (coverage/base.py) — that's per-read unique-vs-multimapper accounting.
    This is a precomputed mappability TRACK (e.g. Umap/GEM-mappability
    bigwig, values ~0..1) read via the same padded per-event windows as
    coverage, purely diagnostic: never affects eligibility/call (same
    "review flag, not a gate" precedent as flank_peakiness/stability).

    map_provider is a BigwigSetProvider (or anything with a matching
    .coverage() method) reading the mappability bigwig; None disables this
    entirely (the common case).

    NB: positions absent from the mappability bigwig collapse to the same
    "0 count" as a genuine low-mappability position in BigwigSetProvider's
    output (see coverage/bigwig.py) — a data gap (wrong chrom name,
    off-contig) reads identically to "confirmed unmappable". Acceptable for
    a diagnostic-only annotation; not resolved here.

    Returns {event_id: {"map_track_mean": float, "map_track_low": bool}}.
    """
    if map_provider is None or ev_chrom.is_empty():
        return {}

    regions = _merge_event_regions(
        str(ev_chrom["chrom"][0]),
        ev_chrom["start"].to_list(),
        ev_chrom["end"].to_list(),
    )
    map_cov = map_provider.coverage(regions)
    if map_cov.is_empty():
        vals: Dict[int, float] = {}
    else:
        pos_val = dict(zip(map_cov["pos"].to_list(), map_cov["count"].to_list()))
        vals = pos_val

    out: Dict[int, dict] = {}
    for eid, start, end in zip(
        ev_chrom["event_id"].to_list(), ev_chrom["start"].to_list(), ev_chrom["end"].to_list()
    ):
        # Each event's OWN padded window (not the merged query regions, which
        # can span several nearby events) — same pad as the coverage fetch,
        # so init/term (single-nt events) get the flank actually read by the
        # step scorers, not just the start/stop codon's own position.
        span = range(max(0, start - _EVENT_FLANK_PAD), end + _EVENT_FLANK_PAD)
        levels = [vals.get(p, 0.0) for p in span]
        mean_val = sum(levels) / len(levels) if levels else 0.0
        out[int(eid)] = {
            "map_track_mean": mean_val,
            "map_track_low": mean_val < thr.mappability_low,
        }
    return out


def _score_events_over_provider(
    events: pl.DataFrame,
    provider,
    *,
    site: str,
    group: str,
    tier: str,
    thr: ScoreThresholds,
    splice_context: Optional[SpliceContext] = None,
    map_provider=None,
) -> pl.DataFrame:
    """Score every event by querying ``provider`` per chromosome.

    Coverage is fetched as padded per-event windows (merged), not one
    chromosome-spanning region: identical scores, but only reads near events are
    pulled — critical on a deep matrix where a whole-chromosome span is enormous.

    ``splice_context``, if given, is forwarded to score_events_vectorised so
    init/term leader/UTR flanks near a splice site are projected across the
    intron instead of read as raw flanking genomic bases.

    ``map_provider``, if given, is a coverage-shaped provider (typically
    BigwigSetProvider reading a mappability track) queried once per
    chromosome via _map_track_for_chrom for a diagnostic map_track_mean/
    map_track_low annotation — see that function's docstring. Never affects
    eligibility/call.
    """
    if events.is_empty():
        from TranslonScorer.scoring.evidence import _RECORD_SCHEMA

        return pl.DataFrame(schema=_RECORD_SCHEMA)

    parts: List[pl.DataFrame] = []
    for chrom in events["chrom"].unique().sort().to_list():
        ev_chrom = events.filter(pl.col("chrom") == chrom)
        regions = _merge_event_regions(
            str(chrom),
            ev_chrom["start"].to_list(),
            ev_chrom["end"].to_list(),
        )
        cov_df = provider.coverage(regions, site=site)
        if cov_df.is_empty():
            continue
        junction_support = _junction_support_for_chrom(provider, ev_chrom)
        map_track = _map_track_for_chrom(map_provider, ev_chrom, thr)
        # Score each strand against its OWN coverage. Ribo-seq is stranded: a
        # + event must see only + reads. Collapsing strands (and building a
        # pos->count dict) would let +/- coverage at the same genomic position
        # overwrite each other — wrong, and order-dependent (nondeterministic).
        has_strand = "strand" in cov_df.columns
        for strand_val in (1, -1):
            ev_s = ev_chrom.filter(pl.col("strand") == strand_val)
            if ev_s.is_empty():
                continue
            cov_s = (
                cov_df.filter(pl.col("strand") == strand_val) if has_strand else cov_df
            ).select(["pos", "count"])
            # collapse any duplicate positions deterministically
            cov_s = cov_s.group_by("pos").agg(pl.col("count").sum()).sort("pos")
            scored = score_events_vectorised(
                ev_s,
                cov_s,
                group=group,
                tier=tier,
                thr=thr,
                splice_context=splice_context,
                junction_support=junction_support,
                map_track=map_track,
            )
            if not scored.is_empty():
                parts.append(scored)

    if not parts:
        from TranslonScorer.scoring.evidence import _RECORD_SCHEMA

        return pl.DataFrame(schema=_RECORD_SCHEMA)
    return pl.concat(parts)


# ---------------------------------------------------------------------------
# score-matrix
# ---------------------------------------------------------------------------


def score_matrix_workflow(
    events_dir: str,
    partition_dirs,
    store_dir: str,
    *,
    data_version: str,
    ref_offset: int = 15,
    psite_index_dir: Optional[str] = None,
    sample_names: Optional[List[str]] = None,
    n_workers: Optional[int] = None,
    site: str = "A",
    group: str = "aggregate",
    tier: str = "aggregate",
    annotation_version: str = "",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
    context_gtf: Optional[str] = None,
    mappability_bigwig: Optional[str] = None,
) -> str:
    """Score events against the sparse annotation-scale matrix; persist results.

    ``psite_index_dir``, if given (a directory produced by ``build-psite-index``),
    routes coverage through per-(sample, length) offset-corrected genomic
    positions instead of the flat ``ref_offset`` — fixes the elong_in_frame
    ~0.33 accuracy floor without leaving this pipeline, unlike the separate
    FrameRollup/``score_matrix_rollup_workflow`` path (elongation-only). This
    is now the recommended way to get accurate matrix scoring: it also gets
    init/term/junction/mappability, which FrameRollup never will. See
    MatrixProvider and psite_index.query_genomic_coverage. ``ref_offset`` is
    ignored when this is set.

    ``context_gtf``, if given, provides the host-transcript exon models used to
    project init/term leader/UTR flanks across nearby introns (see
    _score_events_over_provider). Same GTF you'd pass to extract-events'
    ``context_gtf`` for context junction events.

    ``mappability_bigwig``, if given, is a precomputed mappability track
    (Umap/GEM-mappability style) used to annotate every scored event with a
    diagnostic map_track_mean/map_track_low — see _map_track_for_chrom. Never
    affects eligibility/call.

    Returns the written store path(s) (semicolon-joined), as ``persist_scores``.
    """
    from TranslonScorer.coverage.bigwig import BigwigSetProvider
    from TranslonScorer.matrix.provider import MatrixProvider

    events = read_events(events_dir)
    provider = MatrixProvider(
        partition_dirs,
        ref_offset=ref_offset,
        psite_index_dir=psite_index_dir,
        sample_names=sample_names,
        n_workers=n_workers,
    )
    splice_context = build_splice_context(context_gtf) if context_gtf else None
    map_provider = BigwigSetProvider([mappability_bigwig]) if mappability_bigwig else None
    scored = _score_events_over_provider(
        events,
        provider,
        site=site,
        group=group,
        tier=tier,
        thr=thr,
        splice_context=splice_context,
        map_provider=map_provider,
    )
    return persist_scores(
        scored,
        store_dir,
        data_version=data_version,
        annotation_version=annotation_version,
    )


# ---------------------------------------------------------------------------
# score-matrix-rollup  (FrameRollup path — NFR2, FR2/FR3)
# ---------------------------------------------------------------------------


def score_matrix_rollup_workflow(
    events_dir: str,
    partition_dirs,
    store_dir: str,
    cds_df: pl.DataFrame,
    *,
    data_version: str,
    sample_names: Optional[List[str]] = None,
    n_workers: Optional[int] = None,
    multimap_mode: str = "unique",
    annotation_version: str = "",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
    psite_index_dir: Optional[str] = None,
) -> str:
    """Score elongation events using the FrameRollup (calibrated per-sample, per-length offsets).

    Replaces the flat-offset MatrixProvider path for the elongation aspect.
    init_rise / term_drop scoring is not yet wired here (retained in score_matrix_workflow).

    Parameters
    ----------
    cds_df : output of build_cds_blocks(), used to build the FrameRollup.
    """
    from TranslonScorer.matrix.rollup import (
        build_frame_rollup,
        calibrate_offsets,
        score_frame_rollup,
    )
    from TranslonScorer.scoring.evidence import _RECORD_SCHEMA, event_record
    from TranslonScorer.scoring.run import score_elongation_from_rollup

    events = read_events(events_dir)
    if events.is_empty():
        from TranslonScorer.scoring.evidence import _RECORD_SCHEMA as _S

        return persist_scores(pl.DataFrame(schema=_S), store_dir, data_version=data_version)

    # Build FrameRollup for CDS features referenced by events
    rollup = build_frame_rollup(
        partition_dirs,
        cds_df,
        multimap_mode=multimap_mode,
        sample_names=sample_names,
        n_workers=n_workers,
    )

    # Load or calibrate per-(sample, length) offsets
    if psite_index_dir is not None:
        offsets_path = Path(psite_index_dir) / "offsets.parquet"
        usable_path = Path(psite_index_dir) / "usable_sample_lengths.parquet"
        offsets_df = pl.read_parquet(offsets_path)
        offsets = {
            (str(r["sample_name"]), int(r["length"])): int(r["offset"])
            for r in offsets_df.iter_rows(named=True)
        }
        if usable_path.exists():
            usable_df = (
                pl.read_parquet(usable_path)
                .rename({"sample_id": "sample_name"})
                .select(["sample_name", pl.col("length").cast(pl.Int64)])
                .with_columns(pl.lit(True).alias("_keep"))
            )
            rollup = (
                rollup.join(usable_df, on=["sample_name", "length"], how="left")
                .filter(pl.col("_keep").fill_null(False))
                .drop("_keep")
            )
    else:
        agg = rollup.group_by(["sample_name", "length", "strand", "phase0"]).agg(
            pl.col("count").sum()
        )
        offsets = calibrate_offsets(agg, target_frame=0)

    # Score: aggregate across samples → per-(feature_id, length) elong_in_frame
    scored = score_frame_rollup(rollup, offsets, default_offset=12)

    # Map feature-level scores to events
    elong_ev = score_elongation_from_rollup(events, scored, thr=thr)

    # Serialise to long-form record table
    rows = []
    for r in events.filter(pl.col("type") == "elongation").iter_rows(named=True):
        raw = elong_ev.get(r["event_id"])
        if raw is None:
            continue
        rows.append(event_record(r["event_id"], "elongation", "aggregate", "aggregate", raw, thr.version))

    result = (
        pl.from_dicts(rows, schema=_RECORD_SCHEMA)
        if rows
        else pl.DataFrame(schema=_RECORD_SCHEMA)
    )
    return persist_scores(
        result,
        store_dir,
        data_version=data_version,
        annotation_version=annotation_version,
    )


# ---------------------------------------------------------------------------
# score-bams
# ---------------------------------------------------------------------------


def score_bams_workflow(
    events_dir: str,
    bams: Sequence[Union[str, Path]],
    store_dir: str,
    *,
    data_version: str,
    offsets: OffsetParams = OffsetParams(),
    sample_names: Optional[List[str]] = None,
    multimap: str = "unique",
    site: str = "A",
    group: str = "aggregate",
    tier: str = "aggregate",
    annotation_version: str = "",
    transcriptome: bool = False,
    exon_df: Optional[pl.DataFrame] = None,
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
    context_gtf: Optional[str] = None,
    mappability_bigwig: Optional[str] = None,
) -> str:
    """Score events against a set of genome- (or transcriptome-) aligned BAMs.

    Offsets are calibrated once per BAM/read-length before any locus query
    (handled inside BamSetProvider). When ``transcriptome=True`` reads are
    projected to genome coordinates via ``exon_df`` (required). ``context_gtf``
    projects init/term leader/UTR flanks across nearby introns — see
    score_matrix_workflow. Returns the written store path(s).

    ``offsets.method == "metagene"`` needs genomic start-codon coordinates to
    calibrate against; these are derived here from the extracted ``init``
    events (free — events_dir is already read for scoring) and passed to
    BamSetProvider. Previously score_bams_workflow never supplied them, so
    metagene calibration silently fell back to the fixed global offset every
    time. NB: the metagene histogram is built by fetching directly from the
    BAM by genomic coordinate — it does not go through the transcriptome→
    genome projection, so metagene calibration is only correct for
    ``transcriptome=False`` (genome-aligned) BAMs; with ``transcriptome=True``
    it will still find no reads and fall back to global, same as before.

    ``mappability_bigwig`` — see score_matrix_workflow.
    """
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.coverage.bigwig import BigwigSetProvider

    events = read_events(events_dir)
    start_codons = None
    if offsets.method == "metagene" and not transcriptome:
        init_ev = events.filter(pl.col("type") == "init")
        if not init_ev.is_empty():
            start_codons = list(
                zip(
                    init_ev["chrom"].to_list(),
                    init_ev["start"].to_list(),
                    init_ev["strand"].to_list(),
                )
            )
    provider = BamSetProvider(
        list(bams),
        offsets=offsets,
        multimap=multimap,
        sample_names=sample_names,
        transcriptome=transcriptome,
        exon_df=exon_df,
        start_codons=start_codons,
    )
    splice_context = build_splice_context(context_gtf) if context_gtf else None
    map_provider = BigwigSetProvider([mappability_bigwig]) if mappability_bigwig else None
    scored = _score_events_over_provider(
        events,
        provider,
        site=site,
        group=group,
        tier=tier,
        thr=thr,
        splice_context=splice_context,
        map_provider=map_provider,
    )
    return persist_scores(
        scored,
        store_dir,
        data_version=data_version,
        annotation_version=annotation_version,
    )


# ---------------------------------------------------------------------------
# score-bigwig
# ---------------------------------------------------------------------------


def score_bigwigs_workflow(
    events_dir: str,
    bigwigs: Sequence[Union[str, Path, dict]],
    store_dir: str,
    *,
    data_version: str,
    sample_names: Optional[List[str]] = None,
    stranded: bool = False,
    site: str = "A",
    group: str = "aggregate",
    tier: str = "aggregate",
    annotation_version: str = "",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
    context_gtf: Optional[str] = None,
    mappability_bigwig: Optional[str] = None,
) -> str:
    """Score events against 1-N genomic bigwig coverage tracks.

    Bigwig is a LOSSY coverage source (see coverage/bigwig.py): no P/A-site
    distinction (``site`` is accepted but ignored), no junction spanning, no
    multimapper resolution — init/term/elongation only, junction events fall
    out INSUFFICIENT (BigwigSetProvider has no CIGAR to count spanning reads).

    ``bigwigs``: each entry is a path (unstranded) or, when ``stranded=True``,
    a ``{'forward': path, 'reverse': path}`` dict — REQUIRED for correct
    scoring, since Ribo-seq is stranded and an unstranded bigwig silently
    mixes +/- signal at every position. ``context_gtf`` projects init/term
    leader/UTR flanks across nearby introns — see score_matrix_workflow.
    ``mappability_bigwig`` — see score_matrix_workflow; independent of
    ``bigwigs`` (a separate, always-unstranded track). Returns the written
    store path(s).
    """
    from TranslonScorer.coverage.bigwig import BigwigSetProvider

    events = read_events(events_dir)
    provider = BigwigSetProvider(
        list(bigwigs),
        sample_names=sample_names,
        stranded=stranded,
    )
    splice_context = build_splice_context(context_gtf) if context_gtf else None
    map_provider = BigwigSetProvider([mappability_bigwig]) if mappability_bigwig else None
    scored = _score_events_over_provider(
        events,
        provider,
        site=site,
        group=group,
        tier=tier,
        thr=thr,
        splice_context=splice_context,
        map_provider=map_provider,
    )
    return persist_scores(
        scored,
        store_dir,
        data_version=data_version,
        annotation_version=annotation_version,
    )


# ---------------------------------------------------------------------------
# consequential
# ---------------------------------------------------------------------------


def consequential_workflow(
    report: pl.DataFrame,
    policy: ConsequentialityPolicy = ConsequentialityPolicy(),
) -> pl.DataFrame:
    """Apply a consequentiality policy to a per-translon report (pure pass-through
    over ``consequential.apply_policy``; kept here so the CLI has one import
    surface for orchestration)."""
    return apply_policy(report, policy)


# ---------------------------------------------------------------------------
# report  (scores store + events → per-translon report [+ consequentiality])
# ---------------------------------------------------------------------------


def report_workflow(
    store_dir: str,
    events_dir: str,
    out_path: str,
    *,
    data_version: Optional[str] = None,
    tier: Optional[str] = None,
    policy: Optional[ConsequentialityPolicy] = None,
) -> pl.DataFrame:
    """Compose a per-translon report from the score store and event membership.

    Reads scores (optionally filtered by data_version/tier) and the
    ``feature_event`` mapping, composes one row per translon, applies the
    consequentiality policy, writes Parquet, and returns the report.
    """
    scores = read_scores(store_dir, data_version=data_version, tier=tier)
    feature_event = read_feature_event(events_dir)
    report = compose_report(scores, feature_event)
    report = apply_policy(report, policy or ConsequentialityPolicy())
    report.write_parquet(out_path)
    return report


# ---------------------------------------------------------------------------
# pipeline  (one-shot: extract-events -> score -> report)
# ---------------------------------------------------------------------------


def pipeline_workflow(
    out_dir: str,
    *,
    partition_dirs: Optional[Sequence[Union[str, Path]]] = None,
    bams: Optional[Sequence[Union[str, Path]]] = None,
    data_version: str = "run",
    chroms: Optional[List[str]] = None,
    annotation_version: str = "",
    offsets: OffsetParams = OffsetParams(),
    sample_names: Optional[List[str]] = None,
    multimap: str = "unique",
    site: str = "A",
    transcriptome: bool = False,
    exon_df: Optional[pl.DataFrame] = None,
    policy: Optional[ConsequentialityPolicy] = None,
    cds_df: Optional[pl.DataFrame] = None,
    n_workers: Optional[int] = None,
    context_gtf: Optional[str] = None,
    context_flank: int = 200,
    mappability_bigwig: Optional[str] = None,
    psite_index_dir: Optional[str] = None,
    # feature source for extract-events (exactly one; forwarded verbatim)
    **source_kwargs,
) -> Dict[str, str]:
    """Run the whole event-scoring path in one call.

    extract-events → score (matrix *or* BAMs) → report. Writes
    ``out_dir/{events,scores}`` and ``out_dir/report.parquet``; returns the
    paths produced. Exactly one of ``partition_dirs`` (matrix) or ``bams`` is the
    coverage source; ``source_kwargs`` selects the feature source for
    extract-events (sqlite_path | gtf_path | bed12_path | bigbed_path | fasta_path).

    cds_df: if provided (from build_cds_blocks or build_cds_blocks_from_bigbed),
    matrix mode uses the FrameRollup path with per-(sample,length) calibrated
    P-site offsets — but elongation only (no init/term/junction/mappability).
    Prefer ``psite_index_dir`` instead: same offset accuracy, full event-type
    coverage, one pipeline. Without either, the legacy flat-offset path is used.

    psite_index_dir: directory from ``build-psite-index``; routes matrix
    scoring through per-(sample, length) offset-corrected genomic coverage
    (see score_matrix_workflow, MatrixProvider). Ignored when ``cds_df`` is
    also given (FrameRollup takes precedence, unchanged behaviour).

    context_gtf: host-transcript exon models, forwarded to BOTH extract-events
    (adds context junction events near ORF boundaries) and scoring (projects
    init/term leader/UTR flanks across nearby introns). Previously this was
    only reachable via the standalone extract-events command — pipeline had
    no way to supply it at all.
    """
    if bool(partition_dirs) == bool(bams):
        raise ValueError("provide exactly one of partition_dirs (matrix) or bams")

    out = Path(out_dir)
    events_dir = str(out / "events")
    store_dir = str(out / "scores")
    report_path = str(out / "report.parquet")

    extract_events_workflow(
        events_dir,
        annotation_version=annotation_version,
        chroms=chroms,
        context_gtf=context_gtf,
        context_flank=context_flank,
        **source_kwargs,
    )
    if partition_dirs:
        if cds_df is not None:
            score_matrix_rollup_workflow(
                events_dir,
                list(partition_dirs),
                store_dir,
                cds_df,
                data_version=data_version,
                sample_names=sample_names,
                annotation_version=annotation_version,
                n_workers=n_workers,
            )
        else:
            score_matrix_workflow(
                events_dir,
                list(partition_dirs),
                store_dir,
                data_version=data_version,
                sample_names=sample_names,
                site=site,
                annotation_version=annotation_version,
                context_gtf=context_gtf,
                mappability_bigwig=mappability_bigwig,
                psite_index_dir=psite_index_dir,
            )
    else:
        score_bams_workflow(
            events_dir,
            list(bams or []),
            store_dir,
            data_version=data_version,
            offsets=offsets,
            sample_names=sample_names,
            multimap=multimap,
            site=site,
            annotation_version=annotation_version,
            transcriptome=transcriptome,
            mappability_bigwig=mappability_bigwig,
            context_gtf=context_gtf,
            exon_df=exon_df,
        )
    report_workflow(
        store_dir,
        events_dir,
        report_path,
        data_version=data_version,
        policy=policy,
    )
    return {"events": events_dir, "scores": store_dir, "report": report_path}

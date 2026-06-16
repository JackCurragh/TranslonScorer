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
from TranslonScorer.events import run_extract
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
        blocks, translons, out_dir, annotation_version=annotation_version, chroms=chroms
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


def _score_events_over_provider(
    events: pl.DataFrame,
    provider,
    *,
    site: str,
    group: str,
    tier: str,
    thr: ScoreThresholds,
) -> pl.DataFrame:
    """Score every event by querying ``provider`` per chromosome.

    Coverage is fetched as padded per-event windows (merged), not one
    chromosome-spanning region: identical scores, but only reads near events are
    pulled — critical on a deep matrix where a whole-chromosome span is enormous.
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
            scored = score_events_vectorised(ev_s, cov_s, group=group, tier=tier, thr=thr)
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
    sample_names: Optional[List[str]] = None,
    n_workers: Optional[int] = None,
    site: str = "A",
    group: str = "aggregate",
    tier: str = "aggregate",
    annotation_version: str = "",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
) -> str:
    """Score events against the sparse annotation-scale matrix; persist results.

    Returns the written store path(s) (semicolon-joined), as ``persist_scores``.
    """
    from TranslonScorer.coverage.matrix import MatrixProvider

    events = read_events(events_dir)
    provider = MatrixProvider(
        partition_dirs,
        ref_offset=ref_offset,
        sample_names=sample_names,
        n_workers=n_workers,
    )
    scored = _score_events_over_provider(
        events, provider, site=site, group=group, tier=tier, thr=thr
    )
    return persist_scores(
        scored,
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
) -> str:
    """Score events against a set of genome- (or transcriptome-) aligned BAMs.

    Offsets are calibrated once per BAM/read-length before any locus query
    (handled inside BamSetProvider). When ``transcriptome=True`` reads are
    projected to genome coordinates via ``exon_df`` (required). Returns the
    written store path(s).
    """
    from TranslonScorer.coverage.bam import BamSetProvider

    events = read_events(events_dir)
    provider = BamSetProvider(
        list(bams),
        offsets=offsets,
        multimap=multimap,
        sample_names=sample_names,
        transcriptome=transcriptome,
        exon_df=exon_df,
    )
    scored = _score_events_over_provider(
        events, provider, site=site, group=group, tier=tier, thr=thr
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
    # feature source for extract-events (exactly one; forwarded verbatim)
    **source_kwargs,
) -> Dict[str, str]:
    """Run the whole event-scoring path in one call.

    extract-events → score (matrix *or* BAMs) → report. Writes
    ``out_dir/{events,scores}`` and ``out_dir/report.parquet``; returns the
    paths produced. Exactly one of ``partition_dirs`` (matrix) or ``bams`` is the
    coverage source; ``source_kwargs`` selects the feature source for
    extract-events (sqlite_path | gtf_path | bed12_path | bigbed_path | fasta_path).
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
        **source_kwargs,
    )
    if partition_dirs:
        score_matrix_workflow(
            events_dir,
            list(partition_dirs),
            store_dir,
            data_version=data_version,
            sample_names=sample_names,
            site=site,
            annotation_version=annotation_version,
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

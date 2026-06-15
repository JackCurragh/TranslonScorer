"""Orchestration layer — the imperative shell that wires the functional core.

These functions compose the pure pieces (events extraction, coverage providers,
the event scorer, the consequentiality policy, the score store) into end-to-end
runs. They own *all* the I/O and sequencing; the modules they call stay pure.

Public API
----------
extract_events_workflow  — annotation sqlite → genomic event Parquet store
score_matrix_workflow    — events + sparse matrix partitions → score store
score_bams_workflow      — events + genome BAMs → score store
consequential_workflow   — per-translon report → consequentiality labels

Scoring contract (shared by score_matrix/score_bams, mirrors the golden gate):
a single per-position A-site coverage table per chromosome is fed to
score_events_vectorised, which scores init/term scalar and elongation via prefix
sums. Coverage is queried per chromosome so genomic positions never collide
across chromosomes.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union, cast

import polars as pl

from TranslonScorer.model import (
    ConsequentialityPolicy,
    OffsetParams,
    Region,
    ScoreThresholds,
)
from TranslonScorer.events import run_extract
from TranslonScorer.io.store import (
    persist_scores,
    read_events,
    read_feature_event,
    read_scores,
)
from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS, score_events_vectorised
from TranslonScorer.consequential import apply_policy
from TranslonScorer.report import compose_report


# ---------------------------------------------------------------------------
# extract-events
# ---------------------------------------------------------------------------

def extract_events_workflow(
    sqlite_path: str,
    out_dir: str,
    *,
    annotation_version: str = "",
    chroms: Optional[List[str]] = None,
) -> dict:
    """Extract deduplicated genomic events from an annotation sqlite db.

    Thin shell over ``events.run_extract``; writes events/, feature_event/ and
    event_overlap/ Parquet trees under ``out_dir`` and returns a summary dict.
    ``chroms`` optionally restricts extraction to specific chromosomes.
    """
    return run_extract(
        sqlite_path, out_dir,
        annotation_version=annotation_version, chroms=chroms,
    )


# ---------------------------------------------------------------------------
# shared scoring helper
# ---------------------------------------------------------------------------

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

    Coverage is fetched one chromosome at a time (one spanning Region per chrom)
    so genomic positions are unambiguous; the per-chrom score tables are stacked.
    """
    if events.is_empty():
        from TranslonScorer.scoring.evidence import _RECORD_SCHEMA
        return pl.DataFrame(schema=_RECORD_SCHEMA)

    parts: List[pl.DataFrame] = []
    for chrom in events["chrom"].unique().sort().to_list():
        ev_chrom = events.filter(pl.col("chrom") == chrom)
        start = int(cast(int, ev_chrom["start"].min()))
        end = int(cast(int, ev_chrom["end"].max()))
        region = Region(str(chrom), start, end + 1)
        cov_df = provider.coverage([region], site=site)
        if cov_df.is_empty():
            continue
        cov_df = cov_df.select(["pos", "count"])
        scored = score_events_vectorised(
            ev_chrom, cov_df, group=group, tier=tier, thr=thr
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
        scored, store_dir,
        data_version=data_version, annotation_version=annotation_version,
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
        scored, store_dir,
        data_version=data_version, annotation_version=annotation_version,
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

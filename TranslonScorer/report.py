"""Per-feature report composition — pure transform.

Consumes scored event DataFrames and assembles per-translon / per-feature
summary rows by joining the long-form event scores onto the translon→event
membership table (``feature_event`` from ``events.extract_events``).  No file
I/O; the caller writes parquet/JSON.

Why compose here: events are scored once and shared by many translons (≈78×
dedup). Composition fans the shared evidence back out so two isoforms sharing a
CDS exon get the *same* elongation score, never two contradictory ones.

Public API
----------
compose_report  — scored events [+ feature_event] → per-translon summary
"""
from __future__ import annotations

import polars as pl

# Core translation-chain aspects, in biological order.  Junction is handled
# separately (a translon may have 0..n junctions).
_CHAIN_ASPECTS = ("init", "elongation", "term")

_SUPPORTED = "SUPPORTED"


def compose_report(
    scores: pl.DataFrame,
    feature_event: pl.DataFrame | None = None,
    *,
    group: str = "",
    tier: str = "",
) -> pl.DataFrame:
    """Compose a per-translon summary from a scored-events DataFrame.

    Parameters
    ----------
    scores : output of score_events / score_events_vectorised — must have
             columns event_id, aspect, eligibility, call, group, tier and
             (for composition) n_reads, metric.
    feature_event : translon→event membership (feature_id, event_id, role[,
             rank, phase]) from ``events.extract_events``.  When ``None`` the
             function preserves the legacy filter-only behaviour and returns the
             (optionally group/tier-filtered) scores unchanged.
    group  : filter to this group label (empty string = all groups).
    tier   : filter to this tier label (empty string = all tiers).

    Returns
    -------
    When ``feature_event`` is None: the filtered long-form scores.
    Otherwise: one row per ``feature_id`` with, per chain aspect, the columns
    ``{aspect}_call``, ``{aspect}_metric``, ``{aspect}_n_reads`` and
    ``{aspect}_supported_frac``; plus ``junction_n``, ``junction_supported_frac``;
    plus ``total_reads`` (sum across the translon's events) and ``n_events``.
    """
    df = scores
    if group:
        df = df.filter(pl.col("group") == group)
    if tier:
        df = df.filter(pl.col("tier") == tier)

    if feature_event is None:
        return df

    return _compose_per_translon(df, feature_event)


def _compose_per_translon(
    scores: pl.DataFrame, feature_event: pl.DataFrame
) -> pl.DataFrame:
    """Join scores onto translon membership and aggregate to one row/translon."""
    score_cols = ["event_id", "aspect", "call", "eligibility", "n_reads", "metric"]
    have = [c for c in score_cols if c in scores.columns]
    joined = feature_event.select(["feature_id", "event_id"]).join(
        scores.select(have), on="event_id", how="left"
    )
    if joined.is_empty():
        return pl.DataFrame(schema={"feature_id": pl.Utf8})

    # Per (feature, aspect) aggregates: read-weighted metric, summed reads,
    # supported fraction, event count.
    per_aspect = (
        joined.group_by(["feature_id", "aspect"])
        .agg(
            pl.col("n_reads").sum().alias("n_reads"),
            (
                (pl.col("metric") * pl.col("n_reads")).sum()
                / pl.when(pl.col("n_reads").sum() > 0)
                .then(pl.col("n_reads").sum())
                .otherwise(None)
            ).alias("metric"),
            (pl.col("call") == _SUPPORTED).mean().alias("supported_frac"),
            pl.len().alias("n_events"),
        )
    )

    # Totals across the whole translon (every event, any aspect).
    totals = joined.group_by("feature_id").agg(
        pl.col("n_reads").sum().alias("total_reads"),
        pl.col("event_id").n_unique().alias("n_events"),
    )

    out = totals
    for aspect in _CHAIN_ASPECTS:
        a = per_aspect.filter(pl.col("aspect") == aspect).select(
            "feature_id",
            pl.when(pl.col("supported_frac") >= 0.5)
            .then(pl.lit(_SUPPORTED))
            .otherwise(pl.lit("UNSUPPORTED"))
            .alias(f"{aspect}_call"),
            pl.col("metric").alias(f"{aspect}_metric"),
            pl.col("n_reads").alias(f"{aspect}_n_reads"),
            pl.col("supported_frac").alias(f"{aspect}_supported_frac"),
        )
        out = out.join(a, on="feature_id", how="left")

    junc = per_aspect.filter(pl.col("aspect") == "junction").select(
        "feature_id",
        pl.col("n_events").alias("junction_n"),
        pl.col("supported_frac").alias("junction_supported_frac"),
    )
    out = out.join(junc, on="feature_id", how="left")

    return out.sort("feature_id")

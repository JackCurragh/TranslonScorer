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

_SUPPORTED = "SUPPORTED"

# Display order only — aspects not listed here still appear, sorted after these.
# Membership is discovered from the data, never fixed here, so a new scorer's
# aspect flows into the report with no change to this file.
_ASPECT_ORDER = ("init", "elongation", "term", "junction")


def _ordered_aspects(aspects: list[str]) -> list[str]:
    rank = {a: i for i, a in enumerate(_ASPECT_ORDER)}
    return sorted(aspects, key=lambda a: (rank.get(a, len(_ASPECT_ORDER)), a))


def compose_report(
    scores: pl.DataFrame,
    feature_event: pl.DataFrame | None = None,
    *,
    group: str = "",
    tier: str = "",
    supported_frac_min: float = 0.5,
) -> pl.DataFrame:
    """Compose a per-translon summary from a scored-events DataFrame.

    Parameters
    ----------
    scores : output of score_events — must have
             columns event_id, aspect, eligibility, call, group, tier and
             (for composition) n_reads, metric.
    feature_event : translon→event membership (feature_id, event_id, role[,
             rank, phase]) from ``events.extract_events``.  When ``None`` the
             function preserves the legacy filter-only behaviour and returns the
             (optionally group/tier-filtered) scores unchanged.
    group  : filter to this group label (empty string = all groups).
    tier   : filter to this tier label (empty string = all tiers).
    supported_frac_min : a feature's ``{aspect}_call`` is SUPPORTED when at least
             this fraction of the aspect's events are individually SUPPORTED.

    Returns
    -------
    When ``feature_event`` is None: the filtered long-form scores.
    Otherwise: one row per ``feature_id`` with, for **every** aspect present in
    ``scores``, the columns ``{aspect}_call``, ``{aspect}_metric``,
    ``{aspect}_n_reads`` and ``{aspect}_supported_frac`` (plus
    ``{aspect}_map_track_*`` when a mappability track was used); plus
    ``total_reads`` (sum across the translon's events) and ``n_events``.
    """
    df = scores
    if group:
        df = df.filter(pl.col("group") == group)
    if tier:
        df = df.filter(pl.col("tier") == tier)

    if feature_event is None:
        return df

    return _compose_per_translon(df, feature_event, supported_frac_min)


def _compose_per_translon(
    scores: pl.DataFrame,
    feature_event: pl.DataFrame,
    supported_frac_min: float = 0.5,
) -> pl.DataFrame:
    """Join scores onto translon membership and aggregate to one row/translon."""
    score_cols = [
        "event_id",
        "aspect",
        "call",
        "eligibility",
        "n_reads",
        "metric",
        "map_track_mean",
        "map_track_low",
    ]
    have = [c for c in score_cols if c in scores.columns]
    joined = feature_event.select(["feature_id", "event_id"]).join(
        scores.select(have), on="event_id", how="left"
    )
    if joined.is_empty():
        return pl.DataFrame(schema={"feature_id": pl.Utf8})

    has_map_track = "map_track_mean" in joined.columns
    map_track_aggs = (
        [
            # Unweighted mean (NOT read-weighted like `metric`): read-weighting
            # would wash out the signal for exactly the zero/low-read events
            # this annotation exists to explain.
            pl.col("map_track_mean").mean().alias("map_track_mean"),
            # True if ANY contributing event was individually flagged low —
            # a single spurious/genuine repeat-masked event is worth surfacing.
            pl.col("map_track_low").any().alias("map_track_low"),
        ]
        if has_map_track
        else []
    )

    # Per (feature, aspect) aggregates: read-weighted metric, summed reads,
    # supported fraction, event count, unweighted mappability-track mean.
    per_aspect = joined.group_by(["feature_id", "aspect"]).agg(
        pl.col("n_reads").sum().alias("n_reads"),
        (
            (pl.col("metric") * pl.col("n_reads")).sum()
            / pl.when(pl.col("n_reads").sum() > 0).then(pl.col("n_reads").sum()).otherwise(None)
        ).alias("metric"),
        (pl.col("call") == _SUPPORTED).mean().alias("supported_frac"),
        pl.len().alias("n_events"),
        *map_track_aggs,
    )

    # Totals across the whole translon (every event, any aspect).
    totals = joined.group_by("feature_id").agg(
        pl.col("n_reads").sum().alias("total_reads"),
        pl.col("event_id").n_unique().alias("n_events"),
    )

    out = totals
    for aspect in _ordered_aspects(per_aspect["aspect"].unique().to_list()):
        select_cols = [
            "feature_id",
            pl.when(pl.col("supported_frac") >= supported_frac_min)
            .then(pl.lit(_SUPPORTED))
            .otherwise(pl.lit("UNSUPPORTED"))
            .alias(f"{aspect}_call"),
            pl.col("metric").alias(f"{aspect}_metric"),
            pl.col("n_reads").alias(f"{aspect}_n_reads"),
            pl.col("supported_frac").alias(f"{aspect}_supported_frac"),
        ]
        if has_map_track:
            select_cols += [
                pl.col("map_track_mean").alias(f"{aspect}_map_track_mean"),
                pl.col("map_track_low").alias(f"{aspect}_map_track_low"),
            ]
        a = per_aspect.filter(pl.col("aspect") == aspect).select(*select_cols)
        out = out.join(a, on="feature_id", how="left")

    return out.sort("feature_id")

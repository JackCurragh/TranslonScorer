"""Per-feature report composition — pure transform.

Consumes scored event DataFrames and assembles per-translon / per-feature
summary rows.  No file I/O; the caller writes parquet/JSON.

Public API
----------
compose_report  — scored events + annotation → per-translon summary DataFrame
"""
from __future__ import annotations

import polars as pl


def compose_report(
    scores: pl.DataFrame,
    *,
    group: str = "",
    tier: str = "",
) -> pl.DataFrame:
    """Compose a per-translon summary from a scored-events DataFrame.

    Parameters
    ----------
    scores : output of score_events / score_events_vectorised —
             must have columns event_id, aspect, eligibility, call, group, tier.
    group  : filter to this group label (empty string = all groups).
    tier   : filter to this tier label (empty string = all tiers).

    Returns a tidy DataFrame with one row per (event_id, aspect) summarising
    eligibility, call, and the raw score evidence.
    """
    df = scores
    if group:
        df = df.filter(pl.col("group") == group)
    if tier:
        df = df.filter(pl.col("tier") == tier)
    return df

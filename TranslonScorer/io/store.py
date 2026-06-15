"""Parquet read/write for the fact_event_score store and related tables.

Append-only, partitioned by data_version / tier.  Pure I/O: no computation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import polars as pl


def persist_scores(
    scored: pl.DataFrame,
    store_dir: str,
    *,
    data_version: str,
    annotation_version: str = "",
    shard: str = "part",
) -> str:
    """Append a scored table to the fact_event_score store.

    Append-only, partitioned by data_version / tier so re-scoring never mutates
    prior partitions. Returns the written path(s) as a semicolon-joined string.
    """
    if scored.is_empty():
        return ""
    scored = scored.with_columns(
        [
            pl.lit(data_version).alias("data_version"),
            pl.lit(annotation_version).alias("annotation_version"),
        ]
    )
    out = Path(store_dir) / f"data_version={data_version}"
    written = []
    for tier in scored["tier"].unique().to_list():
        d = out / f"tier={tier}"
        d.mkdir(parents=True, exist_ok=True)
        path = d / f"{shard}.parquet"
        scored.filter(pl.col("tier") == tier).write_parquet(path)
        written.append(str(path))
    return ";".join(written)


def read_scores(
    store_dir: str,
    *,
    data_version: Optional[str] = None,
    tier: Optional[str] = None,
) -> pl.DataFrame:
    """Read event scores from the fact_event_score store, optionally filtered."""
    root = Path(store_dir)
    pattern = "**/*.parquet"
    paths = [str(p) for p in sorted(root.rglob(pattern))]
    if not paths:
        return pl.DataFrame()
    df = pl.read_parquet(paths)
    if data_version is not None:
        df = df.filter(pl.col("data_version") == data_version)
    if tier is not None:
        df = df.filter(pl.col("tier") == tier)
    return df


def _read_event_subdir(base: str, subdir: str) -> pl.DataFrame:
    """Read one event-store sub-tree (one Parquet per chrom).

    Accepts either the sub-tree itself (``.../events``) or the extract-events
    output root (``...`` containing ``events/``), so callers can pass whichever
    path is natural.  Reads only direct ``*.parquet`` children (never descends
    into sibling sub-trees, whose schemas differ).
    """
    root = Path(base)
    if (root / subdir).is_dir():
        root = root / subdir
    paths = sorted(root.glob("*.parquet"))
    if not paths:
        return pl.DataFrame()
    return pl.concat([pl.read_parquet(str(p)) for p in paths])


def read_events(events_dir: str) -> pl.DataFrame:
    """Read events Parquet (events/, one file per chrom)."""
    return _read_event_subdir(events_dir, "events")


def read_feature_event(feature_event_dir: str) -> pl.DataFrame:
    """Read feature_event Parquet (feature_event/, one file per chrom)."""
    return _read_event_subdir(feature_event_dir, "feature_event")

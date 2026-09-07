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
    frames = [pl.read_parquet(str(p)) for p in paths]
    # Polars may infer large integer IDs differently when one chromosome has
    # no rows or only a narrow value range (e.g. Int128 vs UInt64).  Event IDs
    # are semantically one unsigned integer key; normalise later partitions to
    # the first partition's schema before concatenation.
    target = frames[0].schema.copy()
    # Event keys are declared UInt64 throughout the scorer.  Keeping the
    # canonical width here avoids pushing downstream joins onto Polars' much
    # slower Int128 execution path when one partition was inferred as Int128.
    for name in ("event_id", "other_event_id", "junction_id"):
        if name in target:
            target[name] = pl.UInt64
    normalised = []
    for frame in frames:
        casts = [
            pl.col(name).cast(dtype, strict=False).alias(name)
            for name, dtype in target.items()
            if frame.schema.get(name) != dtype
        ]
        normalised.append(frame.with_columns(casts) if casts else frame)
    return pl.concat(normalised)


def read_events(events_dir: str) -> pl.DataFrame:
    """Read events Parquet (events/, one file per chrom)."""
    return _read_event_subdir(events_dir, "events")


def read_feature_event(feature_event_dir: str) -> pl.DataFrame:
    """Read feature_event Parquet (feature_event/, one file per chrom)."""
    return _read_event_subdir(feature_event_dir, "feature_event")


def read_event_overlap(event_overlap_dir: str) -> pl.DataFrame:
    """Read event_overlap Parquet (event_overlap/, one file per chrom).

    Columns: event_id, other_event_id, overlap_start, overlap_end -- NOT
    including comp_phase (see events._contention): the overlap table is
    frame-agnostic on disk, so a caller building ``overlaps_df`` for
    ``scoring.run.score_events`` must join ``other_event_id`` back onto
    ``events``'s own ``phase`` column themselves (see
    ``workflows.py``'s ``_overlaps_df_for_scoring``). Empty (not missing)
    when a run produced no elongation overlaps at all.
    """
    return _read_event_subdir(event_overlap_dir, "event_overlap")

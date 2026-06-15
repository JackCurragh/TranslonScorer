from __future__ import annotations

import os
from typing import Optional

import polars as pl

from ..file_handlers import bigwig as bw_handlers
from ..utils.io import write_csv_safe, write_parquet_safe
from ..utils.logging import log_info
from .panel_manifest import load_panel_manifest, merge_panel_manifest
from .score_schema import (
    add_frame_score_columns,
    compare_score_tables,
    ensure_score_schema,
    load_score_table,
)


def _load_orf_or_exon_table(path: str) -> pl.DataFrame:
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    if path.endswith(".tsv"):
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def compare_raw_frame_scoring(
    *,
    orfs_path: str,
    exons_path: str,
    bigwig_path: str,
    frame_support_path: str,
    out_prefix: str,
    scoring_method: str = "modern",
    sru_range: int = 15,
    profiles_path: Optional[str] = None,
    panel_manifest_path: Optional[str] = None,
    label_column: Optional[str] = None,
    max_workers: Optional[int] = None,
) -> dict[str, str]:
    """Run Gate 1: raw scoring vs frame-weighted scoring on the same ORFs."""
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    old_scoring = scoring_method == "classic"
    orfs = _load_orf_or_exon_table(orfs_path)
    if panel_manifest_path:
        orfs = merge_panel_manifest(orfs, load_panel_manifest(panel_manifest_path))
    exons = _load_orf_or_exon_table(exons_path)
    frame_support = pl.read_parquet(frame_support_path)
    profiles = pl.read_parquet(profiles_path) if profiles_path else None

    log_info("Gate 1: scoring raw ORFs")
    raw = bw_handlers.scoring(
        bigwig_path,
        exons,
        orfs,
        old_scoring,
        sru_range,
        max_workers=max_workers,
    )
    raw = ensure_score_schema(
        raw,
        score_mode="raw",
        input_type="bigwig",
        frame_method="none",
        assignment_model="none",
    )

    log_info("Gate 1: scoring frame-weighted ORFs")
    frame = bw_handlers.scoring(
        bigwig_path,
        exons,
        orfs,
        old_scoring,
        sru_range,
        max_workers=max_workers,
        frame_weighted_scoring=True,
        frame_support_path=frame_support_path,
    )
    frame = ensure_score_schema(
        frame,
        score_mode="frame_weighted",
        input_type="bigwig",
        assignment_model="none",
    )
    frame = add_frame_score_columns(frame, frame_support, profiles=profiles)

    comparison, summary = compare_score_tables(raw, frame, label_column=label_column)

    paths = {
        "raw_scores": f"{out_prefix}.raw_scores.parquet",
        "frame_scores": f"{out_prefix}.frame_scores.parquet",
        "comparison": f"{out_prefix}.comparison.parquet",
        "summary": f"{out_prefix}.summary.csv",
    }
    write_parquet_safe(raw, paths["raw_scores"])
    write_parquet_safe(frame, paths["frame_scores"])
    write_parquet_safe(comparison, paths["comparison"])
    write_csv_safe(summary, paths["summary"])
    return paths


def summarize_existing_score_pair(
    *,
    raw_scores_path: str,
    frame_scores_path: str,
    out_prefix: str,
    label_column: Optional[str] = None,
) -> dict[str, str]:
    """Compare already-produced raw and frame-weighted score tables."""
    raw = load_score_table(raw_scores_path)
    frame = load_score_table(frame_scores_path)
    comparison, summary = compare_score_tables(raw, frame, label_column=label_column)
    paths = {
        "comparison": f"{out_prefix}.comparison.parquet",
        "summary": f"{out_prefix}.summary.csv",
    }
    write_parquet_safe(comparison, paths["comparison"])
    write_csv_safe(summary, paths["summary"])
    return paths

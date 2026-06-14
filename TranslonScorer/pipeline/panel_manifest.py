from __future__ import annotations

import hashlib
from typing import Any

import polars as pl

from .score_schema import ensure_orf_id


RECOMMENDED_PANEL_COLUMNS = [
    "panel_id",
    "category",
    "label",
    "tran_id",
    "start",
    "stop",
    "type",
    "gene_id",
    "evidence",
]


def load_panel_manifest(path: str) -> pl.DataFrame:
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    if path.endswith(".tsv"):
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def manifest_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_panel_manifest(panel: pl.DataFrame) -> pl.DataFrame:
    """Validate the frozen panel manifest and return a one-row report."""
    cols = set(panel.columns)
    has_orf_id = "orf_id" in cols
    has_coord_key = {"tran_id", "start", "stop"}.issubset(cols)
    errors: list[str] = []
    warnings: list[str] = []

    if not has_orf_id and not has_coord_key:
        errors.append("manifest requires either orf_id or tran_id,start,stop")
    if "label" not in cols:
        warnings.append("label column missing; score-gap summaries will be unavailable")
    if "category" not in cols:
        warnings.append("category column missing; panel stratification will be unavailable")
    if panel.is_empty():
        errors.append("manifest is empty")

    keyed = ensure_orf_id(panel) if (has_orf_id or has_coord_key) else panel
    duplicate_orf_ids = 0
    if "orf_id" in keyed.columns:
        duplicate_orf_ids = keyed.height - keyed.select("orf_id").unique().height
        if duplicate_orf_ids:
            errors.append(f"manifest has {duplicate_orf_ids} duplicate orf_id values")

    return pl.DataFrame([{
        "n_rows": panel.height,
        "n_columns": len(panel.columns),
        "has_orf_id": has_orf_id,
        "has_coordinate_key": has_coord_key,
        "duplicate_orf_ids": duplicate_orf_ids,
        "is_valid": len(errors) == 0,
        "errors": "; ".join(errors) if errors else None,
        "warnings": "; ".join(warnings) if warnings else None,
    }])


def merge_panel_manifest(orfs: pl.DataFrame, panel: pl.DataFrame) -> pl.DataFrame:
    """Attach frozen panel metadata to an ORF table by ORF key."""
    panel_report = validate_panel_manifest(panel)
    if not bool(panel_report["is_valid"][0]):
        raise ValueError(panel_report["errors"][0])

    orfs_keyed = ensure_orf_id(orfs)
    panel_keyed = ensure_orf_id(panel)
    panel_cols = [c for c in panel_keyed.columns if c != "orf_id" and c not in orfs_keyed.columns]
    return orfs_keyed.join(panel_keyed.select(["orf_id", *panel_cols]), on="orf_id", how="left")


def panel_freeze_report(path: str) -> pl.DataFrame:
    panel = load_panel_manifest(path)
    report = validate_panel_manifest(panel)
    return report.with_columns(
        pl.lit(path).alias("path"),
        pl.lit(manifest_sha256(path)).alias("sha256"),
    )


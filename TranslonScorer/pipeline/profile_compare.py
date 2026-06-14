from __future__ import annotations

from typing import Any

import numpy as np
import polars as pl

from ..utils.io import write_csv_safe, write_parquet_safe


def _load_profiles(path: str) -> pl.DataFrame:
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    if path.endswith(".tsv"):
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def _normalize_profiles(df: pl.DataFrame, count_col: str = "count") -> pl.DataFrame:
    rename = {}
    if "tran_start" in df.columns and "pos" not in df.columns:
        rename["tran_start"] = "pos"
    if "counts" in df.columns and count_col not in df.columns:
        rename["counts"] = count_col
    if rename:
        df = df.rename(rename)
    required = {"tran_id", "pos", count_col}
    if not required.issubset(set(df.columns)):
        raise ValueError(f"Profile table missing required columns: {sorted(required - set(df.columns))}")
    return (
        df.select(["tran_id", "pos", count_col])
        .with_columns(pl.col("pos").cast(pl.Int64), pl.col(count_col).cast(pl.Float64).alias("count"))
        .group_by(["tran_id", "pos"])
        .agg(pl.col("count").sum().alias("count"))
        .sort(["tran_id", "pos"])
    )


def compare_profiles(
    profile_a: pl.DataFrame,
    profile_b: pl.DataFrame,
    *,
    label_a: str = "a",
    label_b: str = "b",
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Compare two transcript-space profile tables position by position."""
    a = _normalize_profiles(profile_a).rename({"count": f"count_{label_a}"})
    b = _normalize_profiles(profile_b).rename({"count": f"count_{label_b}"})

    joined = (
        a.join(b, on=["tran_id", "pos"], how="full", coalesce=True)
        .with_columns(
            pl.col(f"count_{label_a}").fill_null(0.0),
            pl.col(f"count_{label_b}").fill_null(0.0),
        )
        .with_columns(
            (pl.col(f"count_{label_b}") - pl.col(f"count_{label_a}")).alias("count_delta"),
            (pl.col(f"count_{label_b}") - pl.col(f"count_{label_a}")).abs().alias("abs_delta"),
        )
    )

    x = joined[f"count_{label_a}"].to_numpy()
    y = joined[f"count_{label_b}"].to_numpy()
    if x.size > 1 and float(np.std(x)) > 0 and float(np.std(y)) > 0:
        pearson = float(np.corrcoef(x, y)[0, 1])
    else:
        pearson = None

    total_a = float(np.sum(x)) if x.size else 0.0
    total_b = float(np.sum(y)) if y.size else 0.0
    abs_delta = joined["abs_delta"].to_numpy()
    rmse = float(np.sqrt(np.mean(np.square(joined["count_delta"].to_numpy())))) if joined.height else None

    summary: dict[str, Any] = {
        "label_a": label_a,
        "label_b": label_b,
        "n_positions_union": joined.height,
        "n_positions_a": a.height,
        "n_positions_b": b.height,
        "n_positions_shared_nonzero": int(((x > 0) & (y > 0)).sum()) if x.size else 0,
        "total_count_a": total_a,
        "total_count_b": total_b,
        "total_count_delta": total_b - total_a,
        "total_count_ratio_b_over_a": (total_b / total_a) if total_a > 0 else None,
        "mae": float(np.mean(abs_delta)) if abs_delta.size else None,
        "rmse": rmse,
        "pearson_r": pearson,
        "a_only_mass": float(joined.filter((pl.col(f"count_{label_a}") > 0) & (pl.col(f"count_{label_b}") == 0))[f"count_{label_a}"].sum() or 0.0),
        "b_only_mass": float(joined.filter((pl.col(f"count_{label_b}") > 0) & (pl.col(f"count_{label_a}") == 0))[f"count_{label_b}"].sum() or 0.0),
    }
    return joined, pl.DataFrame([summary])


def compare_profile_files(
    *,
    profile_a_path: str,
    profile_b_path: str,
    out_prefix: str,
    label_a: str = "bigwig",
    label_b: str = "bam",
    write_deltas: bool = False,
) -> dict[str, str]:
    joined, summary = compare_profiles(
        _load_profiles(profile_a_path),
        _load_profiles(profile_b_path),
        label_a=label_a,
        label_b=label_b,
    )
    paths = {"summary": f"{out_prefix}.profile_compare.summary.csv"}
    write_csv_safe(summary, paths["summary"])
    if write_deltas:
        paths["deltas"] = f"{out_prefix}.profile_compare.deltas.parquet"
        write_parquet_safe(joined, paths["deltas"])
    return paths


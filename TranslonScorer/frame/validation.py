"""Validate frame posteriors against annotated CDS frame.

Given frame-support tables built by ``frame_support.build_frame_support``, ask
the only question that can be answered without a second method to compare
against: on CDS interiors, where the true frame is known from the annotation,
does the posterior actually put its mass on that frame?
"""

from __future__ import annotations

import polars as pl


def _argmax_p_expr() -> pl.Expr:
    return (
        pl.when((pl.col("p0") >= pl.col("p1")) & (pl.col("p0") >= pl.col("p2")))
        .then(pl.lit(0))
        .when((pl.col("p1") >= pl.col("p0")) & (pl.col("p1") >= pl.col("p2")))
        .then(pl.lit(1))
        .otherwise(pl.lit(2))
    )


def _true_p_expr() -> pl.Expr:
    return (
        pl.when(pl.col("true_frame") == 0)
        .then(pl.col("p0"))
        .when(pl.col("true_frame") == 1)
        .then(pl.col("p1"))
        .otherwise(pl.col("p2"))
    )


def _cds_interiors(cds: pl.DataFrame, trim_nt: int) -> pl.DataFrame:
    if cds.is_empty():
        return pl.DataFrame(
            schema={
                "tran_id": pl.Utf8,
                "cds_start": pl.Int64,
                "cds_stop": pl.Int64,
                "true_frame": pl.Int64,
            }
        )
    start_col = "tran_start" if "tran_start" in cds.columns else "start"
    stop_col = "tran_stop" if "tran_stop" in cds.columns else "stop"
    required = {"tran_id", start_col, stop_col}
    missing = required - set(cds.columns)
    if missing:
        raise ValueError(f"CDS table missing required columns: {sorted(missing)}")
    return (
        cds.select(
            [
                pl.col("tran_id"),
                pl.col(start_col).cast(pl.Int64).alias("cds_start"),
                pl.col(stop_col).cast(pl.Int64).alias("cds_stop"),
            ]
        )
        .with_columns((pl.col("cds_stop") - pl.col("cds_start")).alias("cds_len"))
        .filter(pl.col("cds_len") > (2 * int(trim_nt) + 3))
        .with_columns(
            [
                (pl.col("cds_start") + int(trim_nt)).alias("interior_start"),
                (pl.col("cds_stop") - int(trim_nt)).alias("interior_stop"),
                (pl.col("cds_start") % 3).cast(pl.Int64).alias("true_frame"),
            ]
        )
        .select(["tran_id", "interior_start", "interior_stop", "true_frame"])
    )


def validate_frame_support_on_cds(
    frame_support: pl.DataFrame,
    cds: pl.DataFrame,
    *,
    method: str,
    trim_nt: int = 30,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Validate frame posteriors against annotated CDS interior frame."""
    interiors = _cds_interiors(cds, trim_nt)
    if frame_support.is_empty() or interiors.is_empty():
        empty_positions = pl.DataFrame()
        summary = pl.DataFrame(
            {
                "method": [method],
                "trim_nt": [trim_nt],
                "n_rows": [0],
                "n_transcripts": [0],
                "total_weight": [0.0],
                "mean_p_true": [None],
                "weighted_mean_p_true": [None],
                "argmax_accuracy": [None],
                "weighted_argmax_accuracy": [None],
                "mean_entropy": [None],
                "weighted_mean_entropy": [None],
                "mean_off_frame_mass": [None],
                "high_conf_wrong_rate": [None],
            }
        )
        return empty_positions, summary

    fs = (
        frame_support.with_columns(
            [
                (pl.col("codon").cast(pl.Int64) * 3).alias("codon_start"),
                (
                    pl.col("total_count").cast(pl.Float64)
                    if "total_count" in frame_support.columns
                    else pl.lit(1.0)
                ).alias("weight"),
            ]
        )
        .join(interiors, on="tran_id", how="inner")
        .filter(
            (pl.col("codon_start") >= pl.col("interior_start"))
            & (pl.col("codon_start") < pl.col("interior_stop"))
        )
        .with_columns(
            [
                _argmax_p_expr().alias("pred_frame"),
                _true_p_expr().alias("p_true"),
            ]
        )
        .with_columns(
            [
                (pl.col("pred_frame") == pl.col("true_frame")).alias("is_correct"),
                pl.max_horizontal([pl.col("p0"), pl.col("p1"), pl.col("p2")]).alias("p_max"),
                (1.0 - pl.col("p_true")).alias("off_frame_mass"),
            ]
        )
        .with_columns(((pl.col("p_max") >= 0.8) & (~pl.col("is_correct"))).alias("high_conf_wrong"))
    )

    if fs.is_empty():
        return validate_frame_support_on_cds(pl.DataFrame(), cds, method=method, trim_nt=trim_nt)

    summary = fs.select(
        [
            pl.lit(method).alias("method"),
            pl.lit(trim_nt).alias("trim_nt"),
            pl.len().alias("n_rows"),
            pl.col("tran_id").n_unique().alias("n_transcripts"),
            pl.col("weight").sum().alias("total_weight"),
            pl.col("p_true").mean().alias("mean_p_true"),
            ((pl.col("p_true") * pl.col("weight")).sum() / pl.col("weight").sum()).alias(
                "weighted_mean_p_true"
            ),
            pl.col("is_correct").cast(pl.Float64).mean().alias("argmax_accuracy"),
            (
                (pl.col("is_correct").cast(pl.Float64) * pl.col("weight")).sum()
                / pl.col("weight").sum()
            ).alias("weighted_argmax_accuracy"),
            (
                pl.col("entropy").mean().alias("mean_entropy")
                if "entropy" in fs.columns
                else pl.lit(None).cast(pl.Float64).alias("mean_entropy")
            ),
            (
                ((pl.col("entropy") * pl.col("weight")).sum() / pl.col("weight").sum()).alias(
                    "weighted_mean_entropy"
                )
                if "entropy" in fs.columns
                else pl.lit(None).cast(pl.Float64).alias("weighted_mean_entropy")
            ),
            pl.col("off_frame_mass").mean().alias("mean_off_frame_mass"),
            pl.col("high_conf_wrong").cast(pl.Float64).mean().alias("high_conf_wrong_rate"),
        ]
    )
    return fs, summary

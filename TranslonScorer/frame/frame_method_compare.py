from __future__ import annotations

from typing import Sequence

import polars as pl

from ..frame_support import build_frame_support
from ..model import FrameSupportParams
from ..utils.io import write_csv_safe, write_parquet_safe

SUPPORTED_FRAME_METHODS = {"linear", "linear+hmm", "deblur+linear+hmm", "latent"}


def _load_table(path: str) -> pl.DataFrame:
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    if path.endswith(".tsv"):
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def _method_slug(method: str) -> str:
    return method.replace("+", "_plus_").replace("/", "_")


def _argmax_expr(prefix: str) -> pl.Expr:
    p0 = pl.col(f"{prefix}_p0")
    p1 = pl.col(f"{prefix}_p1")
    p2 = pl.col(f"{prefix}_p2")
    return (
        pl.when((p0 >= p1) & (p0 >= p2))
        .then(pl.lit(0))
        .when((p1 >= p0) & (p1 >= p2))
        .then(pl.lit(1))
        .otherwise(pl.lit(2))
    )


def _pmax_expr(prefix: str) -> pl.Expr:
    return pl.max_horizontal(
        [pl.col(f"{prefix}_p0"), pl.col(f"{prefix}_p1"), pl.col(f"{prefix}_p2")]
    )


def _pmin_expr(prefix: str) -> pl.Expr:
    return pl.min_horizontal(
        [pl.col(f"{prefix}_p0"), pl.col(f"{prefix}_p1"), pl.col(f"{prefix}_p2")]
    )


def _secondary_mass_expr(prefix: str) -> pl.Expr:
    pmax = _pmax_expr(prefix)
    pmin = _pmin_expr(prefix)
    return pl.col(f"{prefix}_p0") + pl.col(f"{prefix}_p1") + pl.col(f"{prefix}_p2") - pmax - pmin


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


def build_frame_support_for_method(
    profiles: pl.DataFrame,
    cds: pl.DataFrame,
    method: str,
    *,
    frame_by_length: bool = True,
    hmm_lambda: float = 2.0,
    background: str = "flat",
) -> pl.DataFrame:
    method = method.lower()
    if method not in SUPPORTED_FRAME_METHODS:
        raise ValueError(f"Unsupported frame method: {method}")
    params = FrameSupportParams(
        frame_method=method,
        frame_by_length=frame_by_length,
        frame_hmm_lambda=hmm_lambda,
        frame_background=background,
    )
    return build_frame_support(profiles, cds, params)


def compare_frame_support_tables(
    left: pl.DataFrame,
    right: pl.DataFrame,
    *,
    label_left: str = "linear",
    label_right: str = "latent",
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Compare two frame-support tables on shared transcript/codon rows."""
    left_keys = {"tran_id", "codon", "length"} & set(left.columns)
    right_keys = {"tran_id", "codon", "length"} & set(right.columns)
    join_keys = ["tran_id", "codon"] + (
        ["length"] if "length" in left_keys and "length" in right_keys else []
    )

    value_cols = [
        "observed_f0",
        "observed_f1",
        "observed_f2",
        "adjusted_f0",
        "adjusted_f1",
        "adjusted_f2",
        "total_count",
        "p0",
        "p1",
        "p2",
        "entropy",
        "p_max",
        "secondary_frame_mass",
        "frame_periodicity_score",
        "depth_score",
        "support_evidence",
        "method",
    ]
    left_keep = join_keys + [c for c in value_cols if c in left.columns]
    right_keep = join_keys + [c for c in value_cols if c in right.columns]
    left_prepped = left.select(left_keep).rename(
        {c: f"{label_left}_{c}" for c in left_keep if c not in join_keys}
    )
    right_prepped = right.select(right_keep).rename(
        {c: f"{label_right}_{c}" for c in right_keep if c not in join_keys}
    )

    joined = left_prepped.join(right_prepped, on=join_keys, how="inner")
    if joined.is_empty():
        summary = pl.DataFrame(
            {
                "method_left": [label_left],
                "method_right": [label_right],
                "n_joined": [0],
                "argmax_agreement": [None],
                "mean_l1_posterior_delta": [None],
                "mean_max_abs_posterior_delta": [None],
                "mean_entropy_left": [None],
                "mean_entropy_right": [None],
                "mean_entropy_delta": [None],
                "mean_abs_adjusted_count_delta": [None],
                "mean_p_max_left": [None],
                "mean_p_max_right": [None],
                "mean_secondary_frame_mass_left": [None],
                "mean_secondary_frame_mass_right": [None],
                "mean_secondary_frame_mass_delta": [None],
                "hmm_overclean_candidate_rate": [None],
            }
        )
        return joined, summary

    for label in (label_left, label_right):
        total_col = f"{label}_total_count"
        if total_col not in joined.columns:
            joined = joined.with_columns(pl.lit(0.0).alias(total_col))

    joined = (
        joined.with_columns(
            [
                _argmax_expr(label_left).alias(f"{label_left}_argmax"),
                _argmax_expr(label_right).alias(f"{label_right}_argmax"),
                _pmax_expr(label_left).alias(f"{label_left}_p_max"),
                _pmax_expr(label_right).alias(f"{label_right}_p_max"),
                _secondary_mass_expr(label_left).alias(f"{label_left}_secondary_frame_mass"),
                _secondary_mass_expr(label_right).alias(f"{label_right}_secondary_frame_mass"),
            ]
        )
        .with_columns(
            [
                (pl.col(f"{label_right}_p0") - pl.col(f"{label_left}_p0")).alias("p0_delta"),
                (pl.col(f"{label_right}_p1") - pl.col(f"{label_left}_p1")).alias("p1_delta"),
                (pl.col(f"{label_right}_p2") - pl.col(f"{label_left}_p2")).alias("p2_delta"),
                (pl.col(f"{label_right}_entropy") - pl.col(f"{label_left}_entropy")).alias(
                    "entropy_delta"
                ),
                (
                    pl.col(f"{label_right}_secondary_frame_mass")
                    - pl.col(f"{label_left}_secondary_frame_mass")
                ).alias("secondary_frame_mass_delta"),
                (pl.col(f"{label_left}_argmax") == pl.col(f"{label_right}_argmax")).alias(
                    "same_argmax"
                ),
            ]
        )
        .with_columns(
            [
                (
                    pl.col("p0_delta").abs() + pl.col("p1_delta").abs() + pl.col("p2_delta").abs()
                ).alias("l1_posterior_delta"),
                pl.max_horizontal(
                    [
                        pl.col("p0_delta").abs(),
                        pl.col("p1_delta").abs(),
                        pl.col("p2_delta").abs(),
                    ]
                ).alias("max_abs_posterior_delta"),
                (
                    (pl.col(f"{label_left}_total_count").fill_null(0.0) >= 10.0)
                    & (pl.col(f"{label_left}_secondary_frame_mass") >= 0.20)
                    & (pl.col(f"{label_right}_secondary_frame_mass") < 0.10)
                    & (pl.col(f"{label_right}_entropy") < pl.col(f"{label_left}_entropy"))
                ).alias("hmm_overclean_candidate"),
            ]
        )
    )

    if (
        f"{label_left}_adjusted_f0" in joined.columns
        and f"{label_right}_adjusted_f0" in joined.columns
    ):
        joined = joined.with_columns(
            (
                (pl.col(f"{label_right}_adjusted_f0") - pl.col(f"{label_left}_adjusted_f0")).abs()
                + (pl.col(f"{label_right}_adjusted_f1") - pl.col(f"{label_left}_adjusted_f1")).abs()
                + (pl.col(f"{label_right}_adjusted_f2") - pl.col(f"{label_left}_adjusted_f2")).abs()
            ).alias("abs_adjusted_count_delta")
        )
    else:
        joined = joined.with_columns(
            pl.lit(None).cast(pl.Float64).alias("abs_adjusted_count_delta")
        )

    summary = joined.select(
        [
            pl.lit(label_left).alias("method_left"),
            pl.lit(label_right).alias("method_right"),
            pl.len().alias("n_joined"),
            pl.col("same_argmax").cast(pl.Float64).mean().alias("argmax_agreement"),
            pl.col("l1_posterior_delta").mean().alias("mean_l1_posterior_delta"),
            pl.col("max_abs_posterior_delta").mean().alias("mean_max_abs_posterior_delta"),
            pl.col(f"{label_left}_entropy").mean().alias("mean_entropy_left"),
            pl.col(f"{label_right}_entropy").mean().alias("mean_entropy_right"),
            pl.col("entropy_delta").mean().alias("mean_entropy_delta"),
            pl.col("abs_adjusted_count_delta").mean().alias("mean_abs_adjusted_count_delta"),
            pl.col(f"{label_left}_p_max").mean().alias("mean_p_max_left"),
            pl.col(f"{label_right}_p_max").mean().alias("mean_p_max_right"),
            pl.col(f"{label_left}_secondary_frame_mass")
            .mean()
            .alias("mean_secondary_frame_mass_left"),
            pl.col(f"{label_right}_secondary_frame_mass")
            .mean()
            .alias("mean_secondary_frame_mass_right"),
            pl.col("secondary_frame_mass_delta").mean().alias("mean_secondary_frame_mass_delta"),
            pl.col("hmm_overclean_candidate")
            .cast(pl.Float64)
            .mean()
            .alias("hmm_overclean_candidate_rate"),
        ]
    )
    return joined, summary


def compare_frame_methods(
    profiles_path: str,
    cds_path: str,
    out_prefix: str,
    *,
    methods: Sequence[str] = ("linear", "latent"),
    frame_by_length: bool = True,
    hmm_lambda: float = 2.0,
    background: str = "flat",
    trim_nt: int = 30,
    write_validation_rows: bool = False,
) -> dict[str, str]:
    """Build and compare frame-support methods from the same profile/CDS inputs."""
    if len(methods) < 2:
        raise ValueError("At least two methods are required for comparison")

    profiles = _load_table(profiles_path)
    cds = _load_table(cds_path)

    supports: dict[str, pl.DataFrame] = {}
    paths: dict[str, str] = {}
    for method in methods:
        method = method.lower()
        support = build_frame_support_for_method(
            profiles,
            cds,
            method,
            frame_by_length=frame_by_length,
            hmm_lambda=hmm_lambda,
            background=background,
        )
        supports[method] = support
        key = f"{_method_slug(method)}_support"
        paths[key] = write_parquet_safe(
            support, f"{out_prefix}.{_method_slug(method)}.frame_support.parquet"
        )

    validation_summaries: list[pl.DataFrame] = []
    validation_paths: list[str] = []
    for method, support in supports.items():
        validation_rows, validation_summary = validate_frame_support_on_cds(
            support,
            cds,
            method=_method_slug(method),
            trim_nt=trim_nt,
        )
        validation_summaries.append(validation_summary)
        if write_validation_rows:
            path = f"{out_prefix}.{_method_slug(method)}.frame_validation.rows.parquet"
            write_parquet_safe(validation_rows, path)
            validation_paths.append(path)

    left_method = methods[0].lower()
    summaries: list[pl.DataFrame] = []
    comparison_paths: list[str] = []
    for right_method in [m.lower() for m in methods[1:]]:
        comparison, summary = compare_frame_support_tables(
            supports[left_method],
            supports[right_method],
            label_left=_method_slug(left_method),
            label_right=_method_slug(right_method),
        )
        compare_path = f"{out_prefix}.{_method_slug(left_method)}_vs_{_method_slug(right_method)}.frame_compare.parquet"
        write_parquet_safe(comparison, compare_path)
        comparison_paths.append(compare_path)
        summaries.append(summary)

    summary_df = pl.concat(summaries) if summaries else pl.DataFrame()
    paths["summary"] = write_csv_safe(summary_df, f"{out_prefix}.frame_compare.summary.csv")
    validation_df = pl.concat(validation_summaries) if validation_summaries else pl.DataFrame()
    paths["validation_summary"] = write_csv_safe(
        validation_df, f"{out_prefix}.frame_validation.summary.csv"
    )
    if validation_paths:
        paths["validation_rows"] = ",".join(validation_paths)
    paths["comparisons"] = ",".join(comparison_paths)
    return paths

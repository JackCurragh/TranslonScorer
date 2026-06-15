from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
import polars as pl

from ..utils.io import write_csv_safe, write_parquet_safe


SUPPORTED_ASSIGNMENT_METHODS = {
    "unique",
    "fractional",
    "frame_fractional",
    "em",
    "frame_em",
    "rdg_local_frame_em",
    "rdg_gated_frame_em",
}
_EPS = 1e-12


@dataclass(frozen=True)
class AssignmentResult:
    assignments: pl.DataFrame
    abundance: pl.DataFrame
    convergence: pl.DataFrame
    diagnostics: pl.DataFrame | None = None


def _load_table(path: str) -> pl.DataFrame:
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    if path.endswith(".tsv"):
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def _ensure_read_key(df: pl.DataFrame) -> pl.DataFrame:
    if "read_key" in df.columns:
        return df.with_columns(pl.col("read_key").cast(pl.Utf8))
    if "qname" in df.columns:
        return df.with_columns(pl.col("qname").cast(pl.Utf8).alias("read_key"))
    if "read_id" in df.columns:
        return df.with_columns(
            pl.concat_str([pl.lit("read"), pl.col("read_id").cast(pl.Utf8)], separator=":").alias(
                "read_key"
            )
        )
    coord_cols = [c for c in ["chr", "start", "stop", "length", "strand"] if c in df.columns]
    if coord_cols:
        return df.with_columns(
            pl.concat_str([pl.col(c).cast(pl.Utf8) for c in coord_cols], separator=":").alias(
                "read_key"
            )
        )
    return (
        df.with_row_index("_read_row")
        .with_columns(
            pl.concat_str([pl.lit("row"), pl.col("_read_row").cast(pl.Utf8)], separator=":").alias(
                "read_key"
            )
        )
        .drop("_read_row")
    )


def _position_column(df: pl.DataFrame) -> str:
    for col in ("tran_start_bam", "tran_start", "pos", "start_pos_tran"):
        if col in df.columns:
            return col
    raise ValueError("Candidate table needs one transcript-position column")


def _normalise_frame_support(frame_support: pl.DataFrame | None) -> pl.DataFrame | None:
    if frame_support is None or frame_support.is_empty():
        return None
    df = frame_support
    if "transcript_id" in df.columns and "tran_id" not in df.columns:
        df = df.rename({"transcript_id": "tran_id"})
    required = {"tran_id", "codon", "p0", "p1", "p2"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"frame_support missing required columns: {sorted(missing)}")
    return (
        df.select(
            [
                pl.col("tran_id").cast(pl.Utf8),
                pl.col("codon").cast(pl.Int64),
                pl.col("p0").cast(pl.Float64),
                pl.col("p1").cast(pl.Float64),
                pl.col("p2").cast(pl.Float64),
                (
                    pl.col("support_evidence").cast(pl.Float64)
                    if "support_evidence" in df.columns
                    else pl.lit(1.0).alias("support_evidence")
                ),
                (
                    pl.col("total_count").cast(pl.Float64).alias("frame_support_count")
                    if "total_count" in df.columns
                    else pl.lit(None).cast(pl.Float64).alias("frame_support_count")
                ),
            ]
        )
        .group_by(["tran_id", "codon"])
        .agg(
            [
                pl.col("p0").mean(),
                pl.col("p1").mean(),
                pl.col("p2").mean(),
                pl.col("support_evidence").mean(),
                pl.col("frame_support_count").sum(),
            ]
        )
    )


def _signal_frame_expr() -> pl.Expr:
    return (pl.col("pos").cast(pl.Int64) % 3).cast(pl.Int64)


def _cds_relative_frame_expr() -> pl.Expr:
    return (
        pl.when(pl.col("cds_start").is_not_null())
        .then((pl.col("pos").cast(pl.Int64) - pl.col("cds_start").cast(pl.Int64)) % 3)
        .otherwise(pl.col("pos").cast(pl.Int64) % 3)
        .cast(pl.Int64)
    )


def prepare_assignment_candidates(
    candidates: pl.DataFrame,
    *,
    frame_support: pl.DataFrame | None = None,
    cds_tran: pl.DataFrame | None = None,
    abundance_key: str = "tran_id",
    frame_floor: float = 0.02,
    frame_weight: float = 1.0,
    min_frame_support_count: float = 0.0,
    min_frame_support_evidence: float = 0.0,
    require_complete_frame_support: bool = False,
) -> pl.DataFrame:
    """Normalise read-to-origin candidates and attach assignment likelihoods.

    Required candidate fields are a read identity and `tran_id` plus one
    transcript position column. Optional `locus_id`/`gene_id` columns let the
    same machinery operate at genomic-locus level by setting `abundance_key`.
    """
    if candidates.is_empty():
        return pl.DataFrame()
    if not 0.0 <= frame_floor <= 1.0:
        raise ValueError("frame_floor must be between 0 and 1")
    if min_frame_support_count < 0.0:
        raise ValueError("min_frame_support_count must be non-negative")
    if not 0.0 <= min_frame_support_evidence <= 1.0:
        raise ValueError("min_frame_support_evidence must be between 0 and 1")

    pos_col = _position_column(candidates)
    df = (
        _ensure_read_key(candidates)
        .with_row_index("_candidate_row")
        .with_columns(
            [
                pl.col("tran_id").cast(pl.Utf8),
                pl.col(pos_col).cast(pl.Int64).alias("pos"),
                (
                    pl.col("count").cast(pl.Float64)
                    if "count" in candidates.columns
                    else pl.lit(1.0)
                ).alias("read_weight"),
            ]
        )
    )

    target_col = abundance_key
    if "," in abundance_key:
        target_cols = [col.strip() for col in abundance_key.split(",") if col.strip()]
        missing = [col for col in target_cols if col not in df.columns]
        if missing:
            raise ValueError(f"Candidate table missing abundance_key columns: {missing}")
        target_col = "_assignment_target_key"
        df = df.with_columns(
            pl.concat_str([pl.col(col).cast(pl.Utf8) for col in target_cols], separator="|").alias(
                target_col
            )
        )
    elif abundance_key not in df.columns:
        if abundance_key == "locus_id":
            if "gene_id" in df.columns:
                df = df.with_columns(pl.col("gene_id").cast(pl.Utf8).alias("locus_id"))
            else:
                df = df.with_columns(pl.col("tran_id").alias("locus_id"))
        else:
            raise ValueError(f"Candidate table missing abundance_key column: {abundance_key}")

    if "candidate_likelihood" in df.columns:
        df = df.with_columns(
            pl.col("candidate_likelihood")
            .cast(pl.Float64)
            .clip(_EPS, None)
            .alias("alignment_likelihood")
        )
    elif "alignment_likelihood" in df.columns:
        df = df.with_columns(
            pl.col("alignment_likelihood")
            .cast(pl.Float64)
            .clip(_EPS, None)
            .alias("alignment_likelihood")
        )
    elif "mapq" in df.columns:
        df = df.with_columns(
            (1.0 - (((-pl.col("mapq").cast(pl.Float64) / 10.0) * math.log(10.0)).exp()))
            .clip(_EPS, 1.0)
            .alias("alignment_likelihood")
        )
    else:
        df = df.with_columns(pl.lit(1.0).alias("alignment_likelihood"))

    if cds_tran is not None and not cds_tran.is_empty():
        start_col = "tran_start" if "tran_start" in cds_tran.columns else "start"
        cds = cds_tran.select(
            ["tran_id", pl.col(start_col).cast(pl.Int64).alias("cds_start")]
        ).unique(subset=["tran_id"])
        df = df.join(cds, on="tran_id", how="left")
    else:
        df = df.with_columns(pl.lit(None).cast(pl.Int64).alias("cds_start"))

    df = df.with_columns(
        [
            (pl.col("pos") // 3).alias("codon"),
            _signal_frame_expr().alias("signal_frame"),
            _cds_relative_frame_expr().alias("cds_relative_frame"),
            pl.col(target_col).cast(pl.Utf8).alias("assignment_target"),
        ]
    )

    fs = _normalise_frame_support(frame_support)
    if fs is not None:
        df = df.join(fs, on=["tran_id", "codon"], how="left")
    else:
        df = df.with_columns(
            [
                pl.lit(None).cast(pl.Float64).alias("p0"),
                pl.lit(None).cast(pl.Float64).alias("p1"),
                pl.lit(None).cast(pl.Float64).alias("p2"),
                pl.lit(None).cast(pl.Float64).alias("support_evidence"),
                pl.lit(None).cast(pl.Float64).alias("frame_support_count"),
            ]
        )

    df = (
        df.with_columns(
            [
                pl.col("support_evidence")
                .fill_null(0.0)
                .clip(0.0, 1.0)
                .alias("raw_support_evidence"),
                pl.col("frame_support_count")
                .fill_null(0.0)
                .clip(0.0, None)
                .alias("frame_support_count"),
                pl.when(pl.col("signal_frame") == 0)
                .then(pl.col("p0"))
                .when(pl.col("signal_frame") == 1)
                .then(pl.col("p1"))
                .otherwise(pl.col("p2"))
                .fill_null(1.0)
                .clip(frame_floor, 1.0)
                .alias("raw_frame_likelihood"),
            ]
        )
        .with_columns(
            (
                (pl.col("frame_support_count") >= float(min_frame_support_count))
                & (pl.col("raw_support_evidence") >= float(min_frame_support_evidence))
            ).alias("frame_support_passes_gate")
        )
        .with_columns(
            pl.when(pl.col("frame_support_passes_gate"))
            .then(pl.col("raw_support_evidence"))
            .otherwise(0.0)
            .alias("support_evidence")
        )
        .with_columns(
            (
                (1.0 - pl.col("support_evidence"))
                + (pl.col("support_evidence") * pl.col("raw_frame_likelihood"))
            )
            .clip(frame_floor, 1.0)
            .alias("frame_likelihood")
        )
    )

    if require_complete_frame_support:
        target_support = (
            df.group_by(["read_key", "assignment_target"])
            .agg(pl.col("frame_support_passes_gate").any().alias("_target_has_frame_support"))
            .group_by("read_key")
            .agg(
                [
                    pl.len().alias("_n_assignment_targets_for_frame"),
                    pl.col("_target_has_frame_support")
                    .cast(pl.Int64)
                    .sum()
                    .alias("_n_supported_assignment_targets_for_frame"),
                ]
            )
            .with_columns(
                (
                    pl.col("_n_supported_assignment_targets_for_frame")
                    == pl.col("_n_assignment_targets_for_frame")
                ).alias("frame_support_comparable")
            )
        )
        df = (
            df.join(target_support, on="read_key", how="left")
            .with_columns(pl.col("frame_support_comparable").fill_null(False))
            .with_columns(
                [
                    pl.when(pl.col("frame_support_comparable"))
                    .then(pl.col("frame_likelihood"))
                    .otherwise(1.0)
                    .alias("frame_likelihood"),
                    pl.when(pl.col("frame_support_comparable"))
                    .then(pl.col("support_evidence"))
                    .otherwise(0.0)
                    .alias("support_evidence"),
                    (
                        pl.col("frame_support_passes_gate") & pl.col("frame_support_comparable")
                    ).alias("frame_support_passes_gate"),
                ]
            )
        )
    else:
        df = df.with_columns(pl.lit(True).alias("frame_support_comparable"))

    if frame_weight != 1.0:
        df = df.with_columns(
            (pl.col("frame_likelihood") ** float(frame_weight)).alias("frame_likelihood")
        )

    return df.with_columns(
        [
            pl.col("alignment_likelihood").clip(_EPS, None).alias("alignment_likelihood"),
            (pl.col("alignment_likelihood") * pl.col("frame_likelihood"))
            .clip(_EPS, None)
            .alias("frame_aware_likelihood"),
            pl.col("read_weight").fill_null(1.0).clip(0.0, None).alias("read_weight"),
        ]
    )


def _factorize(values: Sequence[object]) -> tuple[np.ndarray, list[str]]:
    ids: dict[str, int] = {}
    labels: list[str] = []
    out = np.empty(len(values), dtype=np.int64)
    for i, value in enumerate(values):
        label = str(value)
        idx = ids.get(label)
        if idx is None:
            idx = len(labels)
            ids[label] = idx
            labels.append(label)
        out[i] = idx
    return out, labels


def _read_weights(read_idx: np.ndarray, weights: np.ndarray, n_reads: int) -> np.ndarray:
    out = np.zeros(n_reads, dtype=np.float64)
    seen = np.zeros(n_reads, dtype=bool)
    for idx, weight in zip(read_idx, weights):
        if not seen[idx]:
            out[idx] = float(weight)
            seen[idx] = True
    return out


def _normalise_within_reads(
    read_idx: np.ndarray, likelihood: np.ndarray, n_reads: int
) -> np.ndarray:
    denom = np.bincount(read_idx, weights=likelihood, minlength=n_reads)
    posterior = np.zeros_like(likelihood, dtype=np.float64)
    valid = denom[read_idx] > _EPS
    posterior[valid] = likelihood[valid] / denom[read_idx][valid]
    return posterior


def _target_likelihood_table(prepared: pl.DataFrame, likelihood_col: str) -> pl.DataFrame:
    return prepared.group_by(["read_key", "assignment_target"], maintain_order=True).agg(
        [
            pl.col("read_weight").first().alias("read_weight"),
            pl.col(likelihood_col).max().clip(_EPS, None).alias("target_likelihood"),
            pl.len().alias("candidate_origin_rows"),
        ]
    )


def _weighted_fraction_expr(flag_col: str, total_weight: float) -> pl.Expr:
    return (
        (pl.col(flag_col).cast(pl.Float64) * pl.col("read_weight")).sum() / total_weight
    ).fill_nan(0.0)


def frame_assignment_gate_diagnostics(
    prepared: pl.DataFrame,
    *,
    frame_gate_min_likelihood_range: float = 0.25,
    frame_gate_min_read_fraction: float = 0.25,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Diagnose when frame support is structurally useful for assignment.

    The gate is RDG-inspired but read-local: a read's frame support can affect
    assignment only when candidate origins have comparable frame support, the
    frame likelihood range is meaningful, and exactly one assignment target has
    the best frame compatibility. A stricter locus/sample consensus flag then
    asks whether enough read-weighted evidence passes the local gate to let EM
    use those frame likelihoods at all.
    """
    if prepared.is_empty():
        empty_gate = pl.DataFrame()
        empty_summary = pl.DataFrame(
            schema={
                "n_reads": pl.Int64,
                "total_read_weight": pl.Float64,
                "rdg_gate_fraction": pl.Float64,
                "frame_contrast_fraction": pl.Float64,
                "alias_group_contrast_fraction": pl.Float64,
                "mean_assignment_targets": pl.Float64,
                "mean_candidate_frames": pl.Float64,
                "rdg_scenario_frame_gate": pl.Boolean,
                "frame_gate_min_likelihood_range": pl.Float64,
                "frame_gate_min_read_fraction": pl.Float64,
            }
        )
        return empty_gate, empty_summary
    if frame_gate_min_likelihood_range < 0.0:
        raise ValueError("frame_gate_min_likelihood_range must be non-negative")
    if not 0.0 <= frame_gate_min_read_fraction <= 1.0:
        raise ValueError("frame_gate_min_read_fraction must be between 0 and 1")

    target = (
        prepared.group_by(["read_key", "assignment_target"])
        .agg(
            [
                pl.col("read_weight").first().alias("read_weight"),
                pl.col("frame_likelihood").max().alias("target_frame_likelihood"),
                pl.col("signal_frame").first().alias("representative_frame"),
            ]
        )
        .with_columns(
            pl.col("target_frame_likelihood")
            .max()
            .over("read_key")
            .alias("_max_target_frame_likelihood")
        )
        .with_columns(
            (pl.col("target_frame_likelihood") == pl.col("_max_target_frame_likelihood")).alias(
                "is_top_frame_target"
            )
        )
    )
    gate = (
        target.group_by("read_key")
        .agg(
            [
                pl.col("read_weight").first().alias("read_weight"),
                pl.len().alias("n_assignment_targets"),
                pl.col("representative_frame").n_unique().alias("n_candidate_frames"),
                pl.col("target_frame_likelihood").max().alias("max_frame_likelihood"),
                pl.col("target_frame_likelihood").min().alias("min_frame_likelihood"),
                pl.col("is_top_frame_target").sum().alias("n_top_frame_targets"),
            ]
        )
        .with_columns(
            [
                (pl.col("max_frame_likelihood") - pl.col("min_frame_likelihood")).alias(
                    "frame_likelihood_range"
                ),
                (
                    (pl.col("max_frame_likelihood") - pl.col("min_frame_likelihood"))
                    >= float(frame_gate_min_likelihood_range)
                ).alias("has_frame_contrast"),
            ]
        )
        .with_columns(
            [
                (pl.col("has_frame_contrast") & (pl.col("n_top_frame_targets") == 1)).alias(
                    "rdg_local_frame_gate"
                ),
                (pl.col("has_frame_contrast") & (pl.col("n_top_frame_targets") > 1)).alias(
                    "frame_separates_only_an_alias_group"
                ),
            ]
        )
    )

    total_weight = float(gate["read_weight"].sum() or 0.0)
    if total_weight <= 0.0:
        total_weight = 1.0
    summary = (
        gate.select(
            [
                pl.len().alias("n_reads"),
                pl.col("read_weight").sum().alias("total_read_weight"),
                _weighted_fraction_expr("rdg_local_frame_gate", total_weight).alias(
                    "rdg_gate_fraction"
                ),
                _weighted_fraction_expr("has_frame_contrast", total_weight).alias(
                    "frame_contrast_fraction"
                ),
                _weighted_fraction_expr("frame_separates_only_an_alias_group", total_weight).alias(
                    "alias_group_contrast_fraction"
                ),
                pl.col("n_assignment_targets").mean().alias("mean_assignment_targets"),
                pl.col("n_candidate_frames").mean().alias("mean_candidate_frames"),
            ]
        )
        .with_columns(
            (
                (pl.col("rdg_gate_fraction") >= float(frame_gate_min_read_fraction))
                & (pl.col("rdg_gate_fraction") > pl.col("alias_group_contrast_fraction"))
            ).alias("rdg_scenario_frame_gate")
        )
        .with_columns(
            [
                pl.lit(float(frame_gate_min_likelihood_range)).alias(
                    "frame_gate_min_likelihood_range"
                ),
                pl.lit(float(frame_gate_min_read_fraction)).alias("frame_gate_min_read_fraction"),
            ]
        )
    )
    scenario_gate = bool(summary["rdg_scenario_frame_gate"][0])
    gate = gate.with_columns(pl.lit(scenario_gate).alias("rdg_scenario_frame_gate"))
    return gate, summary


def _apply_rdg_frame_gate(
    prepared: pl.DataFrame,
    *,
    use_scenario_consensus: bool,
    frame_gate_min_likelihood_range: float,
    frame_gate_min_read_fraction: float,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    gate, summary = frame_assignment_gate_diagnostics(
        prepared,
        frame_gate_min_likelihood_range=frame_gate_min_likelihood_range,
        frame_gate_min_read_fraction=frame_gate_min_read_fraction,
    )
    if prepared.is_empty():
        return prepared, summary

    gated = prepared.join(
        gate.select(["read_key", "rdg_local_frame_gate", "rdg_scenario_frame_gate"]),
        on="read_key",
        how="left",
    ).with_columns(
        [
            pl.col("rdg_local_frame_gate").fill_null(False),
            pl.col("rdg_scenario_frame_gate").fill_null(False),
        ]
    )
    effective_gate = pl.col("rdg_local_frame_gate")
    if use_scenario_consensus:
        effective_gate = effective_gate & pl.col("rdg_scenario_frame_gate")
    gated = gated.with_columns(effective_gate.alias("rdg_frame_gate_applied"))
    return (
        gated.with_columns(
            [
                pl.when(pl.col("rdg_frame_gate_applied"))
                .then(pl.col("frame_likelihood"))
                .otherwise(1.0)
                .alias("rdg_gated_frame_likelihood"),
            ]
        ).with_columns(
            (pl.col("alignment_likelihood") * pl.col("rdg_gated_frame_likelihood"))
            .clip(_EPS, None)
            .alias("rdg_gated_frame_aware_likelihood")
        ),
        summary,
    )


def _candidate_posteriors_from_targets(
    prepared: pl.DataFrame,
    target_posteriors: pl.DataFrame,
    *,
    likelihood_col: str,
) -> pl.DataFrame:
    target_cols = ["read_key", "assignment_target", "target_posterior", "target_likelihood"]
    return (
        prepared.join(
            target_posteriors.select(target_cols), on=["read_key", "assignment_target"], how="left"
        )
        .with_columns(
            pl.col(likelihood_col)
            .sum()
            .over(["read_key", "assignment_target"])
            .clip(_EPS, None)
            .alias("_origin_likelihood_sum")
        )
        .with_columns(
            (
                pl.col(likelihood_col)
                / pl.col("_origin_likelihood_sum")
                * pl.col("target_posterior").fill_null(0.0)
            )
            .fill_nan(0.0)
            .fill_null(0.0)
            .alias("posterior")
        )
    )


def _run_em(
    prepared: pl.DataFrame,
    *,
    likelihood_col: str,
    max_iter: int,
    tol: float,
    abundance_prior: float,
) -> tuple[np.ndarray, np.ndarray, pl.DataFrame]:
    read_idx, _read_labels = _factorize(prepared["read_key"].to_list())
    target_idx, target_labels = _factorize(prepared["assignment_target"].to_list())
    n_reads = int(read_idx.max()) + 1 if read_idx.size else 0
    n_targets = int(target_idx.max()) + 1 if target_idx.size else 0

    likelihood = prepared[likelihood_col].to_numpy().astype(np.float64)
    likelihood[likelihood < _EPS] = _EPS
    weights = _read_weights(
        read_idx, prepared["read_weight"].to_numpy().astype(np.float64), n_reads
    )

    theta = np.bincount(target_idx, weights=likelihood, minlength=n_targets) + float(
        abundance_prior
    )
    theta = theta / max(float(theta.sum()), _EPS)
    posterior = np.zeros_like(likelihood)
    rows: list[dict[str, float | int]] = []

    for iteration in range(int(max_iter)):
        score = theta[target_idx] * likelihood
        posterior = _normalise_within_reads(read_idx, score, n_reads)
        assigned = posterior * weights[read_idx]
        theta_new = np.bincount(target_idx, weights=assigned, minlength=n_targets) + float(
            abundance_prior
        )
        theta_new = theta_new / max(float(theta_new.sum()), _EPS)
        delta = float(np.abs(theta_new - theta).sum())
        rows.append({"iteration": iteration + 1, "theta_l1_delta": delta})
        theta = theta_new
        if delta < tol:
            break

    convergence = (
        pl.DataFrame(rows)
        if rows
        else pl.DataFrame(schema={"iteration": pl.Int64, "theta_l1_delta": pl.Float64})
    )
    abundance = pl.DataFrame({"assignment_target": target_labels, "abundance": theta})
    return posterior, theta, convergence.with_columns(pl.lit("em").alias("estimator"))


def assign_reads(
    candidates: pl.DataFrame,
    *,
    method: str = "frame_em",
    frame_support: pl.DataFrame | None = None,
    cds_tran: pl.DataFrame | None = None,
    abundance_key: str = "tran_id",
    max_iter: int = 100,
    tol: float = 1e-7,
    abundance_prior: float = 1e-3,
    frame_floor: float = 0.02,
    frame_weight: float = 1.0,
    min_frame_support_count: float = 0.0,
    min_frame_support_evidence: float = 0.0,
    require_complete_frame_support: bool = False,
    frame_gate_min_likelihood_range: float = 0.25,
    frame_gate_min_read_fraction: float = 0.25,
) -> AssignmentResult:
    """Assign reads to transcript or locus targets.

    `abundance_key='tran_id'` gives isoform-level assignment. Use
    `abundance_key='locus_id'` or another locus column for genomic-origin
    disambiguation. Reads with one assignment target pass through with posterior
    1.0; only reads with multiple compatible targets can move between targets.
    """
    method = method.lower()
    if method not in SUPPORTED_ASSIGNMENT_METHODS:
        raise ValueError(f"Unsupported read assignment method: {method}")

    prepared = prepare_assignment_candidates(
        candidates,
        frame_support=frame_support,
        cds_tran=cds_tran,
        abundance_key=abundance_key,
        frame_floor=frame_floor,
        frame_weight=frame_weight,
        min_frame_support_count=min_frame_support_count,
        min_frame_support_evidence=min_frame_support_evidence,
        require_complete_frame_support=require_complete_frame_support,
    )
    if prepared.is_empty():
        empty = pl.DataFrame()
        return AssignmentResult(empty, empty, empty)

    diagnostics: pl.DataFrame | None = None
    if method in {"rdg_local_frame_em", "rdg_gated_frame_em"}:
        prepared, diagnostics = _apply_rdg_frame_gate(
            prepared,
            use_scenario_consensus=method == "rdg_gated_frame_em",
            frame_gate_min_likelihood_range=frame_gate_min_likelihood_range,
            frame_gate_min_read_fraction=frame_gate_min_read_fraction,
        )

    if method in {"frame_fractional", "frame_em"}:
        likelihood_col = "frame_aware_likelihood"
    elif method in {"rdg_local_frame_em", "rdg_gated_frame_em"}:
        likelihood_col = "rdg_gated_frame_aware_likelihood"
    else:
        likelihood_col = "alignment_likelihood"
    target_table = _target_likelihood_table(prepared, likelihood_col)
    read_idx, _ = _factorize(target_table["read_key"].to_list())
    n_reads = int(read_idx.max()) + 1 if read_idx.size else 0

    if method == "unique":
        target_counts = np.bincount(read_idx, minlength=n_reads)
        target_posterior = np.where(target_counts[read_idx] == 1, 1.0, 0.0)
        convergence = pl.DataFrame([{"iteration": 0, "theta_l1_delta": 0.0, "estimator": "unique"}])
    elif method in {"fractional", "frame_fractional"}:
        target_posterior = _normalise_within_reads(
            read_idx,
            target_table["target_likelihood"].to_numpy().astype(np.float64),
            n_reads,
        )
        convergence = pl.DataFrame([{"iteration": 0, "theta_l1_delta": 0.0, "estimator": method}])
    else:
        target_posterior, _theta, convergence = _run_em(
            target_table.rename({"target_likelihood": "_em_likelihood"}),
            likelihood_col="_em_likelihood",
            max_iter=max_iter,
            tol=tol,
            abundance_prior=abundance_prior,
        )

    target_posteriors = target_table.with_columns(pl.Series("target_posterior", target_posterior))
    out = _candidate_posteriors_from_targets(
        prepared, target_posteriors, likelihood_col=likelihood_col
    )
    read_idx_out, _ = _factorize(out["read_key"].to_list())
    n_reads_out = int(read_idx_out.max()) + 1 if read_idx_out.size else 0
    weights = _read_weights(
        read_idx_out, out["read_weight"].to_numpy().astype(np.float64), n_reads_out
    )
    out = out.with_columns(
        [
            pl.Series(
                "assigned_count",
                out["posterior"].to_numpy().astype(np.float64) * weights[read_idx_out],
            ),
            pl.lit(method).alias("assignment_method"),
            pl.lit(abundance_key).alias("abundance_key"),
        ]
    )
    abundance = (
        out.group_by("assignment_target")
        .agg(pl.col("assigned_count").sum().alias("assigned_count"))
        .with_columns(
            (pl.col("assigned_count") / pl.col("assigned_count").sum())
            .fill_nan(0.0)
            .alias("abundance")
        )
        .sort("assigned_count", descending=True)
    )
    return AssignmentResult(out, abundance, convergence, diagnostics)


def evaluate_assignments(
    assignments: pl.DataFrame, *, truth_column: str = "is_true"
) -> pl.DataFrame:
    if assignments.is_empty():
        return pl.DataFrame()

    has_truth = truth_column in assignments.columns
    df = assignments
    if has_truth:
        df = df.with_columns(pl.col(truth_column).cast(pl.Boolean).alias("_is_true"))

    target_aggs: list[pl.Expr] = [
        pl.col("read_weight").first().alias("read_weight"),
        pl.col("posterior").sum().alias("target_posterior"),
        pl.len().alias("candidate_origin_rows"),
    ]
    if has_truth:
        target_aggs.append(pl.col("_is_true").any().alias("_target_is_true"))

    per_target = df.group_by(["read_key", "assignment_target"]).agg(target_aggs)
    entropy_terms = (
        per_target.with_columns(
            pl.when(pl.col("target_posterior") > _EPS)
            .then(-(pl.col("target_posterior") * pl.col("target_posterior").log()))
            .otherwise(0.0)
            .alias("_entropy_term")
        )
        .group_by("read_key")
        .agg(pl.col("_entropy_term").sum().alias("assignment_entropy"))
    )
    per_read = (
        per_target.group_by("read_key")
        .agg(
            [
                pl.col("read_weight").first().alias("read_weight"),
                pl.len().alias("n_assignment_targets"),
                pl.col("candidate_origin_rows").sum().alias("n_candidate_origins"),
                pl.col("target_posterior").sum().alias("assigned_posterior_sum"),
                pl.col("target_posterior").max().alias("_max_posterior"),
            ]
        )
        .join(entropy_terms, on="read_key", how="left")
        .with_columns(
            pl.when(pl.col("n_assignment_targets") > 1)
            .then(
                pl.col("assignment_entropy") / pl.col("n_assignment_targets").cast(pl.Float64).log()
            )
            .otherwise(0.0)
            .fill_nan(0.0)
            .alias("assignment_entropy_norm")
        )
    )

    if has_truth:
        true_post = per_target.group_by("read_key").agg(
            pl.when(pl.col("_target_is_true"))
            .then(pl.col("target_posterior"))
            .otherwise(0.0)
            .sum()
            .alias("true_posterior")
        )
        hard = (
            per_target.join(
                per_read.select(["read_key", "_max_posterior"]), on="read_key", how="left"
            )
            .filter(pl.col("target_posterior") == pl.col("_max_posterior"))
            .group_by("read_key")
            .agg(
                [
                    pl.len().alias("_n_tied_max"),
                    pl.col("_target_is_true").cast(pl.Float64).sum().alias("_n_true_tied_max"),
                ]
            )
            .with_columns(
                (pl.col("_n_true_tied_max") / pl.col("_n_tied_max")).alias("hard_correct")
            )
            .select(["read_key", "hard_correct"])
        )
        per_read = (
            per_read.join(true_post, on="read_key", how="left")
            .join(hard, on="read_key", how="left")
            .with_columns(
                [
                    pl.col("true_posterior").fill_null(0.0),
                    pl.when(pl.col("assigned_posterior_sum") > _EPS)
                    .then(pl.col("hard_correct"))
                    .otherwise(0.0)
                    .fill_null(0.0)
                    .alias("hard_correct"),
                ]
            )
        )

    total_weight = float(per_read["read_weight"].sum() or 0.0)
    if total_weight <= 0:
        total_weight = 1.0
    summary = per_read.select(
        [
            pl.len().alias("n_reads"),
            pl.col("read_weight").sum().alias("total_read_weight"),
            pl.col("n_candidate_origins").mean().alias("mean_candidate_origins"),
            pl.col("n_assignment_targets").mean().alias("mean_assignment_targets"),
            (pl.col("n_assignment_targets") > 1)
            .cast(pl.Float64)
            .mean()
            .alias("ambiguous_read_fraction"),
            ((pl.col("assigned_posterior_sum") * pl.col("read_weight")).sum() / total_weight).alias(
                "assigned_fraction"
            ),
            ((pl.col("_max_posterior") * pl.col("read_weight")).sum() / total_weight).alias(
                "mean_max_posterior"
            ),
            (
                (pl.col("assignment_entropy_norm") * pl.col("read_weight")).sum() / total_weight
            ).alias("mean_assignment_entropy_norm"),
            (
                ((pl.col("_max_posterior") >= 0.8).cast(pl.Float64) * pl.col("read_weight")).sum()
                / total_weight
            ).alias("high_confidence_fraction"),
        ]
    )
    if not has_truth:
        return summary
    truth_summary = per_read.select(
        [
            ((pl.col("true_posterior") * pl.col("read_weight")).sum() / total_weight).alias(
                "soft_true_posterior"
            ),
            (
                (pl.col("hard_correct").cast(pl.Float64) * pl.col("read_weight")).sum()
                / total_weight
            ).alias("hard_accuracy"),
        ]
    )
    return pl.concat([summary, truth_summary], how="horizontal")


def assignment_identifiability(
    assignments: pl.DataFrame,
    *,
    confidence_threshold: float = 0.8,
    entropy_resolved_threshold: float = 0.3,
    entropy_ambiguous_threshold: float = 0.7,
) -> pl.DataFrame:
    """Classify per-read assignment uncertainty at the assignment-target level."""
    if assignments.is_empty():
        return pl.DataFrame()
    per_target = assignments.group_by(["read_key", "assignment_target"]).agg(
        [
            pl.col("posterior").sum().alias("target_posterior"),
            pl.col("read_weight").first().alias("read_weight"),
        ]
    )
    entropy = (
        per_target.with_columns(
            pl.when(pl.col("target_posterior") > _EPS)
            .then(-(pl.col("target_posterior") * pl.col("target_posterior").log()))
            .otherwise(0.0)
            .alias("_entropy_term")
        )
        .group_by("read_key")
        .agg(pl.col("_entropy_term").sum().alias("assignment_entropy"))
    )
    return (
        per_target.group_by("read_key")
        .agg(
            [
                pl.col("read_weight").first().alias("read_weight"),
                pl.len().alias("n_assignment_targets"),
                pl.col("target_posterior").sum().alias("assigned_posterior_sum"),
                pl.col("target_posterior").max().alias("max_posterior"),
            ]
        )
        .join(entropy, on="read_key", how="left")
        .with_columns(
            pl.when(pl.col("n_assignment_targets") > 1)
            .then(
                pl.col("assignment_entropy") / pl.col("n_assignment_targets").cast(pl.Float64).log()
            )
            .otherwise(0.0)
            .fill_nan(0.0)
            .alias("assignment_entropy_norm")
        )
        .with_columns(
            pl.when(pl.col("assigned_posterior_sum") <= _EPS)
            .then(pl.lit("unassigned"))
            .when(pl.col("n_assignment_targets") == 1)
            .then(pl.lit("unique"))
            .when(
                (pl.col("max_posterior") >= confidence_threshold)
                & (pl.col("assignment_entropy_norm") <= entropy_resolved_threshold)
            )
            .then(pl.lit("resolved"))
            .when(pl.col("assignment_entropy_norm") >= entropy_ambiguous_threshold)
            .then(pl.lit("ambiguous"))
            .otherwise(pl.lit("data_limited"))
            .alias("identifiability_class")
        )
    )


def summarize_identifiability(classes: pl.DataFrame) -> pl.DataFrame:
    if classes.is_empty():
        return pl.DataFrame()
    total_weight = float(classes["read_weight"].sum() or 1.0)
    return (
        classes.group_by("identifiability_class")
        .agg(
            [
                pl.len().alias("n_reads"),
                pl.col("read_weight").sum().alias("weighted_reads"),
                pl.col("max_posterior").mean().alias("mean_max_posterior"),
                pl.col("assignment_entropy_norm").mean().alias("mean_assignment_entropy_norm"),
            ]
        )
        .with_columns((pl.col("weighted_reads") / total_weight).alias("weighted_fraction"))
        .sort("weighted_reads", descending=True)
    )


def compare_read_assignment_methods(
    *,
    candidates_path: str,
    out_prefix: str,
    frame_support_path: str | None = None,
    cds_path: str | None = None,
    methods: Sequence[str] = ("unique", "fractional", "em", "frame_em"),
    abundance_key: str = "tran_id",
    truth_column: str = "is_true",
    write_assignments: bool = False,
    max_iter: int = 100,
    tol: float = 1e-7,
    abundance_prior: float = 1e-3,
    min_frame_support_count: float = 0.0,
    min_frame_support_evidence: float = 0.0,
    require_complete_frame_support: bool = False,
    frame_gate_min_likelihood_range: float = 0.25,
    frame_gate_min_read_fraction: float = 0.25,
) -> dict[str, str]:
    candidates = _load_table(candidates_path)
    frame_support = _load_table(frame_support_path) if frame_support_path else None
    cds_tran = _load_table(cds_path) if cds_path else None

    summaries: list[pl.DataFrame] = []
    convergence: list[pl.DataFrame] = []
    abundance: list[pl.DataFrame] = []
    identifiability: list[pl.DataFrame] = []
    diagnostics: list[pl.DataFrame] = []
    paths: dict[str, str] = {}
    assignment_paths: list[str] = []

    for method in methods:
        result = assign_reads(
            candidates,
            method=method,
            frame_support=frame_support,
            cds_tran=cds_tran,
            abundance_key=abundance_key,
            max_iter=max_iter,
            tol=tol,
            abundance_prior=abundance_prior,
            min_frame_support_count=min_frame_support_count,
            min_frame_support_evidence=min_frame_support_evidence,
            require_complete_frame_support=require_complete_frame_support,
            frame_gate_min_likelihood_range=frame_gate_min_likelihood_range,
            frame_gate_min_read_fraction=frame_gate_min_read_fraction,
        )
        metric = evaluate_assignments(result.assignments, truth_column=truth_column)
        if metric.is_empty():
            metric = pl.DataFrame(
                [{"n_reads": result.assignments["read_key"].n_unique(), "total_read_weight": None}]
            )
        summaries.append(
            metric.with_columns(
                pl.lit(method).alias("method"), pl.lit(abundance_key).alias("abundance_key")
            )
        )
        convergence.append(result.convergence.with_columns(pl.lit(method).alias("method")))
        abundance.append(
            result.abundance.with_columns(
                pl.lit(method).alias("method"), pl.lit(abundance_key).alias("abundance_key")
            )
        )
        classes = assignment_identifiability(result.assignments)
        identifiability.append(
            summarize_identifiability(classes).with_columns(
                pl.lit(method).alias("method"), pl.lit(abundance_key).alias("abundance_key")
            )
        )
        if result.diagnostics is not None and not result.diagnostics.is_empty():
            diagnostics.append(
                result.diagnostics.with_columns(
                    pl.lit(method).alias("method"), pl.lit(abundance_key).alias("abundance_key")
                )
            )
        if write_assignments:
            path = f"{out_prefix}.{method}.read_assignments.parquet"
            write_parquet_safe(result.assignments, path)
            assignment_paths.append(path)

    paths["summary"] = write_csv_safe(
        pl.concat(summaries), f"{out_prefix}.read_assignment.summary.csv"
    )
    paths["convergence"] = write_csv_safe(
        pl.concat(convergence), f"{out_prefix}.read_assignment.convergence.csv"
    )
    paths["abundance"] = write_csv_safe(
        pl.concat(abundance), f"{out_prefix}.read_assignment.abundance.csv"
    )
    paths["identifiability"] = write_csv_safe(
        pl.concat(identifiability), f"{out_prefix}.read_assignment.identifiability.csv"
    )
    if diagnostics:
        paths["diagnostics"] = write_csv_safe(
            pl.concat(diagnostics), f"{out_prefix}.read_assignment.diagnostics.csv"
        )
    if assignment_paths:
        paths["assignments"] = ",".join(assignment_paths)
    return paths

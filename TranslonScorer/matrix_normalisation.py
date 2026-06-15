"""Normalisation contract for matrix-based locus clustering and scoring.

This module keeps the clustering geometry separate from evidence-scale scoring:
``X_cluster`` is shape-normalised for sample clustering, while ``X_score``
retains raw/depth-normalised counts for cluster aggregates and event scoring.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Dict, List, Mapping, Optional

import numpy as np
import polars as pl


@dataclass(frozen=True)
class MatrixNormalisationResult:
    """Prepared locus matrices and sample QC labels."""

    normalisation_id: str
    X_score: np.ndarray
    X_cluster: np.ndarray
    sample_labels: pl.DataFrame
    params: Dict[str, object]


@dataclass(frozen=True)
class NormalisationResult(MatrixNormalisationResult):
    """Backward-compatible prepared locus result including clusters/aggregates."""

    cluster_aggregates: Dict[int, np.ndarray]
    cluster_summary: pl.DataFrame


def _as_float_matrix(matrix: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, got shape {arr.shape}")
    if np.any(arr < 0):
        raise ValueError("Profile counts must be non-negative")
    return arr


def _size_factor_array(
    sample_names: List[str],
    size_factors: Optional[Mapping[str, float] | np.ndarray],
) -> np.ndarray:
    if size_factors is None:
        factors = np.ones(len(sample_names), dtype=float)
    elif isinstance(size_factors, np.ndarray):
        factors = np.asarray(size_factors, dtype=float)
        if factors.shape != (len(sample_names),):
            raise ValueError(
                "size_factors array must have one value per sample "
                f"({len(sample_names)} expected, got {factors.shape})"
            )
    else:
        missing = [sid for sid in sample_names if sid not in size_factors]
        if missing:
            raise ValueError(f"Missing size factor for sample(s): {', '.join(missing)}")
        factors = np.array([float(size_factors[sid]) for sid in sample_names], dtype=float)

    if np.any(~np.isfinite(factors)) or np.any(factors <= 0):
        raise ValueError("All size factors must be finite and > 0")
    return factors


def depth_normalise_counts(
    matrix: np.ndarray,
    sample_names: List[str],
    *,
    size_factors: Optional[Mapping[str, float] | np.ndarray] = None,
) -> np.ndarray:
    """Apply global per-sample depth correction for score-scale profiles."""
    arr = _as_float_matrix(matrix)
    if arr.shape[0] != len(sample_names):
        raise ValueError(
            f"matrix has {arr.shape[0]} rows but {len(sample_names)} sample names"
        )
    factors = _size_factor_array(sample_names, size_factors)
    return arr / factors[:, None]


def shape_normalise_profiles(
    score_matrix: np.ndarray,
    *,
    method: str = "shifted_clr",
    pseudocount: float = 1e-6,
    profile_scale: float = 1.0,
) -> np.ndarray:
    """Build ``X_cluster`` from score-scale counts.

    Supported methods match ``docs/matrix_normalisation_strategy.md``:
    ``shifted_clr`` (default), ``l1_log``, and ``l2_log``.
    """
    X = _as_float_matrix(score_matrix)
    if pseudocount <= 0:
        raise ValueError("pseudocount must be > 0")
    if profile_scale <= 0:
        raise ValueError("profile_scale must be > 0")

    row_sums = X.sum(axis=1, keepdims=True)
    safe_sums = np.where(row_sums > 0, row_sums, 1.0)

    if method == "shifted_clr":
        u = X / safe_sums
        z = np.log(u + pseudocount)
        return z - z.mean(axis=1, keepdims=True)

    if method == "l1_log":
        u = X / safe_sums
        z = np.log(u + pseudocount)
        z = z - z.mean(axis=1, keepdims=True)
        norms = np.linalg.norm(z, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        return z / norms

    if method == "l2_log":
        z = np.log1p(X / profile_scale)
        norms = np.linalg.norm(z, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        return z / norms

    raise ValueError(
        f"Unknown shape normalisation method: {method!r}. "
        "Choose from: shifted_clr, l1_log, l2_log"
    )


def build_sample_qc_labels(
    raw_matrix: np.ndarray,
    score_matrix: np.ndarray,
    sample_names: List[str],
    size_factors: np.ndarray,
    *,
    min_raw_locus_counts: float = 0.0,
) -> pl.DataFrame:
    """Create the sample label table before clustering."""
    raw_totals = raw_matrix.sum(axis=1)
    score_totals = score_matrix.sum(axis=1)
    included = raw_totals >= min_raw_locus_counts
    reasons = [
        None if ok else f"locus_total_raw < {min_raw_locus_counts:g}"
        for ok in included.tolist()
    ]
    return pl.DataFrame(
        {
            "sample_id": sample_names,
            "locus_total_raw": raw_totals.tolist(),
            "size_factor": size_factors.tolist(),
            "locus_total_depth_norm": score_totals.tolist(),
            "qc_pass": included.tolist(),
            "exclusion_reason": reasons,
            "cluster_id": [-1] * len(sample_names),
        }
    )


def aggregate_score_profiles(
    score_matrix: np.ndarray,
    raw_matrix: np.ndarray,
    labels: np.ndarray,
    *,
    method: str = "sum_depth_norm",
) -> tuple[Dict[int, np.ndarray], pl.DataFrame]:
    """Aggregate clusters from ``X_score``, never from ``X_cluster``."""
    if method not in {"sum_depth_norm", "sum_raw", "mean_depth_norm"}:
        raise ValueError(
            f"Unknown aggregation method: {method!r}. "
            "Choose from: sum_depth_norm, sum_raw, mean_depth_norm"
        )

    aggregates: Dict[int, np.ndarray] = {}
    summary_rows = []
    for cid in sorted(set(labels[labels >= 0].tolist())):
        mask = labels == cid
        if method == "sum_raw":
            profile = raw_matrix[mask].sum(axis=0)
        elif method == "mean_depth_norm":
            profile = score_matrix[mask].mean(axis=0)
        else:
            profile = score_matrix[mask].sum(axis=0)

        aggregates[int(cid)] = profile
        summary_rows.append(
            {
                "cluster_id": int(cid),
                "sample_count": int(mask.sum()),
                "total_raw_counts": float(raw_matrix[mask].sum()),
                "total_depth_norm_counts": float(score_matrix[mask].sum()),
                "aggregation_method": method,
            }
        )

    return aggregates, pl.DataFrame(summary_rows)


def _normalisation_id(params: Mapping[str, object]) -> str:
    payload = json.dumps(params, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def normalise_locus_matrices(
    raw_matrix: np.ndarray,
    sample_names: List[str],
    *,
    size_factors: Optional[Mapping[str, float] | np.ndarray] = None,
    min_raw_locus_counts: float = 0.0,
    shape_method: str = "shifted_clr",
    pseudocount: float = 1e-6,
    profile_scale: float = 1.0,
) -> MatrixNormalisationResult:
    """Prepare score/cluster matrices and sample QC labels.

    Clustering and aggregation are intentionally handled by separate modules.
    Samples failing the raw-count gate remain in ``sample_labels`` with
    ``cluster_id = -1`` for compatibility with downstream joins.
    """
    raw = _as_float_matrix(raw_matrix)
    if raw.shape[0] != len(sample_names):
        raise ValueError(
            f"matrix has {raw.shape[0]} rows but {len(sample_names)} sample names"
        )

    factors = _size_factor_array(sample_names, size_factors)
    X_score = raw / factors[:, None]
    X_cluster = shape_normalise_profiles(
        X_score,
        method=shape_method,
        pseudocount=pseudocount,
        profile_scale=profile_scale,
    )
    sample_labels = build_sample_qc_labels(
        raw,
        X_score,
        sample_names,
        factors,
        min_raw_locus_counts=min_raw_locus_counts,
    )

    params: Dict[str, object] = {
        "shape_method": shape_method,
        "pseudocount": float(pseudocount),
        "profile_scale": float(profile_scale),
        "min_raw_locus_counts": float(min_raw_locus_counts),
        "size_factor_hash": hashlib.sha256(factors.tobytes()).hexdigest()[:16],
    }
    norm_id = _normalisation_id(params)

    return MatrixNormalisationResult(
        normalisation_id=norm_id,
        X_score=X_score,
        X_cluster=X_cluster,
        sample_labels=sample_labels,
        params=params,
    )


def normalise_and_cluster_locus(
    raw_matrix: np.ndarray,
    sample_names: List[str],
    *,
    size_factors: Optional[Mapping[str, float] | np.ndarray] = None,
    min_raw_locus_counts: float = 0.0,
    shape_method: str = "shifted_clr",
    pseudocount: float = 1e-6,
    profile_scale: float = 1.0,
    n_clusters: Optional[int] = 3,
    cluster_method: str = "hierarchical_cosine",
    distance_threshold: Optional[float] = None,
    linkage: str = "average",
    aggregation_method: str = "sum_depth_norm",
    min_cluster_size: int = 1,
    random_state: Optional[int] = 42,
) -> NormalisationResult:
    """Compatibility wrapper: normalise, cluster, and aggregate a locus."""
    from .clustering import aggregate_score_profiles_by_cluster, cluster_locus_profiles

    prepared = normalise_locus_matrices(
        raw_matrix,
        sample_names,
        size_factors=size_factors,
        min_raw_locus_counts=min_raw_locus_counts,
        shape_method=shape_method,
        pseudocount=pseudocount,
        profile_scale=profile_scale,
    )
    clustering = cluster_locus_profiles(
        prepared.X_cluster,
        sample_names,
        prepared.sample_labels,
        method=cluster_method,
        n_clusters=n_clusters,
        distance_threshold=distance_threshold,
        linkage_method=linkage,
        min_cluster_size=min_cluster_size,
        random_state=random_state,
    )
    raw = _as_float_matrix(raw_matrix)
    aggregates, summary = aggregate_score_profiles_by_cluster(
        prepared.X_score,
        raw,
        clustering.sample_labels,
        clustering.cluster_summary,
        method=aggregation_method,
    )
    if not summary.is_empty():
        summary = summary.with_columns(
            pl.lit(prepared.normalisation_id).alias("normalisation_id"),
            pl.lit(shape_method).alias("shape_method"),
        )
    return NormalisationResult(
        normalisation_id=prepared.normalisation_id,
        X_score=prepared.X_score,
        X_cluster=prepared.X_cluster,
        sample_labels=clustering.sample_labels,
        cluster_aggregates=aggregates,
        cluster_summary=summary,
        params={
            **prepared.params,
            "n_clusters": None if n_clusters is None else int(n_clusters),
            "cluster_method": cluster_method,
            "distance_threshold": distance_threshold,
            "linkage": linkage,
            "min_cluster_size": int(min_cluster_size),
            "aggregation_method": aggregation_method,
        },
    )

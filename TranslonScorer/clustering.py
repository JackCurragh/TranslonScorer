"""Profile normalisation and clustering for locus-level Ribo-Seq profiles.

Designed to work on the [n_samples × n_positions] profile matrices built by
locus_profiles.build_locus_profiles_zarr (or the sparse-Parquet analogue).

No heavy dependencies: uses numpy only.  scipy is imported lazily and only
for agglomerative clustering.

Public API
----------
normalise_profiles              — row-normalise a [samples × positions] matrix
cluster_locus_profiles          — full-featured QC-aware clustering
cluster_locus_shape_homogeneous — strict distributional shape clustering
explain_shape_cluster_drivers   — centroid contrast analysis
aggregate_score_profiles_by_cluster — collapse accepted clusters
cluster_profiles                — simple k-means/agglomerative entry point
aggregate_cluster_profiles      — {cluster_id: aggregate_profile}
build_profile_matrix            — tidy profiles DataFrame → dense matrix + pos_vec
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl

from TranslonScorer.utils.logging import log_info, log_warning


@dataclass(frozen=True)
class ProfileClusteringResult:
    """Cluster labels and diagnostics for a locus profile matrix."""

    sample_labels: pl.DataFrame
    cluster_summary: pl.DataFrame
    params: Dict[str, object]


@dataclass(frozen=True)
class ClusterShapeDriverResult:
    """Position-level diagnostics explaining shape-cluster separation."""

    cluster_centroids: pl.DataFrame
    cluster_contrasts: pl.DataFrame
    params: Dict[str, object]


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------


def normalise_profiles(
    matrix: np.ndarray,
    *,
    method: str = "total_count",
) -> np.ndarray:
    """Normalise each row (sample) of a [n_samples × n_positions] matrix.

    Methods
    -------
    total_count  – divide each row by its sum × 1e6 (reads-per-million).
                   Rows summing to zero are left as zero.
    zscore       – centre and scale each row independently.
    minmax       – rescale each row to [0, 1].
    """
    out = matrix.astype(float).copy()
    if method == "total_count":
        totals = out.sum(axis=1, keepdims=True)
        nz = totals > 0
        out[nz.ravel()] = out[nz.ravel()] / totals[nz].reshape(-1, 1) * 1e6
    elif method == "zscore":
        means = out.mean(axis=1, keepdims=True)
        stds = out.std(axis=1, keepdims=True)
        stds[stds == 0] = 1.0
        out = (out - means) / stds
    elif method == "minmax":
        mins = out.min(axis=1, keepdims=True)
        maxs = out.max(axis=1, keepdims=True)
        ranges = maxs - mins
        ranges[ranges == 0] = 1.0
        out = (out - mins) / ranges
    else:
        raise ValueError(
            f"Unknown normalisation method: {method!r}. " "Choose from: total_count, zscore, minmax"
        )
    return out


# ---------------------------------------------------------------------------
# Internal numpy k-means with cosine distance
# ---------------------------------------------------------------------------


def _l2_normalise(X: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return X / norms


def _kmeans_cosine(
    X: np.ndarray,
    k: int,
    *,
    n_init: int = 10,
    max_iter: int = 300,
    random_state: Optional[int] = 42,
) -> np.ndarray:
    """K-means clustering using cosine distance (pure numpy)."""
    rng = np.random.default_rng(random_state)
    X_n = _l2_normalise(X)

    best_labels: Optional[np.ndarray] = None
    best_score = -np.inf

    for _ in range(n_init):
        idx = rng.choice(len(X_n), k, replace=False)
        centers = X_n[idx].copy()

        for _ in range(max_iter):
            sims = X_n @ centers.T
            labels = sims.argmax(axis=1)

            new_centers = np.zeros_like(centers)
            for j in range(k):
                mask = labels == j
                if mask.sum() > 0:
                    c = X_n[mask].mean(axis=0)
                    cn = np.linalg.norm(c)
                    new_centers[j] = c / cn if cn > 0 else c

            if np.allclose(new_centers, centers, atol=1e-6):
                break
            centers = new_centers

        sims = X_n @ centers.T
        score = sims[np.arange(len(sims)), labels].mean()
        if score > best_score:
            best_score = score
            best_labels = labels.copy()

    return best_labels  # type: ignore[return-value]


def _agglomerative_cosine(X: np.ndarray, k: int) -> np.ndarray:
    try:
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import pdist
    except ImportError as exc:
        raise ImportError(
            "scipy is required for agglomerative clustering. " "Install with: pip install scipy"
        ) from exc
    dists = pdist(X, metric="cosine")
    Z = linkage(dists, method="average")
    return fcluster(Z, k, criterion="maxclust") - 1


def _rank_rows(X: np.ndarray) -> np.ndarray:
    try:
        from scipy.stats import rankdata
    except ImportError as exc:
        raise ImportError("scipy is required for Spearman clustering") from exc
    return np.vstack([rankdata(row, method="average") for row in X])


def _distance_space(X: np.ndarray, metric: str) -> np.ndarray:
    if metric == "cosine":
        return X
    if metric in {"correlation", "pearson"}:
        return X - X.mean(axis=1, keepdims=True)
    if metric == "spearman":
        ranks = _rank_rows(X)
        return ranks - ranks.mean(axis=1, keepdims=True)
    raise ValueError("metric must be one of: cosine, correlation, pearson, spearman")


def _pairwise_distances(X: np.ndarray, metric: str = "cosine") -> np.ndarray:
    Y = _distance_space(np.asarray(X, dtype=float), metric)
    Y = _l2_normalise(Y)
    sims = np.clip(Y @ Y.T, -1.0, 1.0)
    d = 1.0 - sims
    np.fill_diagonal(d, 0.0)
    return d


def _probability_profiles(matrix: np.ndarray) -> np.ndarray:
    X = np.clip(np.asarray(matrix, dtype=float), 0.0, None)
    totals = X.sum(axis=1, keepdims=True)
    return np.divide(X, totals, out=np.zeros_like(X), where=totals > 0)


def _smooth_probability_profiles(P: np.ndarray, window: int = 1) -> np.ndarray:
    if window <= 1:
        return P.copy()
    if window % 2 == 0:
        raise ValueError("smooth_window must be odd")
    kernel = np.ones(int(window), dtype=float) / float(window)
    smoothed = np.vstack([np.convolve(row, kernel, mode="same") for row in P])
    totals = smoothed.sum(axis=1, keepdims=True)
    return np.divide(smoothed, totals, out=np.zeros_like(smoothed), where=totals > 0)


def _probability_distance(p: np.ndarray, q: np.ndarray, metric: str) -> float:
    if metric == "jensen_shannon":
        p = np.asarray(p, dtype=float)
        q = np.asarray(q, dtype=float)
        p = p / max(float(p.sum()), 1e-12)
        q = q / max(float(q.sum()), 1e-12)
        m = 0.5 * (p + q)

        def kl(a: np.ndarray, b: np.ndarray) -> float:
            mask = a > 0
            return float(np.sum(a[mask] * np.log2(a[mask] / np.maximum(b[mask], 1e-300))))

        return float(np.sqrt(max(0.0, 0.5 * kl(p, m) + 0.5 * kl(q, m))))
    if metric == "hellinger":
        return float(np.linalg.norm(np.sqrt(p) - np.sqrt(q)) / np.sqrt(2.0))
    raise ValueError("distance_metric must be one of: jensen_shannon, hellinger")


def _probability_distance_contributions(
    p: np.ndarray,
    q: np.ndarray,
    metric: str,
) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / max(float(p.sum()), 1e-12)
    q = q / max(float(q.sum()), 1e-12)
    if metric == "jensen_shannon":
        m = 0.5 * (p + q)
        out = np.zeros_like(p)
        p_mask = p > 0
        q_mask = q > 0
        out[p_mask] += 0.5 * p[p_mask] * np.log2(p[p_mask] / np.maximum(m[p_mask], 1e-300))
        out[q_mask] += 0.5 * q[q_mask] * np.log2(q[q_mask] / np.maximum(m[q_mask], 1e-300))
        return out
    if metric == "hellinger":
        return (np.sqrt(p) - np.sqrt(q)) ** 2 / 2.0
    raise ValueError("distance_metric must be one of: jensen_shannon, hellinger")


def _pairwise_probability_distances(P: np.ndarray, metric: str = "jensen_shannon") -> np.ndarray:
    P = np.asarray(P, dtype=float)
    n = P.shape[0]
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = _probability_distance(P[i], P[j], metric)
            D[i, j] = D[j, i] = d
    return D


def _agglomerative_precomputed(
    distances: np.ndarray,
    *,
    linkage_method: str = "average",
    distance_threshold: float = 0.12,
) -> np.ndarray:
    try:
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform
    except ImportError as exc:
        raise ImportError(
            "scipy is required for homogeneous shape clustering. " "Install with: pip install scipy"
        ) from exc
    if linkage_method not in {"average", "complete", "single"}:
        raise ValueError("linkage_method must be one of: average, complete, single")
    if distances.shape[0] == 1:
        return np.zeros(1, dtype=int)
    condensed = squareform(distances, checks=False)
    Z = linkage(condensed, method=linkage_method)
    return fcluster(Z, float(distance_threshold), criterion="distance").astype(int) - 1


def _agglomerative_distance(
    X: np.ndarray,
    *,
    metric: str = "cosine",
    linkage_method: str = "average",
    n_clusters: Optional[int] = None,
    distance_threshold: Optional[float] = None,
) -> np.ndarray:
    try:
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import pdist
    except ImportError as exc:
        raise ImportError(
            "scipy is required for hierarchical clustering. " "Install with: pip install scipy"
        ) from exc
    if linkage_method not in {"average", "complete", "single"}:
        raise ValueError("linkage_method must be one of: average, complete, single")
    if n_clusters is None and distance_threshold is None:
        distance_threshold = 0.25
    Y = _distance_space(X, metric)
    condensed = pdist(Y, metric="cosine")
    Z = linkage(condensed, method=linkage_method)
    if distance_threshold is not None:
        labels = fcluster(Z, float(distance_threshold), criterion="distance")
    else:
        labels = fcluster(Z, int(n_clusters), criterion="maxclust")
    return labels.astype(int) - 1


def _cluster_diagnostics(
    X: np.ndarray,
    labels: np.ndarray,
    *,
    metric: str = "cosine",
    min_cluster_size: int = 1,
    min_separation: Optional[float] = None,
    max_within_distance: Optional[float] = None,
) -> Tuple[pl.DataFrame, np.ndarray]:
    active = sorted(set(labels[labels >= 0].tolist()))
    if not active:
        return pl.DataFrame(), labels

    distances = _pairwise_distances(X, metric=metric)
    rows = []
    weak = set()
    for cid in active:
        idx = np.where(labels == cid)[0]
        if len(idx) > 1:
            within = distances[np.ix_(idx, idx)]
            tri = within[np.triu_indices(len(idx), k=1)]
            mean_within = float(tri.mean()) if tri.size else 0.0
            max_within = float(tri.max()) if tri.size else 0.0
        else:
            mean_within = 0.0
            max_within = 0.0

        other_dists = []
        for other in active:
            if other == cid:
                continue
            jdx = np.where(labels == other)[0]
            if len(jdx):
                other_dists.append(float(distances[np.ix_(idx, jdx)].mean()))
        nearest = min(other_dists) if other_dists else None

        reasons = []
        if len(idx) < min_cluster_size:
            reasons.append(f"n_samples < {min_cluster_size}")
        if min_separation is not None and nearest is not None and nearest < min_separation:
            reasons.append(f"nearest_other_cluster_distance < {min_separation:g}")
        if max_within_distance is not None and max_within > max_within_distance:
            reasons.append(f"max_within_cluster_distance > {max_within_distance:g}")
        low_conf = bool(reasons)
        if low_conf:
            weak.add(cid)

        rows.append(
            {
                "cluster_id": int(cid),
                "n_samples": int(len(idx)),
                "mean_within_cluster_cosine_distance": mean_within,
                "max_within_cluster_cosine_distance": max_within,
                "nearest_other_cluster_distance": nearest,
                "low_confidence": low_conf,
                "reason": "; ".join(reasons) if reasons else None,
            }
        )

    pruned = labels.copy()
    for cid in weak:
        pruned[labels == cid] = -1
    return pl.DataFrame(rows), pruned


def _sample_assignment_diagnostics(
    X: np.ndarray,
    labels: np.ndarray,
    *,
    metric: str = "cosine",
) -> Tuple[np.ndarray, np.ndarray]:
    active = sorted(set(labels[labels >= 0].tolist()))
    nearest = np.full(len(labels), np.nan, dtype=float)
    margin = np.full(len(labels), np.nan, dtype=float)
    if not active:
        return nearest, margin

    Y = _l2_normalise(_distance_space(X, metric))
    centroids = []
    for cid in active:
        c = Y[labels == cid].mean(axis=0)
        cn = np.linalg.norm(c)
        centroids.append(c / cn if cn > 0 else c)
    C = np.vstack(centroids)
    d = 1.0 - np.clip(Y @ C.T, -1.0, 1.0)
    for i, cid in enumerate(labels):
        if cid < 0:
            continue
        cluster_pos = active.index(int(cid))
        own = float(d[i, cluster_pos])
        nearest[i] = own
        if len(active) > 1:
            other = np.delete(d[i], cluster_pos)
            margin[i] = float(other.min() - own)
    return nearest, margin


def _relabel_active_clusters(labels: np.ndarray) -> Tuple[np.ndarray, Dict[int, int]]:
    active = sorted(set(labels[labels >= 0].tolist()))
    mapping = {old: new for new, old in enumerate(active)}
    relabeled = labels.copy()
    for old, new in mapping.items():
        relabeled[labels == old] = new
    return relabeled, mapping


# ---------------------------------------------------------------------------
# Public clustering functions
# ---------------------------------------------------------------------------


def cluster_locus_profiles(
    X_cluster: np.ndarray,
    sample_names: List[str],
    sample_qc_labels: pl.DataFrame,
    *,
    method: str = "hierarchical_cosine",
    n_clusters: Optional[int] = None,
    distance_threshold: Optional[float] = 0.25,
    linkage_method: str = "average",
    min_cluster_size: int = 1,
    min_separation: Optional[float] = None,
    max_within_distance: Optional[float] = None,
    sample_margin_threshold: Optional[float] = None,
    random_state: Optional[int] = 42,
) -> ProfileClusteringResult:
    """Cluster QC-passing shape profiles and return labels/diagnostics."""
    X = np.asarray(X_cluster, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, got shape {X.shape}")
    if X.shape[0] != len(sample_names):
        raise ValueError(f"matrix has {X.shape[0]} rows but {len(sample_names)} sample names")
    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be >= 1")

    aliases = {
        "agglomerative": "hierarchical_cosine",
        "hierarchical": "hierarchical_cosine",
        "kmeans": "kmeans_cosine",
        "correlation": "hierarchical_correlation",
        "pearson": "hierarchical_correlation",
    }
    method = aliases.get(method, method)
    metric = "cosine"
    if method == "hierarchical_correlation":
        metric = "correlation"
    elif method == "hierarchical_pearson":
        metric = "pearson"
    elif method == "hierarchical_spearman":
        metric = "spearman"
    elif method not in {"hierarchical_cosine", "kmeans_cosine", "dbscan"}:
        raise ValueError(
            f"Unknown clustering method: {method!r}. Choose: hierarchical_cosine, "
            "hierarchical_correlation, hierarchical_spearman, kmeans_cosine, dbscan"
        )

    qc_pass = sample_qc_labels.get_column("qc_pass").to_numpy().astype(bool)
    included_idx = np.where(qc_pass)[0]
    labels = np.full(len(sample_names), -1, dtype=int)
    statuses = np.where(qc_pass, "assigned", "qc_excluded").astype(object)
    summary = pl.DataFrame()

    if len(included_idx) > 0:
        Xi = X[included_idx]
        if method == "kmeans_cosine":
            k = min(int(n_clusters or 1), len(included_idx))
            local = (
                np.zeros(len(included_idx), dtype=int)
                if k == 1
                else _kmeans_cosine(Xi, k, random_state=random_state)
            )
        elif method == "dbscan":
            try:
                from sklearn.cluster import DBSCAN
            except ImportError as exc:
                raise ImportError("scikit-learn is required for DBSCAN clustering") from exc
            eps = 0.25 if distance_threshold is None else float(distance_threshold)
            local = DBSCAN(eps=eps, min_samples=min_cluster_size, metric="cosine").fit_predict(
                _distance_space(Xi, metric)
            )
        else:
            if len(included_idx) == 1 or (n_clusters == 1 and distance_threshold is None):
                local = np.zeros(len(included_idx), dtype=int)
            else:
                local = _agglomerative_distance(
                    Xi,
                    metric=metric,
                    linkage_method=linkage_method,
                    n_clusters=n_clusters,
                    distance_threshold=distance_threshold,
                )
        labels[included_idx] = local.astype(int)
        summary, pruned_local = _cluster_diagnostics(
            Xi,
            local.astype(int),
            metric=metric,
            min_cluster_size=min_cluster_size,
            min_separation=min_separation,
            max_within_distance=max_within_distance,
        )
        if len(pruned_local) == len(included_idx):
            weak_mask = (local >= 0) & (pruned_local < 0)
            labels[included_idx] = pruned_local
            statuses[included_idx[weak_mask]] = "weak_cluster"

        nearest, margin = _sample_assignment_diagnostics(X, labels, metric=metric)
        if sample_margin_threshold is not None:
            ambiguous = (
                (labels >= 0) & np.isfinite(margin) & (margin < float(sample_margin_threshold))
            )
            labels[ambiguous] = -1
            statuses[ambiguous] = "ambiguous"

        labels, relabel_map = _relabel_active_clusters(labels)
        if not summary.is_empty():
            summary = summary.with_columns(
                pl.col("cluster_id").alias("original_cluster_id"),
                pl.col("cluster_id")
                .replace_strict(relabel_map, default=None)
                .cast(pl.Int64)
                .alias("cluster_id"),
            )
    else:
        nearest = np.full(len(sample_names), np.nan, dtype=float)
        margin = np.full(len(sample_names), np.nan, dtype=float)

    out_labels = sample_qc_labels.with_columns(
        pl.Series("cluster_id", labels.tolist()),
        pl.Series("clustering_status", statuses.tolist()),
        pl.Series("nearest_cluster_distance", nearest.tolist()),
        pl.Series("nearest_centroid_margin", margin.tolist()),
    )
    if not summary.is_empty():
        active_after_prune = set(labels[labels >= 0].tolist())
        summary = summary.with_columns(
            pl.col("cluster_id").is_in(active_after_prune).not_().alias("pruned")
        )

    return ProfileClusteringResult(
        sample_labels=out_labels,
        cluster_summary=summary,
        params={
            "method": method,
            "metric": metric,
            "n_clusters": n_clusters,
            "distance_threshold": distance_threshold,
            "linkage_method": linkage_method,
            "min_cluster_size": int(min_cluster_size),
            "min_separation": min_separation,
            "max_within_distance": max_within_distance,
            "sample_margin_threshold": sample_margin_threshold,
        },
    )


def cluster_locus_shape_homogeneous(
    raw_matrix: np.ndarray,
    sample_names: List[str],
    *,
    min_sample_counts: float = 1000.0,
    smooth_window: int = 7,
    distance_metric: str = "jensen_shannon",
    distance_threshold: float = 0.12,
    linkage_method: str = "average",
    min_cluster_size: int = 3,
    max_within_distance: Optional[float] = None,
    max_member_to_centroid_distance: Optional[float] = None,
) -> ProfileClusteringResult:
    """Cluster samples into homogeneous locus-profile shapes (strict, distributional)."""
    X = np.asarray(raw_matrix, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, got shape {X.shape}")
    if X.shape[0] != len(sample_names):
        raise ValueError(f"matrix has {X.shape[0]} rows but {len(sample_names)} sample names")
    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be >= 1")
    if smooth_window < 1 or smooth_window % 2 == 0:
        raise ValueError("smooth_window must be a positive odd integer")
    if max_within_distance is None:
        max_within_distance = distance_threshold
    if max_member_to_centroid_distance is None:
        max_member_to_centroid_distance = distance_threshold

    totals = np.clip(X, 0.0, None).sum(axis=1)
    include_mask = totals >= float(min_sample_counts)
    included_idx = np.where(include_mask)[0]

    labels = np.full(len(sample_names), -1, dtype=int)
    statuses = np.where(include_mask, "shape_outlier", "low_coverage").astype(object)
    member_to_centroid = np.full(len(sample_names), np.nan, dtype=float)
    nearest_cluster_distance = np.full(len(sample_names), np.nan, dtype=float)
    nearest_centroid_margin = np.full(len(sample_names), np.nan, dtype=float)
    summary = pl.DataFrame()

    P_all = _smooth_probability_profiles(_probability_profiles(X), smooth_window)

    if len(included_idx) > 0:
        P = P_all[included_idx]
        D = _pairwise_probability_distances(P, metric=distance_metric)
        local = _agglomerative_precomputed(
            D, linkage_method=linkage_method, distance_threshold=distance_threshold
        )

        accepted_local = local.copy()
        active = sorted(set(local.tolist()))
        cluster_rows = []
        centroids: Dict[int, np.ndarray] = {}
        for cid in active:
            idx = np.where(local == cid)[0]
            if len(idx) > 1:
                within = D[np.ix_(idx, idx)]
                tri = within[np.triu_indices(len(idx), k=1)]
                mean_within = float(tri.mean()) if tri.size else 0.0
                max_within = float(tri.max()) if tri.size else 0.0
            else:
                mean_within = 0.0
                max_within = 0.0

            centroid = P[idx].mean(axis=0)
            centroid = centroid / max(float(centroid.sum()), 1e-12)
            centroid_d = np.array(
                [_probability_distance(P[j], centroid, distance_metric) for j in idx]
            )
            max_centroid = float(centroid_d.max()) if centroid_d.size else 0.0
            mean_centroid = float(centroid_d.mean()) if centroid_d.size else 0.0

            reasons = []
            if len(idx) < min_cluster_size:
                reasons.append(f"n_samples < {min_cluster_size}")
            if max_within > float(max_within_distance):
                reasons.append(f"max_within_distance > {float(max_within_distance):g}")
            if max_centroid > float(max_member_to_centroid_distance):
                reasons.append(
                    f"max_member_to_centroid_distance > {float(max_member_to_centroid_distance):g}"
                )

            low_confidence = bool(reasons)
            if low_confidence:
                accepted_local[idx] = -1
            else:
                centroids[int(cid)] = centroid

            cluster_rows.append(
                {
                    "cluster_id": int(cid),
                    "n_samples": int(len(idx)),
                    "mean_pairwise_distance": mean_within,
                    "max_pairwise_distance": max_within,
                    "mean_member_to_centroid_distance": mean_centroid,
                    "max_member_to_centroid_distance": max_centroid,
                    "low_confidence": low_confidence,
                    "reason": "; ".join(reasons) if reasons else None,
                }
            )

        labels[included_idx] = accepted_local
        accepted_mask = accepted_local >= 0
        statuses[included_idx[accepted_mask]] = "assigned"
        statuses[included_idx[~accepted_mask]] = "shape_outlier"

        if centroids:
            active_centroids = sorted(centroids)
            C_mat = np.vstack([centroids[cid] for cid in active_centroids])
            all_to_centroid = np.zeros((len(sample_names), len(active_centroids)), dtype=float)
            for i in range(len(sample_names)):
                for j, _cid in enumerate(active_centroids):
                    all_to_centroid[i, j] = _probability_distance(
                        P_all[i], C_mat[j], distance_metric
                    )
            for i, cid in enumerate(labels):
                nearest_cluster_distance[i] = float(all_to_centroid[i].min())
                sorted_d = np.sort(all_to_centroid[i])
                if sorted_d.size > 1:
                    nearest_centroid_margin[i] = float(sorted_d[1] - sorted_d[0])
                if cid >= 0:
                    pos = active_centroids.index(int(cid))
                    member_to_centroid[i] = float(all_to_centroid[i, pos])

        labels, relabel_map = _relabel_active_clusters(labels)
        if cluster_rows:
            summary = pl.DataFrame(cluster_rows).with_columns(
                pl.col("cluster_id").alias("original_cluster_id"),
                pl.col("cluster_id")
                .replace_strict(relabel_map, default=None)
                .cast(pl.Int64)
                .alias("cluster_id"),
            )
            active_after_prune = set(labels[labels >= 0].tolist())
            summary = summary.with_columns(
                pl.col("cluster_id").is_in(active_after_prune).not_().alias("pruned")
            )

    out_labels = pl.DataFrame(
        {
            "sample_id": sample_names,
            "qc_pass": include_mask.tolist(),
            "cluster_id": labels.tolist(),
            "clustering_status": statuses.tolist(),
            "locus_total_raw": totals.tolist(),
            "locus_total_depth_norm": totals.tolist(),
            "member_to_centroid_distance": member_to_centroid.tolist(),
            "nearest_cluster_distance": nearest_cluster_distance.tolist(),
            "nearest_centroid_margin": nearest_centroid_margin.tolist(),
        }
    )

    return ProfileClusteringResult(
        sample_labels=out_labels,
        cluster_summary=summary,
        params={
            "method": "homogeneous_shape",
            "profile_normalisation": "probability",
            "distance_metric": distance_metric,
            "smooth_window": int(smooth_window),
            "distance_threshold": float(distance_threshold),
            "linkage_method": linkage_method,
            "min_sample_counts": float(min_sample_counts),
            "min_cluster_size": int(min_cluster_size),
            "max_within_distance": float(max_within_distance),
            "max_member_to_centroid_distance": float(max_member_to_centroid_distance),
        },
    )


def explain_shape_cluster_drivers(
    raw_matrix: np.ndarray,
    sample_names: List[str],
    sample_labels: pl.DataFrame,
    *,
    positions: Optional[np.ndarray] = None,
    smooth_window: int = 7,
    distance_metric: str = "jensen_shannon",
    top_n: int = 25,
) -> ClusterShapeDriverResult:
    """Explain which locus positions separate accepted shape clusters."""
    X = np.asarray(raw_matrix, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, got shape {X.shape}")
    if X.shape[0] != len(sample_names):
        raise ValueError(f"matrix has {X.shape[0]} rows but {len(sample_names)} sample names")
    if smooth_window < 1 or smooth_window % 2 == 0:
        raise ValueError("smooth_window must be a positive odd integer")
    if top_n < 1:
        raise ValueError("top_n must be >= 1")

    if positions is None:
        pos = np.arange(X.shape[1], dtype=int)
    else:
        pos = np.asarray(positions)
        if pos.shape[0] != X.shape[1]:
            raise ValueError(
                f"positions must contain one value per matrix column "
                f"({X.shape[1]} expected, got {pos.shape[0]})"
            )

    label_lookup = {
        str(row["sample_id"]): int(row["cluster_id"])
        for row in sample_labels.select(["sample_id", "cluster_id"]).iter_rows(named=True)
    }
    labels = np.array([label_lookup.get(str(s), -1) for s in sample_names], dtype=int)
    active_clusters = sorted(set(labels[labels >= 0].tolist()))

    P = _smooth_probability_profiles(_probability_profiles(X), smooth_window)
    aggregate = np.clip(X, 0.0, None).sum(axis=0)
    aggregate_fraction = np.divide(
        aggregate,
        max(float(aggregate.sum()), 1e-12),
        out=np.zeros_like(aggregate, dtype=float),
    )

    centroid_rows = []
    centroids: Dict[int, np.ndarray] = {}
    for cluster_id in active_clusters:
        mask = labels == cluster_id
        centroid = P[mask].mean(axis=0)
        centroid = centroid / max(float(centroid.sum()), 1e-12)
        centroids[int(cluster_id)] = centroid
        for i, value in enumerate(centroid):
            centroid_rows.append(
                {
                    "cluster_id": int(cluster_id),
                    "position_index": int(i),
                    "position": int(pos[i]),
                    "centroid_probability": float(value),
                    "aggregate_fraction": float(aggregate_fraction[i]),
                }
            )

    contrast_rows = []
    for a_i, cluster_a in enumerate(active_clusters):
        for cluster_b in active_clusters[a_i + 1 :]:
            p = centroids[int(cluster_a)]
            q = centroids[int(cluster_b)]
            contributions = _probability_distance_contributions(p, q, distance_metric)
            total = float(contributions.sum())
            order = np.argsort(contributions)[::-1]
            cumulative = 0.0
            for rank, idx in enumerate(order[:top_n], start=1):
                contribution = float(contributions[idx])
                cumulative += contribution
                contrast_rows.append(
                    {
                        "cluster_a": int(cluster_a),
                        "cluster_b": int(cluster_b),
                        "rank": int(rank),
                        "position_index": int(idx),
                        "position": int(pos[idx]),
                        "distance_contribution": contribution,
                        "distance_fraction": contribution / total if total > 0 else 0.0,
                        "cumulative_distance_fraction": cumulative / total if total > 0 else 0.0,
                        "centroid_a_probability": float(p[idx]),
                        "centroid_b_probability": float(q[idx]),
                        "probability_difference": float(q[idx] - p[idx]),
                        "aggregate_fraction": float(aggregate_fraction[idx]),
                    }
                )

    return ClusterShapeDriverResult(
        cluster_centroids=pl.DataFrame(centroid_rows),
        cluster_contrasts=pl.DataFrame(contrast_rows),
        params={
            "method": "shape_cluster_drivers",
            "profile_normalisation": "probability",
            "distance_metric": distance_metric,
            "smooth_window": int(smooth_window),
            "top_n": int(top_n),
        },
    )


def aggregate_score_profiles_by_cluster(
    X_score: np.ndarray,
    raw_matrix: np.ndarray,
    sample_labels: pl.DataFrame,
    cluster_summary: Optional[pl.DataFrame] = None,
    *,
    method: str = "sum_depth_norm",
    include_low_confidence: bool = False,
) -> Tuple[Dict[int, np.ndarray], pl.DataFrame]:
    """Aggregate accepted clusters from X_score."""
    if method not in {"sum_depth_norm", "sum_raw", "mean_depth_norm"}:
        raise ValueError(
            f"Unknown aggregation method: {method!r}. "
            "Choose from: sum_depth_norm, sum_raw, mean_depth_norm"
        )
    labels = sample_labels.get_column("cluster_id").to_numpy().astype(int)
    if (
        not include_low_confidence
        and cluster_summary is not None
        and not cluster_summary.is_empty()
    ):
        weak = set(
            cluster_summary.filter(pl.col("low_confidence")).get_column("cluster_id").to_list()
        )
    else:
        weak = set()

    aggregates: Dict[int, np.ndarray] = {}
    rows = []
    for cid in sorted(set(labels[labels >= 0].tolist())):
        if cid in weak:
            continue
        mask = labels == cid
        if method == "sum_raw":
            profile = raw_matrix[mask].sum(axis=0)
        elif method == "mean_depth_norm":
            profile = X_score[mask].mean(axis=0)
        else:
            profile = X_score[mask].sum(axis=0)
        aggregates[int(cid)] = profile
        rows.append(
            {
                "cluster_id": int(cid),
                "sample_count": int(mask.sum()),
                "total_raw_counts": float(raw_matrix[mask].sum()),
                "total_depth_norm_counts": float(X_score[mask].sum()),
                "aggregation_method": method,
            }
        )

    out = pl.DataFrame(rows)
    if cluster_summary is not None and not cluster_summary.is_empty() and not out.is_empty():
        out = out.join(cluster_summary, on="cluster_id", how="left")
    return aggregates, out


# ---------------------------------------------------------------------------
# Simple entry points
# ---------------------------------------------------------------------------


def cluster_profiles(
    matrix: np.ndarray,
    sample_names: List[str],
    *,
    n_clusters: int = 3,
    method: str = "kmeans",
    min_coverage: float = 0.0,
) -> Tuple[np.ndarray, List[str], List[str]]:
    """Cluster samples (rows) by profile shape.

    Returns (labels, included_samples, excluded_samples).
    labels: int array (n_samples,); -1 for excluded samples.
    """
    if len(sample_names) != matrix.shape[0]:
        raise ValueError(f"matrix has {matrix.shape[0]} rows but {len(sample_names)} sample names")

    sample_totals = matrix.sum(axis=1)
    include_mask = sample_totals >= min_coverage
    included_idx = np.where(include_mask)[0]
    excluded_idx = np.where(~include_mask)[0]

    included_samples = [sample_names[i] for i in included_idx]
    excluded_samples = [sample_names[i] for i in excluded_idx]
    labels = np.full(len(sample_names), -1, dtype=int)

    n_avail = len(included_idx)
    if n_avail == 0:
        log_warning("No samples with sufficient coverage to cluster.")
        return labels, included_samples, excluded_samples

    if n_avail < n_clusters:
        log_warning(
            f"Only {n_avail} samples available; reducing n_clusters from {n_clusters} to {n_avail}"
        )
        n_clusters = n_avail

    X = matrix[included_idx]
    if n_clusters == 1:
        cluster_labels = np.zeros(n_avail, dtype=int)
    elif method == "kmeans":
        cluster_labels = _kmeans_cosine(X, n_clusters)
    elif method == "agglomerative":
        cluster_labels = _agglomerative_cosine(X, n_clusters)
    else:
        raise ValueError(f"Unknown clustering method: {method!r}. Choose: kmeans, agglomerative")

    for array_pos, orig_pos in enumerate(included_idx):
        labels[orig_pos] = int(cluster_labels[array_pos])

    log_info(
        f"Clustered {n_avail} samples into {n_clusters} clusters "
        f"(excluded {len(excluded_samples)} below min_coverage={min_coverage})"
    )
    return labels, included_samples, excluded_samples


def aggregate_cluster_profiles(
    matrix: np.ndarray,
    labels: np.ndarray,
    *,
    method: str = "mean",
) -> Dict[int, np.ndarray]:
    """Compute one aggregate profile per cluster. method: 'mean' or 'sum'."""
    unique_clusters = sorted(set(labels[labels >= 0].tolist()))
    result: Dict[int, np.ndarray] = {}
    for cid in unique_clusters:
        mask = labels == cid
        if method == "mean":
            result[cid] = matrix[mask].mean(axis=0)
        elif method == "sum":
            result[cid] = matrix[mask].sum(axis=0)
        else:
            raise ValueError(f"Unknown aggregation method: {method!r}. Choose: mean, sum")
    return result


def build_profile_matrix(
    profiles: pl.DataFrame,
    sample_names: List[str],
    *,
    pos_col: str = "pos",
    count_col: str = "count",
    sample_col: str = "sample_id",
) -> Tuple[np.ndarray, np.ndarray]:
    """Build a [n_samples × n_positions] dense matrix from a tidy profiles DataFrame.

    Returns (matrix float32 [n_samples × n_positions], pos_vec int64 [n_positions]).
    """
    if profiles.is_empty():
        return np.zeros((len(sample_names), 0), dtype=np.float32), np.array([], dtype=np.int64)

    pos_vec = np.sort(profiles.get_column(pos_col).unique().cast(pl.Int64).to_numpy())
    pos_to_idx = {int(p): i for i, p in enumerate(pos_vec)}
    n_pos = len(pos_vec)
    n_samp = len(sample_names)
    samp_to_idx = {s: i for i, s in enumerate(sample_names)}

    matrix = np.zeros((n_samp, n_pos), dtype=np.float32)
    for row in profiles.select([sample_col, pos_col, count_col]).iter_rows():
        sid, pos, count = row
        si = samp_to_idx.get(str(sid))
        pi = pos_to_idx.get(int(pos))
        if si is not None and pi is not None:
            matrix[si, pi] += float(count)

    return matrix, pos_vec

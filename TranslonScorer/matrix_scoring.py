"""Matrix scoring: sparse Parquet → profiles → scored ORFs.

Three scoring modes are exposed:

  aggregate          – collapse all samples into a single pseudo-sample,
                       score once.  Cheapest; good for a first-pass.

  per_sample         – score each sample independently with its own
                       per-sample offsets.  Returns one scored table per
                       sample.

  clustered          – calls profile_clustering, scores the per-cluster
                       aggregate profile, and attaches cluster labels to
                       the output.

All modes share the same downstream scoring code path
(bigwig.score_transcript) so the score column is directly comparable to
BigWig-based scoring runs.
"""
from __future__ import annotations

import gc
from collections import defaultdict
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import numpy as np
import polars as pl

from .utils.logging import log_info, log_warning, log_error
from .file_handlers import bam as bam_handlers
from .file_handlers.bigwig import score_transcript
from .coverage.profiles import _compute_asite_profiles
from .pipeline.score_schema import ensure_score_schema


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _default_offsets(lengths: Iterable[int], default: int = 15) -> Dict[int, int]:
    return {int(L): int(default) for L in set(int(x) for x in lengths)}


def _genomic_to_transcript_profiles(
    genomic_counts: pl.DataFrame,
    exon_df: pl.DataFrame,
    offsets: Dict[int, int],
    *,
    keep_sample: bool = False,
) -> pl.DataFrame:
    """Project a genomic-counts DataFrame to transcript-space A-site profiles.

    genomic_counts columns: chr, start, stop, strand, length, count
      (optionally sample_id when keep_sample=True)

    Returns: tran_id, pos, count [, sample_id]
    """
    if genomic_counts.is_empty():
        schema = {"tran_id": pl.Utf8, "pos": pl.Int64, "count": pl.Float64}
        if keep_sample:
            schema["sample_id"] = pl.Utf8
        return pl.DataFrame(schema=schema)

    mapped = bam_handlers.bamtranscript(genomic_counts, exon_df)
    if mapped.is_empty():
        schema = {"tran_id": pl.Utf8, "pos": pl.Int64, "count": pl.Float64}
        if keep_sample:
            schema["sample_id"] = pl.Utf8
        return pl.DataFrame(schema=schema)

    profiles = _compute_asite_profiles(mapped, offsets, keep_length=False)
    return profiles


def _build_tran_data(prof_sub: pl.DataFrame) -> Dict:
    """Convert a per-transcript profile slice to the sparse dict expected by score_transcript."""
    positions = prof_sub.get_column("pos").cast(pl.Int64).to_list()
    values = prof_sub.get_column("count").cast(pl.Float64).to_list()
    max_pos = int(max(positions)) + 1 if positions else 0
    return {"positions": positions, "values": values, "max_pos": max_pos}


def profiles_to_scored_orfs(
    profiles: pl.DataFrame,
    orf_df: pl.DataFrame,
    *,
    sru_range: int = 15,
    old_scoring: bool = False,
    score_mode: str = "raw",
    input_type: str = "matrix",
) -> pl.DataFrame:
    """Score ORFs directly from A-site profiles (no bigWig required).

    profiles  – tran_id, pos, count
    orf_df    – tran_id, start, stop, type (+ any other columns)

    Returns a scored ORF DataFrame conforming to SCORE_SCHEMA_COLUMNS.
    """
    if profiles.is_empty() or orf_df.is_empty():
        return pl.DataFrame()

    # Build per-transcript lookup
    profiles_by_tran: Dict[str, pl.DataFrame] = {
        str(key[0] if isinstance(key, tuple) else key): sub
        for key, sub in profiles.group_by("tran_id")
    }

    results: List[pl.DataFrame] = []
    tran_ids = orf_df.get_column("tran_id").unique().to_list()

    for tran_id in tran_ids:
        prof = profiles_by_tran.get(str(tran_id))
        if prof is None or prof.is_empty():
            continue

        tran_data = _build_tran_data(prof)
        tran_orfs = orf_df.filter(pl.col("tran_id") == tran_id)
        orf_data = tran_orfs.to_dict(as_series=False)

        scored_dict = score_transcript(
            (str(tran_id), tran_data, orf_data, old_scoring, sru_range)
        )
        if scored_dict is not None:
            results.append(pl.from_dict(scored_dict))

        del tran_data, orf_data
        gc.collect()

    if not results:
        return pl.DataFrame()

    scored = pl.concat(results)
    scored = ensure_score_schema(scored, score_mode=score_mode, input_type=input_type)
    return scored


# ---------------------------------------------------------------------------
# Aggregate mode
# ---------------------------------------------------------------------------

def aggregate_profiles_from_genomic_counts(
    genomic_counts: pl.DataFrame,
    exon_df: pl.DataFrame,
    offsets: Dict[int, int],
) -> pl.DataFrame:
    """Collapse sample-resolved genomic counts into aggregate transcript profiles.

    genomic_counts: sample_id, chr, start, stop, strand, length, count
    Returns: tran_id, pos, count
    """
    if genomic_counts.is_empty():
        return pl.DataFrame(schema={"tran_id": pl.Utf8, "pos": pl.Int64, "count": pl.Float64})

    agg = (
        genomic_counts
        .group_by(["chr", "start", "stop", "strand", "length"])
        .agg(pl.col("count").sum())
    )
    return _genomic_to_transcript_profiles(agg, exon_df, offsets)


def score_aggregate(
    manifest_path: str,
    bam_path: str,
    orf_df: pl.DataFrame,
    exon_df: pl.DataFrame,
    offsets: Dict[int, int],
    *,
    sample_names: Optional[List[str]] = None,
    sru_range: int = 15,
    old_scoring: bool = False,
) -> pl.DataFrame:
    """Full aggregate scoring pipeline: sparse Parquet → aggregate profiles → scored ORFs."""
    from .file_handlers.sparse_parquet import sparse_matrix_genomic_counts

    log_info("Aggregate mode: loading genomic counts from sparse Parquet…")
    genomic_counts = sparse_matrix_genomic_counts(
        bam_path=bam_path,
        manifest_path=manifest_path,
        exon_df=exon_df,
        sample_names=sample_names,
    )

    log_info(f"Loaded {genomic_counts.height:,} genomic-count rows across "
             f"{genomic_counts.get_column('sample_id').n_unique() if not genomic_counts.is_empty() else 0} samples")

    profiles = aggregate_profiles_from_genomic_counts(genomic_counts, exon_df, offsets)
    log_info(f"Aggregate profiles: {profiles.height:,} transcript-position rows")

    return profiles_to_scored_orfs(
        profiles, orf_df, sru_range=sru_range, old_scoring=old_scoring, input_type="sparse_matrix_aggregate"
    )


# ---------------------------------------------------------------------------
# Per-sample mode
# ---------------------------------------------------------------------------

def per_sample_profiles_from_genomic_counts(
    genomic_counts: pl.DataFrame,
    exon_df: pl.DataFrame,
    sample_offsets: pl.DataFrame,
    *,
    default_offset: int = 15,
) -> Iterator[Tuple[str, pl.DataFrame]]:
    """Yield (sample_id, profiles_df) for each sample in genomic_counts.

    sample_offsets: sample_id, length, offset  (per-sample from QC gate)
      If a sample_id is absent, default_offset is used for all lengths.
    """
    if genomic_counts.is_empty():
        return

    # Build per-sample offset lookups
    offsets_by_sample: Dict[str, Dict[int, int]] = {}
    if not sample_offsets.is_empty():
        for key, grp in sample_offsets.group_by("sample_id"):
            sid = str(key[0] if isinstance(key, tuple) else key)
            offsets_by_sample[sid] = {
                int(r[0]): int(r[1])
                for r in grp.select(["length", "offset"]).iter_rows()
            }

    for key, sample_df in genomic_counts.group_by("sample_id"):
        sid = str(key[0] if isinstance(key, tuple) else key)
        offsets = offsets_by_sample.get(
            sid,
            _default_offsets(sample_df.get_column("length").to_list(), default_offset),
        )
        profiles = _genomic_to_transcript_profiles(
            sample_df.drop("sample_id"),
            exon_df,
            offsets,
        )
        yield sid, profiles
        del sample_df, profiles
        gc.collect()


def score_per_sample(
    manifest_path: str,
    bam_path: str,
    orf_df: pl.DataFrame,
    exon_df: pl.DataFrame,
    sample_offsets: pl.DataFrame,
    *,
    sample_names: Optional[List[str]] = None,
    default_offset: int = 15,
    sru_range: int = 15,
    old_scoring: bool = False,
) -> Dict[str, pl.DataFrame]:
    """Full per-sample scoring pipeline.

    sample_offsets: sample_id, length, offset  (empty → default_offset for all)
    Returns: {sample_id → scored_orf_df}
    """
    from .file_handlers.sparse_parquet import sparse_matrix_genomic_counts

    log_info("Per-sample mode: loading genomic counts from sparse Parquet…")
    genomic_counts = sparse_matrix_genomic_counts(
        bam_path=bam_path,
        manifest_path=manifest_path,
        exon_df=exon_df,
        sample_names=sample_names,
    )
    log_info(f"Loaded {genomic_counts.height:,} rows for "
             f"{genomic_counts.get_column('sample_id').n_unique() if not genomic_counts.is_empty() else 0} samples")

    results: Dict[str, pl.DataFrame] = {}
    for sid, profiles in per_sample_profiles_from_genomic_counts(
        genomic_counts, exon_df, sample_offsets, default_offset=default_offset
    ):
        log_info(f"Scoring sample {sid}…")
        scored = profiles_to_scored_orfs(
            profiles, orf_df, sru_range=sru_range, old_scoring=old_scoring,
            input_type="sparse_matrix_per_sample",
        )
        if not scored.is_empty():
            results[sid] = scored

    return results


# ---------------------------------------------------------------------------
# Cluster-aggregate mode  (delegates heavy work to profile_clustering)
# ---------------------------------------------------------------------------

def score_clustered(
    locus_profiles_matrix: np.ndarray,   # [n_samples, n_positions]
    sample_names: List[str],
    pos_vector: np.ndarray,               # genomic/transcript positions (length n_positions)
    tran_id: str,
    orf_df: pl.DataFrame,
    *,
    n_clusters: Optional[int] = None,
    normalise_method: str = "shifted_clr",
    cluster_method: str = "hierarchical_cosine",
    distance_threshold: Optional[float] = 0.25,
    linkage: str = "average",
    size_factors: Optional[Dict[str, float] | np.ndarray] = None,
    min_raw_locus_counts: float = 0.0,
    min_cluster_size: int = 1,
    pseudocount: float = 1e-6,
    aggregation_method: str = "sum_depth_norm",
    sru_range: int = 15,
    old_scoring: bool = False,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Score ORFs for each cluster of a single locus's profile matrix.

    Returns:
        scored  – per-cluster ORF scores with cluster_id column
        labels  – DataFrame: sample_id, cluster_id
    """
    results, labels_df, _summary = score_locus_matrix_levels(
        locus_profiles_matrix,
        sample_names,
        pos_vector,
        tran_id,
        orf_df,
        levels=("cluster",),
        n_clusters=n_clusters,
        normalise_method=normalise_method,
        cluster_method=cluster_method,
        distance_threshold=distance_threshold,
        linkage=linkage,
        size_factors=size_factors,
        min_raw_locus_counts=min_raw_locus_counts,
        min_cluster_size=min_cluster_size,
        pseudocount=pseudocount,
        aggregation_method=aggregation_method,
        sru_range=sru_range,
        old_scoring=old_scoring,
    )
    return results.get("cluster", pl.DataFrame()), labels_df


def _score_profile_vector(
    profile: np.ndarray,
    pos_vector: np.ndarray,
    tran_id: str,
    orf_df: pl.DataFrame,
    *,
    input_type: str,
    sru_range: int,
    old_scoring: bool,
) -> pl.DataFrame:
    profiles = pl.DataFrame(
        {
            "tran_id": [tran_id] * len(pos_vector),
            "pos": pos_vector.astype(np.int64).tolist(),
            "count": np.asarray(profile, dtype=float).tolist(),
        }
    )
    return profiles_to_scored_orfs(
        profiles,
        orf_df,
        sru_range=sru_range,
        old_scoring=old_scoring,
        input_type=input_type,
    )


def score_locus_matrix_levels(
    locus_profiles_matrix: np.ndarray,
    sample_names: List[str],
    pos_vector: np.ndarray,
    tran_id: str,
    orf_df: pl.DataFrame,
    *,
    levels: Tuple[str, ...] = ("sample", "cluster", "aggregate"),
    n_clusters: Optional[int] = None,
    normalise_method: str = "shifted_clr",
    cluster_method: str = "hierarchical_cosine",
    distance_threshold: Optional[float] = 0.25,
    linkage: str = "average",
    size_factors: Optional[Dict[str, float] | np.ndarray] = None,
    min_raw_locus_counts: float = 0.0,
    min_cluster_size: int = 1,
    pseudocount: float = 1e-6,
    aggregation_method: str = "sum_depth_norm",
    include_low_confidence_clusters: bool = False,
    sru_range: int = 15,
    old_scoring: bool = False,
) -> Tuple[Dict[str, pl.DataFrame], pl.DataFrame, pl.DataFrame]:
    """Normalise once, then score sample, cluster, and whole aggregate levels."""
    from .matrix_normalisation import normalise_locus_matrices
    from .clustering import (
        aggregate_score_profiles_by_cluster,
        cluster_locus_profiles,
    )

    wanted = set(levels)
    unknown = wanted - {"sample", "cluster", "aggregate"}
    if unknown:
        raise ValueError(f"Unknown scoring level(s): {', '.join(sorted(unknown))}")

    prepared = normalise_locus_matrices(
        locus_profiles_matrix,
        sample_names,
        size_factors=size_factors,
        min_raw_locus_counts=min_raw_locus_counts,
        shape_method=normalise_method,
        pseudocount=pseudocount,
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
    )
    raw = np.asarray(locus_profiles_matrix, dtype=float)
    cluster_aggs, cluster_summary = aggregate_score_profiles_by_cluster(
        prepared.X_score,
        raw,
        clustering.sample_labels,
        clustering.cluster_summary,
        method=aggregation_method,
        include_low_confidence=include_low_confidence_clusters,
    )
    if not cluster_summary.is_empty():
        cluster_summary = cluster_summary.with_columns(
            pl.lit(prepared.normalisation_id).alias("normalisation_id"),
            pl.lit(normalise_method).alias("shape_method"),
        )

    results: Dict[str, pl.DataFrame] = {}
    if "aggregate" in wanted:
        aggregate_profile = prepared.X_score[
            clustering.sample_labels.get_column("qc_pass").to_numpy().astype(bool)
        ].sum(axis=0)
        scored = _score_profile_vector(
            aggregate_profile,
            pos_vector,
            tran_id,
            orf_df,
            input_type="sparse_matrix_whole_aggregate",
            sru_range=sru_range,
            old_scoring=old_scoring,
        )
        if not scored.is_empty():
            results["aggregate"] = scored.with_columns(
                pl.lit("aggregate").alias("score_level"),
                pl.lit("aggregate").alias("aggregate_id"),
                pl.lit(prepared.normalisation_id).alias("normalisation_id"),
            )

    if "sample" in wanted:
        sample_rows: List[pl.DataFrame] = []
        qc_pass = clustering.sample_labels.get_column("qc_pass").to_numpy().astype(bool)
        for idx, sid in enumerate(sample_names):
            if not qc_pass[idx]:
                continue
            scored = _score_profile_vector(
                prepared.X_score[idx],
                pos_vector,
                tran_id,
                orf_df,
                input_type="sparse_matrix_sample",
                sru_range=sru_range,
                old_scoring=old_scoring,
            )
            if not scored.is_empty():
                sample_rows.append(
                    scored.with_columns(
                        pl.lit("sample").alias("score_level"),
                        pl.lit(sid).alias("sample_id"),
                        pl.lit(prepared.normalisation_id).alias("normalisation_id"),
                    )
                )
        results["sample"] = pl.concat(sample_rows) if sample_rows else pl.DataFrame()

    if "cluster" in wanted:
        cluster_rows: List[pl.DataFrame] = []
        for cid, agg_profile in cluster_aggs.items():
            scored = _score_profile_vector(
                agg_profile,
                pos_vector,
                tran_id,
                orf_df,
                input_type="sparse_matrix_clustered",
                sru_range=sru_range,
                old_scoring=old_scoring,
            )
            if not scored.is_empty():
                cluster_rows.append(
                    scored.with_columns(
                        pl.lit("cluster").alias("score_level"),
                        pl.lit(int(cid)).alias("cluster_id"),
                        pl.lit(prepared.normalisation_id).alias("normalisation_id"),
                    )
                )
        results["cluster"] = pl.concat(cluster_rows) if cluster_rows else pl.DataFrame()

    labels = clustering.sample_labels.with_columns(
        pl.lit(prepared.normalisation_id).alias("normalisation_id")
    )
    return results, labels, cluster_summary

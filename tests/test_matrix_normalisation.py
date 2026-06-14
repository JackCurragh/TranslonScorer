import numpy as np
import polars as pl
import pytest


def test_depth_normalisation_uses_global_size_factors():
    from TranslonScorer.pipeline.matrix_normalisation import depth_normalise_counts

    raw = np.array([[10.0, 30.0], [100.0, 300.0]])
    score = depth_normalise_counts(raw, ["S1", "S2"], size_factors={"S1": 1.0, "S2": 10.0})

    np.testing.assert_allclose(score[0], [10.0, 30.0])
    np.testing.assert_allclose(score[1], [10.0, 30.0])


def test_shifted_clr_removes_magnitude_not_shape():
    from TranslonScorer.pipeline.matrix_normalisation import shape_normalise_profiles

    score = np.array(
        [
            [10.0, 30.0, 0.0],
            [100.0, 300.0, 0.0],
            [30.0, 10.0, 0.0],
        ]
    )
    cluster = shape_normalise_profiles(score, method="shifted_clr", pseudocount=1e-6)

    np.testing.assert_allclose(cluster[0], cluster[1], rtol=1e-7, atol=1e-7)
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(cluster[0], cluster[2], rtol=1e-7, atol=1e-7)
    np.testing.assert_allclose(cluster.mean(axis=1), 0.0, atol=1e-12)


def test_normalise_and_cluster_gates_on_raw_counts_and_sums_score_matrix():
    from TranslonScorer.pipeline.matrix_normalisation import normalise_and_cluster_locus

    raw = np.array(
        [
            [10.0, 30.0, 0.0],
            [100.0, 300.0, 0.0],
            [0.0, 0.0, 5.0],
        ]
    )
    result = normalise_and_cluster_locus(
        raw,
        ["S1", "S2", "low_count"],
        size_factors={"S1": 1.0, "S2": 10.0, "low_count": 1.0},
        min_raw_locus_counts=20.0,
        n_clusters=1,
    )

    labels = {
        row[0]: row[1]
        for row in result.sample_labels.select(["sample_id", "cluster_id"]).iter_rows()
    }
    assert labels["S1"] == 0
    assert labels["S2"] == 0
    assert labels["low_count"] == -1

    np.testing.assert_allclose(result.cluster_aggregates[0], [20.0, 60.0, 0.0])
    assert result.cluster_summary.get_column("total_raw_counts").to_list() == [440.0]
    assert result.cluster_summary.get_column("total_depth_norm_counts").to_list() == [80.0]
    assert result.cluster_summary.get_column("normalisation_id").to_list() == [
        result.normalisation_id
    ]


def test_normalise_locus_matrices_does_not_cluster_or_aggregate():
    from TranslonScorer.pipeline.matrix_normalisation import normalise_locus_matrices

    raw = np.array([[10.0, 30.0, 0.0], [100.0, 300.0, 0.0]])
    result = normalise_locus_matrices(
        raw,
        ["S1", "S2"],
        size_factors={"S1": 1.0, "S2": 10.0},
    )

    np.testing.assert_allclose(result.X_score[0], result.X_score[1])
    assert result.sample_labels.get_column("cluster_id").to_list() == [-1, -1]
    assert not hasattr(result, "cluster_aggregates")


def test_cluster_aggregate_uses_score_matrix_not_cluster_matrix():
    from TranslonScorer.pipeline.matrix_normalisation import normalise_locus_matrices
    from TranslonScorer.pipeline.profile_clustering import (
        aggregate_score_profiles_by_cluster,
        cluster_locus_profiles,
    )

    raw = np.array([[10.0, 30.0, 0.0], [100.0, 300.0, 0.0]])
    prepared = normalise_locus_matrices(
        raw,
        ["S1", "S2"],
        size_factors={"S1": 1.0, "S2": 10.0},
    )
    clustered = cluster_locus_profiles(
        prepared.X_cluster,
        ["S1", "S2"],
        prepared.sample_labels,
        n_clusters=1,
        distance_threshold=None,
    )
    aggregates, summary = aggregate_score_profiles_by_cluster(
        prepared.X_score,
        raw,
        clustered.sample_labels,
        clustered.cluster_summary,
    )

    np.testing.assert_allclose(aggregates[0], [20.0, 60.0, 0.0])
    assert summary.item(0, "total_depth_norm_counts") == pytest.approx(80.0)


def test_homogeneous_shape_clustering_ignores_depth_and_prunes_outliers():
    from TranslonScorer.pipeline.profile_clustering import cluster_locus_shape_homogeneous

    raw = np.array(
        [
            [10, 30, 60, 30, 10, 0, 0, 0],
            [100, 300, 600, 300, 100, 0, 0, 0],
            [8, 24, 48, 24, 8, 0, 0, 0],
            [0, 0, 0, 10, 30, 60, 30, 10],
            [0, 0, 0, 100, 300, 600, 300, 100],
            [0, 0, 0, 8, 24, 48, 24, 8],
            [1, 0, 0, 0, 0, 0, 0, 0],
        ],
        dtype=float,
    )
    result = cluster_locus_shape_homogeneous(
        raw,
        ["A1", "A2", "A3", "B1", "B2", "B3", "low"],
        min_sample_counts=20,
        smooth_window=1,
        distance_threshold=0.20,
        min_cluster_size=3,
        max_within_distance=0.20,
        max_member_to_centroid_distance=0.12,
    )

    labels = {
        row["sample_id"]: row["cluster_id"]
        for row in result.sample_labels.select(["sample_id", "cluster_id"]).iter_rows(named=True)
    }
    statuses = {
        row["sample_id"]: row["clustering_status"]
        for row in result.sample_labels.select(["sample_id", "clustering_status"]).iter_rows(named=True)
    }

    assert labels["A1"] == labels["A2"] == labels["A3"]
    assert labels["B1"] == labels["B2"] == labels["B3"]
    assert labels["A1"] != labels["B1"]
    assert labels["low"] == -1
    assert statuses["low"] == "low_coverage"
    assert result.cluster_summary.filter(~pl.col("pruned")).height == 2
    assert result.params["profile_normalisation"] == "probability"


def test_homogeneous_shape_clustering_rejects_forced_heterogeneous_cluster():
    from TranslonScorer.pipeline.profile_clustering import cluster_locus_shape_homogeneous

    raw = np.array(
        [
            [10, 30, 60, 30, 10, 0, 0, 0],
            [100, 300, 600, 300, 100, 0, 0, 0],
            [0, 0, 0, 100, 300, 600, 300, 100],
        ],
        dtype=float,
    )
    result = cluster_locus_shape_homogeneous(
        raw,
        ["A1", "A2", "B1"],
        min_sample_counts=20,
        smooth_window=1,
        distance_threshold=1.0,
        min_cluster_size=3,
        max_within_distance=0.10,
        max_member_to_centroid_distance=0.10,
    )

    assert result.sample_labels.get_column("cluster_id").to_list() == [-1, -1, -1]
    assert set(result.sample_labels.get_column("clustering_status").to_list()) == {
        "shape_outlier"
    }


def test_shape_cluster_driver_diagnostics_rank_discriminating_positions():
    from TranslonScorer.pipeline.profile_clustering import (
        cluster_locus_shape_homogeneous,
        explain_shape_cluster_drivers,
    )

    raw = np.array(
        [
            [10, 30, 80, 30, 10, 0, 0, 0],
            [100, 300, 800, 300, 100, 0, 0, 0],
            [8, 24, 64, 24, 8, 0, 0, 0],
            [0, 0, 0, 10, 30, 80, 30, 10],
            [0, 0, 0, 100, 300, 800, 300, 100],
            [0, 0, 0, 8, 24, 64, 24, 8],
        ],
        dtype=float,
    )
    clustered = cluster_locus_shape_homogeneous(
        raw,
        ["A1", "A2", "A3", "B1", "B2", "B3"],
        min_sample_counts=20,
        smooth_window=1,
        distance_threshold=0.20,
        min_cluster_size=3,
        max_within_distance=0.20,
        max_member_to_centroid_distance=0.12,
    )
    drivers = explain_shape_cluster_drivers(
        raw,
        ["A1", "A2", "A3", "B1", "B2", "B3"],
        clustered.sample_labels,
        positions=np.arange(raw.shape[1]),
        smooth_window=1,
        top_n=4,
    )

    assert drivers.cluster_centroids.get_column("cluster_id").n_unique() == 2
    top_positions = set(drivers.cluster_contrasts.head(4).get_column("position").to_list())
    assert top_positions & {2, 5}
    assert drivers.cluster_contrasts.get_column("distance_fraction").sum() > 0


def test_score_clustered_exposes_normalisation_labels():
    import polars as pl

    from TranslonScorer.pipeline.matrix_scoring import score_clustered

    tran_id = "TX"
    matrix = np.zeros((3, 60))
    matrix[0, 10:40] = 5.0
    matrix[1, 10:40] = 50.0
    matrix[2, 45:47] = 1.0
    orf_df = pl.DataFrame(
        {
            "tran_id": [tran_id],
            "start": [10],
            "stop": [40],
            "type": ["CDS"],
        }
    )

    scored, labels = score_clustered(
        matrix,
        ["S1", "S2", "low_count"],
        np.arange(60),
        tran_id,
        orf_df,
        n_clusters=1,
        size_factors={"S1": 1.0, "S2": 10.0, "low_count": 1.0},
        min_raw_locus_counts=20.0,
    )

    assert "normalisation_id" in scored.columns
    assert "locus_total_raw" in labels.columns
    assert "exclusion_reason" in labels.columns
    assert labels.filter(pl.col("sample_id") == "low_count").item(0, "cluster_id") == -1

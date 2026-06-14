"""Re-export shim — implementations live in TranslonScorer.clustering.

Keeps the old pipeline.profile_clustering import path working until T14.
"""
from TranslonScorer.clustering import (  # noqa: F401 — re-exported
    ProfileClusteringResult,
    ClusterShapeDriverResult,
    normalise_profiles,
    _l2_normalise,
    _kmeans_cosine,
    _agglomerative_cosine,
    _rank_rows,
    _distance_space,
    _pairwise_distances,
    _probability_profiles,
    _smooth_probability_profiles,
    _probability_distance,
    _probability_distance_contributions,
    _pairwise_probability_distances,
    _agglomerative_precomputed,
    _agglomerative_distance,
    _cluster_diagnostics,
    _sample_assignment_diagnostics,
    _relabel_active_clusters,
    cluster_locus_profiles,
    cluster_locus_shape_homogeneous,
    explain_shape_cluster_drivers,
    aggregate_score_profiles_by_cluster,
    cluster_profiles,
    aggregate_cluster_profiles,
    build_profile_matrix,
)

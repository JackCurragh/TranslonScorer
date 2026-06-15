import json
import numpy as np
import polars as pl

from TranslonScorer.frame.bleed import apply_confusion, apply_confusion_counts, learn_confusion
from TranslonScorer.frame.latent import fit_latent
from TranslonScorer.frame.frame_disambiguation import (
    candidate_frame_table,
    summarize_frame_disambiguation,
)
from TranslonScorer.frame.frame_method_compare import compare_frame_support_tables, validate_frame_support_on_cds
from TranslonScorer.frame_support import build_frame_support
from TranslonScorer.model import FrameSupportParams
from TranslonScorer.pipeline.profile_compare import compare_profiles
from TranslonScorer.pipeline.panel_manifest import merge_panel_manifest, validate_panel_manifest
from TranslonScorer.pipeline.score_schema import (
    add_frame_score_columns,
    compare_score_tables,
    ensure_score_schema,
)
from TranslonScorer.pipeline.rdg_flux_export import export_rdg_flux_v1, rdg_flux_position_table
from TranslonScorer.pipeline.transcript_coords import cds_to_transcript_space


def test_cds_to_transcript_space_maps_minus_strand():
    cds = pl.DataFrame({
        "tran_id": ["tx_plus", "tx_minus"],
        "chr": ["chr1", "chr1"],
        "start": [110, 510],
        "stop": [160, 560],
        "strand": ["+", "-"],
    })
    exons = pl.DataFrame({
        "tran_id": ["tx_plus", "tx_minus"],
        "chr": ["chr1", "chr1"],
        "start": [[100, 200], [500, 600]],
        "stop": [[180, 260], [560, 660]],
        "tran_start": [[0, 80], [0, 60]],
        "tran_stop": [[80, 140], [60, 120]],
        "strand": ["+", "-"],
    })

    mapped = cds_to_transcript_space(cds, exons).sort("tran_id")
    rows = {r["tran_id"]: r for r in mapped.iter_rows(named=True)}
    assert rows["tx_plus"]["start"] == 10
    assert rows["tx_plus"]["stop"] == 60
    assert rows["tx_minus"]["start"] == 0
    assert rows["tx_minus"]["stop"] == 50


def test_ensure_score_schema_and_frame_summary():
    scored = pl.DataFrame({
        "tran_id": ["tx1"],
        "start": [0],
        "stop": [9],
        "type": ["CDS"],
        "score": [3.0],
    })
    frame_support = pl.DataFrame({
        "tran_id": ["tx1", "tx1", "tx1"],
        "codon": [0, 1, 2],
        "p0": [0.9, 0.8, 0.7],
        "p1": [0.05, 0.1, 0.2],
        "p2": [0.05, 0.1, 0.1],
        "entropy": [0.2, 0.3, 0.4],
        "method": ["linear+hmm", "linear+hmm", "linear+hmm"],
    })
    profiles = pl.DataFrame({
        "tran_id": ["tx1", "tx1", "tx1"],
        "pos": [0, 1, 3],
        "count": [10.0, 5.0, 2.0],
    })

    out = add_frame_score_columns(ensure_score_schema(scored), frame_support, profiles=profiles)
    assert out["orf_id"][0] == "tx1:0:9:CDS"
    assert round(out["frame_posterior_mean"][0], 6) == 0.8
    assert round(out["frame_entropy_mean"][0], 6) == 0.3
    assert round(out["frame_weighted_count"][0], 6) == 10.0 * 0.9 + 5.0 * 0.05 + 2.0 * 0.8
    assert out["frame_method"][0] == "linear+hmm"


def test_ensure_score_schema_preserves_empty_input_schema():
    empty = pl.DataFrame(schema={"tran_id": pl.Utf8, "start": pl.Int64, "stop": pl.Int64})
    out = ensure_score_schema(empty)
    assert out.height == 0
    assert "tran_id" in out.columns
    assert "score" in out.columns
    assert "orf_id" in out.columns


def test_compare_score_tables_reports_delta():
    raw = pl.DataFrame({"tran_id": ["tx1"], "start": [0], "stop": [9], "score": [1.0]})
    frame = pl.DataFrame({"tran_id": ["tx1"], "start": [0], "stop": [9], "score": [1.5]})
    comparison, summary = compare_score_tables(raw, frame)
    assert comparison["score_delta"][0] == 0.5
    assert summary["n_joined"][0] == 1
    assert summary["mean_score_delta"][0] == 0.5


def test_panel_manifest_validation_and_merge():
    orfs = pl.DataFrame({"tran_id": ["tx1"], "start": [0], "stop": [9], "type": ["CDS"]})
    panel = pl.DataFrame({
        "panel_id": ["p1"],
        "tran_id": ["tx1"],
        "start": [0],
        "stop": [9],
        "type": ["CDS"],
        "label": ["positive"],
        "category": ["annotated_cds"],
    })
    report = validate_panel_manifest(panel)
    assert bool(report["is_valid"][0])
    merged = merge_panel_manifest(orfs, panel)
    assert merged["label"][0] == "positive"
    assert merged["category"][0] == "annotated_cds"


def test_compare_profiles_summary():
    a = pl.DataFrame({"tran_id": ["tx1", "tx1"], "pos": [0, 1], "count": [1.0, 2.0]})
    b = pl.DataFrame({"tran_id": ["tx1", "tx1"], "pos": [0, 2], "count": [1.0, 3.0]})
    joined, summary = compare_profiles(a, b, label_a="a", label_b="b")
    assert joined.height == 3
    assert summary["n_positions_shared_nonzero"][0] == 1
    assert summary["total_count_a"][0] == 3.0
    assert summary["total_count_b"][0] == 4.0


def test_frame_disambiguation_counts_frame_discordant_candidates():
    candidates = pl.DataFrame({
        "read_key": ["r1", "r1", "r2", "r2", "r3"],
        "tran_id": ["tx1", "tx2", "tx1", "tx2", "tx1"],
        "tran_start_bam": [12, 13, 21, 24, 30],
        "count": [1.0, 1.0, 2.0, 2.0, 1.0],
    })
    cds = pl.DataFrame({
        "tran_id": ["tx1", "tx2"],
        "start": [0, 0],
        "stop": [90, 90],
    })
    frames = candidate_frame_table(candidates, cds)
    classes, summary = summarize_frame_disambiguation(frames)
    r1 = classes.filter(pl.col("read_key") == "r1")
    r2 = classes.filter(pl.col("read_key") == "r2")
    assert bool(r1["is_frame_discordant"][0])
    assert not bool(r2["is_frame_discordant"][0])
    assert summary["n_ambiguous_read_keys"][0] == 2
    assert summary["n_frame_discordant_read_keys"][0] == 1


def test_confusion_learning_preserves_asymmetric_leakage():
    profiles = pl.DataFrame({
        "tran_id": ["tx1", "tx1", "tx1"],
        "pos": [30, 31, 32],
        "count": [70.0, 10.0, 20.0],
    })
    cds = pl.DataFrame({"tran_id": ["tx1"], "start": [0], "stop": [90]})

    C = learn_confusion(profiles, cds, by_length=False)[None]

    assert round(float(C[0, 0]), 3) == 0.715
    assert round(float(C[1, 0]), 3) == 0.095
    assert round(float(C[2, 0]), 3) == 0.190
    assert float(C[1, 0]) != float(C[2, 0])


def test_linear_correction_bumps_true_frame_support_above_observed_fraction():
    leakage = np.array([0.70, 0.10, 0.20])
    C = np.array([[leakage[(observed - true) % 3] for true in range(3)] for observed in range(3)])
    observed = np.array([[80.0, 10.0, 10.0]])

    adjusted = apply_confusion_counts(observed, C, alpha=0.0)
    posterior = apply_confusion(observed, C, alpha=0.0)

    assert round(float(adjusted.sum()), 6) == 100.0
    assert adjusted[0, 0] > observed[0, 0]
    assert posterior[0, 0] > 0.95


def test_latent_em_matches_linear_when_confusion_is_fixed():
    leakage = np.array([0.70, 0.10, 0.20])
    C = np.array([[leakage[(observed - true) % 3] for true in range(3)] for observed in range(3)])
    observed = np.array([[80.0, 10.0, 10.0], [8.0, 80.0, 12.0]])

    linear = apply_confusion(observed, C, alpha=0.0)
    posterior, M_fit, B = fit_latent(observed, C, background="zero", update_M=False)

    assert np.allclose(M_fit, C)
    assert np.allclose(B, np.zeros(3))
    assert posterior.argmax(axis=1).tolist() == linear.argmax(axis=1).tolist()
    assert posterior[0, 0] >= linear[0, 0]
    assert posterior[1, 1] >= linear[1, 1]


def test_latent_em_is_count_weighted_not_sparse_codon_weighted():
    leakage = np.array([0.70, 0.10, 0.20])
    C = np.array([[leakage[(observed - true) % 3] for true in range(3)] for observed in range(3)])
    high_depth = np.array([[8000.0, 1000.0, 1000.0]])
    sparse_noise = np.tile(np.array([[0.0, 1.0, 0.0]]), (100, 1))
    observed = np.vstack([high_depth, sparse_noise])

    posterior, M_fit, _ = fit_latent(observed, C, background="zero", update_M=False)

    assert np.allclose(M_fit, C)
    assert posterior[0, 0] > 0.95


def test_frame_support_outputs_adjusted_counts_for_linear_and_latent():
    profiles = pl.DataFrame({
        "tran_id": ["tx1", "tx1", "tx1", "tx1", "tx1", "tx1"],
        "pos": [30, 31, 32, 33, 34, 35],
        "count": [70.0, 10.0, 20.0, 80.0, 10.0, 10.0],
    })
    cds = pl.DataFrame({"tran_id": ["tx1"], "start": [0], "stop": [90]})

    linear = build_frame_support(profiles, cds, FrameSupportParams(frame_method="linear", frame_by_length=False))
    latent = build_frame_support(profiles, cds, FrameSupportParams(frame_method="latent", frame_by_length=False))

    for support in (linear, latent):
        assert "observed_f0" in support.columns
        assert "adjusted_f0" in support.columns
        assert "total_count" in support.columns
        assert round(float(support["total_count"].sum()), 6) == 200.0

    codon_11 = linear.filter(pl.col("codon") == 11)
    assert codon_11["adjusted_f0"][0] > codon_11["observed_f0"][0]
    assert codon_11["p0"][0] > 0.9

    comparison, summary = compare_frame_support_tables(linear, latent)
    assert comparison.height == 2
    assert summary["n_joined"][0] == 2


def test_frame_validation_uses_cds_start_transcript_frame():
    support = pl.DataFrame({
        "tran_id": ["tx1", "tx1"],
        "codon": [1, 2],
        "p0": [0.05, 0.10],
        "p1": [0.90, 0.80],
        "p2": [0.05, 0.10],
        "entropy": [0.5, 0.7],
        "total_count": [10.0, 5.0],
    })
    cds = pl.DataFrame({"tran_id": ["tx1"], "start": [1], "stop": [30]})

    rows, summary = validate_frame_support_on_cds(support, cds, method="linear", trim_nt=0)

    assert rows["true_frame"].to_list() == [1, 1]
    assert rows["pred_frame"].to_list() == [1, 1]
    assert summary["n_rows"][0] == 2
    assert summary["argmax_accuracy"][0] == 1.0
    assert round(summary["weighted_mean_p_true"][0], 6) == round(((0.90 * 10.0) + (0.80 * 5.0)) / 15.0, 6)


def test_rdg_flux_export_builds_position_level_contract(tmp_path):
    profiles = pl.DataFrame({
        "tran_id": ["tx1", "tx1", "tx1"],
        "pos": [0, 1, 2],
        "count": [80.0, 10.0, 10.0],
    })
    frame_support = pl.DataFrame({
        "tran_id": ["tx1"],
        "codon": [0],
        "adjusted_f0": [98.0],
        "adjusted_f1": [2.0],
        "adjusted_f2": [0.0],
        "total_count": [100.0],
        "p0": [0.98],
        "p1": [0.02],
        "p2": [0.0],
        "entropy": [0.141],
        "method": ["linear"],
    })

    table = rdg_flux_position_table(
        profiles,
        frame_support,
        sample_id="ribocrypt_fwd_full",
        background_probability=0.1,
    )

    assert table.columns == [
        "sample_id",
        "transcript_id",
        "pos",
        "count",
        "p_frame0",
        "p_frame1",
        "p_frame2",
        "p_background",
        "frame_entropy",
        "effective_depth",
        "p_translated",
        "local_periodicity_score",
    ]
    assert table.height == 3
    assert table["sample_id"][0] == "ribocrypt_fwd_full"
    assert round(table["p_frame0"][0] + table["p_frame1"][0] + table["p_frame2"][0] + table["p_background"][0], 6) == 1.0
    assert table["p_frame0"][0] > 0.85
    assert table["effective_depth"][0] == 72.0

    profiles_path = tmp_path / "profiles.parquet"
    frame_path = tmp_path / "frame.parquet"
    out_path = tmp_path / "rdg.parquet"
    profiles.write_parquet(profiles_path)
    frame_support.write_parquet(frame_path)

    paths = export_rdg_flux_v1(
        profiles_path=str(profiles_path),
        frame_support_path=str(frame_path),
        out_path=str(out_path),
        sample_id="ribocrypt_fwd_full",
        transcriptome_fasta="/tmp/hg38.fa",
        annotation_source="/tmp/gencode.gtf",
        psite_offset_model="length_specific_offsets",
    )

    written = pl.read_parquet(paths["parquet"])
    with open(paths["metadata"], "r", encoding="utf-8") as fh:
        metadata = json.load(fh)
    assert written.height == 3
    assert metadata["coordinate_system"] == "transcript_0_based_half_open_positions"
    assert metadata["model_stage"] == "stage1_frame_posterior"
    assert metadata["sample_ids"] == ["ribocrypt_fwd_full"]

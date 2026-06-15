import numpy as np
import polars as pl

from TranslonScorer.frame.read_assignment import (
    assign_reads,
    assignment_identifiability,
    compare_read_assignment_methods,
    evaluate_assignments,
    frame_assignment_gate_diagnostics,
    prepare_assignment_candidates,
    summarize_identifiability,
)


def _metric_value(assignments: pl.DataFrame, name: str) -> float:
    metrics = evaluate_assignments(assignments)
    return float(metrics[name][0])


def test_assignment_frame_likelihood_uses_transcript_signal_frame_not_cds_relative_frame():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1"],
            "tran_id": ["tx_offset"],
            "tran_start_bam": [31],
        }
    )
    cds = pl.DataFrame({"tran_id": ["tx_offset"], "start": [1], "stop": [91]})
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_offset"],
            "codon": [10],
            "p0": [0.02],
            "p1": [0.96],
            "p2": [0.02],
        }
    )

    prepared = prepare_assignment_candidates(candidates, frame_support=frame_support, cds_tran=cds)

    assert prepared["signal_frame"][0] == 1
    assert prepared["cds_relative_frame"][0] == 0
    assert prepared["frame_likelihood"][0] > 0.9


def test_low_support_frame_evidence_is_neutral_for_read_assignment():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1"],
            "tran_id": ["tx_supported", "tx_unsupported"],
            "tran_start_bam": [30, 31],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_supported", "tx_unsupported"],
            "codon": [10, 10],
            "p0": [0.96, 0.96],
            "p1": [0.02, 0.02],
            "p2": [0.02, 0.02],
            "support_evidence": [1.0, 0.0],
        }
    )

    prepared = prepare_assignment_candidates(candidates, frame_support=frame_support).sort(
        "tran_id"
    )
    unsupported = prepared.filter(pl.col("tran_id") == "tx_unsupported")

    assert unsupported["raw_frame_likelihood"][0] == 0.02
    assert unsupported["frame_likelihood"][0] == 1.0


def test_min_frame_support_count_gates_weak_matched_frame_rows():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1"],
            "tran_id": ["tx_supported", "tx_low_count"],
            "tran_start_bam": [30, 31],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_supported", "tx_low_count"],
            "codon": [10, 10],
            "p0": [0.96, 0.96],
            "p1": [0.02, 0.02],
            "p2": [0.02, 0.02],
            "support_evidence": [1.0, 1.0],
            "total_count": [20.0, 2.0],
        }
    )

    prepared = prepare_assignment_candidates(
        candidates,
        frame_support=frame_support,
        min_frame_support_count=10.0,
    ).sort("tran_id")
    low_count = prepared.filter(pl.col("tran_id") == "tx_low_count")
    supported = prepared.filter(pl.col("tran_id") == "tx_supported")

    assert low_count["frame_support_passes_gate"][0] is False
    assert low_count["support_evidence"][0] == 0.0
    assert low_count["frame_likelihood"][0] == 1.0
    assert supported["frame_support_passes_gate"][0] is True
    assert supported["frame_likelihood"][0] > 0.9


def test_complete_frame_support_gate_prevents_no_support_target_bias():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1"],
            "tran_id": ["tx_supported", "tx_missing"],
            "tran_start_bam": [31, 30],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_supported"],
            "codon": [10],
            "p0": [0.96],
            "p1": [0.02],
            "p2": [0.02],
            "support_evidence": [1.0],
            "total_count": [200.0],
        }
    )

    partial = prepare_assignment_candidates(
        candidates,
        frame_support=frame_support,
        min_frame_support_count=10.0,
        require_complete_frame_support=False,
    ).sort("tran_id")
    complete = prepare_assignment_candidates(
        candidates,
        frame_support=frame_support,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    ).sort("tran_id")

    assert partial.filter(pl.col("tran_id") == "tx_supported")["frame_likelihood"][0] == 0.02
    assert partial.filter(pl.col("tran_id") == "tx_missing")["frame_likelihood"][0] == 1.0
    assert complete["frame_support_comparable"].to_list() == [False, False]
    assert complete["frame_likelihood"].to_list() == [1.0, 1.0]
    assert complete["support_evidence"].to_list() == [0.0, 0.0]


def test_single_target_reads_remain_fixed_across_assignment_methods():
    candidates = pl.DataFrame(
        {
            "read_key": ["unique", "ambiguous", "ambiguous"],
            "tran_id": ["tx_fixed", "tx_a", "tx_b"],
            "tran_start_bam": [30, 30, 31],
            "count": [2.0, 1.0, 1.0],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_fixed", "tx_a", "tx_b"],
            "codon": [10, 10, 10],
            "p0": [0.96, 0.96, 0.02],
            "p1": [0.02, 0.02, 0.96],
            "p2": [0.02, 0.02, 0.02],
            "support_evidence": [1.0, 1.0, 1.0],
            "total_count": [100.0, 100.0, 100.0],
        }
    )

    for method in ["fractional", "em", "frame_em", "rdg_gated_frame_em"]:
        result = assign_reads(
            candidates,
            method=method,
            frame_support=frame_support,
            min_frame_support_count=10.0,
            require_complete_frame_support=True,
        )
        fixed = result.assignments.filter(pl.col("read_key") == "unique")

        assert fixed["posterior"].to_list() == [1.0]
        assert fixed["assigned_count"].to_list() == [2.0]


def test_frame_aware_em_beats_frame_blind_em_for_frame_discordant_isoform_candidates():
    rows = []
    for idx in range(25):
        rows.append(
            {
                "read_key": f"amb:{idx}",
                "tran_id": "tx_in_frame",
                "tran_start_bam": 30,
                "is_true": True,
            }
        )
        rows.append(
            {
                "read_key": f"amb:{idx}",
                "tran_id": "tx_shifted",
                "tran_start_bam": 31,
                "is_true": False,
            }
        )

    candidates = pl.DataFrame(rows)
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_in_frame", "tx_shifted"],
            "codon": [10, 10],
            "p0": [0.96, 0.96],
            "p1": [0.02, 0.02],
            "p2": [0.02, 0.02],
        }
    )

    blind = assign_reads(candidates, method="em", frame_support=frame_support)
    aware = assign_reads(candidates, method="frame_em", frame_support=frame_support)

    assert _metric_value(blind.assignments, "soft_true_posterior") == 0.5
    assert _metric_value(aware.assignments, "soft_true_posterior") > 0.95


def test_frame_aware_em_matches_frame_blind_em_when_frame_is_uninformative():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1", "r2", "r2"],
            "tran_id": ["tx_a", "tx_b", "tx_a", "tx_b"],
            "tran_start_bam": [30, 31, 33, 34],
            "is_true": [True, False, True, False],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_a", "tx_a", "tx_b", "tx_b"],
            "codon": [10, 11, 10, 11],
            "p0": [1 / 3, 1 / 3, 1 / 3, 1 / 3],
            "p1": [1 / 3, 1 / 3, 1 / 3, 1 / 3],
            "p2": [1 / 3, 1 / 3, 1 / 3, 1 / 3],
        }
    )

    blind = assign_reads(candidates, method="em", frame_support=frame_support).assignments.sort(
        ["read_key", "tran_id"]
    )
    aware = assign_reads(
        candidates, method="frame_em", frame_support=frame_support
    ).assignments.sort(["read_key", "tran_id"])

    assert np.allclose(blind["posterior"].to_numpy(), aware["posterior"].to_numpy())


def test_frame_aware_em_does_not_resolve_same_frame_isoform_aliases():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1", "r2", "r2"],
            "tran_id": ["tx_true", "tx_same_frame", "tx_true", "tx_same_frame"],
            "tran_start_bam": [30, 60, 33, 63],
            "is_true": [True, False, True, False],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_true", "tx_same_frame", "tx_true", "tx_same_frame"],
            "codon": [10, 20, 11, 21],
            "p0": [0.96, 0.96, 0.96, 0.96],
            "p1": [0.02, 0.02, 0.02, 0.02],
            "p2": [0.02, 0.02, 0.02, 0.02],
            "support_evidence": [1.0, 1.0, 1.0, 1.0],
            "total_count": [100.0, 100.0, 100.0, 100.0],
        }
    )

    blind = assign_reads(candidates, method="em", frame_support=frame_support)
    aware = assign_reads(
        candidates,
        method="frame_em",
        frame_support=frame_support,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )

    assert _metric_value(blind.assignments, "soft_true_posterior") == 0.5
    assert _metric_value(aware.assignments, "soft_true_posterior") == 0.5


def test_dual_frame_support_keeps_frame_discordant_candidates_ambiguous():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1"],
            "tran_id": ["tx_true", "tx_overlap"],
            "tran_start_bam": [30, 31],
            "is_true": [True, False],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_true", "tx_overlap"],
            "codon": [10, 10],
            "p0": [0.49, 0.49],
            "p1": [0.49, 0.49],
            "p2": [0.02, 0.02],
            "support_evidence": [1.0, 1.0],
            "total_count": [100.0, 100.0],
        }
    )

    blind = assign_reads(candidates, method="em", frame_support=frame_support)
    aware = assign_reads(
        candidates,
        method="frame_em",
        frame_support=frame_support,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )

    assert _metric_value(blind.assignments, "soft_true_posterior") == 0.5
    assert _metric_value(aware.assignments, "soft_true_posterior") == 0.5


def test_locus_em_uses_unique_reads_to_resolve_genomic_multimappers():
    rows = []
    for idx in range(20):
        rows.append(
            {
                "read_key": f"uniq_a:{idx}",
                "locus_id": "locus_a",
                "tran_id": "tx_a",
                "tran_start_bam": 30,
                "is_true": True,
            }
        )
    for idx in range(2):
        rows.append(
            {
                "read_key": f"uniq_b:{idx}",
                "locus_id": "locus_b",
                "tran_id": "tx_b",
                "tran_start_bam": 30,
                "is_true": True,
            }
        )
    for idx in range(20):
        rows.append(
            {
                "read_key": f"multi:{idx}",
                "locus_id": "locus_a",
                "tran_id": "tx_a",
                "tran_start_bam": 30,
                "is_true": True,
            }
        )
        rows.append(
            {
                "read_key": f"multi:{idx}",
                "locus_id": "locus_b",
                "tran_id": "tx_b",
                "tran_start_bam": 30,
                "is_true": False,
            }
        )

    candidates = pl.DataFrame(rows)
    fractional = assign_reads(candidates, method="fractional", abundance_key="locus_id")
    em = assign_reads(candidates, method="em", abundance_key="locus_id")

    fractional_soft = _metric_value(fractional.assignments, "soft_true_posterior")
    em_soft = _metric_value(em.assignments, "soft_true_posterior")
    ambiguous_em = em.assignments.filter(
        (pl.col("read_key") == "multi:0") & (pl.col("locus_id") == "locus_a")
    )

    assert fractional_soft < em_soft
    assert ambiguous_em["posterior"][0] > 0.85


def test_composite_abundance_key_supports_joint_locus_transcript_assignment():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1"],
            "locus_id": ["locus_a", "locus_b"],
            "tran_id": ["tx_a", "tx_b"],
            "tran_start_bam": [30, 30],
        }
    )

    result = assign_reads(candidates, method="fractional", abundance_key="locus_id,tran_id")

    assert set(result.assignments["assignment_target"].to_list()) == {
        "locus_a|tx_a",
        "locus_b|tx_b",
    }
    assert np.allclose(result.assignments["posterior"].to_numpy(), [0.5, 0.5])


def test_assignment_collapses_duplicate_targets_before_target_posterior():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1", "r1"],
            "locus_id": ["locus_a", "locus_a", "locus_b"],
            "tran_id": ["tx_a1", "tx_a2", "tx_b1"],
            "tran_start_bam": [30, 30, 30],
        }
    )

    result = assign_reads(candidates, method="fractional", abundance_key="locus_id").assignments
    target_mass = (
        result.group_by("assignment_target")
        .agg(pl.col("posterior").sum().alias("target_posterior"))
        .sort("assignment_target")
    )

    assert np.allclose(target_mass["target_posterior"].to_numpy(), [0.5, 0.5])


def test_assignment_identifiability_classes_reads_by_target_uncertainty():
    candidates = pl.DataFrame(
        {
            "read_key": ["unique", "amb", "amb", "resolved", "resolved"],
            "tran_id": ["tx_a", "tx_a", "tx_b", "tx_a", "tx_b"],
            "tran_start_bam": [30, 30, 30, 30, 31],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_a", "tx_b"],
            "codon": [10, 10],
            "p0": [0.96, 0.96],
            "p1": [0.02, 0.02],
            "p2": [0.02, 0.02],
        }
    )

    result = assign_reads(candidates, method="frame_fractional", frame_support=frame_support)
    classes = assignment_identifiability(result.assignments)
    rows = {row["read_key"]: row["identifiability_class"] for row in classes.iter_rows(named=True)}
    summary = summarize_identifiability(classes)

    assert rows["unique"] == "unique"
    assert rows["amb"] == "ambiguous"
    assert rows["resolved"] == "resolved"
    assert summary["weighted_fraction"].sum() == 1.0


def test_rdg_gated_frame_em_keeps_clean_frame_disambiguation_gain():
    rows = []
    for idx in range(25):
        rows.append(
            {"read_key": f"amb:{idx}", "tran_id": "tx_true", "tran_start_bam": 30, "is_true": True}
        )
        rows.append(
            {
                "read_key": f"amb:{idx}",
                "tran_id": "tx_shifted",
                "tran_start_bam": 31,
                "is_true": False,
            }
        )

    candidates = pl.DataFrame(rows)
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_true", "tx_shifted"],
            "codon": [10, 10],
            "p0": [0.96, 0.96],
            "p1": [0.02, 0.02],
            "p2": [0.02, 0.02],
            "support_evidence": [1.0, 1.0],
            "total_count": [100.0, 100.0],
        }
    )

    blind = assign_reads(candidates, method="em", frame_support=frame_support)
    gated = assign_reads(
        candidates,
        method="rdg_gated_frame_em",
        frame_support=frame_support,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )

    assert _metric_value(blind.assignments, "soft_true_posterior") == 0.5
    assert _metric_value(gated.assignments, "soft_true_posterior") > 0.95
    assert gated.diagnostics is not None
    assert gated.diagnostics["rdg_scenario_frame_gate"][0] is True
    assert gated.diagnostics["rdg_gate_fraction"][0] == 1.0


def test_rdg_gated_frame_em_reports_alias_group_and_falls_back_to_em():
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1", "r1"],
            "tran_id": ["tx_true", "tx_same_frame", "tx_shifted"],
            "tran_start_bam": [30, 60, 31],
            "is_true": [True, False, False],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_true", "tx_same_frame", "tx_shifted"],
            "codon": [10, 20, 10],
            "p0": [0.96, 0.96, 0.96],
            "p1": [0.02, 0.02, 0.02],
            "p2": [0.02, 0.02, 0.02],
            "support_evidence": [1.0, 1.0, 1.0],
            "total_count": [100.0, 100.0, 100.0],
        }
    )
    prepared = prepare_assignment_candidates(
        candidates,
        frame_support=frame_support,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )
    gate, summary = frame_assignment_gate_diagnostics(prepared)
    blind = assign_reads(candidates, method="em", frame_support=frame_support)
    gated = assign_reads(
        candidates,
        method="rdg_gated_frame_em",
        frame_support=frame_support,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )

    assert gate["rdg_local_frame_gate"][0] is False
    assert gate["frame_separates_only_an_alias_group"][0] is True
    assert summary["rdg_scenario_frame_gate"][0] is False
    assert np.allclose(
        blind.assignments.sort("tran_id")["posterior"].to_numpy(),
        gated.assignments.sort("tran_id")["posterior"].to_numpy(),
    )


def test_rdg_consensus_gate_blocks_sparse_misleading_frame_signal():
    rows = []
    for idx in range(99):
        rows.append(
            {
                "read_key": f"alias:{idx}",
                "tran_id": "tx_true",
                "tran_start_bam": 30,
                "is_true": True,
            }
        )
        rows.append(
            {
                "read_key": f"alias:{idx}",
                "tran_id": "tx_wrong_a",
                "tran_start_bam": 31,
                "is_true": False,
            }
        )
        rows.append(
            {
                "read_key": f"alias:{idx}",
                "tran_id": "tx_wrong_b",
                "tran_start_bam": 34,
                "is_true": False,
            }
        )
    rows.append(
        {"read_key": "sparse_unique", "tran_id": "tx_true", "tran_start_bam": 30, "is_true": True}
    )
    rows.append(
        {
            "read_key": "sparse_unique",
            "tran_id": "tx_wrong_unique",
            "tran_start_bam": 31,
            "is_true": False,
        }
    )

    candidates = pl.DataFrame(rows)
    support_rows = (
        candidates.select(
            [
                "tran_id",
                (pl.col("tran_start_bam") // 3).alias("codon"),
            ]
        )
        .unique()
        .with_columns(
            [
                pl.lit(0.02).alias("p0"),
                pl.lit(0.96).alias("p1"),
                pl.lit(0.02).alias("p2"),
                pl.lit(1.0).alias("support_evidence"),
                pl.lit(100.0).alias("total_count"),
            ]
        )
    )

    blind = assign_reads(candidates, method="em", frame_support=support_rows)
    frame = assign_reads(
        candidates,
        method="frame_em",
        frame_support=support_rows,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )
    local = assign_reads(
        candidates,
        method="rdg_local_frame_em",
        frame_support=support_rows,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )
    consensus = assign_reads(
        candidates,
        method="rdg_gated_frame_em",
        frame_support=support_rows,
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )

    blind_truth = _metric_value(blind.assignments, "soft_true_posterior")
    assert _metric_value(frame.assignments, "soft_true_posterior") < blind_truth - 0.2
    assert _metric_value(local.assignments, "soft_true_posterior") < blind_truth
    assert np.isclose(_metric_value(consensus.assignments, "soft_true_posterior"), blind_truth)
    assert consensus.diagnostics is not None
    assert consensus.diagnostics["rdg_scenario_frame_gate"][0] is False
    assert consensus.diagnostics["rdg_gate_fraction"][0] < 0.25
    assert consensus.diagnostics["alias_group_contrast_fraction"][0] > 0.9


def test_compare_read_assignment_writes_rdg_gate_diagnostics(tmp_path):
    candidates = pl.DataFrame(
        {
            "read_key": ["r1", "r1"],
            "tran_id": ["tx_true", "tx_shifted"],
            "tran_start_bam": [30, 31],
            "is_true": [True, False],
        }
    )
    frame_support = pl.DataFrame(
        {
            "tran_id": ["tx_true", "tx_shifted"],
            "codon": [10, 10],
            "p0": [0.96, 0.96],
            "p1": [0.02, 0.02],
            "p2": [0.02, 0.02],
            "support_evidence": [1.0, 1.0],
            "total_count": [100.0, 100.0],
        }
    )
    candidates_path = tmp_path / "candidates.csv"
    frame_support_path = tmp_path / "frame_support.csv"
    out_prefix = tmp_path / "assignment"
    candidates.write_csv(candidates_path)
    frame_support.write_csv(frame_support_path)

    paths = compare_read_assignment_methods(
        candidates_path=str(candidates_path),
        frame_support_path=str(frame_support_path),
        out_prefix=str(out_prefix),
        methods=["em", "rdg_gated_frame_em"],
        min_frame_support_count=10.0,
        require_complete_frame_support=True,
    )
    diagnostics = pl.read_csv(paths["diagnostics"])

    assert "diagnostics" in paths
    assert diagnostics["method"].to_list() == ["rdg_gated_frame_em"]
    assert diagnostics["rdg_scenario_frame_gate"][0] is True

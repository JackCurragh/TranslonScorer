"""Tests for per-translon report composition and dynamic consequentiality.

Covers report.compose_report (scores + feature_event → per-translon rows, with
shared events giving identical aspect scores to both translons) and
consequential.apply_policy (tier-confidence × expression, no hard gates).
"""

from __future__ import annotations

import polars as pl

from TranslonScorer.consequential import apply_policy
from TranslonScorer.model import ConsequentialityPolicy
from TranslonScorer.report import compose_report


def _scores() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "event_id": [1, 2, 3, 4],
            "aspect": ["init", "elongation", "term", "init"],
            "call": ["SUPPORTED", "SUPPORTED", "SUPPORTED", "UNSUPPORTED"],
            "eligibility": ["ELIGIBLE"] * 4,
            "n_reads": [100.0, 200.0, 50.0, 5.0],
            "metric": [2.0, 0.8, 1.5, 0.0],
            "group": ["g"] * 4,
            "tier": ["aggregate"] * 4,
        }
    )


def _feature_event() -> pl.DataFrame:
    # Translon A: events 1(init), 2(elong), 3(term)
    # Translon B: events 4(init), 2(elong), 3(term)  — shares elong+term with A
    return pl.DataFrame(
        {
            "feature_id": ["A", "A", "A", "B", "B", "B"],
            "event_id": [1, 2, 3, 4, 2, 3],
            "role": ["init", "elongation", "term", "init", "elongation", "term"],
        }
    )


def test_compose_report_passthrough_without_feature_event():
    scores = _scores()
    assert compose_report(scores).shape == scores.shape


def test_compose_report_per_translon():
    rep = compose_report(_scores(), _feature_event())
    assert set(rep["feature_id"]) == {"A", "B"}
    a = rep.filter(pl.col("feature_id") == "A").row(0, named=True)
    b = rep.filter(pl.col("feature_id") == "B").row(0, named=True)

    assert a["init_call"] == "SUPPORTED"
    assert b["init_call"] == "UNSUPPORTED"
    # Shared elongation event → identical elong score for both translons.
    assert a["elongation_metric"] == b["elongation_metric"] == 0.8
    assert a["term_call"] == b["term_call"] == "SUPPORTED"
    assert a["total_reads"] == 350.0
    assert b["total_reads"] == 255.0


def test_apply_policy_ranks_by_confidence_and_expression():
    rep = compose_report(_scores(), _feature_event())
    out = apply_policy(rep, ConsequentialityPolicy())
    assert "consequentiality_score" in out.columns
    assert "consequential" in out.columns

    scores = dict(zip(out["feature_id"], out["consequentiality_score"]))
    # A (full chain supported + higher expression) outranks B (init unsupported).
    assert scores["A"] > scores["B"]
    # Default policy has zero floors → nothing is hard-gated out.
    assert out["consequential"].all()


def test_apply_policy_confidence_floor_gates_b_not_a():
    rep = compose_report(_scores(), _feature_event())
    out = apply_policy(rep, ConsequentialityPolicy(min_tier_confidence=0.8))
    flags = dict(zip(out["feature_id"], out["consequential"]))
    assert flags["A"] is True
    assert flags["B"] is False


def test_apply_policy_no_hard_length_gate():
    """A single-event start-stop translon with strong support stays consequential."""
    scores = pl.DataFrame(
        {
            "event_id": [1, 2, 3],
            "aspect": ["init", "elongation", "term"],
            "call": ["SUPPORTED", "SUPPORTED", "SUPPORTED"],
            "eligibility": ["ELIGIBLE"] * 3,
            "n_reads": [30.0, 25.0, 22.0],
            "metric": [1.8, 0.9, 1.4],
            "group": ["g"] * 3,
            "tier": ["aggregate"] * 3,
        }
    )
    fe = pl.DataFrame(
        {"feature_id": ["tiny"] * 3, "event_id": [1, 2, 3], "role": ["init", "elongation", "term"]}
    )
    out = apply_policy(compose_report(scores, fe), ConsequentialityPolicy())
    assert out.filter(pl.col("feature_id") == "tiny")["consequential"][0] is True


def test_apply_policy_empty():
    out = apply_policy(pl.DataFrame(schema={"feature_id": pl.Utf8}))
    assert "consequential" in out.columns
    assert out.is_empty()

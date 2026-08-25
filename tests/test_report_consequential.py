"""Tests for per-translon report composition and dynamic consequentiality.

Covers report.compose_report (scores + feature_event → per-translon rows, with
shared events giving identical aspect scores to both translons) and
consequential.apply_policy (tier-confidence × expression, no hard gates).
"""

from __future__ import annotations

import polars as pl
import pytest

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


# ---------------------------------------------------------------------------
# §6: confidence-weighted tier_confidence (significance_testing_plan.md §6)
# ---------------------------------------------------------------------------


def _scores_with_confidence(init_confidence: float) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "event_id": [1, 2, 3],
            "aspect": ["init", "elongation", "term"],
            "call": ["SUPPORTED", "SUPPORTED", "SUPPORTED"],
            "eligibility": ["ELIGIBLE"] * 3,
            "n_reads": [30.0, 25.0, 22.0],
            "metric": [1.8, 0.9, 1.4],
            "confidence": [init_confidence, 1.0, 1.0],
            "group": ["g"] * 3,
            "tier": ["aggregate"] * 3,
        }
    )


def _fe() -> pl.DataFrame:
    return pl.DataFrame(
        {"feature_id": ["t"] * 3, "event_id": [1, 2, 3], "role": ["init", "elongation", "term"]}
    )


def test_compose_report_carries_per_aspect_confidence():
    out = compose_report(_scores_with_confidence(0.5), _fe())
    assert out["init_confidence"][0] == pytest.approx(0.5)
    assert out["elongation_confidence"][0] == pytest.approx(1.0)


def test_tier_confidence_weights_by_aspect_confidence():
    """All three chain aspects SUPPORTED either way -- a flat mean-of-SUPPORTED
    would give tier_confidence=1.0 regardless of how confident each call was.
    Weighting by confidence should pull it below 1.0 when one aspect's own
    battery barely agreed."""
    low = apply_policy(compose_report(_scores_with_confidence(0.1), _fe()))
    high = apply_policy(compose_report(_scores_with_confidence(1.0), _fe()))
    assert high["tier_confidence"][0] == pytest.approx(1.0)
    assert low["tier_confidence"][0] < high["tier_confidence"][0]


def test_tier_confidence_falls_back_to_equal_weight_without_confidence_column():
    """No confidence column at all (pre-§6 reports) -- exactly the original
    flat mean-of-SUPPORTED behaviour."""
    scores = _scores_with_confidence(0.1).drop("confidence")
    out = apply_policy(compose_report(scores, _fe()))
    assert out["tier_confidence"][0] == pytest.approx(1.0)

"""Translon-level CIF: exact per-block preservation (compose_block_detail)
and the codon-weighted (not read-weighted) approximate composite
(_compose_per_translon's elongation_cif_approx).

Deliberately uses blocks with mismatched codon-count vs read-count so a
read-weighted average and a codon-weighted average disagree — proving the
composite is actually weighting by codons, not silently falling back to the
read-weighting `metric` already uses elsewhere.
"""

from __future__ import annotations

import polars as pl
import pytest

from TranslonScorer.report import compose_block_detail, compose_report


def _scores() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "event_id": pl.Series([1, 2], dtype=pl.UInt64),
            "aspect": ["elongation", "elongation"],
            "call": ["SUPPORTED", "UNSUPPORTED"],
            "eligibility": ["ELIGIBLE", "ELIGIBLE"],
            "n_reads": [1000.0, 100.0],
            "metric": [0.9, 0.2],
            "cif": [1.0, 0.0],
            "n_codons": pl.Series([10, 90], dtype=pl.Int64),
        }
    )


def _feature_event() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "feature_id": ["translonA", "translonA"],
            "event_id": pl.Series([1, 2], dtype=pl.UInt64),
            "rank": pl.Series([1, 2], dtype=pl.Int64),
        }
    )


def test_compose_block_detail_preserves_each_block_exactly():
    detail = compose_block_detail(_scores(), _feature_event())
    assert detail.sort("rank")["cif"].to_list() == [1.0, 0.0]
    assert detail.sort("rank")["n_codons"].to_list() == [10, 90]
    assert detail.sort("rank")["event_id"].to_list() == [1, 2]


def test_compose_block_detail_without_rank_column():
    fe = _feature_event().drop("rank")
    detail = compose_block_detail(_scores(), fe)
    assert "rank" not in detail.columns
    assert set(detail["event_id"].to_list()) == {1, 2}


def test_elongation_cif_approx_is_codon_weighted_not_read_weighted():
    report = compose_report(_scores(), _feature_event())
    row = report.filter(pl.col("feature_id") == "translonA").row(0, named=True)

    codon_weighted = (1.0 * 10 + 0.0 * 90) / (10 + 90)
    read_weighted = (1.0 * 1000 + 0.0 * 100) / (1000 + 100)
    assert codon_weighted == pytest.approx(0.1)
    assert read_weighted == pytest.approx(0.909, abs=1e-3)
    # The two disagree by construction -- this is the point of the fixture.
    assert abs(codon_weighted - read_weighted) > 0.5

    assert row["elongation_cif_approx"] == pytest.approx(codon_weighted)
    # elongation_metric (n_reads-weighted, unchanged) should NOT match the
    # cif composite -- they are different weightings on purpose.
    assert row["elongation_cif_approx"] != pytest.approx(row["elongation_metric"])


def test_cif_approx_absent_for_non_elongation_aspects():
    scores = pl.DataFrame(
        {
            "event_id": pl.Series([1, 2, 3], dtype=pl.UInt64),
            "aspect": ["init", "elongation", "term"],
            "call": ["SUPPORTED", "SUPPORTED", "SUPPORTED"],
            "eligibility": ["ELIGIBLE"] * 3,
            "n_reads": [50.0, 1000.0, 50.0],
            "metric": [2.0, 0.9, 2.0],
            "cif": [None, 1.0, None],
            "n_codons": pl.Series([None, 10, None], dtype=pl.Int64),
        }
    )
    fe = pl.DataFrame(
        {
            "feature_id": ["translonA"] * 3,
            "event_id": pl.Series([1, 2, 3], dtype=pl.UInt64),
            "rank": pl.Series([0, 1, 2], dtype=pl.Int64),
        }
    )
    report = compose_report(scores, fe)
    assert "elongation_cif_approx" in report.columns
    assert "init_cif_approx" not in report.columns
    assert "term_cif_approx" not in report.columns

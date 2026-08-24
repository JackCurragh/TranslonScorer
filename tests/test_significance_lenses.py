"""New significance lenses from docs/significance_testing_plan.md that don't
fit naturally alongside the boundary-axes tests (elongation, CIF,
attribution). Termination dropoff significance/continuity lives in
test_boundary_axes.py instead, next to the periodicity test it reuses.
"""

from __future__ import annotations

import numpy as np
import pytest

from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.aspects import _elong_frame_chisq, score_elongation_event
from TranslonScorer.scoring.signature import cif_codon_contiguity, cif_codon_significance

_THR = ScoreThresholds()


# ---------------------------------------------------------------------------
# Elongation level lens: chi-square goodness-of-fit (significance_testing_plan §2)
# ---------------------------------------------------------------------------


def test_elong_frame_chisq_detects_real_periodicity():
    p = _elong_frame_chisq([300.0, 20.0, 15.0])
    assert p is not None
    assert p < 0.01


def test_elong_frame_chisq_none_for_uniform_tally():
    p = _elong_frame_chisq([100.0, 100.0, 100.0])
    assert p is None or p > 0.05


def test_elong_frame_chisq_none_with_no_signal():
    assert _elong_frame_chisq([0.0, 0.0, 0.0]) is None


def test_elong_frame_chisq_catches_split_off_frame_signal():
    # Frame 0 is NOT dominant, but frames 1/2 are unevenly split -- still
    # non-uniform overall, which a frame-0-vs-rest binomial test would miss.
    p = _elong_frame_chisq([34.0, 33.0, 100.0])
    assert p is not None
    assert p < 0.01


def test_score_elongation_event_carries_frame_chisq_p():
    coverage = {p: (10.0 if p % 3 == 0 else 1.0) for p in range(0, 300)}
    s = score_elongation_event(0, 300, 0, 1, coverage, thr=_THR)
    assert "frame_chisq_p" in s
    assert s["frame_chisq_p"] is not None
    assert s["frame_chisq_p"] < 0.01


# ---------------------------------------------------------------------------
# CIF per-codon significance (significance_testing_plan §3)
# ---------------------------------------------------------------------------


def test_cif_codon_significance_none_for_wrong_length():
    assert cif_codon_significance(np.zeros(7)) is None


def test_cif_codon_significance_no_testable_codons_below_floor():
    # 30 codons, 5 reads each -- below the default min_reads=10 floor.
    signal = np.tile([3.0, 1.0, 1.0], 30)
    result = cif_codon_significance(signal, min_reads=10)
    assert result is not None
    assert result["n_testable"] == 0
    assert result["frac_significant"] is None
    assert result["n_codons"] == 30


def test_cif_codon_significance_detects_dominant_codons():
    # 20 codons, each with 30 reads, 27 of them frame-0 (90%) -- clearly
    # significant vs the 1/3 null at n=30.
    signal = np.tile([27.0, 2.0, 1.0], 20)
    result = cif_codon_significance(signal, min_reads=10, alpha=0.05)
    assert result["n_testable"] == 20
    assert result["frac_significant"] == pytest.approx(1.0)
    assert all(result["sig_mask"])


def test_cif_codon_significance_low_for_uniform_codons():
    # Depth clears the floor, but frame-0 share is exactly the chance rate.
    signal = np.tile([10.0, 10.0, 10.0], 20)
    result = cif_codon_significance(signal, min_reads=10, alpha=0.05)
    assert result["n_testable"] == 20
    assert result["frac_significant"] == pytest.approx(0.0)


def test_cif_codon_contiguity_none_with_no_significant_codons():
    assert cif_codon_contiguity([False] * 10) is None


def test_cif_codon_contiguity_clustered_vs_scattered():
    clustered = cif_codon_contiguity([False, True, True, True, True, False, False, False])
    scattered = cif_codon_contiguity([True, False, True, False, True, False, True, False])
    assert clustered["max_run_frac"] == pytest.approx(1.0)
    assert clustered["n_runs"] == 1
    assert scattered["max_run_frac"] == pytest.approx(0.25)
    assert scattered["n_runs"] == 4

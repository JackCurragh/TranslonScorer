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

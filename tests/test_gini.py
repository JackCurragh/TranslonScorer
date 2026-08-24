"""Gini coefficient over a per-bin signal vector.

Not a Chothani/Sonia reference score — a locally-added uniformity metric for
the boundary-detection axes, checked against hand-computed values rather than
an external reference (see ``test_signature_scores.py`` for the R-reference
checked scores).
"""

from __future__ import annotations

import numpy as np
import pytest

from TranslonScorer.scoring.signature import gini


def test_gini_uniform_signal_is_zero():
    assert gini(np.array([4.0, 4.0, 4.0, 4.0])) == pytest.approx(0.0, abs=1e-9)


def test_gini_single_spike_approaches_one():
    # All mass in one bin among many: (n+1-2)/n -> 1 - 1/n, approaching 1 as
    # n grows, not exactly 1 (the standard discrete Gini formula's known bias).
    x = np.array([0.0, 0.0, 0.0, 10.0])
    n = len(x)
    assert gini(x) == pytest.approx(1 - 1 / n)


def test_gini_matches_hand_computed_two_point_example():
    # Mean absolute difference / (2 * n * mean), hand-computed directly from
    # the definition rather than the cumulative-sum formula under test.
    x = np.array([1.0, 3.0])
    mad = abs(1.0 - 3.0) * 2 / (len(x) ** 2)  # mean absolute difference over all pairs incl. self
    expected = mad / (2 * np.mean(x))
    assert gini(x) == pytest.approx(expected)


def test_gini_is_none_when_no_signal():
    assert gini(np.zeros(5)) is None


def test_gini_is_none_on_empty_input():
    assert gini(np.array([])) is None


def test_gini_ignores_nan_entries():
    assert gini(np.array([4.0, 4.0, np.nan, 4.0])) == pytest.approx(0.0, abs=1e-9)

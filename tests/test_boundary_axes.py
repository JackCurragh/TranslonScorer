"""Multi-axis init/term boundary detection: periodicity, uniformity
(peakiness + Gini), breadth, and the periodicity significance test.

These are evidence, not headline metrics (see scoring_model.md) -- riding in
`boundary_axes_by_flank`, per flank length, both sides of the boundary. The
one place they *can* change behaviour is `_decide_step`'s borderline band,
covered separately in this file and in test_splice_aware_flank.py.
"""

from __future__ import annotations

import numpy as np
import pytest

from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.aspects import (
    _boundary_axes_for_flank,
    _breadth,
    _codon_bins,
    _peakiness,
    _periodicity,
    _periodicity_significance,
    score_initiation_event,
    score_termination_event,
)

_THR = ScoreThresholds()


def test_codon_bins_frame0_is_first_of_each_triple():
    coverage = {0: 10.0, 1: 1.0, 2: 1.0, 3: 20.0, 4: 2.0, 5: 2.0}
    total, frame0 = _codon_bins(coverage, list(range(6)))
    assert total.tolist() == [12.0, 24.0]
    assert frame0.tolist() == [10.0, 20.0]


def test_codon_bins_drops_trailing_partial_codon():
    coverage = {p: 1.0 for p in range(8)}
    total, _ = _codon_bins(coverage, list(range(8)))
    assert len(total) == 2  # 8 // 3, trailing 2 positions dropped


def test_periodicity_matches_pif_style_ratio():
    total = np.array([4.0, 4.0])
    frame0 = np.array([3.0, 1.0])
    assert _periodicity(total, frame0) == pytest.approx(4 / 8)


def test_periodicity_is_none_with_no_signal():
    assert _periodicity(np.zeros(3), np.zeros(3)) is None


def test_peakiness_flat_signal_is_one():
    assert _peakiness(np.array([4.0, 4.0, 4.0])) == pytest.approx(1.0)


def test_peakiness_empty_is_zero():
    assert _peakiness(np.array([])) == 0.0


def test_breadth_counts_nonzero_bins():
    assert _breadth(np.array([1.0, 0.0, 2.0, 0.0])) == pytest.approx(0.5)


def test_breadth_empty_is_zero():
    assert _breadth(np.array([])) == 0.0


def test_periodicity_significance_none_below_codon_floor():
    small = np.array([1.0, 1.0])
    p = _periodicity_significance(small, small, small, small, min_codons=5)
    assert p is None


def test_periodicity_significance_detects_a_real_change():
    # Body: clean triplet periodicity. Outer: flat, non-periodic.
    body_total = np.array([40.0] * 10)
    body_frame0 = np.array([30.0] * 10)
    out_total = np.array([24.0] * 10)
    out_frame0 = np.array([8.0] * 10)  # exactly 1/3 -> no frame-0 preference
    p = _periodicity_significance(body_total, body_frame0, out_total, out_frame0, min_codons=5)
    assert p is not None
    assert p < 0.01


def test_periodicity_significance_none_when_shares_identical():
    total = np.array([9.0] * 10)
    frame0 = np.array([3.0] * 10)  # exactly 1/3 both sides -- no change to detect
    p = _periodicity_significance(total, frame0, total, frame0, min_codons=5)
    # Not necessarily None (ties still produce a p-value), but should not be
    # a small/"significant" p-value for a genuinely identical distribution.
    assert p is None or p > 0.05


def test_boundary_axes_for_flank_shape():
    coverage = {p: (30.0 if p % 3 == 0 else 5.0) for p in range(0, 60)}
    coverage.update({p: 8.0 for p in range(-60, 0)})
    axes = _boundary_axes_for_flank(coverage, list(range(0, 60)), list(range(-60, 0)), min_codons=5)
    expected_keys = {
        "periodicity_body",
        "periodicity_outer",
        "periodicity_delta",
        "peakiness_body",
        "peakiness_outer",
        "gini_body",
        "gini_outer",
        "breadth_body",
        "breadth_outer",
        "periodicity_p",
    }
    assert set(axes) == expected_keys
    assert axes["periodicity_delta"] == pytest.approx(axes["periodicity_body"] - axes["periodicity_outer"])


def test_periodicity_resolves_ambiguous_on_minus_strand_too():
    """Mirror of the plus-strand case in test_splice_aware_flank.py: a real
    triplet-periodic body against a flat, non-periodic leader should resolve
    an AMBIGUOUS call the same way regardless of strand. Transcript distance
    from the start on minus strand is `(start_pos - p)`, not `p - start_pos`
    -- this is exactly the kind of off-by-one a strand-naive fixture gets
    wrong, so it's worth pinning down explicitly rather than only checking
    plus strand."""
    start_pos = 199
    coverage = {p: (30.0 if (start_pos - p) % 3 == 0 else 5.0) for p in range(100, 200)}
    coverage.update({p: 8.0 for p in range(200, 260)})  # leader: flat, non-periodic
    s = score_initiation_event(start_pos, -1, coverage, thr=_THR)
    assert s["metric"] < 2.0  # depth step alone is honestly weak
    assert s["call"] == "SUPPORTED"
    assert s["periodicity_resolved_ambiguous"] is True


def test_score_initiation_event_carries_boundary_axes_by_flank():
    coverage = {p: (30.0 if p % 3 == 0 else 5.0) for p in range(100, 200)}
    s = score_initiation_event(100, 1, coverage, thr=_THR)
    assert set(s["boundary_axes_by_flank"]) == {9, 18, 30, 60}
    for axes in s["boundary_axes_by_flank"].values():
        assert "periodicity_body" in axes


def test_score_termination_event_minus_strand_spliced_does_not_crash():
    """Termination, minus strand, AND a splice-projected flank together --
    the combination least likely to be exercised by any single existing
    test. Not asserting a specific call (that's covered elsewhere); this is
    a "the new axes survive the hardest combination of existing features"
    check."""
    coverage = {p: (30.0 if p % 3 == 0 else 5.0) for p in range(100, 200)}
    coverage.update({p: 6.0 for p in range(-40, 20)})
    splice_ctx = {("chr1", "-"): [(20, 100)]}
    s = score_termination_event(
        100, -1, coverage, thr=_THR, chrom="chr1", splice_context=splice_ctx
    )
    assert s["flank_spliced"] is True
    assert set(s["boundary_axes_by_flank"]) == {9, 18, 30, 60}

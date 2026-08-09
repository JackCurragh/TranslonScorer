"""Translation signature scores checked against the R reference.

Expected values are hand-computed from the formulas in
``Sonia_translation signature scores.txt`` (Chothani et al., Molecular Cell
2022), not from a previous run of this code, so these tests can actually
disagree with the implementation.

R uses 1-based indexing; the offsets below are the 0-based equivalents:

  PIF        sum(psites[seq(1,n,3)]) / sum(psites)          -> signal[0::3] / signal
  CIF        which(p1*100/tot > 33.33) / (n/3)              -> denominator is ALL codons
  drop-off   sum(w[seq(1,18,3)]) / sum(w[c(seq(1,18,3), seq(19,33,3))])
             -> numerator 0,3,6,9,12,15; denominator adds 18,21,24,27,30
"""

from __future__ import annotations

import numpy as np
import pytest

from TranslonScorer.scoring.signature import (
    DROPOFF_WINDOW_NT,
    cif,
    dropoff,
    dropoff_window_indices,
    pif,
    signal_from_intervals,
)

# ---------------------------------------------------------------------------
# PIF
# ---------------------------------------------------------------------------


def test_pif_perfect_periodicity_is_one():
    assert pif(np.array([3.0, 0, 0, 3, 0, 0])) == 1.0


def test_pif_uniform_signal_is_one_third():
    assert pif(np.array([1.0] * 6)) == pytest.approx(2 / 6)


def test_pif_denominator_is_all_positions_not_just_covered():
    # 2 in frame 0, 6 total -> 1/3, NOT 2/2 from the nonzero positions alone.
    assert pif(np.array([1.0, 2, 0, 1, 2, 0])) == pytest.approx(2 / 6)


def test_pif_is_none_when_there_is_no_signal():
    # None, not 0.0: "no signal" and "no in-frame signal" are different claims.
    assert pif(np.zeros(6)) is None


# ---------------------------------------------------------------------------
# CIF
# ---------------------------------------------------------------------------


def test_cif_perfect_periodicity_is_one():
    assert cif(np.array([3.0, 0, 0, 3, 0, 0])) == 1.0


def test_cif_counts_perfectly_uniform_codons_as_dominant():
    """The R threshold is the literal 33.33, not 100/3.

    An equal-thirds codon scores 33.333... which IS > 33.33, so it counts.
    Very likely unintended upstream, but it is what the published numbers
    reflect, so this test pins the quirk rather than silently fixing it.
    """
    assert cif(np.array([1.0, 1, 1])) == 1.0
    # Just below a third fails, confirming the boundary is real and not slop.
    assert cif(np.array([0.999, 1, 1])) == 0.0


def test_cif_denominator_includes_zero_signal_codons():
    """Uncovered codons penalise CIF; they are not excluded from the denominator.

    R: numerator drops NaN via which(); denominator stays nrow(psites)/3.
    """
    # codon 1 fully in frame, codon 2 empty -> 1 of 2, not 1 of 1.
    assert cif(np.array([3.0, 0, 0, 0, 0, 0])) == 0.5


def test_cif_is_none_for_non_codon_length():
    # R would recycle mismatched seq() vectors and return a wrong number.
    assert cif(np.array([1.0, 2, 3, 4])) is None
    assert cif(np.array([])) is None


# ---------------------------------------------------------------------------
# drop-off
# ---------------------------------------------------------------------------


def _window(before: float = 0.0, after: float = 0.0, other: float = 0.0) -> np.ndarray:
    """33-nt window with distinct values on the in-frame before/after offsets."""
    w = np.full(DROPOFF_WINDOW_NT, other, dtype=float)
    for i in range(0, 18, 3):
        w[i] = before
    for i in range(18, 33, 3):
        w[i] = after
    return w


def test_dropoff_all_signal_before_stop_is_one():
    assert dropoff(_window(before=10.0, after=0.0)) == 1.0


def test_dropoff_all_signal_after_stop_is_zero():
    assert dropoff(_window(before=0.0, after=10.0)) == 0.0


def test_dropoff_uses_six_before_and_five_after_offsets():
    # Equal per-position signal -> 6 positions before, 5 after.
    assert dropoff(_window(before=1.0, after=1.0)) == pytest.approx(6 / 11)


def test_dropoff_ignores_out_of_frame_positions():
    """Only offsets 0,3,...,15 and 18,...,30 contribute; the rest are noise."""
    w = _window(before=1.0, after=1.0, other=1000.0)
    assert dropoff(w) == pytest.approx(6 / 11)


def test_dropoff_is_none_when_denominator_is_zero():
    assert dropoff(np.zeros(DROPOFF_WINDOW_NT)) is None


def test_dropoff_rejects_a_wrong_sized_window():
    with pytest.raises(ValueError, match="33 nt"):
        dropoff(np.zeros(32))


# ---------------------------------------------------------------------------
# window placement
# ---------------------------------------------------------------------------


def test_dropoff_window_needs_17_before_and_15_after():
    assert dropoff_window_indices(17, 100) == range(0, 33)
    # One short upstream: R flags "15bp no upstream flank" at stop_coord <= 17
    # (1-based), i.e. 0-based index 16.
    assert dropoff_window_indices(16, 100) is None
    # One short downstream.
    assert dropoff_window_indices(17, 32) is None
    assert dropoff_window_indices(17, 33) == range(0, 33)


# ---------------------------------------------------------------------------
# interval lookup, and the R-compat divergence
# ---------------------------------------------------------------------------


def _track():
    # One interval covering 1-based positions 11..20 (BED start=10, end=20).
    return np.array([10]), np.array([20]), np.array([5.0])


def test_interval_lookup_is_half_open_start_exclusive_end_inclusive():
    starts, ends, values = _track()
    got = signal_from_intervals([10, 11, 20, 21], starts, ends, values)
    # R: bedgraph[,2] < x & bedgraph[,3] >= x
    assert got.tolist() == [0.0, 5.0, 5.0, 0.0]


def test_uncovered_positions_are_zero_and_keep_the_vector_aligned():
    starts, ends, values = _track()
    got = signal_from_intervals([11, 50, 12], starts, ends, values)
    assert got.tolist() == [5.0, 0.0, 5.0]


def test_strict_r_compat_drops_uncovered_positions_and_shifts_frame():
    """The reference's unlist() silently discards uncovered positions.

    This is the divergence that makes "emulate exactly" and "preserve zero
    signal" incompatible: on a sparse track the vector comes back short and
    every downstream position changes frame, so PIF and CIF differ.
    """
    starts, ends, values = _track()
    positions = [11, 50, 12]  # middle position uncovered
    corrected = signal_from_intervals(positions, starts, ends, values)
    strict = signal_from_intervals(positions, starts, ends, values, strict_r_compat=True)
    assert len(corrected) == 3
    assert len(strict) == 2, "strict mode must drop the uncovered position"
    assert strict.tolist() == [5.0, 5.0]


def test_strict_and_corrected_agree_on_a_dense_track():
    """On a track with no gaps the divergence never fires -- worth knowing,
    because it means the two modes only disagree on sparse input."""
    starts = np.array([0, 10, 20])
    ends = np.array([10, 20, 30])
    values = np.array([1.0, 2.0, 3.0])
    positions = list(range(1, 31))
    corrected = signal_from_intervals(positions, starts, ends, values)
    strict = signal_from_intervals(positions, starts, ends, values, strict_r_compat=True)
    assert corrected.tolist() == strict.tolist()


# ---------------------------------------------------------------------------
# agreement with the literal R transliteration (independent oracle)
# ---------------------------------------------------------------------------


def _random_signal(rng, n_codons: int, sparsity: float) -> np.ndarray:
    """Signal with a realistic mix of zeros, small counts and spikes."""
    x = rng.integers(0, 50, size=n_codons * 3).astype(float)
    x[rng.random(len(x)) < sparsity] = 0.0
    return x


@pytest.mark.parametrize("sparsity", [0.0, 0.3, 0.7, 0.95])
def test_pif_and_cif_match_the_r_transliteration(sparsity):
    from tests.reference_signature import r_cif, r_pif

    rng = np.random.default_rng(20260806)
    for _ in range(200):
        x = _random_signal(rng, int(rng.integers(1, 40)), sparsity)
        mine_pif, ref_pif = pif(x), r_pif(x.tolist())
        if ref_pif is None:
            assert mine_pif is None
        else:
            assert mine_pif == pytest.approx(ref_pif)
        mine_cif, ref_cif = cif(x), r_cif(x.tolist())
        if ref_cif is None:
            assert mine_cif is None
        else:
            assert mine_cif == pytest.approx(ref_cif)


def test_dropoff_matches_the_r_transliteration():
    from tests.reference_signature import r_dropoff

    rng = np.random.default_rng(7)
    for _ in range(500):
        w = rng.integers(0, 30, size=DROPOFF_WINDOW_NT).astype(float)
        mine, ref = dropoff(w), r_dropoff(w.tolist())
        if ref is None:
            assert mine is None
        else:
            assert mine == pytest.approx(ref)

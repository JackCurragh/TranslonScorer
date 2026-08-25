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


def test_a_dropped_whole_codon_preserves_frame_but_still_moves_cif():
    """The divergence is not uniform: what it costs depends on the gap length.

    A gap that is a multiple of 3 keeps seq(1,n,3) on frame 0, so PIF survives
    -- but CIF still shifts, because the dropped codons leave the denominator
    rather than penalising the score. Only a non-multiple-of-3 gap corrupts PIF.
    """
    # 30 nt, periodic (9,1,1 per codon), with codon 5 (positions 13-15) absent.
    cov = {}
    for c in range(10):
        if c == 4:
            continue
        cov[1 + 3 * c] = 9.0
        cov[2 + 3 * c] = 1.0
        cov[3 + 3 * c] = 1.0
    starts = np.array([p - 1 for p in cov])
    ends = np.array(list(cov))
    values = np.array(list(cov.values()), dtype=float)
    positions = list(range(1, 31))

    corrected = signal_from_intervals(positions, starts, ends, values)
    strict = signal_from_intervals(positions, starts, ends, values, strict_r_compat=True)
    assert len(strict) == 27, "a whole codon's worth of positions was dropped"

    # Frame survives a multiple-of-3 gap, so PIF agrees.
    assert pif(strict) == pytest.approx(pif(corrected))
    # CIF does not: the empty codon stops counting against the score.
    assert cif(corrected) == pytest.approx(0.9)
    assert cif(strict) == pytest.approx(1.0)


def test_a_dropped_partial_codon_shifts_frame_and_breaks_pif():
    """One uncovered position is enough to put everything downstream out of frame."""
    cov = {}
    for c in range(10):
        cov[1 + 3 * c] = 9.0
        cov[2 + 3 * c] = 1.0
        cov[3 + 3 * c] = 1.0
    del cov[14]  # a single interior position
    starts = np.array([p - 1 for p in cov])
    ends = np.array(list(cov))
    values = np.array(list(cov.values()), dtype=float)
    positions = list(range(1, 31))

    corrected = signal_from_intervals(positions, starts, ends, values)
    strict = signal_from_intervals(positions, starts, ends, values, strict_r_compat=True)
    assert len(strict) == 29
    assert pif(corrected) != pytest.approx(pif(strict))
    # 29 is not a whole number of codons, so CIF cannot be computed at all.
    assert cif(strict) is None


# ---------------------------------------------------------------------------
# transcript-order vector construction (the strand trap)
# ---------------------------------------------------------------------------


def _sparse(d):
    pos = np.array(sorted(d), dtype=np.int64)
    return pos, np.array([float(d[p]) for p in pos])


def test_plus_strand_vector_starts_on_a_codon_boundary():
    from TranslonScorer.scoring.signature import orf_signal_vector

    # Span [10,19), frame-0 residue 1 -> frame-0 genomic positions 10,13,16.
    cov = {10: 9.0, 11: 1.0, 12: 1.0, 13: 9.0, 14: 1.0, 15: 1.0, 16: 9.0, 17: 1.0, 18: 1.0}
    pos, cnt = _sparse(cov)
    v = orf_signal_vector(pos, cnt, 10, 19, expected_frame=10 % 3, strand=1)
    assert v.tolist() == [9, 1, 1, 9, 1, 1, 9, 1, 1]
    assert pif(v) == pytest.approx(27 / 33)
    assert cif(v) == 1.0


def test_minus_strand_vector_is_reversed_so_codons_group_correctly():
    """On the minus strand the codon's first base is the HIGHEST coordinate.

    Same signal shape as the plus-strand case but laid out 3'->5' genomically:
    the strong base of each codon sits at 18, 15, 12. Read in transcript order
    the vector must come back 9,1,1 repeated -- if the reversal or the offset
    were wrong, the 9s would land at index 1 or 2 and CIF would collapse.
    """
    from TranslonScorer.scoring.signature import orf_signal_vector

    cov = {18: 9.0, 17: 1.0, 16: 1.0, 15: 9.0, 14: 1.0, 13: 1.0, 12: 9.0, 11: 1.0, 10: 1.0}
    pos, cnt = _sparse(cov)
    v = orf_signal_vector(pos, cnt, 10, 19, expected_frame=18 % 3, strand=-1)
    assert v.tolist() == [9, 1, 1, 9, 1, 1, 9, 1, 1]
    assert cif(v) == 1.0


def test_vector_is_trimmed_to_whole_codons_at_both_ends():
    from TranslonScorer.scoring.signature import orf_signal_vector

    # Span starts one base before the first frame-0 position and ends two late.
    cov = {p: 1.0 for p in range(10, 22)}
    pos, cnt = _sparse(cov)
    v = orf_signal_vector(pos, cnt, 10, 22, expected_frame=11 % 3, strand=1)
    assert len(v) % 3 == 0
    assert len(v) == 9, "one leading base trimmed, two trailing"


def test_uncovered_positions_inside_the_orf_become_zero():
    from TranslonScorer.scoring.signature import orf_signal_vector

    pos, cnt = _sparse({10: 5.0, 16: 5.0})
    v = orf_signal_vector(pos, cnt, 10, 19, expected_frame=10 % 3, strand=1)
    assert v.tolist() == [5, 0, 0, 0, 0, 0, 5, 0, 0]
    # The empty middle codon counts against CIF rather than vanishing.
    assert cif(v) == pytest.approx(2 / 3)


def test_empty_span_returns_empty_vector():
    from TranslonScorer.scoring.signature import orf_signal_vector

    pos, cnt = _sparse({10: 1.0})
    assert orf_signal_vector(pos, cnt, 10, 10, expected_frame=1, strand=1).size == 0


# ---------------------------------------------------------------------------
# drop-off wired onto the termination aspect
# ---------------------------------------------------------------------------


def test_dropoff_window_is_frame_locked_to_term_pos_on_both_strands():
    """window[17] must be term_pos itself, or every offset is out of register."""
    from TranslonScorer.scoring.aspects import _dropoff_at

    # All signal on the six in-frame positions up to and including the stop.
    for strand in (1, -1):
        term = 1000
        cov = {}
        for k in range(0, 18, 3):
            cov[term - (17 - k) * strand] = 5.0
        assert _dropoff_at(term, strand, cov) == 1.0, f"strand {strand}"

    # All signal after the stop -> 0.0.
    for strand in (1, -1):
        term = 1000
        cov = {term + (k - 17) * strand: 5.0 for k in range(18, 33, 3)}
        assert _dropoff_at(term, strand, cov) == 0.0, f"strand {strand}"


def test_dropoff_is_none_off_contig():
    from TranslonScorer.scoring.aspects import _dropoff_at

    # A stop 5 nt from position 0 cannot carry a 17 nt upstream flank.
    assert _dropoff_at(5, 1, {5: 1.0}) is None


def test_dropoff_is_none_without_signal_either_side():
    from TranslonScorer.scoring.aspects import _dropoff_at

    assert _dropoff_at(1000, 1, {}) is None


def test_term_evidence_carries_dropoff():
    from TranslonScorer.scoring.aspects import score_termination_event

    cov = {1000 - (17 - k): 5.0 for k in range(0, 18, 3)}
    s = score_termination_event(1000, 1, cov)
    assert "dropoff" in s
    assert s["dropoff"] == 1.0
    # It must not have displaced the headline metric.
    assert s["metric_name"] == "term_drop"


def test_flank_fallback_and_splice_aware_agree_on_order():
    """The fallback must be the splice-aware path's order, not its reverse.

    They diverged on the minus strand: fallback ascending genomic, splice-aware
    descending (transcript order). No score moved, because _codon_levels bins
    in 3s and every flank length is a multiple of 3, so the bin set was
    identical -- but any frame-locked consumer would have read the window
    backwards.
    """
    from TranslonScorer.scoring.aspects import _project_flank

    far_intron = {("chr1", "+"): [(50_000, 50_100)], ("chr1", "-"): [(50_000, 50_100)]}
    for strand in (1, -1):
        for body_side in (True, False):
            fallback, _ = _project_flank(1000, strand, 9, body_side=body_side)
            aware, _ = _project_flank(
                1000, strand, 9, body_side=body_side, chrom="chr1", splice_context=far_intron
            )
            assert fallback == aware, f"strand={strand} body_side={body_side}"

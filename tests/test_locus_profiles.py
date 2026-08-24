"""Strand-aware A-site placement for locus-profile export.

Not previously covered by any test — `build_locus_profiles_zarr` has no test
file at all. This targets the extracted `_a_site_positions` helper directly
rather than the full zarr/parquet-backed export function, since the actual
fix is entirely in that one calculation and the rest of the export pipeline
is unchanged.
"""

from __future__ import annotations

import numpy as np

from TranslonScorer.coverage.locus_profiles import _a_site_positions


def test_plus_strand_matches_naive_start_plus_offset():
    """Plus strand is unaffected by the fix -- 5' end is genomically `start`."""
    starts = np.array([1000, 2000])
    stops = np.array([1030, 2029])
    strands = ["+", "+"]
    offsets = np.array([15, 12])
    result = _a_site_positions(starts, stops, strands, offsets)
    assert result.tolist() == [1015, 2012]


def test_minus_strand_anchors_at_stop_minus_one_not_start():
    """Minus strand: 5' end is genomically `stop - 1`, offset retreats from
    there. The pre-fix code used `start + offset` for both strands, which
    anchors a minus-strand read at its 3' end instead of its 5' end -- wrong
    by (length - 1) before the offset is even applied, easily enough to
    place the "A-site" outside the read entirely and destroy any
    frame-based signal."""
    start, stop, offset = 1000, 1030, 15  # length 30
    strands = ["-"]
    correct = _a_site_positions(np.array([start]), np.array([stop]), strands, np.array([offset]))
    naive_buggy = start + offset  # the pre-fix formula, for contrast
    assert correct.tolist() == [stop - 1 - offset]
    assert correct[0] != naive_buggy
    assert correct[0] == 1014
    assert naive_buggy == 1015


def test_mixed_strand_array_handled_elementwise():
    starts = np.array([1000, 5000, 9000])
    stops = np.array([1030, 5029, 9030])
    strands = ["+", "-", "+"]
    offsets = np.array([15, 12, 20])
    result = _a_site_positions(starts, stops, strands, offsets)
    assert result.tolist() == [
        1000 + 15,  # plus: start + offset
        5029 - 1 - 12,  # minus: stop - 1 - offset
        9000 + 20,  # plus: start + offset
    ]


def test_strand_as_numpy_str_array_not_just_python_list():
    """idx.get_column('strand').to_list() returns Python str objects in
    practice, but guard against the numpy-string-scalar case too since the
    comparison is `str(s) == '-'`, not `s == '-'`."""
    starts = np.array([1000])
    stops = np.array([1030])
    strands = np.array(["-"], dtype="<U1")
    offsets = np.array([15])
    result = _a_site_positions(starts, stops, list(strands), offsets)
    assert result.tolist() == [1030 - 1 - 15]

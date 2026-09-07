"""BedgraphSetProvider: coordinates, summing, and capability boundaries.

The provider exists so the normal scorer can eat the R reference's exact input,
so the coordinate convention is the thing most worth pinning: bedGraph is
0-based half-open, and coverage() emits 0-based positions to match
BigwigSetProvider.
"""

from __future__ import annotations

import pytest

from TranslonScorer.coverage.bedgraph import BedgraphSetProvider
from TranslonScorer.model import Region


def _write(path, rows, header: bool = False):
    lines = ["track type=bedGraph\n"] if header else []
    lines += [f"{c}\t{s}\t{e}\t{v}\n" for c, s, e, v in rows]
    path.write_text("".join(lines))
    return path


def test_bedgraph_is_zero_based_half_open(tmp_path):
    bg = _write(tmp_path / "a.bedgraph", [("chr1", 10, 13, 4.0)])
    p = BedgraphSetProvider([str(bg)])
    df = p.coverage([Region("chr1", 0, 20)])
    # [10,13) covers 0-based 10, 11, 12 -- not 13.
    assert df["pos"].to_list() == [10, 11, 12]
    assert df["count"].to_list() == [4.0, 4.0, 4.0]


def test_region_clips_the_interval(tmp_path):
    bg = _write(tmp_path / "a.bedgraph", [("chr1", 0, 100, 1.0)])
    p = BedgraphSetProvider([str(bg)])
    df = p.coverage([Region("chr1", 5, 8)])
    assert df["pos"].to_list() == [5, 6, 7]


def test_files_are_summed_position_wise(tmp_path):
    """The five atomic runs must be summed BEFORE scoring, not averaged after."""
    a = _write(tmp_path / "a.bedgraph", [("chr1", 10, 12, 3.0)])
    b = _write(tmp_path / "b.bedgraph", [("chr1", 11, 13, 5.0)])
    p = BedgraphSetProvider([str(a), str(b)])
    df = p.coverage([Region("chr1", 0, 20)]).sort("pos")
    assert df["pos"].to_list() == [10, 11, 12]
    assert df["count"].to_list() == [3.0, 8.0, 5.0]  # 11 is in both


def test_by_sample_keeps_files_apart(tmp_path):
    a = _write(tmp_path / "a.bedgraph", [("chr1", 10, 11, 3.0)])
    b = _write(tmp_path / "b.bedgraph", [("chr1", 10, 11, 5.0)])
    p = BedgraphSetProvider([str(a), str(b)], sample_names=["s1", "s2"])
    df = p.coverage([Region("chr1", 0, 20)], by_sample=True).sort("sample_id")
    assert df["sample_id"].to_list() == ["s1", "s2"]
    assert df["count"].to_list() == [3.0, 5.0]


def test_gaps_are_absent_not_zero_rows(tmp_path):
    """A bedGraph gap means zero signal; the frame stays sparse."""
    bg = _write(tmp_path / "a.bedgraph", [("chr1", 10, 11, 1.0), ("chr1", 20, 21, 1.0)])
    p = BedgraphSetProvider([str(bg)])
    df = p.coverage([Region("chr1", 0, 30)])
    assert df["pos"].to_list() == [10, 20]


def test_zero_valued_intervals_are_dropped(tmp_path):
    bg = _write(tmp_path / "a.bedgraph", [("chr1", 10, 12, 0.0), ("chr1", 12, 13, 2.0)])
    p = BedgraphSetProvider([str(bg)])
    df = p.coverage([Region("chr1", 0, 20)])
    assert df["pos"].to_list() == [12]


def test_missing_chromosome_yields_empty_frame_with_schema(tmp_path):
    bg = _write(tmp_path / "a.bedgraph", [("chr1", 10, 12, 1.0)])
    p = BedgraphSetProvider([str(bg)])
    df = p.coverage([Region("chrX", 0, 20)])
    assert df.height == 0
    assert set(df.columns) == {"pos", "count"}


def test_track_and_comment_lines_are_skipped(tmp_path):
    bg = _write(tmp_path / "a.bedgraph", [("chr1", 10, 11, 7.0)], header=True)
    p = BedgraphSetProvider([str(bg)])
    assert p.coverage([Region("chr1", 0, 20)])["count"].to_list() == [7.0]


def test_stranded_pairs_emit_a_strand_column(tmp_path):
    fwd = _write(tmp_path / "f.bedgraph", [("chr1", 10, 11, 1.0)])
    rev = _write(tmp_path / "r.bedgraph", [("chr1", 20, 21, 2.0)])
    p = BedgraphSetProvider([{"forward": str(fwd), "reverse": str(rev)}], stranded=True)
    df = p.coverage([Region("chr1", 0, 30)]).sort("pos")
    assert df["strand"].to_list() == [1, -1]
    assert df["count"].to_list() == [1.0, 2.0]


def test_stranded_requires_forward_reverse_dicts(tmp_path):
    bg = _write(tmp_path / "a.bedgraph", [("chr1", 10, 11, 1.0)])
    p = BedgraphSetProvider([str(bg)], stranded=True)
    with pytest.raises(ValueError, match="forward"):
        p.coverage([Region("chr1", 0, 20)])


def test_capabilities_it_does_not_have_fail_loudly(tmp_path):
    bg = _write(tmp_path / "a.bedgraph", [("chr1", 10, 11, 1.0)])
    p = BedgraphSetProvider([str(bg)])
    with pytest.raises(NotImplementedError, match="CIGAR"):
        p.junction_support()
    with pytest.raises(NotImplementedError, match="multimapper"):
        p.mappability_ledger(None)
    assert p.size_factors() == {"a": 1.0}


def test_empty_and_malformed_files_are_rejected(tmp_path):
    empty = tmp_path / "empty.bedgraph"
    empty.write_text("")
    with pytest.raises(ValueError, match="empty"):
        BedgraphSetProvider([str(empty)]).coverage([Region("chr1", 0, 1)])

    bad = tmp_path / "bad.bedgraph"
    bad.write_text("chr1\t10\t12\n")
    with pytest.raises(ValueError, match="malformed"):
        BedgraphSetProvider([str(bad)]).coverage([Region("chr1", 0, 1)])


def test_agrees_with_the_bigwig_provider_on_identical_data(tmp_path):
    """Same signal, two container formats, same coverage frame.

    This is what makes a bedGraph-vs-bigwig score comparison meaningful: any
    difference in downstream numbers is then the data, not the reader.
    """
    pyBigWig = pytest.importorskip("pyBigWig")

    intervals = [("chr1", 5, 8, 2.0), ("chr1", 12, 15, 7.0), ("chr1", 30, 31, 1.0)]
    bg = _write(tmp_path / "a.bedgraph", intervals)

    bw_path = tmp_path / "a.bw"
    bw = pyBigWig.open(str(bw_path), "w")
    bw.addHeader([("chr1", 100)])
    bw.addEntries(
        [c for c, _s, _e, _v in intervals],
        [s for _c, s, _e, _v in intervals],
        ends=[e for _c, _s, e, _v in intervals],
        values=[v for _c, _s, _e, v in intervals],
    )
    bw.close()

    from TranslonScorer.coverage.bigwig import BigwigSetProvider

    region = [Region("chr1", 0, 100)]
    from_bg = BedgraphSetProvider([str(bg)]).coverage(region).sort("pos")
    from_bw = BigwigSetProvider([str(bw_path)]).coverage(region).sort("pos")

    assert from_bg["pos"].to_list() == from_bw["pos"].to_list()
    assert from_bg["count"].to_list() == pytest.approx(from_bw["count"].to_list())

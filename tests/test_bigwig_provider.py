"""BigwigSetProvider stranded-coverage fix.

Previously `stranded` was a stub: forward/reverse bigwig values were always
sum-merged into one `count` per position regardless of the flag, so genuine
per-strand scoring never worked even when the provider was wired in. These
tests write small real bigwig files (pyBigWig) and assert coverage() emits a
correct `strand` column, and that unstranded/legacy behaviour is unchanged.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

pyBigWig = pytest.importorskip("pyBigWig")

from TranslonScorer.coverage.bigwig import BigwigSetProvider
from TranslonScorer.model import Region


def _write_bw(path: Path, chrom: str, chrom_len: int, intervals: list) -> None:
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(chrom, chrom_len)])
    starts = [s for s, e, v in intervals]
    ends = [e for s, e, v in intervals]
    values = [float(v) for s, e, v in intervals]
    bw.addEntries([chrom] * len(intervals), starts, ends=ends, values=values)
    bw.close()


@pytest.fixture
def stranded_bigwigs(tmp_path: Path):
    fwd = tmp_path / "fwd.bw"
    rev = tmp_path / "rev.bw"
    _write_bw(fwd, "chr1", 1000, [(100, 110, 5.0)])
    _write_bw(rev, "chr1", 1000, [(100, 110, 9.0)])
    return {"forward": str(fwd), "reverse": str(rev)}


def test_stranded_coverage_emits_strand_column(stranded_bigwigs):
    provider = BigwigSetProvider([stranded_bigwigs], stranded=True)
    cov = provider.coverage([Region("chr1", 100, 110)])

    assert "strand" in cov.columns
    fwd_rows = cov.filter(pl.col("strand") == 1)
    rev_rows = cov.filter(pl.col("strand") == -1)
    assert fwd_rows.height == 10
    assert rev_rows.height == 10
    assert set(fwd_rows["count"].to_list()) == {5.0}
    assert set(rev_rows["count"].to_list()) == {9.0}
    # forward and reverse must NOT be merged into one row per position
    assert cov.height == 20


def test_unstranded_still_merges_forward_and_reverse(stranded_bigwigs):
    """stranded=False (default): legacy merge behaviour, unchanged."""
    provider = BigwigSetProvider([stranded_bigwigs], stranded=False)
    cov = provider.coverage([Region("chr1", 100, 110)])

    assert "strand" not in cov.columns
    assert cov.height == 10  # one row per position, forward+reverse summed
    assert set(cov["count"].to_list()) == {14.0}  # 5.0 + 9.0


def test_stranded_true_requires_forward_reverse_dict(tmp_path: Path):
    bw = tmp_path / "plain.bw"
    _write_bw(bw, "chr1", 1000, [(100, 110, 3.0)])
    provider = BigwigSetProvider([str(bw)], stranded=True)
    with pytest.raises(ValueError, match="forward.*reverse"):
        provider.coverage([Region("chr1", 100, 110)])


def test_plain_paths_unstranded_unaffected(tmp_path: Path):
    bw = tmp_path / "plain.bw"
    _write_bw(bw, "chr1", 1000, [(200, 205, 4.0)])
    provider = BigwigSetProvider([str(bw)])
    cov = provider.coverage([Region("chr1", 200, 205)])
    assert "strand" not in cov.columns
    assert cov.height == 5
    assert set(cov["count"].to_list()) == {4.0}

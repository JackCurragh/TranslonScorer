"""Unit tests for the P-site index build + query round-trip.

No BAM or matrix files required — synthetic DataFrames exercise the
query path (_assign_psite_to_features, query_frame_rollup, query_coverage_index).
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from TranslonScorer.psite_index import (
    _apply_offsets_vectorised,
    _assign_psite_to_features,
    _exon_intervals_for_chrom,
    available_chroms,
    query_coverage_index,
    query_frame_rollup,
)


# ---------------------------------------------------------------------------
# Helpers to build synthetic index shards on disk
# ---------------------------------------------------------------------------


def _write_shard(out_dir: Path, chrom: str, rows: list[dict]) -> None:
    shard_dir = out_dir / f"chrom={chrom}"
    shard_dir.mkdir(parents=True, exist_ok=True)
    df = pl.DataFrame(
        {
            "pos5": pl.Series([r["pos5"] for r in rows], dtype=pl.Int32),
            "strand": pl.Series([r["strand"] for r in rows], dtype=pl.Boolean),
            "length": pl.Series([r["length"] for r in rows], dtype=pl.UInt8),
            "sample_id": pl.Series([r["sample_id"] for r in rows], dtype=pl.UInt16),
            "count": pl.Series([r["count"] for r in rows], dtype=pl.Float32),
        }
    ).sort("pos5")
    df.write_parquet(shard_dir / "data.parquet", statistics=True)

    pl.DataFrame({"sample_id": pl.Series([0, 1], dtype=pl.UInt16),
                  "sample_name": pl.Series(["s0", "s1"], dtype=pl.Utf8)}).write_parquet(
        out_dir / "samples.parquet"
    )


def _simple_cds_df(chrom: str = "chr1", strand: str = "+") -> pl.DataFrame:
    """Single-exon CDS: chrom, pos 1000–1300, tran_id='feat1'."""
    return pl.DataFrame({
        "tran_id": ["feat1"],
        "gene_id": ["feat1"],
        "chr": [chrom],
        "strand": [strand],
        "start": [[1000]],
        "stop": [[1300]],
        "tran_start": [[0]],
    })


def _two_exon_cds_df() -> pl.DataFrame:
    """Two-exon CDS: exon1 1000-1100 (tran_start=0), exon2 2000-2200 (tran_start=100)."""
    return pl.DataFrame({
        "tran_id": ["feat2"],
        "gene_id": ["feat2"],
        "chr": ["chr1"],
        "strand": ["+"],
        "start": [[1000, 2000]],
        "stop": [[1100, 2200]],
        "tran_start": [[0, 100]],
    })


# ---------------------------------------------------------------------------
# _apply_offsets_vectorised
# ---------------------------------------------------------------------------


def test_apply_offsets_plus_strand():
    pos5 = np.array([100, 200], dtype=np.int32)
    strand = np.array([True, True])
    length = np.array([28, 29], dtype=np.uint8)
    sample_id = np.array([0, 0], dtype=np.uint16)
    offsets = {(0, 28): 12, (0, 29): 13}
    result = _apply_offsets_vectorised(pos5, strand, length, sample_id, offsets, 15)
    assert result[0] == 112
    assert result[1] == 213


def test_apply_offsets_minus_strand():
    pos5 = np.array([200], dtype=np.int32)
    strand = np.array([False])
    length = np.array([28], dtype=np.uint8)
    sample_id = np.array([0], dtype=np.uint16)
    offsets = {(0, 28): 12}
    result = _apply_offsets_vectorised(pos5, strand, length, sample_id, offsets, 15)
    assert result[0] == 188  # 200 - 12


def test_apply_offsets_default_fallback():
    pos5 = np.array([100], dtype=np.int32)
    strand = np.array([True])
    length = np.array([30], dtype=np.uint8)
    sample_id = np.array([0], dtype=np.uint16)
    result = _apply_offsets_vectorised(pos5, strand, length, sample_id, {}, 15)
    assert result[0] == 115  # 100 + 15


# ---------------------------------------------------------------------------
# _exon_intervals_for_chrom
# ---------------------------------------------------------------------------


def test_exon_intervals_single_exon():
    cds = _simple_cds_df()
    ivs = _exon_intervals_for_chrom(cds, "chr1")
    assert len(ivs) == 1
    fid, strand, g_min, g_max, n_exons, exons = ivs[0]
    assert fid == "feat1"
    assert strand == "+"
    assert g_min == 1000
    assert g_max == 1300
    assert n_exons == 1
    assert exons[0] == (1000, 1300, 0)


def test_exon_intervals_two_exons():
    cds = _two_exon_cds_df()
    ivs = _exon_intervals_for_chrom(cds, "chr1")
    assert len(ivs) == 1
    _, _, g_min, g_max, n_exons, exons = ivs[0]
    assert g_min == 1000
    assert g_max == 2200
    assert n_exons == 2
    assert exons[1] == (2000, 2200, 100)


def test_exon_intervals_wrong_chrom():
    cds = _simple_cds_df(chrom="chr1")
    ivs = _exon_intervals_for_chrom(cds, "chr22")
    assert ivs == []


# ---------------------------------------------------------------------------
# _assign_psite_to_features
# ---------------------------------------------------------------------------


def test_assign_in_frame_read():
    cds = _simple_cds_df()
    features = _exon_intervals_for_chrom(cds, "chr1")
    # read at pos5=985, offset=15 → p_site=1000 → tx_pos=0 → frame 0
    result = _assign_psite_to_features(
        np.array([985], dtype=np.int32),
        np.array([True]),
        features,
        {(0, 28): 15},
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.height == 1
    assert result["feature_id"][0] == "feat1"
    assert result["tx_pos"][0] == 0


def test_assign_out_of_span():
    cds = _simple_cds_df()
    features = _exon_intervals_for_chrom(cds, "chr1")
    # p_site = 500 → not in [1000, 1300)
    result = _assign_psite_to_features(
        np.array([485], dtype=np.int32),
        np.array([True]),
        features,
        {(0, 28): 15},
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.is_empty()


def test_assign_intron_read():
    cds = _two_exon_cds_df()
    features = _exon_intervals_for_chrom(cds, "chr1")
    # p_site = 1500 → in genomic span but not in any exon (1100-2000 is intron)
    result = _assign_psite_to_features(
        np.array([1485], dtype=np.int32),
        np.array([True]),
        features,
        {(0, 28): 15},
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.is_empty()


def test_assign_second_exon():
    cds = _two_exon_cds_df()
    features = _exon_intervals_for_chrom(cds, "chr1")
    # p_site = 2000 → exon2, tx_pos = 100
    result = _assign_psite_to_features(
        np.array([1985], dtype=np.int32),
        np.array([True]),
        features,
        {(0, 28): 15},
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.height == 1
    assert result["tx_pos"][0] == 100


# ---------------------------------------------------------------------------
# query_frame_rollup end-to-end with synthetic index on disk
# ---------------------------------------------------------------------------


def test_query_frame_rollup_basic(tmp_path):
    # 3 in-frame reads (p_site at tx_pos 0, 3, 6 → frame 0) + 1 out-of-frame
    # feature: pos 1000–1300 on chr1 +
    # reads: pos5=985 (offset 15 → p_site=1000, tx_pos=0, frame=0)
    #        pos5=988 (→ p_site=1003, tx_pos=3, frame=0)
    #        pos5=991 (→ p_site=1006, tx_pos=6, frame=0)
    #        pos5=986 (→ p_site=1001, tx_pos=1, frame=1)
    rows = [
        {"pos5": 985, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"pos5": 988, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"pos5": 991, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"pos5": 986, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
    ]
    _write_shard(tmp_path, "chr1", rows)

    cds = _simple_cds_df()
    offsets = {(0, 28): 15}
    result = query_frame_rollup(tmp_path, cds, offsets, chroms=["chr1"])

    assert result.height == 1
    row = result.row(0, named=True)
    assert row["feature_id"] == "feat1"
    assert row["n_reads"] == pytest.approx(4.0)
    assert row["frame0"] == pytest.approx(3.0)
    assert row["frame1"] == pytest.approx(1.0)
    assert row["frame2"] == pytest.approx(0.0)


def test_query_frame_rollup_strand_filter(tmp_path):
    # minus-strand read should NOT match a plus-strand feature
    rows = [
        {"pos5": 985, "strand": False, "length": 28, "sample_id": 0, "count": 1.0},
    ]
    _write_shard(tmp_path, "chr1", rows)
    cds = _simple_cds_df(strand="+")
    result = query_frame_rollup(tmp_path, cds, {(0, 28): 15}, chroms=["chr1"])
    assert result.is_empty()


def test_query_frame_rollup_missing_chrom(tmp_path):
    rows = [{"pos5": 985, "strand": True, "length": 28, "sample_id": 0, "count": 1.0}]
    _write_shard(tmp_path, "chr22", rows)
    cds = _simple_cds_df(chrom="chr1")
    result = query_frame_rollup(tmp_path, cds, {(0, 28): 15}, chroms=["chr22"])
    assert result.is_empty()


def test_query_frame_rollup_multiple_samples(tmp_path):
    rows = [
        {"pos5": 985, "strand": True, "length": 28, "sample_id": 0, "count": 2.0},
        {"pos5": 985, "strand": True, "length": 28, "sample_id": 1, "count": 3.0},
    ]
    _write_shard(tmp_path, "chr1", rows)
    cds = _simple_cds_df()
    offsets = {(0, 28): 15, (1, 28): 15}
    result = query_frame_rollup(tmp_path, cds, offsets, chroms=["chr1"])
    assert result.height == 2
    totals = dict(zip(result["sample_id"].to_list(), result["n_reads"].to_list()))
    assert totals[0] == pytest.approx(2.0)
    assert totals[1] == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# query_coverage_index
# ---------------------------------------------------------------------------


def test_query_coverage_index_basic(tmp_path):
    rows = [
        {"pos5": 985, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"pos5": 988, "strand": True, "length": 28, "sample_id": 0, "count": 2.0},
    ]
    _write_shard(tmp_path, "chr1", rows)
    cds = _simple_cds_df()
    result = query_coverage_index(tmp_path, cds, {(0, 28): 15}, chroms=["chr1"])
    assert result.height == 2
    pos_to_count = dict(zip(result["tx_pos"].to_list(), result["count"].to_list()))
    assert pos_to_count[0] == pytest.approx(1.0)
    assert pos_to_count[3] == pytest.approx(2.0)


def test_query_coverage_index_feature_filter(tmp_path):
    # Two features; filter to feat1 only
    cds = pl.concat([
        _simple_cds_df(),
        pl.DataFrame({
            "tran_id": ["feat2"], "gene_id": ["feat2"], "chr": ["chr1"], "strand": ["+"],
            "start": [[2000]], "stop": [[2300]], "tran_start": [[0]],
        }),
    ])
    rows = [
        {"pos5": 985, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"pos5": 1985, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
    ]
    _write_shard(tmp_path, "chr1", rows)
    result = query_coverage_index(
        tmp_path, cds, {(0, 28): 15}, feature_ids=["feat1"], chroms=["chr1"]
    )
    assert all(r == "feat1" for r in result["feature_id"].to_list())


# ---------------------------------------------------------------------------
# available_chroms
# ---------------------------------------------------------------------------


def test_available_chroms(tmp_path):
    for c in ["chr1", "chr22", "chrX"]:
        (tmp_path / f"chrom={c}").mkdir()
    (tmp_path / "samples.parquet").touch()
    chroms = available_chroms(tmp_path)
    assert set(chroms) == {"chr1", "chr22", "chrX"}

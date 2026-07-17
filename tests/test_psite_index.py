"""Unit tests for the P-site index build + query round-trip.

No BAM or matrix files required.  Synthetic shards store pre-offset-corrected
p_site values (as the real index does after Phase 2 bakes offsets in).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from TranslonScorer.matrix.psite_index import (
    _assign_psite_to_features,
    _exon_intervals_for_chrom,
    available_chroms,
    query_coverage_index,
    query_frame_rollup,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_shard(out_dir: Path, chrom: str, rows: list[dict]) -> None:
    """Write a synthetic p_site shard (offsets already baked in)."""
    shard_dir = out_dir / f"chrom={chrom}"
    shard_dir.mkdir(parents=True, exist_ok=True)
    df = pl.DataFrame(
        {
            "p_site": pl.Series([r["p_site"] for r in rows], dtype=pl.Int32),
            "strand": pl.Series([r["strand"] for r in rows], dtype=pl.Boolean),
            "length": pl.Series([r["length"] for r in rows], dtype=pl.UInt8),
            "sample_id": pl.Series([r["sample_id"] for r in rows], dtype=pl.UInt16),
            "count": pl.Series([r["count"] for r in rows], dtype=pl.Float32),
        }
    ).sort("p_site")
    df.write_parquet(shard_dir / "data.parquet", statistics=True)

    pl.DataFrame(
        {
            "sample_id": pl.Series([0, 1], dtype=pl.UInt16),
            "sample_name": pl.Series(["s0", "s1"], dtype=pl.Utf8),
        }
    ).write_parquet(out_dir / "samples.parquet")


def _simple_cds_df(chrom: str = "chr1", strand: str = "+") -> pl.DataFrame:
    """Single-exon CDS: chr1 pos 1000–1300."""
    return pl.DataFrame(
        {
            "tran_id": ["feat1"],
            "gene_id": ["feat1"],
            "chr": [chrom],
            "strand": [strand],
            "start": [[1000]],
            "stop": [[1300]],
            "tran_start": [[0]],
        }
    )


def _two_exon_cds_df() -> pl.DataFrame:
    """Two-exon: exon1 1000–1100 (tran_start=0), exon2 2000–2200 (tran_start=100)."""
    return pl.DataFrame(
        {
            "tran_id": ["feat2"],
            "gene_id": ["feat2"],
            "chr": ["chr1"],
            "strand": ["+"],
            "start": [[1000, 2000]],
            "stop": [[1100, 2200]],
            "tran_start": [[0, 100]],
        }
    )


# ---------------------------------------------------------------------------
# _exon_intervals_for_chrom
# ---------------------------------------------------------------------------


def test_exon_intervals_single_exon():
    ivs = _exon_intervals_for_chrom(_simple_cds_df(), "chr1")
    assert len(ivs) == 1
    fid, strand, g_min, g_max, n_exons, exons = ivs[0]
    assert fid == "feat1"
    assert strand == "+"
    assert g_min == 1000
    assert g_max == 1300
    assert n_exons == 1
    assert exons[0] == (1000, 1300, 0)


def test_exon_intervals_two_exons():
    ivs = _exon_intervals_for_chrom(_two_exon_cds_df(), "chr1")
    _, _, g_min, g_max, n_exons, exons = ivs[0]
    assert g_min == 1000
    assert g_max == 2200
    assert n_exons == 2
    assert exons[1] == (2000, 2200, 100)


def test_exon_intervals_wrong_chrom():
    assert _exon_intervals_for_chrom(_simple_cds_df(chrom="chr1"), "chr22") == []


# ---------------------------------------------------------------------------
# _assign_psite_to_features (p_site already offset-corrected)
# ---------------------------------------------------------------------------


def test_assign_in_frame_read():
    features = _exon_intervals_for_chrom(_simple_cds_df(), "chr1")
    # p_site=1000 → tx_pos=0, frame=0
    result = _assign_psite_to_features(
        np.array([1000], dtype=np.int32),
        np.array([True]),
        features,
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.height == 1
    assert result["feature_id"][0] == "feat1"
    assert result["tx_pos"][0] == 0


def test_assign_out_of_span():
    features = _exon_intervals_for_chrom(_simple_cds_df(), "chr1")
    result = _assign_psite_to_features(
        np.array([500], dtype=np.int32),
        np.array([True]),
        features,
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.is_empty()


def test_assign_intron_read():
    features = _exon_intervals_for_chrom(_two_exon_cds_df(), "chr1")
    # p_site=1500 → in genomic span but between exons
    result = _assign_psite_to_features(
        np.array([1500], dtype=np.int32),
        np.array([True]),
        features,
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.is_empty()


def test_assign_second_exon():
    features = _exon_intervals_for_chrom(_two_exon_cds_df(), "chr1")
    # p_site=2000 → exon2, tx_pos=100
    result = _assign_psite_to_features(
        np.array([2000], dtype=np.int32),
        np.array([True]),
        features,
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.height == 1
    assert result["tx_pos"][0] == 100


def test_assign_strand_filter():
    features = _exon_intervals_for_chrom(_simple_cds_df(strand="+"), "chr1")
    # minus-strand read should not match plus-strand feature
    result = _assign_psite_to_features(
        np.array([1000], dtype=np.int32),
        np.array([False]),
        features,
        np.array([28], dtype=np.uint8),
        np.array([0], dtype=np.uint16),
        np.array([1.0], dtype=np.float32),
    )
    assert result.is_empty()


# ---------------------------------------------------------------------------
# query_frame_rollup (end-to-end with synthetic disk shard)
# ---------------------------------------------------------------------------


def test_query_frame_rollup_basic(tmp_path):
    # 3 in-frame reads (tx_pos 0, 3, 6 → frame 0) + 1 out-of-frame (tx_pos 1 → frame 1)
    # p_site values are already offset-corrected
    rows = [
        {"p_site": 1000, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"p_site": 1003, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"p_site": 1006, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"p_site": 1001, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
    ]
    _write_shard(tmp_path, "chr1", rows)

    result = query_frame_rollup(tmp_path, _simple_cds_df(), chroms=["chr1"])
    assert result.height == 1
    row = result.row(0, named=True)
    assert row["feature_id"] == "feat1"
    assert row["n_reads"] == pytest.approx(4.0)
    assert row["frame0"] == pytest.approx(3.0)
    assert row["frame1"] == pytest.approx(1.0)
    assert row["frame2"] == pytest.approx(0.0)


def test_query_frame_rollup_strand_filter(tmp_path):
    rows = [{"p_site": 1000, "strand": False, "length": 28, "sample_id": 0, "count": 1.0}]
    _write_shard(tmp_path, "chr1", rows)
    result = query_frame_rollup(tmp_path, _simple_cds_df(strand="+"), chroms=["chr1"])
    assert result.is_empty()


def test_query_frame_rollup_missing_chrom(tmp_path):
    rows = [{"p_site": 1000, "strand": True, "length": 28, "sample_id": 0, "count": 1.0}]
    _write_shard(tmp_path, "chr22", rows)
    result = query_frame_rollup(tmp_path, _simple_cds_df(chrom="chr1"), chroms=["chr22"])
    assert result.is_empty()


def test_query_frame_rollup_multiple_samples(tmp_path):
    rows = [
        {"p_site": 1000, "strand": True, "length": 28, "sample_id": 0, "count": 2.0},
        {"p_site": 1000, "strand": True, "length": 28, "sample_id": 1, "count": 3.0},
    ]
    _write_shard(tmp_path, "chr1", rows)
    result = query_frame_rollup(tmp_path, _simple_cds_df(), chroms=["chr1"])
    assert result.height == 2
    totals = dict(zip(result["sample_id"].to_list(), result["n_reads"].to_list()))
    assert totals[0] == pytest.approx(2.0)
    assert totals[1] == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# query_coverage_index
# ---------------------------------------------------------------------------


def test_query_coverage_index_basic(tmp_path):
    rows = [
        {"p_site": 1000, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"p_site": 1003, "strand": True, "length": 28, "sample_id": 0, "count": 2.0},
    ]
    _write_shard(tmp_path, "chr1", rows)
    result = query_coverage_index(tmp_path, _simple_cds_df(), chroms=["chr1"])
    assert result.height == 2
    pos_to_count = dict(zip(result["tx_pos"].to_list(), result["count"].to_list()))
    assert pos_to_count[0] == pytest.approx(1.0)
    assert pos_to_count[3] == pytest.approx(2.0)


def test_query_coverage_index_feature_filter(tmp_path):
    cds = pl.concat(
        [
            _simple_cds_df(),
            pl.DataFrame(
                {
                    "tran_id": ["feat2"],
                    "gene_id": ["feat2"],
                    "chr": ["chr1"],
                    "strand": ["+"],
                    "start": [[2000]],
                    "stop": [[2300]],
                    "tran_start": [[0]],
                }
            ),
        ]
    )
    rows = [
        {"p_site": 1000, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
        {"p_site": 2000, "strand": True, "length": 28, "sample_id": 0, "count": 1.0},
    ]
    _write_shard(tmp_path, "chr1", rows)
    result = query_coverage_index(tmp_path, cds, feature_ids=["feat1"], chroms=["chr1"])
    assert all(r == "feat1" for r in result["feature_id"].to_list())


# ---------------------------------------------------------------------------
# available_chroms
# ---------------------------------------------------------------------------


def test_available_chroms(tmp_path):
    for c in ["chr1", "chr22", "chrX"]:
        (tmp_path / f"chrom={c}").mkdir()
    (tmp_path / "samples.parquet").touch()
    assert set(available_chroms(tmp_path)) == {"chr1", "chr22", "chrX"}

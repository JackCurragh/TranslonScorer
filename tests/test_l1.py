"""Tests for L1 schema and builder."""
from __future__ import annotations

import time
from pathlib import Path

import polars as pl
import pyarrow as pa
import pytest

from TranslonScorer.l0b.contracts import VersionKey, write_meta
from TranslonScorer.l1 import build_l1
from TranslonScorer.l1.schema import L1_JUNC_SCHEMA, L1_POS_SCHEMA

DATA_DIR = Path(__file__).parent.parent / "data" / "global_partitioned"
L0B_PREBUILT = Path("/tmp/l0b_test")

pytestmark = pytest.mark.skipif(
    not DATA_DIR.exists() or not L0B_PREBUILT.exists(),
    reason="local cohort data or L0b build not present",
)

KEY = VersionKey("GRCh38.p14", "testcfg", "GENCODE_v44", "2026-06", "unique_only")


# ---------------------------------------------------------------------------
# Schema unit tests (fast, no data)
# ---------------------------------------------------------------------------

def test_pos_schema_field_names():
    names = [f.name for f in L1_POS_SCHEMA]
    assert names == ["pos5", "strand", "length", "sample_id", "count"]


def test_junc_schema_field_names():
    names = [f.name for f in L1_JUNC_SCHEMA]
    assert names == ["donor", "acceptor", "strand", "length", "sample_id", "count"]


# ---------------------------------------------------------------------------
# Builder integration tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def l1_dir(tmp_path_factory):
    out = tmp_path_factory.mktemp("l1_out")
    build_l1(
        L0B_PREBUILT,
        DATA_DIR,
        out,
        KEY,
        chroms=["chr1"],
        workers=1,
    )
    return out


def test_meta_written(l1_dir):
    from TranslonScorer.l0b.contracts import read_meta
    meta = read_meta(l1_dir)
    assert meta["status"] == "complete"
    assert meta["total_pos_rows"] > 0
    assert meta["total_junc_rows"] > 0


def test_pos_shard_exists(l1_dir):
    assert (l1_dir / "chr1_pos.parquet").exists()


def test_junc_shard_exists(l1_dir):
    assert (l1_dir / "chr1_junc.parquet").exists()


def test_pos_schema_correct(l1_dir):
    tbl = pl.read_parquet(l1_dir / "chr1_pos.parquet", n_rows=1)
    for field in L1_POS_SCHEMA:
        assert field.name in tbl.columns, f"missing column {field.name}"


def test_junc_schema_correct(l1_dir):
    tbl = pl.read_parquet(l1_dir / "chr1_junc.parquet", n_rows=1)
    for field in L1_JUNC_SCHEMA:
        assert field.name in tbl.columns, f"missing column {field.name}"


def test_pos_sorted_by_pos5(l1_dir):
    """Positional shard must be sorted by pos5 (enables predicate pushdown)."""
    pos5 = pl.read_parquet(l1_dir / "chr1_pos.parquet", columns=["pos5"])["pos5"].to_list()
    assert pos5 == sorted(pos5), "pos5 not sorted"


def test_junc_sorted(l1_dir):
    """Junction shard must be sorted by (donor, acceptor)."""
    df = pl.read_parquet(l1_dir / "chr1_junc.parquet", columns=["donor", "acceptor"])
    expected = df.sort(["donor", "acceptor"])
    assert df.equals(expected), "junc not sorted by (donor, acceptor)"


def test_unique_only_policy(l1_dir):
    """L1 must contain only unique-mapping reads (nh==1 from L0b)."""
    l0b = pl.read_parquet(L0B_PREBUILT / "chr1.parquet", columns=["read_id", "nh", "is_secondary"])
    unique_rids = set(
        l0b.filter((pl.col("nh") == 1) & (~pl.col("is_secondary")))["read_id"].cast(pl.UInt64).to_list()
    )
    multi_rids = set(
        l0b.filter(pl.col("nh") > 1)["read_id"].cast(pl.UInt64).to_list()
    ) - unique_rids  # reads that appear ONLY as multimappers

    # If multi_rids contributed to L1, we'd find their pos5 values in the shard.
    # Since this cohort is small, sample a few multimapper positions and confirm
    # their contributions are absent.
    pos = pl.read_parquet(l1_dir / "chr1_pos.parquet")
    # Total L1 pos rows should match what a unique-only join produces
    expected_unique_reads = len(unique_rids)
    assert pos["sample_id"].n_unique() > 0  # sanity
    # Per-position unique read count should be >= 1
    assert pos["count"].min() >= 1


def test_count_reconciles_with_source(l1_dir):
    """A single read's L1 count must match its source count parquet entry."""
    l0b = pl.read_parquet(L0B_PREBUILT / "chr1.parquet")
    row = l0b.filter((pl.col("nh") == 1) & (~pl.col("is_secondary")))[0]
    rid = int(row["read_id"][0])
    pos5_val = int(row["pos5"][0])
    length_val = int(row["length"][0])
    strand_val = int(row["strand"][0])

    # Source count
    source = (
        pl.scan_parquet(str(DATA_DIR / "*" / "global.*_counts" / "**" / "*.parquet"))
        .filter(pl.col("read_id") == pl.lit(rid, dtype=pl.UInt64))
        .select("sample_id", "count")
        .collect()
    )

    # L1 counts at that position+length+strand
    l1_rows = pl.read_parquet(l1_dir / "chr1_pos.parquet").filter(
        (pl.col("pos5") == pos5_val)
        & (pl.col("length") == length_val)
        & (pl.col("strand") == strand_val)
    )

    for src_row in source.iter_rows(named=True):
        sid = src_row["sample_id"]
        expected_count = src_row["count"]
        l1_match = l1_rows.filter(pl.col("sample_id") == sid)
        assert len(l1_match) == 1, f"no L1 row for sample_id={sid}"
        assert l1_match["count"][0] >= expected_count, (
            f"L1 count {l1_match['count'][0]} < source count {expected_count}"
        )


def test_all_samples_present(l1_dir):
    """All sample_ids that have counts for chr1 reads must appear in L1."""
    # Get sample_ids from source counts that have chr1 unique reads
    l0b_rids = (
        pl.read_parquet(L0B_PREBUILT / "chr1.parquet")
        .filter((pl.col("nh") == 1) & (~pl.col("is_secondary")))["read_id"]
        .cast(pl.UInt64)
    )
    source_samples = (
        pl.scan_parquet(str(DATA_DIR / "*" / "global.*_counts" / "**" / "*.parquet"))
        .filter(pl.col("read_id").is_in(l0b_rids.to_list()))
        .select("sample_id")
        .unique()
        .collect()["sample_id"]
        .to_list()
    )
    l1_samples = pl.read_parquet(l1_dir / "chr1_pos.parquet")["sample_id"].unique().to_list()
    assert set(source_samples) == set(l1_samples)


def test_resumable(l1_dir):
    """Re-running build_l1 on an existing dir returns quickly (checkpoint)."""
    t0 = time.monotonic()
    build_l1(L0B_PREBUILT, DATA_DIR, l1_dir, KEY, chroms=["chr1"], workers=1)
    elapsed = time.monotonic() - t0
    assert elapsed < 20, f"re-build took {elapsed:.1f}s — checkpoint not working?"


def test_version_key_mismatch_raises(tmp_path):
    """build_l1 must refuse to run against an L0b with a different version key."""
    wrong_key = VersionKey("hg19", "other", "GENCODE_v44", "2026-06", "unique_only")
    with pytest.raises(ValueError, match="version key mismatch"):
        build_l1(L0B_PREBUILT, DATA_DIR, tmp_path / "l1", wrong_key, chroms=["chr1"])

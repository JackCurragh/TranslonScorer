"""Tests for the alignment schema, provenance and builder.

Correctness checks against the local global_partitioned cohort.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq
import pytest

from TranslonScorer.alignments.builder import build_alignments
from TranslonScorer.alignments.provenance import VersionKey, read_meta, write_meta
from TranslonScorer.alignments.schema import ALIGNMENT_SCHEMA

DATA_DIR = Path(__file__).parent.parent / "data" / "global_partitioned"
pytestmark = pytest.mark.skipif(not DATA_DIR.exists(), reason="local cohort data not present")


# ---------------------------------------------------------------------------
# Contract tests (fast, no I/O beyond schema inspection)
# ---------------------------------------------------------------------------


def test_schema_field_names():
    names = [f.name for f in ALIGNMENT_SCHEMA]
    required = [
        "read_id",
        "chrom",
        "pos5",
        "end",
        "strand",
        "length",
        "cigar",
        "mapq",
        "nh",
        "is_secondary",
        "aln_score",
        "mismatches",
        "junctions_crossed",
        "weight",
    ]
    assert names == required


def test_version_key_fingerprint_stable():
    k = VersionKey("GRCh38", "aabbccdd1234", "GENCODE_v44", "2026-06", "unique_only")
    fp1 = k.fingerprint()
    fp2 = k.fingerprint()
    assert fp1 == fp2
    assert len(fp1) == 12


def test_version_key_different_field_gives_different_fingerprint():
    k1 = VersionKey("GRCh38", "aabbccdd1234", "GENCODE_v44", "2026-06", "unique_only")
    k2 = VersionKey("GRCh38", "aabbccdd1234", "GENCODE_v45", "2026-06", "unique_only")
    assert k1.fingerprint() != k2.fingerprint()


def test_write_read_meta(tmp_path):
    k = VersionKey("GRCh38", "abc", "v44", "2026-06", "unique_only")
    write_meta(tmp_path, k, {"foo": "bar"})
    meta = read_meta(tmp_path)
    assert meta["version_key"]["genome_id"] == "GRCh38"
    assert meta["foo"] == "bar"
    assert "fingerprint" in meta


def test_read_meta_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_meta(tmp_path)


# ---------------------------------------------------------------------------
# Builder tests (need local cohort data + samtools)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def alignment_dir(tmp_path_factory):
    """Build alignments for chr1 only (fast)."""
    out = tmp_path_factory.mktemp("alignments_out")
    key = VersionKey("GRCh38.p14", "testcfg", "GENCODE_v44", "2026-06", "unique_only")
    build_alignments(
        DATA_DIR,
        out,
        key,
        chroms=["chr1"],
        workers=1,
        samtools_threads=1,
    )
    return out


def test_meta_written(alignment_dir):
    meta = read_meta(alignment_dir)
    assert meta["status"] == "complete"
    assert meta["version_key"]["genome_id"] == "GRCh38.p14"


def test_shard_exists(alignment_dir):
    shard = alignment_dir / "chr1.parquet"
    assert shard.exists(), "chr1.parquet not found"


def test_shard_schema(alignment_dir):
    shard = alignment_dir / "chr1.parquet"
    tbl = pq.read_table(str(shard))
    assert tbl.schema.equals(ALIGNMENT_SCHEMA, check_metadata=False)


def test_shard_pos5_sorted_within_row_groups(alignment_dir):
    """Each Parquet row group must have pos5 sorted internally (enables predicate pushdown).

    Global sort is NOT required: reverse-strand reads use pos5=ref_end-1 so
    they can be out-of-order relative to adjacent batches from a ref_start-sorted
    BAM.  Within each batch the sort is applied before writing.
    """
    import pyarrow.parquet as pq

    shard = alignment_dir / "chr1.parquet"
    pf = pq.ParquetFile(str(shard))
    for rg_idx in range(pf.metadata.num_row_groups):
        tbl = pf.read_row_group(rg_idx, columns=["pos5"])
        pos = tbl.column("pos5").to_pylist()
        assert pos == sorted(pos), f"row group {rg_idx} pos5 not sorted"


def test_shard_pos5_rowgroup_statistics(alignment_dir):
    """Row-group min/max statistics must be written for pos5 (Parquet predicate pushdown)."""
    shard = alignment_dir / "chr1.parquet"
    pf = pq.ParquetFile(str(shard))
    md = pf.metadata
    pos5_col_idx = pq.read_schema(str(shard)).get_field_index("pos5")
    for rg_idx in range(md.num_row_groups):
        rg = md.row_group(rg_idx)
        col_md = rg.column(pos5_col_idx)
        stats = col_md.statistics
        assert stats is not None, f"row group {rg_idx} has no statistics"
        assert stats.has_min_max, f"row group {rg_idx} missing min/max for pos5"


def test_shard_row_count_vs_samtools(alignment_dir):
    """Row count in Parquet should equal number of chr1 alignments across all BAMs."""
    import subprocess

    import polars as pl

    bam_paths = sorted(
        p
        for subdir in DATA_DIR.iterdir()
        if subdir.is_dir()
        for p in subdir.glob("unique_reads.*.bam")
    )
    total_from_samtools = 0
    for bam in bam_paths:
        r = subprocess.run(
            ["samtools", "view", "-c", str(bam), "chr1"],
            capture_output=True,
            text=True,
            check=True,
        )
        total_from_samtools += int(r.stdout.strip())

    shard = alignment_dir / "chr1.parquet"
    n_parquet = len(pl.read_parquet(shard))
    assert (
        n_parquet == total_from_samtools
    ), f"row count mismatch: parquet={n_parquet}, samtools={total_from_samtools}"


def test_read_id_uniqueness_not_enforced_at_alignment_level(alignment_dir):
    """read_id is NOT unique in the alignment table — multimappers appear once per alignment."""
    import polars as pl

    shard = alignment_dir / "chr1.parquet"
    df = pl.read_parquet(shard)
    n_rows = len(df)
    n_unique_ids = df["read_id"].n_unique()
    # For a multimapper-retaining table, n_rows >= n_unique_ids
    assert n_rows >= n_unique_ids


def test_nh_tag_used_not_count(alignment_dir):
    """NH values > 1 must appear (not just 1s from bad fallback)."""
    import polars as pl

    shard = alignment_dir / "chr1.parquet"
    df = pl.read_parquet(shard)
    assert df.filter(pl.col("nh") > 1).height > 0, "expected multimapper rows with NH>1"


def test_chrom_column_correct(alignment_dir):
    import polars as pl

    shard = alignment_dir / "chr1.parquet"
    df = pl.read_parquet(shard)
    assert df["chrom"].unique().to_list() == ["chr1"]


def test_resumable_no_rebuild(alignment_dir):
    """Re-running build_alignments with the same output dir returns immediately (checkpoint)."""
    import time

    key = VersionKey("GRCh38.p14", "testcfg", "GENCODE_v44", "2026-06", "unique_only")
    t0 = time.monotonic()
    build_alignments(DATA_DIR, alignment_dir, key, chroms=["chr1"], workers=1, samtools_threads=1)
    elapsed = time.monotonic() - t0
    # A checkpoint hit should be much faster than a full build
    assert elapsed < 30, f"Re-build took {elapsed:.1f}s — checkpoint not working?"

"""psite_index.query_genomic_coverage + MatrixProvider(psite_index_dir=...).

Closes the "two matrix pipelines" gap: build-psite-index's Phase 2 already
writes genomic, per-(sample, length)-offset-corrected P-site positions
(chrom-sharded Parquet, `_SHARD_SCHEMA`), but that data was previously only
ever re-projected into transcript space to feed FrameRollup's elongation-only
scorer. query_genomic_coverage reads it as ordinary genomic coverage instead,
so MatrixProvider can feed it straight into the same
_score_events_over_provider pipeline used for BAM/bigwig sources — accurate
offsets AND init/term/junction/mappability, no separate pipeline.

No real BAMs/pysam needed: the index is just Parquet, built directly with
polars, same fixture style as tests/test_bigwig_provider.py.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from TranslonScorer.matrix.provider import MatrixProvider
from TranslonScorer.model import Region, ScoreThresholds
from TranslonScorer.matrix.psite_index import query_genomic_coverage
from TranslonScorer.workflows import _score_events_over_provider

_THR = ScoreThresholds()


def _write_index(
    index_dir: Path,
    *,
    chrom: str = "chr1",
    rows: list | None = None,
    samples: dict | None = None,
    usable: list | None = None,
) -> None:
    rows = rows or []
    samples = samples or {1: "S1"}
    (index_dir / f"chrom={chrom}").mkdir(parents=True)
    pl.DataFrame(
        {
            "sample_id": pl.Series([r[0] for r in rows], dtype=pl.UInt16),
            "p_site": pl.Series([r[1] for r in rows], dtype=pl.Int32),
            "strand": pl.Series([r[2] for r in rows], dtype=pl.Boolean),
            "length": pl.Series([r[3] for r in rows], dtype=pl.UInt8),
            "count": pl.Series([r[4] for r in rows], dtype=pl.Float32),
        }
    ).write_parquet(index_dir / f"chrom={chrom}" / "data.parquet")
    pl.DataFrame(
        {
            "sample_id": pl.Series(list(samples.keys()), dtype=pl.UInt16),
            "sample_name": list(samples.values()),
        }
    ).write_parquet(index_dir / "samples.parquet")
    if usable is not None:
        pl.DataFrame(
            {
                "sample_id": [u[0] for u in usable],  # sample NAME, despite the column name
                "length": [u[1] for u in usable],
            }
        ).write_parquet(index_dir / "usable_sample_lengths.parquet")


# ---------------------------------------------------------------------------
# query_genomic_coverage
# ---------------------------------------------------------------------------


def test_query_genomic_coverage_basic(tmp_path: Path):
    idx = tmp_path / "idx"
    _write_index(
        idx,
        rows=[
            (1, 100, True, 29, 5.0),
            (1, 150, True, 30, 3.0),
            (2, 100, False, 29, 2.0),
            (1, 900, True, 29, 9.0),  # outside the query window
        ],
        samples={1: "S1", 2: "S2"},
    )
    out = query_genomic_coverage(str(idx), [("chr1", 90, 200)], site="P")
    assert out["pos"].dtype == pl.Int64
    assert set(zip(out["pos"].to_list(), out["strand"].to_list(), out["count"].to_list())) == {
        (100, 1, 5.0),
        (100, -1, 2.0),
        (150, 1, 3.0),
    }


def test_query_genomic_coverage_a_site_shift(tmp_path: Path):
    idx = tmp_path / "idx"
    _write_index(idx, rows=[(1, 100, True, 29, 5.0)])
    out = query_genomic_coverage(str(idx), [("chr1", 90, 200)], site="A")
    assert out["pos"].to_list() == [103]


def test_query_genomic_coverage_sample_filter(tmp_path: Path):
    idx = tmp_path / "idx"
    _write_index(
        idx,
        rows=[(1, 100, True, 29, 5.0), (2, 100, True, 29, 2.0)],
        samples={1: "S1", 2: "S2"},
    )
    out = query_genomic_coverage(str(idx), [("chr1", 90, 200)], site="P", sample_names=["S1"])
    assert out["count"].to_list() == [5.0]


def test_query_genomic_coverage_by_sample(tmp_path: Path):
    idx = tmp_path / "idx"
    _write_index(
        idx,
        rows=[(1, 100, True, 29, 5.0), (2, 100, True, 29, 2.0)],
        samples={1: "S1", 2: "S2"},
    )
    out = query_genomic_coverage(str(idx), [("chr1", 90, 200)], site="P", group_level="sample")
    assert set(zip(out["sample_name"].to_list(), out["count"].to_list())) == {
        ("S1", 5.0),
        ("S2", 2.0),
    }


def test_query_genomic_coverage_usable_sample_lengths_gate(tmp_path: Path):
    idx = tmp_path / "idx"
    _write_index(
        idx,
        rows=[
            (1, 100, True, 29, 5.0),  # usable
            (1, 100, True, 30, 3.0),  # NOT usable (length 30 excluded for S1)
            (2, 100, True, 29, 2.0),  # usable
        ],
        samples={1: "S1", 2: "S2"},
        usable=[("S1", 29), ("S2", 29)],
    )
    out = query_genomic_coverage(str(idx), [("chr1", 90, 200)], site="P")
    assert out["count"].sum() == 7.0  # 5.0 + 2.0, the length-30 row excluded


def test_query_genomic_coverage_empty_region_returns_empty_schema(tmp_path: Path):
    idx = tmp_path / "idx"
    _write_index(idx, rows=[(1, 100, True, 29, 5.0)])
    out = query_genomic_coverage(str(idx), [("chr1", 5000, 6000)], site="P")
    assert out.is_empty()
    assert out.columns == ["pos", "count", "strand"]


def test_query_genomic_coverage_missing_chrom_shard_skipped(tmp_path: Path):
    idx = tmp_path / "idx"
    _write_index(idx, rows=[(1, 100, True, 29, 5.0)])
    out = query_genomic_coverage(str(idx), [("chr2", 0, 1000)], site="P")
    assert out.is_empty()


# ---------------------------------------------------------------------------
# MatrixProvider(psite_index_dir=...)
# ---------------------------------------------------------------------------


def test_matrix_provider_psite_index_mode(tmp_path: Path):
    idx = tmp_path / "idx"
    _write_index(idx, rows=[(1, 100, True, 29, 5.0), (1, 150, True, 30, 3.0)])
    provider = MatrixProvider([str(tmp_path / "fake_partition")], psite_index_dir=str(idx))
    cov = provider.coverage([Region("chr1", 90, 200)], site="P")
    assert cov["count"].sum() == 8.0


def test_matrix_provider_default_mode_unaffected(monkeypatch, tmp_path: Path):
    """Backward compatibility: without psite_index_dir, coverage() still goes
    through the original region_coverage/ref_offset path unchanged."""
    import TranslonScorer.matrix.provider as matrix_mod

    called = {}

    def _fake_region_coverage(dirs, regions, *, ref_offset, group_level, sample_names, n_workers):
        called["ref_offset"] = ref_offset
        return pl.DataFrame({"pos": [42], "count": [1.0]})

    monkeypatch.setattr("TranslonScorer.matrix.rollup.region_coverage", _fake_region_coverage)
    provider = matrix_mod.MatrixProvider([str(tmp_path)], ref_offset=15)
    cov = provider.coverage([Region("chr1", 0, 100)], site="P")
    assert called["ref_offset"] == 15
    assert cov["count"].to_list() == [1.0]


# ---------------------------------------------------------------------------
# End-to-end through _score_events_over_provider
# ---------------------------------------------------------------------------


def test_end_to_end_psite_index_scores_gapdh_shape(tmp_path: Path):
    """Same GAPDH-shaped fixture as tests/test_golden.py, sourced from a
    synthetic P-site index instead of a dict — proves the full pipeline
    (splice-flank/junction-wiring/mappability machinery, all unmodified)
    works against this new coverage source, not just query_genomic_coverage
    in isolation."""
    idx = tmp_path / "idx"
    rows = []
    for p in range(40, 99):
        rows.append((1, p, True, 29, 2.0))
    for p in range(99, 191):
        rows.append((1, p, True, 29, 30.0 if p % 3 == 0 else 5.0))
    for p in range(191, 260):
        rows.append((1, p, True, 29, 1.0))
    _write_index(idx, rows=rows)

    events = pl.DataFrame(
        [
            {
                "event_id": 1,
                "type": "init",
                "chrom": "chr1",
                "strand": 1,
                "start": 99,
                "end": 100,
                "phase": None,
            },
            {
                "event_id": 2,
                "type": "elongation",
                "chrom": "chr1",
                "strand": 1,
                "start": 100,
                "end": 190,
                "phase": 0,
            },
            {
                "event_id": 3,
                "type": "term",
                "chrom": "chr1",
                "strand": 1,
                "start": 190,
                "end": 191,
                "phase": None,
            },
        ],
        schema={
            "event_id": pl.UInt64,
            "type": pl.Utf8,
            "chrom": pl.Utf8,
            "strand": pl.Int64,
            "start": pl.Int64,
            "end": pl.Int64,
            "phase": pl.Int64,
        },
    )
    # site="P": query positions match the synthetic p_site values directly
    # (no +3 shift), so the body/flank boundaries land exactly where the
    # fixture intends — mirrors how tests/test_golden.py's coverage dict is
    # pre-positioned with no shift applied at the scoring layer.
    provider = MatrixProvider([str(tmp_path / "fake_partition")], psite_index_dir=str(idx))
    scored = _score_events_over_provider(
        events, provider, site="P", group="g", tier="t", thr=_THR
    ).sort("event_id")
    assert scored["call"].to_list() == ["SUPPORTED", "SUPPORTED", "SUPPORTED"]

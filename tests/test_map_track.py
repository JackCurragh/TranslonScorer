"""Mappability-track diagnostic annotation (map_track_mean / map_track_low).

Not the BAM-NH-tag SupportsMappability/mappability_ledger concept
(coverage/base.py) — that's per-read unique-vs-multimapper accounting. This
is a precomputed mappability TRACK (Umap/GEM-mappability style bigwig) read
via BigwigSetProvider and used purely as a diagnostic annotation: it must
never affect eligibility/call (same "review flag, not a gate" precedent as
flank_peakiness/stability), and it must actually reach the composed
per-translon report (not just the raw score store) to be useful.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

pyBigWig = pytest.importorskip("pyBigWig")

from TranslonScorer.coverage.bigwig import BigwigSetProvider
from TranslonScorer.model import ScoreThresholds
from TranslonScorer.report import compose_report
from TranslonScorer.workflows import _map_track_for_chrom, _score_events_over_provider

_THR = ScoreThresholds()


def _write_bw(path: Path, chrom: str, length: int, intervals: list) -> None:
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(chrom, length)])
    starts = [s for s, e, v in intervals]
    ends = [e for s, e, v in intervals]
    values = [float(v) for s, e, v in intervals]
    bw.addEntries([chrom] * len(intervals), starts, ends=ends, values=values)
    bw.close()


def _events() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "event_id": pl.Series([1, 2], dtype=pl.UInt64),
            "type": ["init", "elongation"],
            "chrom": ["chr1", "chr1"],
            "strand": [1, 1],
            "start": [950, 10000],
            "end": [951, 10100],
            "phase": [None, 0],
        }
    )


def test_map_track_for_chrom_real_bigwig(tmp_path: Path):
    bw = tmp_path / "map.bw"
    _write_bw(bw, "chr1", 20000, [(0, 5000, 0.1), (5000, 15000, 0.95), (15000, 20000, 0.1)])
    provider = BigwigSetProvider([str(bw)])

    out = _map_track_for_chrom(provider, _events(), _THR)
    assert out[1]["map_track_low"] is True  # window around 950 stays in the 0.1 band
    assert out[2]["map_track_low"] is False  # window around 10000-10100 in the 0.95 band


def test_map_track_disabled_when_no_provider():
    assert _map_track_for_chrom(None, _events(), _THR) == {}


def test_map_track_gap_reads_as_unknown_not_zero(tmp_path: Path):
    """A position with no interval in the mappability bigwig (wrong chrom,
    off-contig, assembly gap) must be distinguishable from a position with a
    confirmed value of 0.0 — collapsing both to the same default silently
    turns a data gap into evidence of low mappability."""
    bw = tmp_path / "map.bw"
    # Interval only covers [0, 500); the event's window (around 950) falls
    # in a true gap with no interval at all, not a confirmed 0.0.
    _write_bw(bw, "chr1", 20000, [(0, 500, 0.9)])
    provider = BigwigSetProvider([str(bw)])

    out = _map_track_for_chrom(provider, _events(), _THR)
    assert out[1]["map_track_mean"] is None
    assert out[1]["map_track_low"] is None


def test_map_track_confirmed_zero_still_reads_as_low(tmp_path: Path):
    """A position with an explicit 0.0 value (genuinely unmappable) must
    still be counted — only a true gap should read as unknown."""
    bw = tmp_path / "map.bw"
    _write_bw(bw, "chr1", 20000, [(0, 20000, 0.0)])
    provider = BigwigSetProvider([str(bw)])

    out = _map_track_for_chrom(provider, _events(), _THR)
    assert out[1]["map_track_mean"] == 0.0
    assert out[1]["map_track_low"] is True


def test_map_track_never_changes_call_or_eligibility(tmp_path: Path):
    """Same coverage, scored with vs without a mappability track: call/eligibility
    must be byte-identical — the annotation is diagnostic only."""
    bw = tmp_path / "map.bw"
    _write_bw(bw, "chr1", 300, [(0, 300, 0.05)])  # aggressively "low" everywhere
    map_provider = BigwigSetProvider([str(bw)])

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
    cov = {p: 2.0 for p in range(40, 99)}
    cov.update({p: (30.0 if p % 3 == 0 else 5.0) for p in range(99, 191)})
    cov.update({p: 1.0 for p in range(191, 260)})
    cov_rows = [{"pos": p, "count": v} for p, v in cov.items()]

    class _CovProvider:
        def coverage(self, regions, *, site="A"):
            return pl.DataFrame(cov_rows, schema={"pos": pl.Int64, "count": pl.Float64})

    without = _score_events_over_provider(
        events, _CovProvider(), site="A", group="g", tier="t", thr=_THR
    )
    withit = _score_events_over_provider(
        events, _CovProvider(), site="A", group="g", tier="t", thr=_THR, map_provider=map_provider
    )

    without = without.sort("event_id")
    withit = withit.sort("event_id")
    assert without["call"].to_list() == withit["call"].to_list()
    assert without["eligibility"].to_list() == withit["eligibility"].to_list()
    assert without["metric"].to_list() == withit["metric"].to_list()

    # map_track_* is null without a track, populated (and flagged low) with one.
    assert without["map_track_mean"].is_null().all()
    assert withit["map_track_mean"].is_not_null().all()
    assert withit["map_track_low"].to_list() == [True, True, True]


def test_map_track_survives_into_composed_report(tmp_path: Path):
    """The exact gap found this session: compose_report only ever selected a
    fixed column subset from the score store, silently dropping anything left
    in the `evidence` JSON blob. map_track_mean/map_track_low must be
    first-class columns that actually reach the per-translon report."""
    bw = tmp_path / "map.bw"
    _write_bw(bw, "chr1", 300, [(0, 300, 0.05)])
    map_provider = BigwigSetProvider([str(bw)])

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
    cov = {p: 2.0 for p in range(40, 99)}
    cov.update({p: (30.0 if p % 3 == 0 else 5.0) for p in range(99, 191)})
    cov.update({p: 1.0 for p in range(191, 260)})
    cov_rows = [{"pos": p, "count": v} for p, v in cov.items()]

    class _CovProvider:
        def coverage(self, regions, *, site="A"):
            return pl.DataFrame(cov_rows, schema={"pos": pl.Int64, "count": pl.Float64})

    scored = _score_events_over_provider(
        events, _CovProvider(), site="A", group="g", tier="t", thr=_THR, map_provider=map_provider
    )
    feature_event = pl.DataFrame(
        {
            "feature_id": ["translonA", "translonA", "translonA"],
            "event_id": pl.Series([1, 2, 3], dtype=pl.UInt64),
            "role": ["init", "elongation", "term"],
        }
    )
    report = compose_report(scored, feature_event)

    assert "init_map_track_mean" in report.columns
    assert "init_map_track_low" in report.columns
    assert "term_map_track_mean" in report.columns
    row = report.filter(pl.col("feature_id") == "translonA").row(0, named=True)
    assert row["init_map_track_low"] is True
    assert row["term_map_track_low"] is True
    assert row["init_map_track_mean"] < 0.5

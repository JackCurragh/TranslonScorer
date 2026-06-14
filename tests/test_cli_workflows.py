"""Smoke tests for the Phase 4 event-scoring CLI + workflows (T13).

Covers the new subcommands (extract-events, score-matrix, score-bams,
consequential), the shared per-chrom scoring helper, and the deprecation
notices on the legacy ORF-composite commands. Does not require external data:
the scoring helper is exercised through a tiny in-process coverage provider.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from TranslonScorer.cli import cli


# ---------------------------------------------------------------------------
# Command registration + help
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["extract-events", "score-matrix", "score-bams", "consequential"])
def test_new_subcommands_registered(name):
    runner = CliRunner()
    result = runner.invoke(cli, [name, "--help"])
    assert result.exit_code == 0, result.output
    assert name in cli.commands


@pytest.mark.parametrize("name", ["score-orfs", "feature-metrics", "orf-composite"])
def test_legacy_commands_marked_deprecated(name):
    runner = CliRunner()
    result = runner.invoke(cli, [name, "--help"])
    assert result.exit_code == 0, result.output
    assert "DEPRECATED" in result.output


# ---------------------------------------------------------------------------
# Shared per-chrom scoring helper (the core of score-matrix/score-bams)
# ---------------------------------------------------------------------------

def test_score_events_over_provider_matches_direct():
    """_score_events_over_provider reproduces a direct vectorised score."""
    from tests.test_golden import _gapdh_events, _gapdh_coverage
    from TranslonScorer.workflows import _score_events_over_provider
    from TranslonScorer.scoring.run import score_events_vectorised, DEFAULT_THRESHOLDS

    events = _gapdh_events().with_columns(pl.lit("chr12").alias("chrom"))
    cov = _gapdh_coverage()
    cov_df = pl.DataFrame(
        {"pos": list(cov.keys()), "count": [float(v) for v in cov.values()]}
    )

    class _Provider:
        def coverage(self, regions, *, site="A"):
            return cov_df

    via_helper = _score_events_over_provider(
        events, _Provider(), site="A", group="gapdh", tier="aggregate",
        thr=DEFAULT_THRESHOLDS,
    )
    direct = score_events_vectorised(
        events, cov_df, group="gapdh", tier="aggregate", thr=DEFAULT_THRESHOLDS
    )
    assert via_helper.sort("event_id").equals(direct.sort("event_id"))
    assert via_helper.height == events.height


def test_score_events_over_provider_empty():
    from TranslonScorer.workflows import _score_events_over_provider
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS

    class _Provider:
        def coverage(self, regions, *, site="A"):
            return pl.DataFrame(schema={"pos": pl.Int64, "count": pl.Float64})

    out = _score_events_over_provider(
        pl.DataFrame(schema={"event_id": pl.UInt64, "type": pl.Utf8, "chrom": pl.Utf8,
                             "strand": pl.Int64, "start": pl.Int64, "end": pl.Int64,
                             "phase": pl.Int64}),
        _Provider(), site="A", group="g", tier="aggregate", thr=DEFAULT_THRESHOLDS,
    )
    assert out.is_empty()


# ---------------------------------------------------------------------------
# consequential end-to-end via CliRunner
# ---------------------------------------------------------------------------

def test_consequential_cmd_roundtrip(tmp_path: Path):
    report = pl.DataFrame({"translon_id": ["a", "b"], "score": [0.9, 0.1]})
    report_path = tmp_path / "report.parquet"
    out_path = tmp_path / "labelled.parquet"
    report.write_parquet(report_path)

    runner = CliRunner()
    result = runner.invoke(
        cli, ["consequential", "--report", str(report_path), "--out", str(out_path)]
    )
    assert result.exit_code == 0, result.output
    out = pl.read_parquet(out_path)
    assert "consequential" in out.columns
    assert out.height == 2

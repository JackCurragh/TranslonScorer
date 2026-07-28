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


def test_pipeline_workflow_requires_exactly_one_mode(tmp_path: Path):
    from TranslonScorer.workflows import pipeline_workflow

    out = str(tmp_path / "out")
    with pytest.raises(ValueError):
        pipeline_workflow(out)  # neither matrix nor bams
    with pytest.raises(ValueError):
        pipeline_workflow(out, partition_dirs=["p"], bams=["b"])  # both coverage modes


@pytest.mark.parametrize(
    "name",
    [
        "extract-events",
        "score-matrix",
        "score-bams",
        "score-bigwig",
        "report",
        "consequential",
        "pipeline",
    ],
)
def test_new_subcommands_registered(name):
    runner = CliRunner()
    result = runner.invoke(cli, [name, "--help"])
    assert result.exit_code == 0, result.output
    assert name in cli.commands


@pytest.mark.parametrize(
    "name", ["score-orfs", "feature-metrics", "orf-composite", "all", "find-orfs"]
)
def test_deprecated_commands_removed(name):
    """The deprecated ORF-composite/legacy-pipeline commands are gone."""
    assert name not in cli.commands


# ---------------------------------------------------------------------------
# Shared per-chrom scoring helper (the core of score-matrix/score-bams)
# ---------------------------------------------------------------------------


def test_score_events_over_provider_matches_direct():
    """_score_events_over_provider reproduces a direct vectorised score."""
    from tests.test_golden import _gapdh_coverage, _gapdh_events
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS, score_events_vectorised
    from TranslonScorer.workflows import _score_events_over_provider

    events = _gapdh_events().with_columns(pl.lit("chr12").alias("chrom"))
    cov = _gapdh_coverage()
    cov_df = pl.DataFrame({"pos": list(cov.keys()), "count": [float(v) for v in cov.values()]})

    class _Provider:
        def coverage(self, regions, *, site="A"):
            return cov_df

    via_helper = _score_events_over_provider(
        events,
        _Provider(),
        site="A",
        group="gapdh",
        tier="aggregate",
        thr=DEFAULT_THRESHOLDS,
    )
    direct = score_events_vectorised(
        events, cov_df, group="gapdh", tier="aggregate", thr=DEFAULT_THRESHOLDS
    )
    assert via_helper.sort("event_id").equals(direct.sort("event_id"))
    assert via_helper.height == events.height


def test_per_event_regions_match_whole_span():
    """Padded per-event windows give IDENTICAL scores to a whole-span fetch.

    A region-respecting provider returns coverage only within the requested
    regions; scoring through it must equal scoring against the full coverage,
    proving the flank pad covers every position the scorers read.
    """
    from tests.test_golden import _gapdh_coverage, _gapdh_events
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    events = _gapdh_events().with_columns(pl.lit("chr12").alias("chrom"))
    cov = _gapdh_coverage()
    full = pl.DataFrame({"pos": list(cov.keys()), "count": [float(v) for v in cov.values()]})

    class _WholeSpan:
        def coverage(self, regions, *, site="A"):
            return full

    class _RegionRespecting:
        # only returns coverage at positions inside the requested regions
        def coverage(self, regions, *, site="A"):
            keep = full
            mask = pl.lit(False)
            for r in regions:
                mask = mask | ((pl.col("pos") >= r.start) & (pl.col("pos") < r.end))
            return keep.filter(mask)

    kw = dict(site="A", group="gapdh", tier="aggregate", thr=DEFAULT_THRESHOLDS)
    whole = _score_events_over_provider(events, _WholeSpan(), **kw)
    windowed = _score_events_over_provider(events, _RegionRespecting(), **kw)
    assert windowed.sort("event_id").equals(whole.sort("event_id"))


def test_score_events_strand_isolation_is_deterministic():
    """+ and - events must be scored against their OWN strand coverage — never
    cross-contaminated (the cause of the matrix nondeterminism bug)."""
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    # two init events at the SAME genomic position, opposite strands
    events = pl.DataFrame(
        {
            "event_id": pl.Series([1, 2], dtype=pl.UInt64),
            "type": ["init", "init"],
            "chrom": ["chr1", "chr1"],
            "strand": [1, -1],
            "start": [1000, 1000],
            "end": [1001, 1001],
            "phase": [None, None],
        }
    )
    # strand-aware coverage: heavy on +, empty on - at the body positions
    pos = list(range(900, 1100))
    cov = pl.DataFrame(
        {
            "strand": [1] * len(pos) + [-1] * len(pos),
            "pos": pos + pos,
            "count": [100.0] * len(pos) + [0.0] * len(pos),
        }
    )

    class _StrandProvider:
        def coverage(self, regions, *, site="A"):
            return cov

    out1 = _score_events_over_provider(
        events, _StrandProvider(), site="A", group="g", tier="aggregate", thr=DEFAULT_THRESHOLDS
    )
    out2 = _score_events_over_provider(
        events, _StrandProvider(), site="A", group="g", tier="aggregate", thr=DEFAULT_THRESHOLDS
    )
    # deterministic, and the two strands get DIFFERENT n_reads (not mixed)
    assert out1.sort("event_id").equals(out2.sort("event_id"))
    nreads = dict(zip(out1["event_id"].to_list(), out1["n_reads"].to_list()))
    assert nreads[1] != nreads[2]


def test_merge_event_regions_merges_and_pads():
    from TranslonScorer.workflows import _EVENT_FLANK_PAD, _merge_event_regions

    # two close events merge into one padded interval; a far one stays separate
    regs = _merge_event_regions("chr1", [1000, 1010, 50000], [1001, 1011, 50001])
    assert len(regs) == 2
    assert regs[0].start == 1000 - _EVENT_FLANK_PAD
    assert regs[0].end == 1011 + _EVENT_FLANK_PAD
    assert regs[1].start == 50000 - _EVENT_FLANK_PAD


def test_score_events_over_provider_empty():
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    class _Provider:
        def coverage(self, regions, *, site="A"):
            return pl.DataFrame(schema={"pos": pl.Int64, "count": pl.Float64})

    out = _score_events_over_provider(
        pl.DataFrame(
            schema={
                "event_id": pl.UInt64,
                "type": pl.Utf8,
                "chrom": pl.Utf8,
                "strand": pl.Int64,
                "start": pl.Int64,
                "end": pl.Int64,
                "phase": pl.Int64,
            }
        ),
        _Provider(),
        site="A",
        group="g",
        tier="aggregate",
        thr=DEFAULT_THRESHOLDS,
    )
    assert out.is_empty()


# ---------------------------------------------------------------------------
# consequential end-to-end via CliRunner
# ---------------------------------------------------------------------------


def test_report_workflow_end_to_end(tmp_path: Path):
    """persist_scores + feature_event on disk -> report_workflow -> consequential."""
    from TranslonScorer.io.store import persist_scores
    from TranslonScorer.workflows import report_workflow

    scores = pl.DataFrame(
        {
            "event_id": [1, 2, 3],
            "aspect": ["init", "elongation", "term"],
            "group": ["aggregate"] * 3,
            "tier": ["aggregate"] * 3,
            "n_reads": [100.0, 200.0, 50.0],
            "metric": [2.0, 0.8, 1.5],
            "metric_name": ["rise", "elong_in_frame", "drop"],
            "eligibility": ["ELIGIBLE"] * 3,
            "call": ["SUPPORTED", "SUPPORTED", "SUPPORTED"],
            "evidence": ["{}"] * 3,
            "thresholds_version": ["v1"] * 3,
        }
    )
    store = tmp_path / "store"
    persist_scores(scores, str(store), data_version="d1")

    events_dir = tmp_path / "events"
    fe_dir = events_dir / "feature_event"
    fe_dir.mkdir(parents=True)
    pl.DataFrame(
        {
            "feature_id": ["A", "A", "A"],
            "event_id": [1, 2, 3],
            "role": ["init", "elongation", "term"],
        }
    ).write_parquet(fe_dir / "chr1.parquet")

    out_path = tmp_path / "report.parquet"
    rep = report_workflow(str(store), str(events_dir), str(out_path), data_version="d1")
    assert out_path.exists()
    assert rep.filter(pl.col("feature_id") == "A")["consequential"][0] is True
    assert "consequentiality_score" in rep.columns


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


# ---------------------------------------------------------------------------
# Junction support wiring (previously never called — every junction event
# scored against an empty support dict regardless of provider capability)
# ---------------------------------------------------------------------------


def test_junction_events_scored_via_provider_junction_support():
    """_score_events_over_provider must call provider.junction_support() and
    feed it into scoring — junction events should not be silently INSUFFICIENT
    just because the provider is capable of answering."""
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    events = pl.DataFrame(
        {
            "event_id": pl.Series([42], dtype=pl.UInt64),
            "type": ["junction"],
            "chrom": ["chr1"],
            "strand": [1],
            "start": [1000],  # donor
            "end": [1200],  # acceptor
            "phase": [None],
        }
    )

    class _JunctionCapableProvider:
        def coverage(self, regions, *, site="A"):
            # non-empty (but irrelevant to this junction event) so the chrom
            # isn't short-circuited by the empty-coverage skip.
            return pl.DataFrame({"pos": [1], "count": [1.0]})

        def junction_support(self, junctions, *, by_sample=False):
            assert junctions == [("chr1", 1000, 1200, 1, 42)]
            return pl.DataFrame(
                {
                    "junction_id": pl.Series([42], dtype=pl.UInt64),
                    "kind": ["span_conf"],
                    "count": [25.0],
                }
            )

    scored = _score_events_over_provider(
        events,
        _JunctionCapableProvider(),
        site="A",
        group="g",
        tier="t",
        thr=DEFAULT_THRESHOLDS,
    )
    row = scored.filter(pl.col("event_id") == 42).row(0, named=True)
    assert row["call"] == "SUPPORTED"
    assert row["n_reads"] == 25.0


def test_junction_events_degrade_gracefully_without_junction_support():
    """A provider that doesn't expose junction_support() (a coverage-only
    source) must not crash — the hasattr gate skips it and junction events
    just fall out INSUFFICIENT, as before."""
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    events = pl.DataFrame(
        {
            "event_id": pl.Series([42], dtype=pl.UInt64),
            "type": ["junction"],
            "chrom": ["chr1"],
            "strand": [1],
            "start": [1000],
            "end": [1200],
            "phase": [None],
        }
    )

    class _CoverageOnlyProvider:
        # non-empty (but irrelevant to this junction event) so the chrom
        # isn't short-circuited by the empty-coverage skip, and the
        # junction_support code path is actually exercised.
        def coverage(self, regions, *, site="A"):
            return pl.DataFrame({"pos": [1], "count": [1.0]})

    scored = _score_events_over_provider(
        events,
        _CoverageOnlyProvider(),
        site="A",
        group="g",
        tier="t",
        thr=DEFAULT_THRESHOLDS,
    )
    row = scored.filter(pl.col("event_id") == 42).row(0, named=True)
    assert row["eligibility"] == "INSUFFICIENT"
    assert row["n_reads"] == 0.0


def test_junction_events_dont_crash_provider_that_raises_not_implemented():
    """A provider (e.g. BigwigSetProvider) that DEFINES junction_support() but
    always raises NotImplementedError must not crash scoring — the hasattr
    gate only checks the method exists, not that it works, so
    NotImplementedError is caught explicitly too."""
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    events = pl.DataFrame(
        {
            "event_id": pl.Series([42], dtype=pl.UInt64),
            "type": ["junction"],
            "chrom": ["chr1"],
            "strand": [1],
            "start": [1000],
            "end": [1200],
            "phase": [None],
        }
    )

    class _RaisingJunctionProvider:
        def coverage(self, regions, *, site="A"):
            return pl.DataFrame({"pos": [1], "count": [1.0]})

        def junction_support(self, junctions, *, by_sample=False):
            raise NotImplementedError("bigwigs have no CIGAR")

    scored = _score_events_over_provider(
        events,
        _RaisingJunctionProvider(),
        site="A",
        group="g",
        tier="t",
        thr=DEFAULT_THRESHOLDS,
    )
    row = scored.filter(pl.col("event_id") == 42).row(0, named=True)
    assert row["eligibility"] == "INSUFFICIENT"


# ---------------------------------------------------------------------------
# Metagene offset calibration wiring (previously score_bams_workflow never
# supplied genomic start-codon coordinates, so --offset-method metagene
# silently fell back to the fixed global offset every time)
# ---------------------------------------------------------------------------


def test_score_bams_workflow_supplies_start_codons_for_metagene(tmp_path: Path, monkeypatch):
    import TranslonScorer.coverage.bam as bam_mod
    import TranslonScorer.workflows as workflows_mod
    from TranslonScorer.model import OffsetParams

    events = pl.DataFrame(
        {
            "event_id": pl.Series([1, 2], dtype=pl.UInt64),
            "type": ["init", "elongation"],
            "chrom": ["chr1", "chr1"],
            "strand": [1, 1],
            "start": [500, 600],
            "end": [501, 700],
            "phase": [None, 0],
        }
    )
    monkeypatch.setattr(workflows_mod, "read_events", lambda events_dir: events)

    captured: dict = {}

    class _FakeProvider:
        def __init__(self, bams, **kwargs):
            captured.update(kwargs)

        def coverage(self, regions, *, site="A"):
            return pl.DataFrame(schema={"pos": pl.Int64, "count": pl.Float64})

    monkeypatch.setattr(bam_mod, "BamSetProvider", _FakeProvider)

    workflows_mod.score_bams_workflow(
        "unused_events_dir",
        ["fake.bam"],
        str(tmp_path / "scores"),
        data_version="d1",
        offsets=OffsetParams(method="metagene"),
    )
    assert captured["start_codons"] == [("chr1", 500, 1)]


def test_score_bams_workflow_omits_start_codons_for_global_method(tmp_path: Path, monkeypatch):
    """No behaviour change for the (default) global/file methods."""
    import TranslonScorer.coverage.bam as bam_mod
    import TranslonScorer.workflows as workflows_mod
    from TranslonScorer.model import OffsetParams

    events = pl.DataFrame(
        {
            "event_id": pl.Series([1], dtype=pl.UInt64),
            "type": ["init"],
            "chrom": ["chr1"],
            "strand": [1],
            "start": [500],
            "end": [501],
            "phase": [None],
        }
    )
    monkeypatch.setattr(workflows_mod, "read_events", lambda events_dir: events)

    captured: dict = {}

    class _FakeProvider:
        def __init__(self, bams, **kwargs):
            captured.update(kwargs)

        def coverage(self, regions, *, site="A"):
            return pl.DataFrame(schema={"pos": pl.Int64, "count": pl.Float64})

    monkeypatch.setattr(bam_mod, "BamSetProvider", _FakeProvider)

    workflows_mod.score_bams_workflow(
        "unused_events_dir",
        ["fake.bam"],
        str(tmp_path / "scores"),
        data_version="d1",
        offsets=OffsetParams(method="global"),
    )
    assert captured["start_codons"] is None


# ---------------------------------------------------------------------------
# score-bigwig CLI validation
# ---------------------------------------------------------------------------


def test_score_bigwig_cmd_rejects_both_plain_and_stranded(tmp_path: Path):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "score-bigwig",
            "--events-dir",
            str(tmp_path / "events"),
            "--store-dir",
            str(tmp_path / "scores"),
            "--data-version",
            "d1",
            "--bigwig",
            "a.bw",
            "--forward-bigwig",
            "f.bw",
            "--reverse-bigwig",
            "r.bw",
        ],
    )
    assert result.exit_code != 0
    assert "not both, not neither" in result.output


def test_score_bigwig_cmd_rejects_neither_mode(tmp_path: Path):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "score-bigwig",
            "--events-dir",
            str(tmp_path / "events"),
            "--store-dir",
            str(tmp_path / "scores"),
            "--data-version",
            "d1",
        ],
    )
    assert result.exit_code != 0
    assert "not both, not neither" in result.output


def test_score_bigwig_cmd_rejects_mismatched_strand_pair_counts(tmp_path: Path):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "score-bigwig",
            "--events-dir",
            str(tmp_path / "events"),
            "--store-dir",
            str(tmp_path / "scores"),
            "--data-version",
            "d1",
            "--forward-bigwig",
            "f1.bw",
            "--forward-bigwig",
            "f2.bw",
            "--reverse-bigwig",
            "r1.bw",
        ],
    )
    assert result.exit_code != 0
    assert "counts must match" in result.output


# ---------------------------------------------------------------------------
# --psite-index is the only matrix coverage strategy: required, and forwarded
# ---------------------------------------------------------------------------


def test_score_matrix_cmd_requires_psite_index(tmp_path: Path):
    """Omitting --psite-index is a usage error, not a silent fall back to a
    flat offset (which would pin elong_in_frame to the ~0.33 random floor)."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "score-matrix",
            "--events-dir",
            str(tmp_path / "events"),
            "--matrix-dir",
            str(tmp_path / "matrix"),
            "--store-dir",
            str(tmp_path / "scores"),
            "--data-version",
            "d1",
        ],
    )
    assert result.exit_code != 0
    assert "psite-index" in result.output


def test_score_matrix_cmd_forwards_psite_index(monkeypatch, tmp_path: Path):
    import TranslonScorer.io.matrix as io_matrix
    import TranslonScorer.workflows as workflows_mod

    monkeypatch.setattr(io_matrix, "discover_partitions", lambda matrix_dir: ["p1"])
    captured = {}

    def _fake_score_matrix_workflow(*args, **kwargs):
        captured.update(kwargs)
        return "fake_store"

    monkeypatch.setattr(workflows_mod, "score_matrix_workflow", _fake_score_matrix_workflow)

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "score-matrix",
            "--events-dir",
            str(tmp_path / "events"),
            "--matrix-dir",
            str(tmp_path / "matrix"),
            "--store-dir",
            str(tmp_path / "scores"),
            "--data-version",
            "d1",
            "--psite-index",
            str(tmp_path / "psite"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["psite_index_dir"] == str(tmp_path / "psite")


def test_pipeline_cmd_forwards_psite_index_dir(monkeypatch, tmp_path: Path):
    import TranslonScorer.io.matrix as io_matrix
    import TranslonScorer.workflows as workflows_mod

    monkeypatch.setattr(io_matrix, "discover_partitions", lambda matrix_dir: ["p1"])
    captured = {}

    def _fake_pipeline_workflow(*args, **kwargs):
        captured.update(kwargs)
        return {"events": "e", "scores": "s", "report": "r"}

    monkeypatch.setattr(workflows_mod, "pipeline_workflow", _fake_pipeline_workflow)

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "pipeline",
            "--out-dir",
            str(tmp_path / "out"),
            "--sqlite",
            str(tmp_path / "annot.sqlite"),
            "--matrix-dir",
            str(tmp_path / "matrix"),
            "--psite-index",
            str(tmp_path / "psite"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["psite_index_dir"] == str(tmp_path / "psite")

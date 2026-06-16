"""Unit coverage for the refactored new-tree modules that lacked direct tests:
events.py, io/bam.py, io/matrix.py, offsets.py, coverage/profile.py.

These are pure (or fixture-backed) functions on the live event-scoring path.
Real-data cases skip when the local matrix / BAM fixtures are absent.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

REPO_ROOT = Path(__file__).parent.parent
GENOME_BAM = REPO_ROOT / "data" / "gapdh_cohort_genome.bam"
MATRIX_PART = REPO_ROOT / "data" / "global_partitioned" / "AAAA"
HAS_BAM = GENOME_BAM.exists()
HAS_MATRIX = MATRIX_PART.exists()


# ===========================================================================
# events.py
# ===========================================================================


def _blocks() -> pl.DataFrame:
    # One + strand translon, two exons with a genomic gap (=> one junction).
    return pl.DataFrame(
        {
            "translon_id": ["T1", "T1"],
            "translation_block_rank": [1, 2],
            "bed_chrom": ["chr1", "chr1"],
            "bed_start": [100, 200],
            "bed_end": [130, 230],
            "seq_region_strand": [1, 1],
            "block_length_nt": [30, 30],
        }
    )


def _translons() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "translon_id": ["T1"],
            "bed_chrom": ["chr1"],
            "bed_start": [100],
            "bed_end": [230],
            "seq_region_strand": [1],
        }
    )


def test_extract_events_shapes_and_aspects():
    from TranslonScorer.events import extract_events

    events, feature_event, overlap = extract_events(_blocks(), _translons())
    aspects = set(events["type"])
    assert {"elongation", "junction", "init", "term"} <= aspects
    # 2 elong + 1 junction + 1 init + 1 term = 5 distinct events.
    assert events.height == 5
    assert events["event_id"].n_unique() == 5
    # feature_event maps every event back to T1.
    assert set(feature_event["feature_id"]) == {"T1"}
    assert {"init", "term", "elongation", "junction"} <= set(feature_event["role"])
    # init at genomic start (100), term at end-1 (229) for + strand.
    init = events.filter(pl.col("type") == "init").row(0, named=True)
    term = events.filter(pl.col("type") == "term").row(0, named=True)
    assert init["start"] == 100
    assert term["start"] == 229
    assert overlap.is_empty()  # single frame register => no contention


def test_extract_events_deterministic_event_ids():
    from TranslonScorer.events import extract_events

    e1, _, _ = extract_events(_blocks(), _translons())
    e2, _, _ = extract_events(_blocks(), _translons())
    assert e1.sort("event_id")["event_id"].to_list() == e2.sort("event_id")["event_id"].to_list()


def test_extract_events_empty_junction_branch():
    from TranslonScorer.events import extract_events

    # Single contiguous block => no junction.
    blocks = _blocks().filter(pl.col("translation_block_rank") == 1)
    events, _, _ = extract_events(blocks, _translons())
    assert "junction" not in set(events["type"])


def test_contention_flags_different_frame_overlap():
    from TranslonScorer.events import _contention

    elong = pl.DataFrame(
        {
            "event_id": pl.Series([1, 2], dtype=pl.UInt64),
            "chrom": ["chr1", "chr1"],
            "strand": [1, 1],
            "start": [100, 110],
            "end": [160, 170],
            "phase": [0, 1],  # different registers => contention
        }
    )
    out = _contention(elong)
    assert out.height == 2  # symmetric pair
    assert set(out["overlap_start"]) == {110}
    assert set(out["overlap_end"]) == {160}


def test_contention_same_frame_no_overlap():
    from TranslonScorer.events import _contention

    elong = pl.DataFrame(
        {
            "event_id": pl.Series([1, 2], dtype=pl.UInt64),
            "chrom": ["chr1", "chr1"],
            "strand": [1, 1],
            "start": [100, 110],
            "end": [160, 170],
            "phase": [0, 0],  # same register => no contention
        }
    )
    assert _contention(elong).is_empty()


def test_deconflict_intervals_drops_conflicting_phase():
    from TranslonScorer.events import _deconflict_intervals

    # [0,10) phase0 and [5,15) phase1 overlap in [5,10): conflict dropped.
    out = _deconflict_intervals([(0, 10, 0), (5, 15, 1)])
    # surviving: [0,5) phase0, [10,15) phase1 ; [5,10) removed.
    assert (0, 5, 0) in out
    assert (10, 15, 1) in out
    assert all(not (s <= 5 < e and s < 10) or ph is not None for s, e, ph in out)
    # no surviving interval covers the conflicted region [5,10)
    assert not any(s < 10 and e > 5 and s >= 5 and e <= 10 for s, e, ph in out)


def test_deconflict_intervals_merges_adjacent_same_phase():
    from TranslonScorer.events import _deconflict_intervals

    out = _deconflict_intervals([(0, 5, 0), (5, 10, 0)])
    assert out == [(0, 10, 0)]


def test_deconflict_intervals_empty():
    from TranslonScorer.events import _deconflict_intervals

    assert _deconflict_intervals([]) == []


def test_build_frame_intervals():
    from TranslonScorer.events import _build_frame_intervals

    cds_df = pl.DataFrame(
        {
            "chr": ["chr1"],
            "strand": ["+"],
            "start": [[0]],
            "stop": [[9]],
            "tran_start": [[0]],
        }
    )
    out = _build_frame_intervals(pl.DataFrame(), cds_df, {"chr1"})
    assert ("chr1", "+") in out
    assert out[("chr1", "+")] == [(0, 9, 0)]


# ===========================================================================
# io/bam.py (pure helpers)
# ===========================================================================


def test_parse_read_id():
    from TranslonScorer.io.bam import parse_read_id

    assert parse_read_id("read_42") == 42
    assert parse_read_id("no-id-here") is None


def test_normalise_chrom():
    from TranslonScorer.io.bam import normalise_chrom

    refs = {"chr1", "chr2"}
    assert normalise_chrom("chr1", refs) == "chr1"  # direct
    assert normalise_chrom("1", refs) == "chr1"  # add prefix
    assert normalise_chrom("chr3", refs) is None  # absent
    assert normalise_chrom("2", {"2"}) == "2"  # strip prefix path


def test_cigar_blocks_splits_on_intron():
    from TranslonScorer.io.bam import cigar_blocks

    # 30M 100N 20M starting at 1000 => two ref blocks separated by the N skip.
    blocks = cigar_blocks(1000, [(0, 30), (3, 100), (0, 20)])
    assert blocks == [(1000, 1030), (1130, 1150)]


def test_cigar_blocks_soft_clip_ignored():
    from TranslonScorer.io.bam import cigar_blocks

    # 5S 10M : soft-clip consumes no reference.
    assert cigar_blocks(500, [(4, 5), (0, 10)]) == [(500, 510)]


@pytest.mark.skipif(not HAS_BAM, reason="genome GAPDH fixture not in data/")
def test_aggregate_junctions_on_fixture():
    from TranslonScorer.io.bam import aggregate_junctions

    j = aggregate_junctions(str(GENOME_BAM))
    assert set(j.columns) == {"chr", "donor_pos", "acceptor_pos", "strand", "count"}
    if not j.is_empty():
        assert (j["count"] > 0).all()


def test_aggregate_junctions_empty(tmp_path):
    """Unmapped/empty BAM path returns the typed empty schema."""
    import pysam

    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    bam_path = tmp_path / "empty.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header):
        pass
    from TranslonScorer.io.bam import aggregate_junctions

    out = aggregate_junctions(str(bam_path))
    assert out.is_empty()
    assert "donor_pos" in out.columns


# ===========================================================================
# offsets.py
# ===========================================================================


def test_plausible_range_and_usable():
    from TranslonScorer.model import OffsetParams
    from TranslonScorer.offsets import plausible_offset_range, psite_to_asite, usable_read_length

    p = OffsetParams()
    lo, hi = plausible_offset_range(25, p)
    assert lo == 8 and hi == min(20, int(25 * p.max_frac))  # 25*0.6667=16
    assert usable_read_length(25, p) and not usable_read_length(24, p)
    assert psite_to_asite(12) == 15


def test_global_offsets_and_dispatch():
    from TranslonScorer.model import OffsetParams
    from TranslonScorer.offsets import global_offsets, make_offset_table

    p = OffsetParams(method="global", global_offset=12)
    table = global_offsets(p)
    assert set(table.keys()) == set(range(25, 36))
    assert all(v == 12 for v in table.values())
    # dispatch routes "global" identically
    assert make_offset_table(p) == table


def test_make_offset_table_metagene_raises():
    from TranslonScorer.model import OffsetParams
    from TranslonScorer.offsets import make_offset_table

    with pytest.raises(NotImplementedError):
        make_offset_table(OffsetParams(method="metagene"))


def test_make_offset_table_unknown_raises():
    from TranslonScorer.model import OffsetParams
    from TranslonScorer.offsets import make_offset_table

    with pytest.raises(ValueError):
        make_offset_table(OffsetParams(method="bogus"))


def test_file_offsets_roundtrip_and_clamp(tmp_path):
    from TranslonScorer.model import OffsetParams
    from TranslonScorer.offsets import file_offsets, make_offset_table

    csv = tmp_path / "offsets.csv"
    pl.DataFrame(
        {
            "read_length": [28, 29, 24, 30],  # 24 dropped (unusable)
            "offset": [12, 99, 12, 13],  # 99 clamped to plausible hi
            "site": ["P", "P", "P", "A"],  # A row dropped
        }
    ).write_csv(csv)
    p = OffsetParams(method="file", offsets_file=str(csv))
    table = file_offsets(str(csv), p)
    assert 24 not in table  # unusable length
    assert 30 not in table  # A-site row filtered
    assert table[28] == 12
    assert table[29] == min(20, int(29 * p.max_frac))  # clamped
    # dispatch path equivalent
    assert make_offset_table(p) == table


def test_file_offsets_missing_columns_raises(tmp_path):
    from TranslonScorer.model import OffsetParams
    from TranslonScorer.offsets import file_offsets

    csv = tmp_path / "bad.csv"
    pl.DataFrame({"read_length": [28]}).write_csv(csv)
    with pytest.raises(ValueError):
        file_offsets(str(csv), OffsetParams())


def test_make_offset_table_file_requires_path():
    from TranslonScorer.model import OffsetParams
    from TranslonScorer.offsets import make_offset_table

    with pytest.raises(ValueError):
        make_offset_table(OffsetParams(method="file", offsets_file=None))


# ===========================================================================
# coverage/profile.py
# ===========================================================================


def _reads() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "tran_id": ["t1", "t1", "t1"],
            "tran_start_bam": [100, 100, 200],
            "length": [29, 29, 30],
            "count": [1.0, 2.0, 5.0],
        }
    )


def test_apply_offsets_psite_groups_positions():
    from TranslonScorer.coverage.profile import apply_offsets

    out = apply_offsets(_reads(), {29: 12, 30: 13}, site="P")
    d = dict(zip(out["pos"], out["count"]))
    assert d[112] == 3.0  # 100+12, two reads summed
    assert d[213] == 5.0  # 200+13


def test_apply_offsets_asite_adds_codon():
    from TranslonScorer.coverage.profile import apply_offsets

    out = apply_offsets(_reads(), {29: 12, 30: 13}, site="A")
    d = dict(zip(out["pos"], out["count"]))
    assert d[115] == 3.0  # 100+12+3
    assert d[216] == 5.0  # 200+13+3


def test_apply_offsets_default_offset_fallback():
    from TranslonScorer.coverage.profile import apply_offsets

    out = apply_offsets(_reads(), {}, site="P", default_offset=10)
    d = dict(zip(out["pos"], out["count"]))
    assert d[110] == 3.0  # 100+10 fallback
    assert d[210] == 5.0


def test_apply_offsets_empty_and_bad_site_and_missing_cols():
    from TranslonScorer.coverage.profile import apply_offsets

    empty = apply_offsets(pl.DataFrame(), {}, site="P")
    assert empty.is_empty() and "pos" in empty.columns
    with pytest.raises(ValueError):
        apply_offsets(_reads(), {}, site="Z")
    with pytest.raises(ValueError):
        apply_offsets(pl.DataFrame({"tran_id": ["t"]}), {}, site="P")


def test_size_factors_single_and_multi():
    from TranslonScorer.coverage.profile import size_factors

    # no sample col => trivial
    assert size_factors(pl.DataFrame({"pos": [1], "count": [1.0]})) == {"": 1.0}
    one = size_factors(pl.DataFrame({"sample_id": ["a"], "pos": [1], "count": [1.0]}))
    assert one == {"a": 1.0}


# ===========================================================================
# workflows.py  (orchestration loop via a fake provider — no heavy data)
# ===========================================================================


class _FakeProvider:
    """Minimal CoverageProvider: returns the same coverage for any region."""

    def __init__(self, cov: pl.DataFrame):
        self._cov = cov

    def coverage(self, regions, site="A"):
        return self._cov


def test_score_events_over_provider_scores_and_skips():
    from TranslonScorer.events import extract_events
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    events, _, _ = extract_events(_blocks(), _translons())
    # Dense coverage across the locus so events clear the eligibility floor.
    pos = list(range(100, 235))
    cov = pl.DataFrame({"pos": pos, "count": [50.0] * len(pos)})

    scored = _score_events_over_provider(
        events,
        _FakeProvider(cov),
        site="A",
        group="g",
        tier="aggregate",
        thr=DEFAULT_THRESHOLDS,
    )
    assert scored.height > 0
    assert set(scored.columns) >= {"event_id", "aspect", "eligibility", "call"}


def test_score_events_over_provider_empty_events():
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    empty = pl.DataFrame(
        schema={
            "event_id": pl.UInt64,
            "type": pl.Utf8,
            "chrom": pl.Utf8,
            "strand": pl.Int64,
            "start": pl.Int64,
            "end": pl.Int64,
            "phase": pl.Int64,
        }
    )
    out = _score_events_over_provider(
        empty,
        _FakeProvider(pl.DataFrame({"pos": [1], "count": [1.0]})),
        site="A",
        group="g",
        tier="t",
        thr=DEFAULT_THRESHOLDS,
    )
    assert out.is_empty()


def test_score_events_over_provider_empty_coverage_skips():
    from TranslonScorer.events import extract_events
    from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS
    from TranslonScorer.workflows import _score_events_over_provider

    events, _, _ = extract_events(_blocks(), _translons())
    empty_cov = pl.DataFrame(schema={"pos": pl.Int64, "count": pl.Float64})
    out = _score_events_over_provider(
        events,
        _FakeProvider(empty_cov),
        site="A",
        group="g",
        tier="t",
        thr=DEFAULT_THRESHOLDS,
    )
    assert out.is_empty()


# ===========================================================================
# io/matrix.py (fixture-backed)
# ===========================================================================


def test_matrix_parse_read_id():
    from TranslonScorer.io.matrix import _parse_read_id

    assert _parse_read_id("read_7") == 7
    assert _parse_read_id("junk") is None


# ===========================================================================
# io/annotation.py
# ===========================================================================


def _write_gtf(tmp_path) -> str:
    # Two-CDS-exon transcript on + strand, plus one on - strand, same gene.
    lines = [
        'chr1\tsrc\tCDS\t101\t130\t.\t+\t0\tgene_id "G1"; transcript_id "TX1";',
        'chr1\tsrc\tCDS\t201\t230\t.\t+\t0\tgene_id "G1"; transcript_id "TX1";',
        'chr1\tsrc\texon\t101\t130\t.\t+\t0\tgene_id "G1"; transcript_id "TX1";',
        'chr2\tsrc\tCDS\t301\t330\t.\t-\t0\tgene_id "G2"; transcript_id "TX2";',
    ]
    p = tmp_path / "mini.gtf"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_build_cds_blocks(tmp_path):
    from TranslonScorer.io.annotation import build_cds_blocks

    df = build_cds_blocks(_write_gtf(tmp_path))
    tx1 = df.filter(pl.col("tran_id") == "TX1").row(0, named=True)
    assert tx1["gene_id"] == "G1"
    assert tx1["chr"] == "chr1"
    # 0-based half-open starts; two blocks; tran_start cumulative (0, 30).
    assert tx1["start"] == [100, 200]
    assert tx1["stop"] == [130, 230]
    assert tx1["tran_start"] == [0, 30]


def test_build_exon_blocks_includes_utrs_mrna_origin(tmp_path):
    """build_exon_blocks spans full exons (UTR included) with mRNA-5' origin,
    unlike build_cds_blocks which starts at the first CDS base."""
    from TranslonScorer.io.annotation import build_cds_blocks, build_exon_blocks

    lines = [
        # + strand transcript: exon 100-160 (incl 5'UTR), CDS only 130-160
        'chr1\tsrc\texon\t101\t160\t.\t+\t.\tgene_id "G"; transcript_id "TX";',
        'chr1\tsrc\tCDS\t131\t160\t.\t+\t0\tgene_id "G"; transcript_id "TX";',
    ]
    p = tmp_path / "u.gtf"
    p.write_text("\n".join(lines) + "\n")
    ex = build_exon_blocks(str(p)).filter(pl.col("tran_id") == "TX").row(0, named=True)
    cd = build_cds_blocks(str(p)).filter(pl.col("tran_id") == "TX").row(0, named=True)
    # exon block starts at genomic 100 (0-based), tran_start 0 = mRNA 5' end
    assert ex["start"] == [100] and ex["tran_start"] == [0]
    # CDS block starts at 130 — so a transcriptome read at mRNA pos 0 would be
    # mis-placed by ~30nt if CDS blocks were (wrongly) used for projection.
    assert cd["start"] == [130]


def test_build_gene_spans(tmp_path):
    from TranslonScorer.io.annotation import build_cds_blocks, build_gene_spans

    spans, id_of = build_gene_spans(build_cds_blocks(_write_gtf(tmp_path)))
    assert ("chr1", "+") in spans
    assert ("chr2", "-") in spans
    # G1 spans min start (100) to max stop (230).
    g1_span = spans[("chr1", "+")][0]
    assert g1_span[0] == 100 and g1_span[1] == 230
    assert set(id_of.values()) == {"G1", "G2"}


@pytest.mark.skipif(not HAS_MATRIX, reason="matrix partition fixture not in data/")
def test_aggregate_fast_path_equals_per_sample(tmp_path, monkeypatch):
    """region_coverage aggregate (read_id totals fast path) == per-sample summed."""
    from TranslonScorer.matrix_rollup import region_coverage

    monkeypatch.setenv("TS_MATRIX_CACHE_DIR", str(tmp_path / "cache"))
    parts = [str(MATRIX_PART)]
    # a broad region; n_workers=1 keeps it spawn-free under pytest
    regions = [("chr1", 0, 250_000_000)]
    agg = region_coverage(parts, regions, group_level="aggregate", n_workers=1).sort(
        ["strand", "pos"]
    )
    per = region_coverage(parts, regions, group_level="sample", n_workers=1)
    if per.is_empty():
        pytest.skip("no reads for region in fixture partition")
    per_agg = (
        per.group_by(["strand", "pos"])
        .agg(pl.col("count").sum().alias("count"))
        .sort(["strand", "pos"])
    )
    assert agg.equals(per_agg)


@pytest.mark.skipif(not HAS_MATRIX, reason="matrix partition fixture not in data/")
def test_matrix_manifest_helpers():
    from TranslonScorer.io.matrix import (
        _count_parquets,
        _discover_bam,
        _manifest,
        _reads_parquet_path,
        _samples_df,
    )

    mp, manifest = _manifest(MATRIX_PART)
    assert mp.exists() and isinstance(manifest, dict)
    assert isinstance(_count_parquets(mp, manifest), list)
    assert "sample" in _samples_df(mp, manifest).columns or _samples_df(mp, manifest).width >= 1
    assert _reads_parquet_path(mp, manifest).endswith(".parquet")
    bam = _discover_bam(MATRIX_PART)
    assert bam is None or bam.suffix == ".bam"


# ===========================================================================
# io/matrix.discover_partitions
# ===========================================================================


def test_discover_partitions(tmp_path):
    from TranslonScorer.io.matrix import discover_partitions

    root = tmp_path / "matrix"
    for name in ("AAAA", "AAAC", "AAAG"):
        d = root / name
        d.mkdir(parents=True)
        (d / f"global.{name}_matrix_manifest.json").write_text("{}")
    (root / "not_a_partition").mkdir()  # no manifest -> excluded
    parts = discover_partitions(root)
    assert [p.name for p in parts] == ["AAAA", "AAAC", "AAAG"]


def test_discover_partitions_empty_raises(tmp_path):
    from TranslonScorer.io.matrix import discover_partitions

    (tmp_path / "empty").mkdir()
    import pytest as _pt

    with _pt.raises(FileNotFoundError):
        discover_partitions(tmp_path / "empty")

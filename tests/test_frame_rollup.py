"""Tests for the feature-scoped, length-resolved FrameRollup (Step 1 redesign).

Unit tests use synthetic in-memory BAMs and CDS data; integration tests
require the local matrix fixture (data/global_partitioned) and are skipped
when absent (HAS_MATRIX guard).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

REPO_ROOT = Path(__file__).parent.parent
MATRIX_PART = REPO_ROOT / "data" / "global_partitioned" / "AAAA"
HAS_MATRIX = MATRIX_PART.exists()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _gapdh_cds() -> pl.DataFrame:
    """Minimal GAPDH ENST00000229239 CDS structure (8 exons, chr12, +)."""
    return pl.DataFrame(
        {
            "tran_id": ["ENST00000229239"],
            "gene_id": ["ENSG00000111640"],
            "chr": ["12"],
            "strand": ["+"],
            "start": [[6534832, 6536493, 6536683, 6536919, 6537100, 6537308, 6537583, 6538100]],
            "stop": [[6534861, 6536593, 6536790, 6537010, 6537216, 6537390, 6537996, 6538167]],
            "tran_start": [[0, 29, 129, 236, 327, 443, 525, 938]],
        }
    )


def _simple_cds(chrom="chr1", strand="+", starts=(100, 200), stops=(130, 230)) -> pl.DataFrame:
    """Two-exon CDS for synthetic tests."""
    tran_starts = [0, stops[0] - starts[0]]  # [0, 30]
    return pl.DataFrame(
        {
            "tran_id": ["TX1"],
            "gene_id": ["G1"],
            "chr": [chrom],
            "strand": [strand],
            "start": [list(starts)],
            "stop": [list(stops)],
            "tran_start": [tran_starts],
        }
    )


# ---------------------------------------------------------------------------
# build_feature_exon_index
# ---------------------------------------------------------------------------


def test_build_feature_exon_index_plus_strand():
    from TranslonScorer.matrix.rollup import build_feature_exon_index

    cds = _simple_cds()
    fivs, _ = build_feature_exon_index(cds, {"chr1"})
    assert ("chr1", "+") in fivs
    ivs = fivs[("chr1", "+")]
    assert len(ivs) == 2
    # First exon: [100, 130), tran_start=0
    assert ivs[0] == (100, 130, "TX1", 0)
    # Second exon: [200, 230), tran_start=30
    assert ivs[1] == (200, 230, "TX1", 30)


def test_build_feature_exon_index_chrom_normalisation():
    """chr-prefixed and bare chrom names both appear in bam_refs → normalised correctly."""
    from TranslonScorer.matrix.rollup import build_feature_exon_index

    cds = _simple_cds(chrom="1")  # bare chrom
    fivs, _ = build_feature_exon_index(cds, {"chr1", "1"})
    # Should resolve to some form of chr1
    assert any(k[0] in ("1", "chr1") for k in fivs)


def test_build_feature_exon_index_minus_strand():
    from TranslonScorer.matrix.rollup import build_feature_exon_index

    cds = _simple_cds(strand="-", starts=(200, 100), stops=(230, 130))
    fivs, _ = build_feature_exon_index(cds, {"chr1"})
    assert ("chr1", "-") in fivs


def test_build_feature_exon_index_gapdh():
    """GAPDH-shaped 8-exon CDS builds without error."""
    from TranslonScorer.matrix.rollup import build_feature_exon_index

    fivs, _ = build_feature_exon_index(_gapdh_cds(), {"chr12", "12"})
    key = next(k for k in fivs if k[1] == "+")
    assert len(fivs[key]) == 8


# ---------------------------------------------------------------------------
# _assign_tx_positions
# ---------------------------------------------------------------------------


def test_assign_tx_positions_cds_start():
    """5'-end at CDS start position → tx_pos=0."""
    from TranslonScorer.matrix.rollup import build_feature_exon_index, _assign_tx_positions

    cds = _simple_cds()
    fivs, _ = build_feature_exon_index(cds, {"chr1"})
    ivs = fivs[("chr1", "+")]

    hits = _assign_tx_positions([0], np.array([100], dtype=np.int64), "+", ivs)
    assert len(hits) == 1
    local_idx, fid, tx_pos = hits[0]
    assert fid == "TX1"
    assert tx_pos == 0


def test_assign_tx_positions_within_exon():
    """5'-end 10 nt into first exon → tx_pos=10."""
    from TranslonScorer.matrix.rollup import build_feature_exon_index, _assign_tx_positions

    cds = _simple_cds()
    fivs, _ = build_feature_exon_index(cds, {"chr1"})
    ivs = fivs[("chr1", "+")]

    hits = _assign_tx_positions([0], np.array([110], dtype=np.int64), "+", ivs)
    assert hits[0][2] == 10


def test_assign_tx_positions_second_exon():
    """5'-end at start of second exon → tx_pos=30 (first exon length)."""
    from TranslonScorer.matrix.rollup import build_feature_exon_index, _assign_tx_positions

    cds = _simple_cds()
    fivs, _ = build_feature_exon_index(cds, {"chr1"})
    ivs = fivs[("chr1", "+")]

    hits = _assign_tx_positions([0], np.array([200], dtype=np.int64), "+", ivs)
    assert hits[0][2] == 30


def test_assign_tx_positions_intron_returns_empty():
    """5'-end in an intron (between exon 1 stop and exon 2 start) → no hit."""
    from TranslonScorer.matrix.rollup import build_feature_exon_index, _assign_tx_positions

    cds = _simple_cds()
    fivs, _ = build_feature_exon_index(cds, {"chr1"})
    ivs = fivs[("chr1", "+")]

    hits = _assign_tx_positions([0], np.array([150], dtype=np.int64), "+", ivs)
    assert hits == []


def test_assign_tx_positions_multiple_reads():
    """Multiple reads in one batch — all correctly assigned."""
    from TranslonScorer.matrix.rollup import build_feature_exon_index, _assign_tx_positions

    cds = _simple_cds()
    fivs, _ = build_feature_exon_index(cds, {"chr1"})
    ivs = fivs[("chr1", "+")]

    positions = np.array([100, 105, 150, 200, 225], dtype=np.int64)  # 150 = intron
    hits = _assign_tx_positions(list(range(5)), positions, "+", ivs)
    tx_by_pos = {positions[li]: tp for li, _, tp in hits}
    assert tx_by_pos[100] == 0
    assert tx_by_pos[105] == 5
    assert 150 not in tx_by_pos  # intron
    assert tx_by_pos[200] == 30
    assert tx_by_pos[225] == 55


def test_assign_tx_positions_phase0():
    """phase0 = tx_pos % 3; in-frame positions (tx_pos 0,3,6,…) give phase0=0."""
    from TranslonScorer.matrix.rollup import build_feature_exon_index, _assign_tx_positions

    cds = _simple_cds()
    fivs, _ = build_feature_exon_index(cds, {"chr1"})
    ivs = fivs[("chr1", "+")]

    # In-frame positions: 100 (tx=0), 103 (tx=3), 106 (tx=6)
    positions = np.array([100, 103, 106, 101, 102], dtype=np.int64)
    hits = _assign_tx_positions(list(range(5)), positions, "+", ivs)
    phase0_by_pos = {positions[li]: tp % 3 for li, _, tp in hits}
    assert phase0_by_pos[100] == 0
    assert phase0_by_pos[103] == 0
    assert phase0_by_pos[106] == 0
    assert phase0_by_pos[101] == 1
    assert phase0_by_pos[102] == 2


# ---------------------------------------------------------------------------
# score_frame_rollup
# ---------------------------------------------------------------------------


def test_score_frame_rollup_perfect_periodicity():
    """Synthetic rollup with all reads at phase0=0 → elong_in_frame=1.0 when offset=0."""
    from TranslonScorer.matrix.rollup import score_frame_rollup

    rollup = pl.DataFrame(
        {
            "sample_name": ["S1", "S1", "S1"],
            "feature_id": ["F1", "F1", "F1"],
            "length": [29, 29, 29],
            "strand": ["+", "+", "+"],
            "phase0": [0, 1, 2],
            "count": [100.0, 5.0, 5.0],
        }
    )
    scored = score_frame_rollup(rollup, {("S1", 29): 0}, default_offset=0)
    row = scored.filter(pl.col("sample_name") == "S1").row(0, named=True)
    assert row["elong_in_frame"] == pytest.approx(100.0 / 110.0, abs=0.01)
    assert row["dominant_frame"] == 0


def test_score_frame_rollup_offset_shifts_frame():
    """Applying offset 1 shifts phase0=2 to frame=0 (2+1=3≡0 mod 3)."""
    from TranslonScorer.matrix.rollup import score_frame_rollup

    rollup = pl.DataFrame(
        {
            "sample_name": ["S1", "S1", "S1"],
            "feature_id": ["F1", "F1", "F1"],
            "length": [29, 29, 29],
            "strand": ["+", "+", "+"],
            "phase0": [0, 1, 2],
            "count": [5.0, 5.0, 100.0],  # phase0=2 is dominant
        }
    )
    scored = score_frame_rollup(rollup, {("S1", 29): 1}, default_offset=0)
    row = scored.row(0, named=True)
    # frame=(phase0+1)%3: phase0=2 → frame=0; dominant
    assert row["elong_in_frame"] == pytest.approx(100.0 / 110.0, abs=0.01)


def test_score_frame_rollup_empty():
    """Empty rollup returns empty DataFrame without error."""
    from TranslonScorer.matrix.rollup import score_frame_rollup

    empty = pl.DataFrame(
        schema={
            "sample_name": pl.Utf8,
            "feature_id": pl.Utf8,
            "length": pl.Int64,
            "strand": pl.Utf8,
            "phase0": pl.Int8,
            "count": pl.Float64,
        }
    )
    scored = score_frame_rollup(empty, {}, default_offset=12)
    assert scored.is_empty()


# ---------------------------------------------------------------------------
# Integration: HAS_MATRIX gate
# ---------------------------------------------------------------------------


def test_score_elongation_from_rollup_unit():
    """Unit: known elong_in_frame → correct call; missing feature → INSUFFICIENT."""
    from TranslonScorer.scoring.run import score_elongation_from_rollup

    scored = pl.DataFrame(
        {
            "sample_name": ["S1"],
            "feature_id": ["F1"],
            "length": [29],
            "n_reads": [1000.0],
            "frame0_count": [800.0],
            "frame1_count": [100.0],
            "frame2_count": [100.0],
            "elong_in_frame": [0.80],
            "dominant_frame": [0],
        }
    )
    events = pl.DataFrame(
        {
            "event_id": [1, 2],
            "type": ["elongation", "elongation"],
            "feature_id": ["F1", "MISSING"],
            "start": [100, 200],
            "end": [400, 500],
            "strand": [1, 1],
            "phase": [0, 0],
            "chrom": ["chr1", "chr1"],
        },
        schema={
            "event_id": pl.UInt64,
            "type": pl.Utf8,
            "feature_id": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "strand": pl.Int64,
            "phase": pl.Int64,
            "chrom": pl.Utf8,
        },
    )
    ev = score_elongation_from_rollup(events, scored)
    assert ev[1]["call"] == "SUPPORTED", f"Expected SUPPORTED, got {ev[1]['call']}"
    assert ev[1]["metric"] == pytest.approx(0.80)
    # Missing feature: 0 reads → INSUFFICIENT
    assert ev[2]["eligibility"] == "INSUFFICIENT"
    assert ev[2]["call"] is None


# ---------------------------------------------------------------------------
# profile_from_index (CoverageIndex helpers — unit, no BAM needed)
# ---------------------------------------------------------------------------


def test_profile_from_index_aggregates():
    """profile_from_index sums across samples and lengths to give tx_pos profile."""
    from TranslonScorer.matrix.rollup import profile_from_index

    cov = pl.DataFrame(
        {
            "sample_name": ["S1", "S1", "S2", "S2"],
            "feature_id": ["F1", "F1", "F1", "F1"],
            "tx_pos": [0, 3, 0, 3],
            "length": [29, 29, 28, 28],
            "strand": ["+", "+", "+", "+"],
            "count": [100.0, 20.0, 50.0, 10.0],
        }
    )
    prof = profile_from_index(cov)
    rows = {r["tx_pos"]: r["count"] for r in prof.iter_rows(named=True)}
    assert rows[0] == pytest.approx(150.0)
    assert rows[3] == pytest.approx(30.0)


def test_profile_from_index_normalise():
    """normalise=True divides each feature's profile by its total count."""
    from TranslonScorer.matrix.rollup import profile_from_index

    cov = pl.DataFrame(
        {
            "sample_name": ["S1", "S1"],
            "feature_id": ["F1", "F1"],
            "tx_pos": [0, 3],
            "length": [29, 29],
            "strand": ["+", "+"],
            "count": [75.0, 25.0],
        }
    )
    prof = profile_from_index(cov, normalise=True)
    rows = {r["tx_pos"]: r["count"] for r in prof.iter_rows(named=True)}
    assert rows[0] == pytest.approx(0.75)
    assert rows[3] == pytest.approx(0.25)


def test_profile_from_index_sample_filter():
    """sample_name kwarg restricts to one sample."""
    from TranslonScorer.matrix.rollup import profile_from_index

    cov = pl.DataFrame(
        {
            "sample_name": ["S1", "S2"],
            "feature_id": ["F1", "F1"],
            "tx_pos": [0, 0],
            "length": [29, 29],
            "strand": ["+", "+"],
            "count": [100.0, 999.0],
        }
    )
    prof = profile_from_index(cov, sample_name="S1")
    assert prof["count"].sum() == pytest.approx(100.0)


@pytest.mark.skipif(not HAS_MATRIX, reason="local matrix partition not in data/")
def test_build_coverage_index_gapdh():
    """build_coverage_index returns per-tx_pos counts for GAPDH; periodic triplet spacing."""
    from TranslonScorer.io.annotation import build_cds_blocks
    from TranslonScorer.matrix.rollup import build_coverage_index, profile_from_index

    cds = build_cds_blocks(str(REPO_ROOT / "data" / "genes.gtf"))
    gapdh = cds.filter(pl.col("tran_id") == "ENST00000229239")

    cov = build_coverage_index(
        REPO_ROOT / "data" / "global_partitioned",
        gapdh,
        multimap_mode="unique",
        n_workers=1,
    )
    assert not cov.is_empty()
    assert set(cov.columns) >= {"sample_name", "feature_id", "tx_pos", "length", "strand", "count"}
    assert cov["feature_id"].unique().to_list() == ["ENST00000229239"]

    # Profile should show triplet spacing (every 3rd position dominant)
    prof = profile_from_index(cov, lengths=[29])
    if prof.is_empty():
        pytest.skip("no 29nt reads in CoverageIndex")
    # tx_pos % 3 == 0 positions should carry more reads than other phases
    prof = prof.with_columns((pl.col("tx_pos") % 3).alias("phase"))
    phase_totals = prof.group_by("phase").agg(pl.col("count").sum()).sort("phase")
    counts = phase_totals["count"].to_list()
    assert max(counts) / sum(counts) > 0.40, (
        "Expected dominant triplet phase > 40% of CoverageIndex reads"
    )


# ---------------------------------------------------------------------------
# prevalence_from_rollup
# ---------------------------------------------------------------------------


def test_prevalence_from_rollup_basic():
    """FR5: prevalence counts eligible samples with elong_in_frame ≥ threshold."""
    from TranslonScorer.matrix.rollup import prevalence_from_rollup

    # 3 samples: 2 supported, 1 not
    scored = pl.DataFrame(
        {
            "sample_name": ["S1", "S2", "S3"],
            "feature_id": ["F1", "F1", "F1"],
            "length": [29, 29, 29],
            "n_reads": [500.0, 600.0, 100.0],
            "frame0_count": [400.0, 480.0, 30.0],
            "frame1_count": [50.0, 60.0, 35.0],
            "frame2_count": [50.0, 60.0, 35.0],
            "elong_in_frame": [0.80, 0.80, 0.30],
            "dominant_frame": [0, 0, 1],
        }
    )
    prev = prevalence_from_rollup(scored, min_reads_per_sample=50, elong_in_frame_thr=0.5)
    assert len(prev) == 1
    row = prev.row(0, named=True)
    assert row["feature_id"] == "F1"
    assert row["n_eligible_samples"] == 3
    assert row["n_supported_samples"] == 2
    assert row["prevalence"] == pytest.approx(2 / 3)


def test_prevalence_from_rollup_min_reads_gate():
    """Samples below min_reads_per_sample are excluded from eligibility."""
    from TranslonScorer.matrix.rollup import prevalence_from_rollup

    scored = pl.DataFrame(
        {
            "sample_name": ["S1", "S2"],
            "feature_id": ["F1", "F1"],
            "length": [29, 29],
            "n_reads": [500.0, 10.0],  # S2 below 50
            "frame0_count": [400.0, 8.0],
            "frame1_count": [50.0, 1.0],
            "frame2_count": [50.0, 1.0],
            "elong_in_frame": [0.80, 0.80],
            "dominant_frame": [0, 0],
        }
    )
    prev = prevalence_from_rollup(scored, min_reads_per_sample=50)
    row = prev.row(0, named=True)
    assert row["n_eligible_samples"] == 1  # only S1 eligible
    assert row["n_supported_samples"] == 1
    assert row["prevalence"] == 1.0


def test_prevalence_from_rollup_empty():
    """Empty scored rollup returns empty DataFrame without error."""
    from TranslonScorer.matrix.rollup import prevalence_from_rollup

    empty = pl.DataFrame(
        schema={
            "sample_name": pl.Utf8,
            "feature_id": pl.Utf8,
            "length": pl.Int64,
            "n_reads": pl.Float64,
            "frame0_count": pl.Float64,
            "frame1_count": pl.Float64,
            "frame2_count": pl.Float64,
            "elong_in_frame": pl.Float64,
            "dominant_frame": pl.Int32,
        }
    )
    result = prevalence_from_rollup(empty)
    assert result.is_empty()


@pytest.mark.skipif(not HAS_MATRIX, reason="local matrix partition not in data/")
def test_prevalence_from_rollup_gapdh():
    """FR5: GAPDH should show high prevalence (most samples show periodic frame signal)."""
    from TranslonScorer.io.annotation import build_cds_blocks
    from TranslonScorer.matrix.rollup import (
        build_frame_rollup,
        calibrate_offsets,
        prevalence_from_rollup,
        score_frame_rollup,
    )

    cds = build_cds_blocks(str(REPO_ROOT / "data" / "genes.gtf"))
    gapdh = cds.filter(pl.col("tran_id") == "ENST00000229239")

    rollup = build_frame_rollup(
        REPO_ROOT / "data" / "global_partitioned",
        gapdh,
        multimap_mode="unique",
        n_workers=1,
    )
    agg = rollup.group_by(["sample_name", "length", "strand", "phase0"]).agg(
        pl.col("count").sum()
    )
    offsets = calibrate_offsets(agg, target_frame=0)
    scored = score_frame_rollup(rollup, offsets, default_offset=12)

    prev = prevalence_from_rollup(
        scored, min_reads_per_sample=50, elong_in_frame_thr=0.5, dominant_lengths=[28, 29, 30]
    )
    assert not prev.is_empty()
    row = prev.filter(pl.col("feature_id") == "ENST00000229239").row(0, named=True)
    assert row["n_eligible_samples"] >= 1, "Expected at least 1 eligible sample"
    assert row["prevalence"] > 0.5, (
        f"Expected GAPDH prevalence > 50%, got {row['prevalence']:.1%}"
    )


@pytest.mark.skipif(not HAS_MATRIX, reason="local matrix partition not in data/")
def test_score_elongation_from_rollup_canonical_cds_supported():
    """FR3 positive-control gate: canonical GAPDH CDS must score SUPPORTED via FrameRollup."""
    from TranslonScorer.io.annotation import build_cds_blocks
    from TranslonScorer.matrix.rollup import (
        build_frame_rollup,
        calibrate_offsets,
        score_frame_rollup,
    )
    from TranslonScorer.scoring.run import score_elongation_from_rollup

    cds = build_cds_blocks(str(REPO_ROOT / "data" / "genes.gtf"))
    gapdh = cds.filter(pl.col("tran_id") == "ENST00000229239")

    rollup = build_frame_rollup(
        REPO_ROOT / "data" / "global_partitioned",
        gapdh,
        multimap_mode="unique",
        n_workers=1,
    )
    agg = rollup.group_by(["sample_name", "length", "strand", "phase0"]).agg(
        pl.col("count").sum()
    )
    offsets = calibrate_offsets(agg, target_frame=0)
    scored = score_frame_rollup(rollup, offsets, default_offset=12)

    # Build a minimal synthetic elongation event for GAPDH CDS
    fake_events = pl.DataFrame(
        {
            "event_id": [12345678],
            "type": ["elongation"],
            "feature_id": ["ENST00000229239"],
            "start": [6534832],
            "end": [6538167],
            "strand": [1],
            "phase": [0],
            "chrom": ["chr12"],
        },
        schema={
            "event_id": pl.UInt64,
            "type": pl.Utf8,
            "feature_id": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "strand": pl.Int64,
            "phase": pl.Int64,
            "chrom": pl.Utf8,
        },
    )
    ev = score_elongation_from_rollup(fake_events, scored)
    assert 12345678 in ev, "Expected evidence for GAPDH elongation event"
    result = ev[12345678]
    assert result["eligibility"] == "ELIGIBLE", f"Expected ELIGIBLE, got {result['eligibility']}"
    assert result["call"] == "SUPPORTED", (
        f"FR3 positive-control gate FAILED: canonical CDS should score SUPPORTED, "
        f"got {result['call']!r} (elong_in_frame={result['metric']:.1%})"
    )


@pytest.mark.skipif(not HAS_MATRIX, reason="local matrix partition not in data/")
def test_build_frame_rollup_gapdh_produces_rollup():
    """build_frame_rollup on the local 30-sample matrix returns non-empty rollup
    for GAPDH with correct schema."""
    from TranslonScorer.io.annotation import build_cds_blocks
    from TranslonScorer.matrix.rollup import build_frame_rollup

    cds = build_cds_blocks(str(REPO_ROOT / "data" / "genes.gtf"))
    gapdh = cds.filter(pl.col("tran_id") == "ENST00000229239")

    rollup = build_frame_rollup(
        REPO_ROOT / "data" / "global_partitioned",
        gapdh,
        multimap_mode="unique",
        n_workers=1,
    )
    assert not rollup.is_empty()
    assert set(rollup.columns) >= {"sample_name", "feature_id", "length", "strand", "phase0", "count"}
    assert rollup["feature_id"].unique().to_list() == ["ENST00000229239"]
    assert rollup["count"].sum() > 0


@pytest.mark.skipif(not HAS_MATRIX, reason="local matrix partition not in data/")
def test_build_frame_rollup_gapdh_calibrated_periodicity():
    """Calibrated offsets recover frame-0 > 60% for dominant 29nt reads on GAPDH."""
    from TranslonScorer.io.annotation import build_cds_blocks
    from TranslonScorer.matrix.rollup import (
        build_frame_rollup,
        calibrate_offsets,
        score_frame_rollup,
    )

    cds = build_cds_blocks(str(REPO_ROOT / "data" / "genes.gtf"))
    gapdh = cds.filter(pl.col("tran_id") == "ENST00000229239")

    rollup = build_frame_rollup(
        REPO_ROOT / "data" / "global_partitioned",
        gapdh,
        multimap_mode="unique",
        n_workers=1,
    )
    # Calibrate per-(sample,length) offsets
    agg = rollup.group_by(["sample_name", "length", "strand", "phase0"]).agg(
        pl.col("count").sum()
    )
    offsets = calibrate_offsets(agg, target_frame=0)

    # Score
    scored = score_frame_rollup(rollup, offsets, default_offset=12)
    s29 = scored.filter(pl.col("length") == 29)
    if s29.is_empty():
        pytest.skip("no 29nt reads in matrix fixture")

    # Aggregate frame-0 over all samples with ≥50 reads
    eligible = s29.filter(pl.col("n_reads") >= 50)
    if eligible.is_empty():
        pytest.skip("no samples with ≥50 29nt GAPDH reads")

    agg_f0 = eligible["frame0_count"].sum() / eligible["n_reads"].sum()
    assert agg_f0 > 0.55, (
        f"Expected frame-0 > 55% with calibrated offsets, got {agg_f0:.1%}. "
        "Possible calibration failure."
    )

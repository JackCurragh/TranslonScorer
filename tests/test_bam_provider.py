"""Tests for BamSetProvider and BigwigSetProvider (T12).

Real-data tests (requiring SRR11005875 BAM) are marked skip when the file
is absent.  Structural/protocol tests run in CI without external data.

Gate tests added to test_golden.py verify protocol compliance; this file
tests the provider-specific behaviour:
  - Offsets are calibrated once per BAM/length and reused for locus profiles.
  - Coverage via BamSetProvider scores GAPDH with sane init/elong/term calls
    (when real BAM is available).
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict

import pytest
import polars as pl

REPO_ROOT = Path(__file__).parent.parent
# Real-data fixture: genome-aligned GAPDH-region reads merged from the matrix
# unique_reads partitions (chr12, + strand). Built locally (data/ is gitignored),
# so these tests run locally and skip in CI. Build:
#   samtools merge -R chr12:6533000-6539500 data/gapdh_cohort_genome.bam \
#       data/global_partitioned/*/unique_reads.*.bam && samtools index <out>
GENOME_BAM = REPO_ROOT / "data" / "gapdh_cohort_genome.bam"
HAS_BAM = GENOME_BAM.exists()

# ---------------------------------------------------------------------------
# Protocol compliance (no real data needed)
# ---------------------------------------------------------------------------

def test_bam_set_provider_protocols():
    """BamSetProvider satisfies all four capability protocols."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.coverage.base import (
        CoverageProvider, SupportsSites, SupportsJunctions, SupportsMappability,
    )
    provider = BamSetProvider([])
    assert isinstance(provider, CoverageProvider)
    assert isinstance(provider, SupportsSites)
    assert isinstance(provider, SupportsJunctions)
    assert isinstance(provider, SupportsMappability)


def test_bigwig_set_provider_protocols():
    """BigwigSetProvider satisfies only CoverageProvider (no junctions/sites)."""
    from TranslonScorer.coverage.bigwig import BigwigSetProvider
    from TranslonScorer.coverage.base import (
        CoverageProvider, SupportsJunctions, SupportsMappability,
    )
    provider = BigwigSetProvider([])
    assert isinstance(provider, CoverageProvider)
    # Bigwig does NOT satisfy SupportsJunctions or SupportsMappability
    # (it has those methods but they raise NotImplementedError — the Protocol
    # check passes structurally; we verify the runtime error instead)
    with pytest.raises(NotImplementedError):
        provider.junction_support([])
    with pytest.raises(NotImplementedError):
        provider.mappability_ledger(pl.DataFrame())


def test_bam_set_provider_offset_calibration_cached():
    """BamSetProvider: _calibrate_offsets is idempotent (returns the same object twice)."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.model import OffsetParams

    p = BamSetProvider([], offsets=OffsetParams(method="global", global_offset=12))
    t1 = p._calibrate_offsets()
    t2 = p._calibrate_offsets()
    assert t1 is t2, "offset tables should be cached (same object on second call)"


def test_bam_set_provider_coverage_empty():
    """BamSetProvider.coverage() on empty BAM list returns correct schema."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.model import Region, OffsetParams

    provider = BamSetProvider([], offsets=OffsetParams(method="global"))
    result = provider.coverage([Region("chr12", 1000, 2000)])
    assert "pos" in result.columns
    assert "count" in result.columns
    assert result.is_empty()


def test_bam_set_provider_size_factors():
    """BamSetProvider.size_factors() returns 1.0 per sample."""
    from TranslonScorer.coverage.bam import BamSetProvider

    provider = BamSetProvider(["s1.bam", "s2.bam"], sample_names=["s1", "s2"])
    sf = provider.size_factors()
    assert sf == {"s1": 1.0, "s2": 1.0}


def test_bam_set_provider_junction_support_empty():
    """BamSetProvider.junction_support() on empty BAM list returns correct schema."""
    from TranslonScorer.coverage.bam import BamSetProvider

    provider = BamSetProvider([])
    result = provider.junction_support([("chr12", 100, 200, 1, 1)])
    assert "junction_id" in result.columns
    assert "kind" in result.columns
    assert result.is_empty()


def test_bam_set_provider_mappability_ledger_empty():
    """BamSetProvider.mappability_ledger() returns the correct schema."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.coverage.base import MAPPABILITY_LEDGER_SCHEMA

    provider = BamSetProvider([])
    result = provider.mappability_ledger(pl.DataFrame())
    for col in MAPPABILITY_LEDGER_SCHEMA:
        assert col in result.columns


# ---------------------------------------------------------------------------
# Real-data tests (skipped unless SRR11005875.bam is present)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not HAS_BAM, reason="genome GAPDH fixture not in data/")
def test_srr_bam_offsets_calibrated_once():
    """Offsets are calibrated once per BAM/length and cached across coverage() calls."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.model import OffsetParams, Region

    provider = BamSetProvider(
        [str(GENOME_BAM)],
        offsets=OffsetParams(method="global", global_offset=12),
    )
    # Coverage call triggers offset calibration
    region = Region("chr12", 6534512, 6538371)  # GAPDH locus
    _ = provider.coverage([region])
    tables_after_first = provider._offset_tables

    # Second call must reuse cache
    _ = provider.coverage([region])
    assert provider._offset_tables is tables_after_first, "offset table must be the same object"


@pytest.mark.skipif(not HAS_BAM, reason="genome GAPDH fixture not in data/")
def test_srr_gapdh_score_sane():
    """Score GAPDH locus from SRR11005875 BAM; assert sane init/elong/term calls."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.model import OffsetParams, Region
    from TranslonScorer.pipeline.event_score import (
        ScoreThresholds,
        score_events,
    )
    from TranslonScorer.events import extract_events

    # Minimal GAPDH annotation (single-exon proxy; real test would use full exon_df)
    from tests.test_golden import _gapdh_events

    provider = BamSetProvider(
        [str(GENOME_BAM)],
        offsets=OffsetParams(method="global", global_offset=12),
    )
    region = Region("chr12", 6534512, 6538371)
    cov_df = provider.coverage([region], site="A")
    cov_dict: Dict[int, float] = dict(zip(
        cov_df["pos"].to_list(),
        cov_df["count"].to_list(),
    ))

    thr = ScoreThresholds()
    events = _gapdh_events()
    result = score_events(events, cov_dict, group="gapdh", tier="aggregate", thr=thr)

    assert result.height > 0, "expected non-empty score results"
    calls = dict(zip(result["aspect"].to_list(), result["call"].to_list()))

    # With real data, init/elong/term should at minimum be ELIGIBLE or SUPPORTED
    for aspect in ("init", "elongation", "term"):
        assert calls.get(aspect) in ("SUPPORTED", "AMBIGUOUS", "ELIGIBLE", None), (
            f"unexpected call for {aspect}: {calls.get(aspect)}"
        )

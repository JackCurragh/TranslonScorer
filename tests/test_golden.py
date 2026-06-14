"""Golden regression gate for the TranslonScorer refactor.

Gate command:  make gate   (= pytest tests/test_golden.py -q)

Assertions
----------
1. GAPDH golden: score_events on the synthetic GAPDH-locus fixture reproduces
   outputs/reference_scores_gapdh.parquet with max|Δmetric|==0 and 0
   call/eligibility mismatches.
2. Scalar ≡ vectorised: score_events and score_events_vectorised agree (Δ=0)
   on the same fixture, including a contended elongation pair.

To regenerate the golden (only needed when intentionally changing the scorer):
    python3 tests/test_golden.py
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from TranslonScorer.events import extract_events
from TranslonScorer.pipeline.event_score import (
    ScoreThresholds,
    score_events,
    score_events_vectorised,
)

REPO_ROOT = Path(__file__).parent.parent
GOLDEN_PATH = REPO_ROOT / "outputs" / "reference_scores_gapdh.parquet"

_THR = ScoreThresholds()

# ---------------------------------------------------------------------------
# Synthetic GAPDH-locus fixture
# ---------------------------------------------------------------------------
# Single-exon translon on + strand.  Positions chosen so that frame-0
# in-frame signal falls at p % 3 == 0.
#
# Events
#   init   : pos 99   (just before CDS body)
#   elong  : [100, 190)  phase 0  a_e = 0  (in-frame at p%3==0)
#   term   : pos 190  (just after CDS body)
#
# Coverage
#   outer flank [40, 99) : 2 reads/nt  (low – scorer sees step UP into body)
#   CDS body [99, 191)   : 30 at p%3==0, 5 elsewhere (strong frame-0 signal)
#   UTR after [191, 260) : 1 read/nt  (low – scorer sees step DOWN out of body)
#
# Expected calls  init→SUPPORTED, elong→SUPPORTED, term→SUPPORTED

_GAPDH_SCHEMA = {
    "event_id": pl.UInt64,
    "type": pl.Utf8,
    "start": pl.Int64,
    "end": pl.Int64,
    "strand": pl.Int64,
    "phase": pl.Int64,
}


def _gapdh_events() -> pl.DataFrame:
    return pl.DataFrame(
        [
            {"event_id": 1, "type": "init",      "start": 99,  "end": 100, "strand": 1, "phase": None},
            {"event_id": 2, "type": "elongation", "start": 100, "end": 190, "strand": 1, "phase": 0},
            {"event_id": 3, "type": "term",       "start": 190, "end": 191, "strand": 1, "phase": None},
        ],
        schema=_GAPDH_SCHEMA,
    )


def _gapdh_coverage() -> dict:
    cov: dict = {}
    for p in range(40, 99):    # outer flank for init
        cov[p] = 2.0
    for p in range(99, 191):   # CDS body + init pos + term pos
        cov[p] = 30.0 if p % 3 == 0 else 5.0
    for p in range(191, 260):  # UTR after term
        cov[p] = 1.0
    return cov


def _run_gapdh() -> pl.DataFrame:
    return score_events(
        _gapdh_events(),
        _gapdh_coverage(),
        group="gapdh",
        tier="aggregate",
        thr=_THR,
    )


# ---------------------------------------------------------------------------
# Contention fixture for scalar≡vectorised comparison
# ---------------------------------------------------------------------------
# Two elongation events that overlap in different frames:
#   elong A: [200, 290)  phase 0  a_e=0  (in-frame at p%3==0)
#   elong B: [250, 330)  phase 1  a_e=2  (in-frame at p%3==2)
# Overlap: [250, 290).  Coverage has mixed frame signal in the overlap.

_CONTEND_SCHEMA = {
    "event_id": pl.UInt64,
    "type": pl.Utf8,
    "start": pl.Int64,
    "end": pl.Int64,
    "strand": pl.Int64,
    "phase": pl.Int64,
}


def _contend_events() -> pl.DataFrame:
    return pl.DataFrame(
        [
            {"event_id": 10, "type": "elongation", "start": 200, "end": 290, "strand": 1, "phase": 0},
            {"event_id": 11, "type": "elongation", "start": 250, "end": 330, "strand": 1, "phase": 1},
        ],
        schema=_CONTEND_SCHEMA,
    )


def _contend_coverage() -> dict:
    cov: dict = {}
    for p in range(200, 330):
        # Mix of frame-0 and frame-2 signal
        if p % 3 == 0:
            cov[p] = 20.0
        elif p % 3 == 1:
            cov[p] = 15.0
        else:
            cov[p] = 3.0
    return cov


def _contend_overlaps_dict() -> dict:
    return {
        10: [(250, 290, 1, 11)],   # A: overlap with B (comp_phase=1)
        11: [(250, 290, 0, 10)],   # B: overlap with A (comp_phase=0)
    }


def _contend_overlaps_df() -> pl.DataFrame:
    return pl.DataFrame(
        [
            {"event_id": 10, "other_event_id": 11, "overlap_start": 250, "overlap_end": 290, "comp_phase": 1},
            {"event_id": 11, "other_event_id": 10, "overlap_start": 250, "overlap_end": 290, "comp_phase": 0},
        ],
        schema={
            "event_id": pl.UInt64,
            "other_event_id": pl.UInt64,
            "overlap_start": pl.Int64,
            "overlap_end": pl.Int64,
            "comp_phase": pl.Int64,
        },
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _assert_identical(a: pl.DataFrame, b: pl.DataFrame, label: str) -> None:
    key = ["event_id", "aspect"]
    a = a.sort(key)
    b = b.sort(key)
    assert a.shape == b.shape, f"[{label}] shape mismatch: {a.shape} vs {b.shape}"

    m_a = a["metric"].fill_null(0.0)
    m_b = b["metric"].fill_null(0.0)
    delta = (m_a - m_b).abs().max()
    assert delta == 0.0, f"[{label}] max|Δmetric|={delta}"

    call_mm = (a["call"] != b["call"]).sum()
    elig_mm = (a["eligibility"] != b["eligibility"]).sum()
    assert call_mm == 0, f"[{label}] {call_mm} call mismatch(es)"
    assert elig_mm == 0, f"[{label}] {elig_mm} eligibility mismatch(es)"


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_gapdh_golden():
    """Score synthetic GAPDH fixture; must reproduce committed golden exactly."""
    assert GOLDEN_PATH.exists(), (
        f"Golden not found: {GOLDEN_PATH}\n"
        "Run  python3 tests/test_golden.py  to create it."
    )
    result = _run_gapdh()
    golden = pl.read_parquet(GOLDEN_PATH)
    _assert_identical(result, golden, "gapdh_golden")


def test_scalar_equals_vectorised_gapdh():
    """score_events and score_events_vectorised agree on the GAPDH fixture."""
    events = _gapdh_events()
    cov_dict = _gapdh_coverage()
    cov_df = pl.DataFrame(
        {"pos": list(cov_dict.keys()), "count": list(cov_dict.values())},
        schema={"pos": pl.Int64, "count": pl.Float64},
    )
    scalar = score_events(events, cov_dict, group="gapdh", tier="aggregate", thr=_THR)
    vec = score_events_vectorised(events, cov_df, group="gapdh", tier="aggregate", thr=_THR)
    _assert_identical(scalar, vec, "scalar_vs_vec_gapdh")


def test_scalar_equals_vectorised_contended():
    """score_events and score_events_vectorised agree on a contended-elong fixture."""
    events = _contend_events()
    cov_dict = _contend_coverage()
    cov_df = pl.DataFrame(
        {"pos": list(cov_dict.keys()), "count": list(cov_dict.values())},
        schema={"pos": pl.Int64, "count": pl.Float64},
    )
    overlaps_dict = _contend_overlaps_dict()
    overlaps_df = _contend_overlaps_df()

    scalar = score_events(
        events, cov_dict, overlaps=overlaps_dict, group="test", tier="aggregate", thr=_THR
    )
    vec = score_events_vectorised(
        events, cov_df, overlaps_df=overlaps_df, group="test", tier="aggregate", thr=_THR
    )
    _assert_identical(scalar, vec, "scalar_vs_vec_contended")


# ---------------------------------------------------------------------------
# T6 — chr22 event-count reproduction (synthetic)
# ---------------------------------------------------------------------------
# Two non-overlapping single-exon translons on chr22, + strand.
#   T1: [100, 200)   → 1 elongation, 1 init, 1 term
#   T2: [300, 400)   → 1 elongation, 1 init, 1 term
# Distinct events because neither position nor phase overlaps.
# Expected totals: 2 elongation, 2 init, 2 term = 6 events total, 0 junctions.

_BLOCKS_SCHEMA = {
    "translon_id": pl.Utf8,
    "translation_block_rank": pl.Int64,
    "bed_chrom": pl.Utf8,
    "bed_start": pl.Int64,
    "bed_end": pl.Int64,
    "seq_region_strand": pl.Int64,
    "block_length_nt": pl.Int64,
}
_TRANSLONS_SCHEMA = {
    "translon_id": pl.Utf8,
    "bed_chrom": pl.Utf8,
    "bed_start": pl.Int64,
    "bed_end": pl.Int64,
    "seq_region_strand": pl.Int64,
}


def _chr22_blocks() -> pl.DataFrame:
    return pl.DataFrame(
        [
            {"translon_id": "t1", "translation_block_rank": 1,
             "bed_chrom": "chr22", "bed_start": 100, "bed_end": 200,
             "seq_region_strand": 1, "block_length_nt": 100},
            {"translon_id": "t2", "translation_block_rank": 1,
             "bed_chrom": "chr22", "bed_start": 300, "bed_end": 400,
             "seq_region_strand": 1, "block_length_nt": 100},
        ],
        schema=_BLOCKS_SCHEMA,
    )


def _chr22_translons() -> pl.DataFrame:
    return pl.DataFrame(
        [
            {"translon_id": "t1", "bed_chrom": "chr22",
             "bed_start": 100, "bed_end": 200, "seq_region_strand": 1},
            {"translon_id": "t2", "bed_chrom": "chr22",
             "bed_start": 300, "bed_end": 400, "seq_region_strand": 1},
        ],
        schema=_TRANSLONS_SCHEMA,
    )


def test_chr22_event_counts():
    """Synthetic chr22 extract_events: 2 translons → 6 events, 0 junctions."""
    events, fe, overlap = extract_events(_chr22_blocks(), _chr22_translons())
    by_type = dict(events.group_by("type").len().iter_rows())
    assert by_type.get("elongation", 0) == 2, f"expected 2 elongation, got {by_type}"
    assert by_type.get("init", 0) == 2, f"expected 2 init, got {by_type}"
    assert by_type.get("term", 0) == 2, f"expected 2 term, got {by_type}"
    assert by_type.get("junction", 0) == 0, f"expected 0 junction, got {by_type}"
    assert events.height == 6, f"expected 6 total events, got {events.height}"
    assert fe.height == 6, f"expected 6 feature_event rows, got {fe.height}"
    assert overlap.is_empty(), "expected no contention (events don't overlap)"


# ---------------------------------------------------------------------------
# Golden creation  (python3 tests/test_golden.py)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    GOLDEN_PATH.parent.mkdir(parents=True, exist_ok=True)
    result = _run_gapdh()
    result.write_parquet(GOLDEN_PATH)
    print(f"Written {result.height} rows → {GOLDEN_PATH}")
    for r in result.iter_rows(named=True):
        print(f"  event={r['event_id']} type={r['aspect']} "
              f"elig={r['eligibility']} call={r['call']} metric={r['metric']:.4f}")

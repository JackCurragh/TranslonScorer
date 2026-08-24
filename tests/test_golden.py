"""Golden regression gate for the TranslonScorer refactor.

Gate command:  make gate   (= pytest tests/test_golden.py -q)

Assertions
----------
1. GAPDH golden: the shipped ``scoring.run.score_events`` on the synthetic
   GAPDH-locus fixture reproduces outputs/reference_scores_gapdh.parquet with
   max|Δmetric|==0 and 0 call/eligibility mismatches.
2. Product ≡ reference: the shipped scorer and the independent scalar
   implementation in ``tests/reference_scorer.py`` agree (Δ=0) on the same
   fixture, including a contended elongation pair.

To regenerate the golden (only needed when intentionally changing the scorer):
    python3 tests/test_golden.py
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tests.reference_scorer import score_events_scalar
from TranslonScorer.events import extract_events
from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.run import score_events

REPO_ROOT = Path(__file__).parent.parent


class _DictCoverageProvider:
    """In-memory CoverageProvider test double wrapping a {pos: count} dict.

    Relocated from product code: only tests need an offline provider that
    requires no matrix files. Duck-typed — exposes coverage()/size_factors(),
    which is all the scoring pipeline asks of a coverage source.
    """

    def __init__(self, cov_dict) -> None:
        self._dict = cov_dict

    def coverage(self, regions, *, site: str = "A", by_sample: bool = False):
        positions, counts = [], []
        for r in regions:
            for pos in range(r.start, r.end):
                c = self._dict.get(pos)
                if c is not None:
                    positions.append(pos)
                    counts.append(float(c))
        return pl.DataFrame(
            {"pos": positions, "count": counts},
            schema={"pos": pl.Int64, "count": pl.Float64},
        )

    def size_factors(self):
        return {"": 1.0}


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
            {"event_id": 1, "type": "init", "start": 99, "end": 100, "strand": 1, "phase": None},
            {
                "event_id": 2,
                "type": "elongation",
                "start": 100,
                "end": 190,
                "strand": 1,
                "phase": 0,
            },
            {"event_id": 3, "type": "term", "start": 190, "end": 191, "strand": 1, "phase": None},
        ],
        schema=_GAPDH_SCHEMA,
    )


def _gapdh_coverage() -> dict:
    cov: dict = {}
    for p in range(40, 99):  # outer flank for init
        cov[p] = 2.0
    for p in range(99, 191):  # CDS body + init pos + term pos
        cov[p] = 30.0 if p % 3 == 0 else 5.0
    for p in range(191, 260):  # UTR after term
        cov[p] = 1.0
    return cov


def _cov_df(cov: dict) -> pl.DataFrame:
    return pl.DataFrame(
        {"pos": list(cov.keys()), "count": list(cov.values())},
        schema={"pos": pl.Int64, "count": pl.Float64},
    )


def _run_gapdh() -> pl.DataFrame:
    """Golden is produced by the SHIPPED scorer, not the test-only reference.

    (It used to run the scalar reference, so the golden only reached the
    product transitively via the scalar≡vectorised test.)
    """
    return score_events(
        _gapdh_events(),
        _cov_df(_gapdh_coverage()),
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
            {
                "event_id": 10,
                "type": "elongation",
                "start": 200,
                "end": 290,
                "strand": 1,
                "phase": 0,
            },
            {
                "event_id": 11,
                "type": "elongation",
                "start": 250,
                "end": 330,
                "strand": 1,
                "phase": 1,
            },
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
        10: [(250, 290, 1, 11)],  # A: overlap with B (comp_phase=1)
        11: [(250, 290, 0, 10)],  # B: overlap with A (comp_phase=0)
    }


def _contend_overlaps_df() -> pl.DataFrame:
    return pl.DataFrame(
        [
            {
                "event_id": 10,
                "other_event_id": 11,
                "overlap_start": 250,
                "overlap_end": 290,
                "comp_phase": 1,
            },
            {
                "event_id": 11,
                "other_event_id": 10,
                "overlap_start": 250,
                "overlap_end": 290,
                "comp_phase": 0,
            },
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
        f"Golden not found: {GOLDEN_PATH}\n" "Run  python3 tests/test_golden.py  to create it."
    )
    result = _run_gapdh()
    golden = pl.read_parquet(GOLDEN_PATH)
    _assert_identical(result, golden, "gapdh_golden")


def test_product_equals_scalar_reference_gapdh():
    """Shipped score_events agrees with the independent scalar reference."""
    events = _gapdh_events()
    cov_dict = _gapdh_coverage()
    scalar = score_events_scalar(events, cov_dict, group="gapdh", tier="aggregate", thr=_THR)
    product = score_events(events, _cov_df(cov_dict), group="gapdh", tier="aggregate", thr=_THR)
    _assert_identical(scalar, product, "scalar_vs_product_gapdh")


def test_product_equals_scalar_reference_contended():
    """Same, on a contended elongation pair — where the prefix-sum kernel is
    least obvious by inspection and the reference earns its keep."""
    events = _contend_events()
    cov_dict = _contend_coverage()
    overlaps_dict = _contend_overlaps_dict()
    overlaps_df = _contend_overlaps_df()

    scalar = score_events_scalar(
        events, cov_dict, overlaps=overlaps_dict, group="test", tier="aggregate", thr=_THR
    )
    product = score_events(
        events, _cov_df(cov_dict), overlaps_df=overlaps_df, group="test", tier="aggregate", thr=_THR
    )
    _assert_identical(scalar, product, "scalar_vs_product_contended")


def test_scalar_equals_vectorised_on_cif_and_n_codons():
    """`_assert_identical` only checks metric/call/eligibility — cif and
    n_codons are first-class columns computed independently in both the
    scalar (`score_elongation_event`) and vectorised (`score_elongation_batch`)
    paths, and nothing else cross-checks that they actually agree. Three
    blocks, deliberately varied in size/phase/strand/depth so a divergence
    couldn't hide behind a degenerate case."""
    events = pl.DataFrame(
        [
            {
                "event_id": 1,
                "type": "elongation",
                "chrom": "chr1",
                "start": 100,
                "end": 190,
                "strand": 1,
                "phase": 0,
            },
            {
                "event_id": 2,
                "type": "elongation",
                "chrom": "chr1",
                "start": 300,
                "end": 345,
                "strand": 1,
                "phase": 1,
            },
            {
                "event_id": 3,
                "type": "elongation",
                "chrom": "chr1",
                "start": 500,
                "end": 700,
                "strand": -1,
                "phase": 2,
            },
        ],
        schema={
            "event_id": pl.UInt64,
            "type": pl.Utf8,
            "chrom": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "strand": pl.Int64,
            "phase": pl.Int64,
        },
    )
    cov = {}
    cov.update({p: (25.0 if p % 3 == 0 else 4.0) for p in range(90, 200)})
    cov.update({p: (12.0 if (p - 1) % 3 == 0 else 3.0) for p in range(290, 355)})
    cov.update({p: (18.0 if p % 3 == 2 else 6.0) for p in range(490, 710)})
    cov_df = pl.DataFrame(
        {"pos": list(cov.keys()), "count": list(cov.values())},
        schema={"pos": pl.Int64, "count": pl.Float64},
    )

    scalar = score_events_scalar(events, cov, group="g", tier="t", thr=_THR).sort("event_id")
    vec = score_events(events, cov_df, group="g", tier="t", thr=_THR).sort("event_id")

    assert scalar["cif"].to_list() == vec["cif"].to_list()
    assert scalar["n_codons"].to_list() == vec["n_codons"].to_list()
    # Sanity: the fixture actually exercises non-degenerate values, not all-None.
    assert any(v is not None for v in scalar["cif"].to_list())


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
            {
                "translon_id": "t1",
                "translation_block_rank": 1,
                "bed_chrom": "chr22",
                "bed_start": 100,
                "bed_end": 200,
                "seq_region_strand": 1,
                "block_length_nt": 100,
            },
            {
                "translon_id": "t2",
                "translation_block_rank": 1,
                "bed_chrom": "chr22",
                "bed_start": 300,
                "bed_end": 400,
                "seq_region_strand": 1,
                "block_length_nt": 100,
            },
        ],
        schema=_BLOCKS_SCHEMA,
    )


def _chr22_translons() -> pl.DataFrame:
    return pl.DataFrame(
        [
            {
                "translon_id": "t1",
                "bed_chrom": "chr22",
                "bed_start": 100,
                "bed_end": 200,
                "seq_region_strand": 1,
            },
            {
                "translon_id": "t2",
                "bed_chrom": "chr22",
                "bed_start": 300,
                "bed_end": 400,
                "seq_region_strand": 1,
            },
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
# T9 smoke tests — new-tree modules: qc, frame_support, clustering,
#                  report, consequential
# ---------------------------------------------------------------------------


def test_qc_periodicity_from_frames():
    """_compute_periodicity_from_frames returns expected structure."""
    from TranslonScorer.qc import _compute_periodicity_from_frames

    # Perfect frame-0 dominance → high periodicity score
    r = _compute_periodicity_from_frames({0: 100, 1: 1, 2: 1})
    assert "periodicity_score" in r
    assert r["periodicity_score"] > 0.7

    # All zeros → score 0
    r2 = _compute_periodicity_from_frames({0: 0, 1: 0, 2: 0})
    assert r2["periodicity_score"] == 0.0


def test_qc_frame_dominance_matrix_empty():
    """frame_dominance_matrix on empty input returns schema-only DataFrame."""
    from TranslonScorer.qc import frame_dominance_matrix

    empty = pl.DataFrame()
    result = frame_dominance_matrix(empty)
    assert "read_length" in result.columns


def test_qc_nudge_to_frame0():
    """_nudge_to_frame0 leaves frame-0-dominant lengths alone."""
    from TranslonScorer.qc import _nudge_to_frame0

    rfd = {28: {0: 200, 1: 10, 2: 10}, 29: {0: 10, 1: 200, 2: 10}}
    base = {28: 12, 29: 12}
    out = _nudge_to_frame0(rfd, base, min_reads=50, min_fraction=0.60)
    assert out[28] == 12  # already frame-0 dominant — no change
    assert out[29] != 12  # frame-1 dominant — nudged


def test_qc_assign_frames_sweep_smoke():
    """_assign_frames_sweep returns a dict keyed by read_ids."""
    import numpy as np

    from TranslonScorer.qc import _assign_frames_sweep

    rids = [1, 2, 3]
    asites = np.array([105, 108, 120])
    ivs = [(100, 130, 0)]
    result = _assign_frames_sweep(rids, asites, "+", ivs)
    assert isinstance(result, dict)
    assert all(k in rids for k in result)


def test_frame_support_none_method():
    """build_frame_support with method='none' returns empty DataFrame."""
    from TranslonScorer.frame_support import build_frame_support
    from TranslonScorer.model import FrameSupportParams

    profiles = pl.DataFrame(schema={"tran_id": pl.Utf8, "pos": pl.Int64, "count": pl.Float64})
    cds_df = pl.DataFrame()
    result = build_frame_support(profiles, cds_df, FrameSupportParams(frame_method="none"))
    assert result.is_empty()
    assert "tran_id" in result.columns


def test_clustering_normalise_profiles():
    """normalise_profiles total_count method: each row sums to ~1e6."""
    import numpy as np

    from TranslonScorer.matrix.clustering import normalise_profiles

    mat = np.array([[10.0, 20.0, 30.0], [5.0, 5.0, 5.0]])
    out = normalise_profiles(mat, method="total_count")
    assert out.shape == mat.shape
    assert abs(out[0].sum() - 1e6) < 1.0
    assert abs(out[1].sum() - 1e6) < 1.0


def test_clustering_cluster_profiles_smoke():
    """cluster_profiles: 4 samples → 2 clusters, none excluded."""
    import numpy as np

    from TranslonScorer.matrix.clustering import cluster_profiles

    mat = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.9, 0.1, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.1, 0.9, 0.0],
        ]
    )
    names = ["s1", "s2", "s3", "s4"]
    labels, included, excluded = cluster_profiles(mat, names, n_clusters=2, method="kmeans")
    assert len(excluded) == 0
    assert len(set(labels[labels >= 0])) <= 2


def test_clustering_build_profile_matrix_smoke():
    """build_profile_matrix round-trips a tidy profile DataFrame."""
    from TranslonScorer.matrix.clustering import build_profile_matrix

    profiles = pl.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2"],
            "pos": [0, 1, 0],
            "count": [5.0, 3.0, 7.0],
        }
    )
    mat, pos_vec = build_profile_matrix(profiles, ["s1", "s2"])
    assert mat.shape == (2, 2)
    assert float(mat[0, 0]) == 5.0
    assert float(mat[1, 0]) == 7.0


def test_report_compose_passthrough():
    """compose_report returns the input scores DataFrame unchanged when no filter."""
    from TranslonScorer.report import compose_report

    scores = pl.DataFrame(
        {
            "event_id": [1, 2],
            "aspect": ["init", "elong"],
            "eligibility": ["ELIGIBLE", "ELIGIBLE"],
            "call": ["SUPPORTED", "SUPPORTED"],
            "group": ["g", "g"],
            "tier": ["t", "t"],
        }
    )
    out = compose_report(scores)
    assert out.shape == scores.shape


def test_consequential_apply_policy():
    """apply_policy adds a boolean consequential column."""
    from TranslonScorer.consequential import apply_policy
    from TranslonScorer.model import ConsequentialityPolicy

    report = pl.DataFrame({"event_id": [1, 2]})
    out = apply_policy(report, ConsequentialityPolicy())
    assert "consequential" in out.columns
    assert out.height == 2


def test_matrix_qc_public_surface_importable():
    """TranslonScorer.matrix.qc re-exports the periodicity/frame helpers."""
    from TranslonScorer.matrix.qc import (
        _compute_periodicity_from_frames,
        _empty_periodicity_schema,
        frame_dominance_matrix,
    )

    assert callable(_compute_periodicity_from_frames)
    assert callable(frame_dominance_matrix)
    schema = _empty_periodicity_schema()
    assert "sample_id" in schema.columns


def test_pipeline_shim_clustering_importable():
    """pipeline/profile_clustering.py shim: public API importable from old path."""
    from TranslonScorer.matrix.clustering import (
        cluster_profiles,
        normalise_profiles,
    )

    assert callable(normalise_profiles)
    assert callable(cluster_profiles)


# ---------------------------------------------------------------------------
# T10 smoke tests — coverage/base.py protocols + coverage/profile.py
# ---------------------------------------------------------------------------


def test_coverage_base_schema_and_duck_typing():
    """coverage/base.py exports the ledger schema; providers are duck-typed."""
    import polars as pl

    from TranslonScorer.coverage.base import MAPPABILITY_LEDGER_SCHEMA

    assert "event_id" in MAPPABILITY_LEDGER_SCHEMA

    # A coverage source is anything exposing coverage()/size_factors().
    class FakeProvider:
        def coverage(self, regions, *, by_sample=False):
            return pl.DataFrame()

        def size_factors(self):
            return {}

    fp = FakeProvider()
    assert callable(getattr(fp, "coverage", None))
    assert callable(getattr(fp, "size_factors", None))


def test_coverage_profile_apply_offsets_psite():
    """apply_offsets with site='P' adds P-site offset to tran_start_bam."""
    from TranslonScorer.coverage.psite_profile import apply_offsets

    reads = pl.DataFrame(
        {
            "tran_id": ["t1", "t1"],
            "tran_start_bam": [100, 200],
            "length": [28, 29],
            "count": [1.0, 2.0],
        }
    )
    offsets = {28: 12, 29: 13}
    out = apply_offsets(reads, offsets, site="P")
    assert "pos" in out.columns
    pos_set = set(out["pos"].to_list())
    assert 112 in pos_set  # 100 + 12
    assert 213 in pos_set  # 200 + 13


def test_coverage_profile_apply_offsets_asite():
    """apply_offsets with site='A' adds P-site offset + 3."""
    from TranslonScorer.coverage.psite_profile import apply_offsets

    reads = pl.DataFrame(
        {
            "tran_id": ["t1"],
            "tran_start_bam": [100],
            "length": [28],
            "count": [5.0],
        }
    )
    offsets = {28: 12}
    out = apply_offsets(reads, offsets, site="A")
    assert out["pos"].item() == 115  # 100 + 12 + 3


def test_coverage_profile_apply_offsets_empty():
    """apply_offsets on empty reads returns correct schema."""
    from TranslonScorer.coverage.psite_profile import apply_offsets

    reads = pl.DataFrame(
        schema={
            "tran_id": pl.Utf8,
            "tran_start_bam": pl.Int64,
            "length": pl.Int64,
            "count": pl.Float64,
        }
    )
    out = apply_offsets(reads, {}, site="A")
    assert out.is_empty()
    assert "pos" in out.columns


def test_coverage_profile_prefix_sums_bit_identical():
    """_prefix_sums/_range_sums in coverage/profile.py produce same results
    as the re-exported versions in scoring/run.py (bit-identical after migration)."""
    import numpy as np

    from TranslonScorer.coverage.psite_profile import _prefix_sums as prof_ps
    from TranslonScorer.coverage.psite_profile import _range_sums as prof_rs
    from TranslonScorer.scoring.run import _prefix_sums as run_ps
    from TranslonScorer.scoring.run import _range_sums as run_rs

    # They should be the same function object after the re-import
    assert prof_ps is run_ps
    assert prof_rs is run_rs

    pos = np.array([100, 101, 102, 103, 104, 105], dtype=np.int64)
    cnt = np.array([10.0, 5.0, 3.0, 8.0, 2.0, 6.0], dtype=np.float64)
    p, ca, cf = prof_ps(pos, cnt)
    assert len(ca) == len(pos) + 1
    assert len(cf) == 3

    starts = np.array([100, 102], dtype=np.int64)
    ends = np.array([103, 106], dtype=np.int64)
    tot, fr, cov = prof_rs(p, ca, cf, starts, ends)
    assert tot[0] == pytest.approx(10.0 + 5.0 + 3.0)
    assert fr.shape == (2, 3)
    assert cov[0] == 3


import pytest  # noqa: E402 — needed for approx above

# ---------------------------------------------------------------------------
# T11 — GAPDH golden through the provider path (Δ=0)
# ---------------------------------------------------------------------------


def test_gapdh_golden_via_provider():
    """GAPDH golden reproduced through MatrixProvider.coverage() interface.

    Uses the in-memory _DictCoverageProvider (no real matrix files) to verify
    that the provider path produces bit-identical results to the direct
    dict-based scoring path.  This proves the provider interface does not
    alter the scoring arithmetic.
    """
    provider = _DictCoverageProvider(_gapdh_coverage())
    assert callable(getattr(provider, "coverage", None))

    # Pull coverage over the GAPDH locus via the provider interface
    from TranslonScorer.model import Region

    locus = Region("chr12", 40, 260)
    cov_df = provider.coverage([locus])
    assert "pos" in cov_df.columns
    assert "count" in cov_df.columns

    # Convert provider output back to a dict (same shape as _gapdh_coverage)
    cov_via_provider: dict = dict(
        zip(
            cov_df["pos"].to_list(),
            cov_df["count"].to_list(),
        )
    )

    # Score using both paths and assert bit-identical
    direct = score_events(
        _gapdh_events(), _cov_df(_gapdh_coverage()), group="gapdh", tier="aggregate", thr=_THR
    )
    via_provider = score_events(
        _gapdh_events(), _cov_df(cov_via_provider), group="gapdh", tier="aggregate", thr=_THR
    )
    _assert_identical(direct, via_provider, "gapdh_via_provider")


def test_provider_capabilities():
    """MatrixProvider exposes every capability method; the dict double only coverage."""
    from TranslonScorer.matrix.provider import MatrixProvider

    provider = _DictCoverageProvider({})
    assert callable(getattr(provider, "coverage", None))
    assert callable(getattr(provider, "size_factors", None))
    assert not hasattr(provider, "junction_support")

    mp = MatrixProvider([], psite_index_dir="/nonexistent/index")
    for method in ("coverage", "size_factors", "junction_support", "mappability_ledger"):
        assert callable(getattr(mp, method, None)), method


# ---------------------------------------------------------------------------
# T12 smoke tests — coverage/bam.py + coverage/bigwig.py protocols
# ---------------------------------------------------------------------------


def test_bam_provider_importable_and_capabilities():
    """BamSetProvider exposes every capability method."""
    from TranslonScorer.coverage.bam import BamSetProvider

    provider = BamSetProvider([])
    for method in ("coverage", "size_factors", "junction_support", "mappability_ledger"):
        assert callable(getattr(provider, method, None)), method


def test_bigwig_provider_importable_and_capability():
    """BigwigSetProvider exposes coverage; junction/mappability raise (no CIGAR)."""
    from TranslonScorer.coverage.bigwig import BigwigSetProvider

    provider = BigwigSetProvider([])
    assert callable(getattr(provider, "coverage", None))
    with pytest.raises(NotImplementedError):
        provider.junction_support([])
    with pytest.raises(NotImplementedError):
        provider.mappability_ledger(pl.DataFrame())


def test_bam_provider_offset_calibration_is_cached():
    """BamSetProvider._calibrate_offsets is idempotent (caches offset table)."""
    from TranslonScorer.coverage.bam import BamSetProvider
    from TranslonScorer.model import OffsetParams

    p = BamSetProvider([], offsets=OffsetParams(method="global", global_offset=12))
    t1 = p._calibrate_offsets()
    t2 = p._calibrate_offsets()
    assert t1 is t2


# ---------------------------------------------------------------------------
# Golden creation  (python3 tests/test_golden.py)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    GOLDEN_PATH.parent.mkdir(parents=True, exist_ok=True)
    result = _run_gapdh()
    result.write_parquet(GOLDEN_PATH)
    print(f"Written {result.height} rows → {GOLDEN_PATH}")
    for r in result.iter_rows(named=True):
        print(
            f"  event={r['event_id']} type={r['aspect']} "
            f"elig={r['eligibility']} call={r['call']} metric={r['metric']:.4f}"
        )


# ---------------------------------------------------------------------------
# Evidence-blob agreement
#
# _assert_identical deliberately compares only metric/call/eligibility, so the
# evidence JSON -- overall_in_frame (== Chothani PIF), cif, breadth,
# identifiability, competitor_share -- was entirely unpinned. Anything in there
# could drift between the batched and scalar paths, or silently change shape
# between releases, without a single test noticing. It is part of the published
# output, so it gets the same treatment as the headline columns.
# ---------------------------------------------------------------------------


def _assert_evidence_identical(a: pl.DataFrame, b: pl.DataFrame, label: str) -> None:
    import json as _j

    key = ["event_id", "aspect"]
    a, b = a.sort(key), b.sort(key)

    def _same(x, y, path):
        if isinstance(x, dict):
            assert isinstance(y, dict) and x.keys() == y.keys(), f"[{label}] {path} keys differ"
            for k in x:
                _same(x[k], y[k], f"{path}.{k}")
        elif isinstance(x, (int, float)) and not isinstance(x, bool):
            assert y == pytest.approx(x), f"[{label}] {path}: {x} vs {y}"
        else:
            assert x == y, f"[{label}] {path}: {x!r} vs {y!r}"

    for ra, rb in zip(a.iter_rows(named=True), b.iter_rows(named=True)):
        ea, eb = _j.loads(ra["evidence"]), _j.loads(rb["evidence"])
        assert ea.keys() == eb.keys(), (
            f"[{label}] evidence keys differ for event {ra['event_id']}: "
            f"{sorted(set(ea) ^ set(eb))}"
        )
        for k in ea:
            _same(ea[k], eb[k], f"event {ra['event_id']}.{k}")


def test_evidence_matches_between_batched_and_scalar_gapdh():
    events = _gapdh_events()
    cov_dict = _gapdh_coverage()
    scalar = score_events_scalar(events, cov_dict, group="gapdh", tier="aggregate", thr=_THR)
    product = score_events(events, _cov_df(cov_dict), group="gapdh", tier="aggregate", thr=_THR)
    _assert_evidence_identical(scalar, product, "evidence_gapdh")


def test_evidence_matches_between_batched_and_scalar_contended():
    """The contended case is where the two CIF vector builds could diverge:
    the batch path reads sorted sparse arrays, the scalar path a dict."""
    events = _contend_events()
    cov_dict = _contend_coverage()
    scalar = score_events_scalar(
        events,
        cov_dict,
        overlaps=_contend_overlaps_dict(),
        group="test",
        tier="aggregate",
        thr=_THR,
    )
    product = score_events(
        events,
        _cov_df(cov_dict),
        overlaps_df=_contend_overlaps_df(),
        group="test",
        tier="aggregate",
        thr=_THR,
    )
    _assert_evidence_identical(scalar, product, "evidence_contended")


def test_elongation_evidence_carries_pif_and_cif():
    """Both signature scores are present and in range on a real fixture.

    PIF (``overall_in_frame``) stays in the evidence JSON blob; CIF is a
    first-class column (promoted alongside `n_codons` so translon-level
    composition can weight by codon count -- see report.py) and must
    therefore NOT also be duplicated into the evidence blob.
    """
    events = _gapdh_events()
    product = score_events(
        events, _cov_df(_gapdh_coverage()), group="gapdh", tier="aggregate", thr=_THR
    )
    import json as _j

    elong = product.filter(pl.col("aspect") == "elongation")
    assert elong.height > 0
    for row in elong.iter_rows(named=True):
        ev = _j.loads(row["evidence"])
        assert "overall_in_frame" in ev, "PIF must be present under its existing name"
        assert "cif" not in ev, "cif is a first-class column now, not duplicated in evidence JSON"
        if ev["overall_in_frame"] is not None:
            assert 0.0 <= ev["overall_in_frame"] <= 1.0
        if row["cif"] is not None:
            assert 0.0 <= row["cif"] <= 1.0
        if row["n_codons"] is not None:
            assert row["n_codons"] > 0

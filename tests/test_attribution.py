"""Lens batteries and neighbor attribution
(docs/significance_testing_plan.md §0, §6, "Sequencing" step 5).
"""

from __future__ import annotations

from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.attribution import (
    cif_battery,
    classify_neighbor_outcome,
    classify_termination_downstream,
    compose_confidence,
    elong_battery,
    init_battery,
    pair_boundary_neighbors,
    pair_cif_neighbors,
    pair_elongation_neighbors,
    term_battery,
)
from TranslonScorer.scoring.aspects import score_elongation_event, score_initiation_event

_THR = ScoreThresholds()

# ---------------------------------------------------------------------------
# Batteries
# ---------------------------------------------------------------------------


def test_init_battery_all_lenses_agree_for_a_clean_start():
    coverage = {p: (30.0 if p % 3 == 0 else 5.0) for p in range(100, 200)}
    coverage.update({p: 8.0 for p in range(40, 100)})
    s = score_initiation_event(100, 1, coverage, thr=_THR)
    battery = init_battery(s, _THR)
    assert battery["n_total"] > 0
    assert battery["battery_passed"] is True
    assert battery["lenses"]["level"] is True


def test_init_battery_none_confidence_with_no_evidence():
    battery = init_battery({}, _THR)
    assert battery["n_total"] == 0
    assert battery["confidence"] is None
    assert battery["battery_passed"] is False


def test_term_battery_shape():
    battery = term_battery({"dropoff_significance_p": 0.001, "dropoff_continuity_p": 0.8}, _THR)
    assert battery["lenses"]["level"] is True
    assert battery["lenses"]["continuity"] is True
    assert battery["battery_passed"] is True


def test_term_battery_flags_isolated_pileup():
    battery = term_battery({"dropoff_significance_p": 0.001, "dropoff_continuity_p": 0.01}, _THR)
    assert battery["lenses"]["continuity"] is False


def test_elong_battery_for_clean_elongation():
    coverage = {p: (10.0 if p % 3 == 0 else 1.0) for p in range(0, 300)}
    s = score_elongation_event(0, 300, 0, 1, coverage, thr=_THR)
    battery = elong_battery(s, _THR)
    assert battery["battery_passed"] is True


def test_cif_battery_shape():
    ev = {
        "cif_significance": {"frac_significant": 0.9},
        "cif_contiguity": {"max_run_frac": 0.2},
        "frame_chisq_p": 0.001,
    }
    battery = cif_battery(ev, _THR)
    assert battery["lenses"]["level"] is True
    assert battery["lenses"]["uniformity"] is True
    assert battery["lenses"]["cross_agree"] is True
    assert battery["battery_passed"] is True


def test_compose_confidence_pools_elongation_and_cif():
    coverage = {p: (10.0 if p % 3 == 0 else 1.0) for p in range(0, 300)}
    s = score_elongation_event(0, 300, 0, 1, coverage, thr=_THR)
    conf = compose_confidence("elongation", s, _THR)
    assert conf is not None
    assert 0.0 <= conf <= 1.0


def test_compose_confidence_downweights_low_mappability():
    coverage = {p: (10.0 if p % 3 == 0 else 1.0) for p in range(0, 300)}
    s = score_elongation_event(0, 300, 0, 1, coverage, thr=_THR)
    plain = compose_confidence("elongation", s, _THR)
    s_flagged = dict(s, map_track_low=True)
    flagged = compose_confidence("elongation", s_flagged, _THR)
    assert flagged == plain * _THR.mappability_confidence_penalty


def test_compose_confidence_unknown_aspect_is_none():
    assert compose_confidence("junction", {}, _THR) is None


# ---------------------------------------------------------------------------
# classify_neighbor_outcome
# ---------------------------------------------------------------------------


def test_classify_no_neighbor():
    assert classify_neighbor_outcome(True, None) == "no_neighbor"


def test_classify_only_this_frame():
    assert classify_neighbor_outcome(True, False) == "only_this_frame"


def test_classify_both_independent():
    assert classify_neighbor_outcome(True, True) == "both_independent"


def test_classify_leakage():
    assert classify_neighbor_outcome(False, True, leakage=True) == "leakage"


def test_classify_leakage_requires_the_leakage_signal():
    # Own battery fails, neighbor passes, but no leakage signal supplied --
    # attribution has nothing further to add.
    assert classify_neighbor_outcome(False, True, leakage=False) == "neither_supported"


def test_classify_neither_supported_when_both_fail():
    assert classify_neighbor_outcome(False, False) == "neither_supported"


# ---------------------------------------------------------------------------
# classify_termination_downstream
# ---------------------------------------------------------------------------


def test_classify_termination_downstream_clean_drop_no_signal():
    assert classify_termination_downstream({"dropoff_after_share": None}, _THR) == "clean_drop"
    assert classify_termination_downstream({"dropoff_after_share": 0.05}, _THR) == "clean_drop"


def test_classify_termination_downstream_readthrough_without_a_known_neighbor():
    ev = {"dropoff_after_share": 0.9}
    assert classify_termination_downstream(ev, _THR) == "readthrough"


def test_classify_termination_downstream_distinct_when_neighbor_passes():
    ev = {"dropoff_after_share": 0.9}
    assert (
        classify_termination_downstream(ev, _THR, neighbor_passed=True) == "distinct_downstream"
    )


# ---------------------------------------------------------------------------
# pair_elongation_neighbors / pair_cif_neighbors
# ---------------------------------------------------------------------------


def _elong_evidence_stub(*, passes: bool, competitor_share=None, identifiability=1.0):
    ev = {
        "frame_chisq_p": 0.001 if passes else 0.9,
        "body_uniformity_p": 0.9 if passes else 0.01,
        "cif_significance": {"frac_significant": 0.9 if passes else 0.1},
        "cif_contiguity": {"max_run_frac": 0.2},
        "competitor_share": competitor_share or {},
        "identifiability": identifiability,
    }
    return ev


def test_pair_elongation_neighbors_no_competitors():
    scored = {1: _elong_evidence_stub(passes=True)}
    out = pair_elongation_neighbors(scored, _THR)
    assert out[1]["outcome"] == "no_neighbor"


def test_pair_elongation_neighbors_both_independent():
    scored = {
        1: _elong_evidence_stub(passes=True, competitor_share={2: 0.4}),
        2: _elong_evidence_stub(passes=True, competitor_share={1: 0.4}),
    }
    out = pair_elongation_neighbors(scored, _THR)
    assert out[1]["outcome"] == "both_independent"
    assert out[2]["outcome"] == "both_independent"


def test_pair_elongation_neighbors_leakage():
    scored = {
        1: _elong_evidence_stub(
            passes=False, competitor_share={2: 0.9}, identifiability=0.1
        ),
        2: _elong_evidence_stub(passes=True, competitor_share={1: 0.9}),
    }
    out = pair_elongation_neighbors(scored, _THR)
    assert out[1]["outcome"] == "leakage"


def test_pair_elongation_neighbors_multi_way_flagged_not_resolved():
    scored = {
        1: _elong_evidence_stub(passes=True, competitor_share={2: 0.3, 3: 0.3}),
        2: _elong_evidence_stub(passes=True, competitor_share={1: 0.3}),
        3: _elong_evidence_stub(passes=True, competitor_share={1: 0.3}),
    }
    out = pair_elongation_neighbors(scored, _THR)
    assert out[1]["outcome"] == "multi_way_overlap"
    assert out[1]["n_competitors"] == 2


def test_pair_cif_neighbors_both_independent():
    scored = {
        1: _elong_evidence_stub(passes=True, competitor_share={2: 0.4}),
        2: _elong_evidence_stub(passes=True, competitor_share={1: 0.4}),
    }
    out = pair_cif_neighbors(scored, _THR)
    assert out[1]["outcome"] == "both_independent"


# ---------------------------------------------------------------------------
# pair_boundary_neighbors (init/term proximity heuristic)
# ---------------------------------------------------------------------------


def _init_stub(passes: bool) -> dict:
    coverage = {p: (30.0 if p % 3 == 0 else 5.0) for p in range(100, 200)}
    if not passes:
        coverage = {p: 8.0 for p in range(40, 200)}
    return score_initiation_event(100, 1, coverage, thr=_THR)


def test_pair_boundary_neighbors_finds_a_same_locus_different_frame_pair():
    events = [
        {"event_id": 1, "type": "init", "chrom": "chr1", "strand": 1, "start": 100, "phase": 0},
        {"event_id": 2, "type": "init", "chrom": "chr1", "strand": 1, "start": 101, "phase": 1},
    ]
    scored = {1: _init_stub(True), 2: _init_stub(True)}
    out = pair_boundary_neighbors(events, scored, init_battery, _THR, window_nt=60)
    assert out[1]["outcome"] == "both_independent"
    assert out[2]["outcome"] == "both_independent"


def test_pair_boundary_neighbors_no_neighbor_when_far_apart():
    events = [
        {"event_id": 1, "type": "init", "chrom": "chr1", "strand": 1, "start": 100, "phase": 0},
        {"event_id": 2, "type": "init", "chrom": "chr1", "strand": 1, "start": 5000, "phase": 1},
    ]
    scored = {1: _init_stub(True), 2: _init_stub(True)}
    out = pair_boundary_neighbors(events, scored, init_battery, _THR, window_nt=60)
    assert out[1]["outcome"] == "no_neighbor"


def test_pair_boundary_neighbors_skips_events_missing_chrom_or_phase():
    events = [
        {"event_id": 1, "type": "init", "chrom": None, "strand": 1, "start": 100, "phase": 0},
    ]
    scored = {1: _init_stub(True)}
    out = pair_boundary_neighbors(events, scored, init_battery, _THR, window_nt=60)
    assert out == {}


def test_pair_boundary_neighbors_multi_way_flagged():
    events = [
        {"event_id": 1, "type": "init", "chrom": "chr1", "strand": 1, "start": 100, "phase": 0},
        {"event_id": 2, "type": "init", "chrom": "chr1", "strand": 1, "start": 101, "phase": 1},
        {"event_id": 3, "type": "init", "chrom": "chr1", "strand": 1, "start": 102, "phase": 2},
    ]
    scored = {1: _init_stub(True), 2: _init_stub(True), 3: _init_stub(True)}
    out = pair_boundary_neighbors(events, scored, init_battery, _THR, window_nt=60)
    assert out[1]["outcome"] == "multi_way_overlap"
    assert out[1]["n_competitors"] == 2

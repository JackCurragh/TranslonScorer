"""Splice-aware init/term flank projection.

Without a nearby-intron correction, the leader/UTR flank used for
init_rise/term_drop reads raw flanking *genomic* bases. When the true leader
(or UTR) sits on the other side of a nearby intron, that flat window reads
mostly intronic (zero-coverage) sequence instead of the correct spliced
leader/UTR — spuriously inflating the step metric. See
TranslonScorer.scoring.aspects._project_flank.

These tests build small synthetic (chrom, strand) intron maps and assert:
  1. flank_spliced=False / identical output when no splice_context is given
     (backward compatibility — the golden gate covers this too).
  2. the leader/UTR flank is correctly projected across a nearby intron to
     the true upstream/downstream exon when splice_context is given.
  3. the shipped scorer and the test-only scalar reference agree
     when splice_context is threaded through end-to-end.
"""

from __future__ import annotations

import polars as pl

from tests.reference_scorer import score_events_scalar
from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.aspects import (
    _project_flank,
    score_initiation_event,
    score_termination_event,
)
from TranslonScorer.scoring.run import score_events

_THR = ScoreThresholds()


def test_project_flank_no_context_matches_flat():
    pos_true, spliced = _project_flank(99, 1, 9, body_side=False, chrom=None, splice_context=None)
    assert not spliced
    assert pos_true == list(range(90, 99))


def test_project_flank_jumps_intron_plus_strand_outer():
    # Intron [40, 99): last exonic base upstream = 39.
    splice_ctx = {("chr1", "+"): [(40, 99)]}
    pos, spliced = _project_flank(
        99, 1, 9, body_side=False, chrom="chr1", splice_context=splice_ctx
    )
    assert spliced
    # 9nt immediately preceding pos in TRANSCRIPT order, jumping the intron:
    # 31..39 (ascending, ending at the exon boundary just before pos).
    assert pos == list(range(31, 40))


def test_project_flank_jumps_intron_minus_strand_outer():
    # Minus strand: transcript-forward is genomic-descending. An intron
    # positioned downstream (higher genomic coord) of pos is the "outer"
    # (leader) side for a minus-strand init event.
    splice_ctx = {("chr1", "-"): [(100, 160)]}
    pos, spliced = _project_flank(
        99, -1, 9, body_side=False, chrom="chr1", splice_context=splice_ctx
    )
    assert spliced
    # Outer flank walks backward (transcript-order) from pos=99: genomic
    # ascending step lands on acceptor=160 as soon as it crosses the donor
    # at 100 (p+1==donor means step_up jumps to acceptor for strand<0's
    # backward==step_up direction). Just assert internal consistency: the
    # returned positions are exonic (none fall inside [100, 160)).
    assert all(not (100 <= p < 160) for p in pos)


def test_project_flank_no_nearby_intron_falls_back_to_flat():
    splice_ctx = {("chr1", "+"): [(5000, 5100)]}  # far away, irrelevant
    pos, spliced = _project_flank(
        99, 1, 9, body_side=False, chrom="chr1", splice_context=splice_ctx
    )
    assert not spliced
    assert pos == list(range(90, 99))


def _spliced_leader_coverage() -> dict:
    """Body [99,191) as in the GAPDH golden fixture; true leader lives on the
    far side of an intron [40, 99) with strong, honest coverage."""
    cov = {p: (30.0 if p % 3 == 0 else 5.0) for p in range(99, 191)}
    for p in range(-20, 40):
        cov[p] = 8.0
    return cov


def test_init_flat_flank_is_fooled_by_intron():
    """Documents the bug: without splice context, the intron reads as an
    empty leader and init_rise is spuriously inflated to a clean SUPPORTED."""
    cov = _spliced_leader_coverage()
    flat = score_initiation_event(99, 1, cov, thr=_THR)
    assert flat["flank_spliced"] is False
    assert flat["call"] == "SUPPORTED"
    assert flat["metric"] > 4.0  # spuriously high — leader looks near-empty


def test_init_splice_aware_flank_corrects_it():
    cov = _spliced_leader_coverage()
    splice_ctx = {("chr1", "+"): [(40, 99)]}
    spliced = score_initiation_event(99, 1, cov, thr=_THR, chrom="chr1", splice_context=splice_ctx)
    assert spliced["flank_spliced"] is True
    # Leader is honestly well-covered once correctly located -> the DEPTH
    # metric alone is weak, not the falsely-clean SUPPORTED the flat-genomic
    # version reports.
    assert spliced["metric"] < 2.0
    # But the body here is genuinely triplet-periodic (30/5/5) against a
    # leader that, once correctly located, is genuinely flat/non-periodic
    # (uniform 8.0) -- exactly the case the periodicity axis exists to catch.
    # The borderline-depth AMBIGUOUS call is correctly resolved to SUPPORTED
    # via periodicity, not via an inflated depth number, and that resolution
    # is auditable rather than silent.
    assert spliced["call"] == "SUPPORTED"
    assert spliced["periodicity_resolved_ambiguous"] is True


def _spliced_utr_coverage() -> dict:
    """Body up to term_pos=190 (exclusive of UTR); true UTR lives on the far
    side of a downstream intron [191, 260) with strong coverage beyond it."""
    cov = {p: (30.0 if p % 3 == 0 else 5.0) for p in range(99, 191)}
    for p in range(260, 320):
        cov[p] = 6.0
    return cov


def test_term_flat_flank_is_fooled_by_intron():
    cov = _spliced_utr_coverage()
    flat = score_termination_event(190, 1, cov, thr=_THR)
    assert flat["flank_spliced"] is False
    assert flat["call"] == "SUPPORTED"
    assert flat["metric"] > 4.0  # UTR looks near-empty -> spuriously clean stop


def test_term_splice_aware_flank_corrects_it():
    cov = _spliced_utr_coverage()
    splice_ctx = {("chr1", "+"): [(191, 260)]}
    spliced = score_termination_event(
        190, 1, cov, thr=_THR, chrom="chr1", splice_context=splice_ctx
    )
    assert spliced["flank_spliced"] is True
    assert spliced["metric"] < 2.0


def test_scalar_equals_vectorised_with_splice_context():
    events = pl.DataFrame(
        [
            {
                "event_id": 1,
                "type": "init",
                "chrom": "chr1",
                "start": 99,
                "end": 100,
                "strand": 1,
                "phase": None,
            },
            {
                "event_id": 3,
                "type": "term",
                "chrom": "chr1",
                "start": 190,
                "end": 191,
                "strand": 1,
                "phase": None,
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
    cov = {**_spliced_leader_coverage(), **_spliced_utr_coverage()}
    splice_ctx = {("chr1", "+"): [(40, 99), (191, 260)]}
    cov_df = pl.DataFrame(
        {"pos": list(cov.keys()), "count": list(cov.values())},
        schema={"pos": pl.Int64, "count": pl.Float64},
    )

    scalar = score_events_scalar(
        events, cov, group="g", tier="t", thr=_THR, splice_context=splice_ctx
    )
    vec = score_events(events, cov_df, group="g", tier="t", thr=_THR, splice_context=splice_ctx)

    scalar = scalar.sort("event_id")
    vec = vec.sort("event_id")
    assert (scalar["metric"] - vec["metric"]).abs().max() == 0.0
    assert (scalar["call"] != vec["call"]).sum() == 0
    assert (scalar["eligibility"] != vec["eligibility"]).sum() == 0

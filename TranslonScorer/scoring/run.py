"""End-to-end event scoring: scalar reference and vectorised batch orchestration.

score_events           — scalar reference scorer (one event at a time)
score_events_vectorised — vectorised scorer (elongation via prefix sums)
score_elongation_batch — vectorised elongation kernel (reusable)

Prefix-sum helpers (_prefix_sums, _range_sums) will migrate to
coverage/profile.py in T10 when the coverage provider layer is built.

Public API
----------
DEFAULT_THRESHOLDS      — default ScoreThresholds instance
score_elongation_batch  — vectorised elongation: prefix sums → {event_id: dict}
score_events            — scalar reference scorer → pl.DataFrame
score_events_vectorised — vectorised scorer → pl.DataFrame
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as _np
import polars as pl

from TranslonScorer.events import SpliceContext
from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.aspects import (
    score_elongation_event,
    score_initiation_event,
    score_junction_event,
    score_termination_event,
)
from TranslonScorer.scoring.evidence import (
    _RECORD_SCHEMA,
    _elong_evidence,
    event_record,
)

DEFAULT_THRESHOLDS = ScoreThresholds()


# ---------------------------------------------------------------------------
# Prefix-sum helpers — canonical implementation in coverage/profile.py
# ---------------------------------------------------------------------------
from TranslonScorer.coverage.profile import (  # noqa: E402, F401
    _prefix_sums,
    _range_sums,
)

# ---------------------------------------------------------------------------
# Vectorised elongation kernel
# ---------------------------------------------------------------------------


def _accumulate_contended_coverage(
    overlaps_df: Optional[pl.DataFrame],
    elong: pl.DataFrame,
    row_of_event: Dict[int, int],
    strand: "_np.ndarray",
    pos_sorted: "_np.ndarray",
    cum_total: "_np.ndarray",
    cum_by_frame: List,
    n_events: int,
):
    """Sum competitor (overlapping-event) coverage per elongation event.

    For each elongation event, the coverage contributed by other events that
    overlap it (its "competitors") must be subtracted before scoring, so
    contended reads aren't double-counted. Returns four aligned outputs:

      contended_total[n_events]        competitor read total per event
      contended_frame[n_events, 3]     competitor reads per frame per event
      contended_nt[n_events]           competitor-covered nt per event
      competitor_frames{row: {comp_id: comp_a_site_frame}}
    """
    contended_total = _np.zeros(n_events)
    contended_frame = _np.zeros((n_events, 3))
    contended_nt = _np.zeros(n_events)
    competitor_frames: Dict[int, dict] = {}
    if overlaps_df is None or overlaps_df.is_empty():
        return contended_total, contended_frame, contended_nt, competitor_frames

    overlaps = (
        overlaps_df.with_columns(
            [
                pl.col("event_id").cast(pl.UInt64),
                pl.col("other_event_id").cast(pl.UInt64),
            ]
        )
        .filter(pl.col("event_id").is_in(elong["event_id"].cast(pl.UInt64)))
        .sort(["event_id", "overlap_start"])
    )
    if overlaps.is_empty():
        return contended_total, contended_frame, contended_nt, competitor_frames

    # Merge each event's overlap intervals, collecting flat (start, end, row)
    # segments to sum in one vectorised _range_sums call.
    seg_starts: List[int] = []
    seg_ends: List[int] = []
    seg_rows: List[int] = []
    for event_id, grp in overlaps.group_by("event_id", maintain_order=True):
        row = row_of_event[int(event_id[0] if isinstance(event_id, tuple) else event_id)]
        intervals = sorted(zip(grp["overlap_start"].to_list(), grp["overlap_end"].to_list()))
        cur_start, cur_end = intervals[0]
        merged = []
        for start, end in intervals[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((cur_start, cur_end))
        for start, end in merged:
            seg_starts.append(start)
            seg_ends.append(end)
            seg_rows.append(row)
        for comp_id, comp_phase in zip(
            grp["other_event_id"].to_list(), grp["comp_phase"].to_list()
        ):
            comp_frame = (-comp_phase) % 3 if strand[row] > 0 else comp_phase % 3
            competitor_frames.setdefault(row, {})[int(comp_id)] = int(comp_frame)

    seg_start = _np.array(seg_starts, dtype=_np.int64)
    seg_end = _np.array(seg_ends, dtype=_np.int64)
    seg_row = _np.array(seg_rows)
    seg_total, seg_frame, _ = _range_sums(pos_sorted, cum_total, cum_by_frame, seg_start, seg_end)
    _np.add.at(contended_total, seg_row, seg_total)
    _np.add.at(contended_frame, seg_row, seg_frame)
    _np.add.at(contended_nt, seg_row, (seg_end - seg_start))
    return contended_total, contended_frame, contended_nt, competitor_frames


def score_elongation_batch(
    elong: pl.DataFrame,
    overlaps_df: Optional[pl.DataFrame],
    cov_pos: "_np.ndarray",
    cov_cnt: "_np.ndarray",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
) -> Dict[int, dict]:
    """Vectorised elongation: per-event per-frame sums via coverage prefix sums,
    then shared _elong_evidence. Identical decision logic to the scalar path;
    only the summation is vectorised."""
    if elong.is_empty():
        return {}
    pos_sorted, cum_total, cum_by_frame = _prefix_sums(cov_pos, cov_cnt)
    event_ids = elong["event_id"].to_numpy()
    starts = elong["start"].to_numpy().astype(_np.int64)
    ends = elong["end"].to_numpy().astype(_np.int64)
    phase = elong["phase"].to_numpy().astype(_np.int64)
    strand = elong["strand"].to_numpy().astype(_np.int64)
    expected_frame = _np.where(strand > 0, _np.mod(-phase, 3), _np.mod(phase, 3))

    total_reads, frame_reads, covered_positions = _range_sums(
        pos_sorted, cum_total, cum_by_frame, starts, ends
    )
    row_of_event = {int(e): i for i, e in enumerate(event_ids)}

    (
        contended_total,
        contended_frame,
        contended_nt,
        competitor_frames,
    ) = _accumulate_contended_coverage(
        overlaps_df,
        elong,
        row_of_event,
        strand,
        pos_sorted,
        cum_total,
        cum_by_frame,
        len(event_ids),
    )

    out: Dict[int, dict] = {}
    for row, event_id in enumerate(event_ids):
        clean_total = float(total_reads[row] - contended_total[row])
        clean_in_frame = float(
            frame_reads[row][expected_frame[row]] - contended_frame[row][expected_frame[row]]
        )
        out[int(event_id)] = _elong_evidence(
            float(total_reads[row]),
            int(covered_positions[row]),
            clean_total,
            clean_in_frame,
            float(contended_total[row]),
            list(map(float, contended_frame[row])),
            int(expected_frame[row]),
            competitor_frames.get(row, {}),
            int(contended_nt[row]),
            int(ends[row] - starts[row]),
            thr,
        )
    return out


# ---------------------------------------------------------------------------
# Main entry points
# ---------------------------------------------------------------------------


def score_events(
    events: pl.DataFrame,
    coverage: Dict[int, float],
    *,
    overlaps: Optional[Dict[int, list]] = None,
    junction_support: Optional[Dict[int, dict]] = None,
    group: str = "aggregate",
    tier: str = "aggregate",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
    splice_context: Optional[SpliceContext] = None,
    map_track: Optional[Dict[int, dict]] = None,
) -> pl.DataFrame:
    """Reference scorer: score every event in `events` → granular long-form table.

    events columns: event_id, type, start, end, strand, phase[, chrom].
    overlaps: {elongation event_id → [(o_start, o_end, competitor_phase, competitor_id)]}.
    junction_support: {junction event_id → spanning count}.
    splice_context: optional {(chrom, strand) → [(donor, acceptor), ...]} (see
        events.build_splice_context). When given (and events carries a
        `chrom` column), init/term leader/UTR flanks are projected across
        nearby introns instead of read as raw flanking genomic bases.
    map_track: optional {event_id → {"map_track_mean": float, "map_track_low":
        bool}} (see workflows._map_track_for_chrom) — a region-context
        mappability annotation, merged into every event type uniformly. Never
        affects eligibility/call.

    Scalar reference implementation; score_events_vectorised must reproduce this.
    """
    overlaps = overlaps or {}
    junction_support = junction_support or {}
    map_track = map_track or {}
    rows: List[dict] = []
    for r in events.iter_rows(named=True):
        t, eid = r["type"], r["event_id"]
        if t == "init":
            raw = score_initiation_event(
                r["start"],
                r["strand"],
                coverage,
                thr=thr,
                chrom=r.get("chrom"),
                splice_context=splice_context,
            )
        elif t == "term":
            raw = score_termination_event(
                r["start"],
                r["strand"],
                coverage,
                thr=thr,
                chrom=r.get("chrom"),
                splice_context=splice_context,
            )
        elif t == "elongation":
            raw = score_elongation_event(
                r["start"],
                r["end"],
                r["phase"],
                r["strand"],
                coverage,
                overlaps.get(eid),
                thr=thr,
            )
        elif t == "junction":
            raw = score_junction_event(junction_support.get(eid, {}), thr=thr)
        else:
            continue
        if eid in map_track:
            raw = {**raw, **map_track[eid]}
        rows.append(event_record(eid, t, group, tier, raw, thr.version))
    return (
        pl.from_dicts(rows, schema=_RECORD_SCHEMA) if rows else pl.DataFrame(schema=_RECORD_SCHEMA)
    )


def score_events_vectorised(
    events: pl.DataFrame,
    cov_df: pl.DataFrame,
    *,
    overlaps_df: Optional[pl.DataFrame] = None,
    junction_support: Optional[Dict[int, dict]] = None,
    group: str = "aggregate",
    tier: str = "aggregate",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
    splice_context: Optional[SpliceContext] = None,
    map_track: Optional[Dict[int, dict]] = None,
) -> pl.DataFrame:
    """Vectorised scorer: elongation via prefix sums, init/term scalar.

    cov_df columns: pos (Int64), count (Float64).
    overlaps_df columns: event_id, other_event_id, overlap_start, overlap_end, comp_phase.
    splice_context: see score_events — projects init/term leader/UTR flanks
        across nearby introns when given (and events carries `chrom`).
    map_track: see score_events — region-context mappability annotation.

    Reproduces score_events exactly; only elongation summation is vectorised.
    """
    junction_support = junction_support or {}
    map_track = map_track or {}
    coverage = dict(zip(cov_df["pos"].to_list(), cov_df["count"].to_list()))
    cov_pos = cov_df["pos"].to_numpy()
    cov_cnt = cov_df["count"].to_numpy()

    elong = events.filter(pl.col("type") == "elongation")
    elong_raw = score_elongation_batch(elong, overlaps_df, cov_pos, cov_cnt, thr)

    rows: List[dict] = []
    for r in events.iter_rows(named=True):
        t, eid = r["type"], r["event_id"]
        if t == "elongation":
            raw = elong_raw.get(eid)
            if raw is None:
                continue
        elif t == "init":
            raw = score_initiation_event(
                r["start"],
                r["strand"],
                coverage,
                thr=thr,
                chrom=r.get("chrom"),
                splice_context=splice_context,
            )
        elif t == "term":
            raw = score_termination_event(
                r["start"],
                r["strand"],
                coverage,
                thr=thr,
                chrom=r.get("chrom"),
                splice_context=splice_context,
            )
        elif t == "junction":
            raw = score_junction_event(junction_support.get(eid, {}), thr=thr)
        else:
            continue
        if eid in map_track:
            raw = {**raw, **map_track[eid]}
        rows.append(event_record(eid, t, group, tier, raw, thr.version))
    return (
        pl.from_dicts(rows, schema=_RECORD_SCHEMA) if rows else pl.DataFrame(schema=_RECORD_SCHEMA)
    )

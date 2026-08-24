"""End-to-end event scoring.

One scorer for every coverage source. Elongation is batched through prefix
sums; init/term/junction are a per-event loop. An independent scalar
implementation is kept in ``tests/reference_scorer.py`` and diffed against
this one — deliberately outside product code, so there is only ever one
scoring path shipping.

Public API
----------
DEFAULT_THRESHOLDS      — default ScoreThresholds instance
score_elongation_batch  — batched elongation: prefix sums → {event_id: dict}
score_events            — score an event table against coverage → pl.DataFrame
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as _np
import polars as pl

from TranslonScorer.events import SpliceContext
from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.aspects import (
    _elong_body_uniformity,
    _elong_frame_chisq,
    score_initiation_event,
    score_junction_event,
    score_termination_event,
)
from TranslonScorer.scoring.evidence import (
    _RECORD_SCHEMA,
    _elong_evidence,
    event_record,
)
from TranslonScorer.scoring.signature import (
    cif,
    cif_codon_contiguity,
    cif_codon_significance,
    orf_signal_vector,
)

DEFAULT_THRESHOLDS = ScoreThresholds()


# ---------------------------------------------------------------------------
# Prefix-sum helpers — canonical implementation in coverage/profile.py
# ---------------------------------------------------------------------------
from TranslonScorer.coverage.psite_profile import (  # noqa: E402, F401
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
    # cum_total is a prefix array (len n+1, [0] = 0), so this recovers the
    # sorted counts without a second argsort. Needed because CIF is per-codon
    # and cannot be derived from event-level frame sums.
    cnt_sorted = _np.diff(cum_total)
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
        vec = orf_signal_vector(
            pos_sorted,
            cnt_sorted,
            int(starts[row]),
            int(ends[row]),
            int(expected_frame[row]),
            int(strand[row]),
        )
        cif_sig = cif_codon_significance(
            vec, min_reads=thr.cif_codon_min_reads, alpha=thr.periodicity_significance_alpha
        )
        cif_contig = cif_codon_contiguity(cif_sig["sig_mask"]) if cif_sig else None
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
            cif_value=cif(vec),
            n_codons=(len(vec) // 3) if len(vec) else None,
            frame_chisq_p=_elong_frame_chisq(frame_reads[row].tolist()),
            cif_significance=cif_sig,
            cif_contiguity=cif_contig,
            body_uniformity_p=_elong_body_uniformity(vec, min_codons=thr.elong_uniformity_min_codons),
        )
    return out


# ---------------------------------------------------------------------------
# Main entry points
# ---------------------------------------------------------------------------


def score_events(
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
    """Score every event in `events` → granular long-form table.

    The one scorer. Elongation is batched through prefix sums (where the volume
    is); init/term/junction are a plain per-event loop (far fewer point events,
    and clearer than a fiddly vectorisation).

    events columns: event_id, type, start, end, strand, phase[, chrom].
    cov_df: per-position coverage — pos (Int64), count (Float64).
    overlaps_df: event_id, other_event_id, overlap_start, overlap_end, comp_phase.
    junction_support: {junction event_id → spanning count}.
    splice_context: optional {(chrom, strand) → [(donor, acceptor), ...]} (see
        events.build_splice_context). When given (and events carries a
        `chrom` column), init/term leader/UTR flanks are projected across
        nearby introns instead of read as raw flanking genomic bases.
    map_track: optional {event_id → {"map_track_mean": float, "map_track_low":
        bool}} (see workflows._map_track_for_chrom) — a region-context
        mappability annotation, merged into every event type uniformly. Never
        affects eligibility/call.

    An independent scalar implementation lives in ``tests/reference_scorer.py``
    and is diffed against this one on every run; it is not product code.
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

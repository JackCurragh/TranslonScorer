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

from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.aspects import (
    score_initiation_event,
    score_termination_event,
    score_elongation_event,
    score_junction_event,
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


def score_elongation_batch(
    elong: pl.DataFrame,
    overlaps_df: pl.DataFrame,
    cov_pos: "_np.ndarray",
    cov_cnt: "_np.ndarray",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
) -> Dict[int, dict]:
    """Vectorised elongation: per-event per-frame sums via coverage prefix sums,
    then shared _elong_evidence. Identical decision logic to the scalar path;
    only the summation is vectorised."""
    if elong.is_empty():
        return {}
    p, cum_all, cum_f = _prefix_sums(cov_pos, cov_cnt)
    eids = elong["event_id"].to_numpy()
    starts = elong["start"].to_numpy().astype(_np.int64)
    ends = elong["end"].to_numpy().astype(_np.int64)
    phase = elong["phase"].to_numpy().astype(_np.int64)
    strand = elong["strand"].to_numpy().astype(_np.int64)
    a_e = _np.where(strand > 0, _np.mod(-phase, 3), _np.mod(phase, 3))

    tot_w, fr_w, ncov_w = _range_sums(p, cum_all, cum_f, starts, ends)
    idx_of = {int(e): i for i, e in enumerate(eids)}

    cont_tot = _np.zeros(len(eids))
    cont_fr = _np.zeros((len(eids), 3))
    cont_nt = _np.zeros(len(eids))
    comp_frames: Dict[int, dict] = {}
    if overlaps_df is not None and not overlaps_df.is_empty():
        od = (
            overlaps_df.with_columns(
                [
                    pl.col("event_id").cast(pl.UInt64),
                    pl.col("other_event_id").cast(pl.UInt64),
                ]
            )
            .filter(pl.col("event_id").is_in(elong["event_id"].cast(pl.UInt64)))
            .sort(["event_id", "overlap_start"])
        )
        if not od.is_empty():
            m_start, m_end, m_idx = [], [], []
            for eid, grp in od.group_by("event_id", maintain_order=True):
                i = idx_of[int(eid[0] if isinstance(eid, tuple) else eid)]
                ivs = sorted(zip(grp["overlap_start"].to_list(), grp["overlap_end"].to_list()))
                cs, ce = ivs[0]
                merged = []
                for s, e in ivs[1:]:
                    if s <= ce:
                        ce = max(ce, e)
                    else:
                        merged.append((cs, ce))
                        cs, ce = s, e
                merged.append((cs, ce))
                for s, e in merged:
                    m_start.append(s)
                    m_end.append(e)
                    m_idx.append(i)
                for cid, cph in zip(grp["other_event_id"].to_list(), grp["comp_phase"].to_list()):
                    af = (-cph) % 3 if strand[i] > 0 else cph % 3
                    comp_frames.setdefault(i, {})[int(cid)] = int(af)
            m_start = _np.array(m_start, dtype=_np.int64)
            m_end = _np.array(m_end, dtype=_np.int64)
            m_idx = _np.array(m_idx)
            o_tot, o_fr, _ = _range_sums(p, cum_all, cum_f, m_start, m_end)
            _np.add.at(cont_tot, m_idx, o_tot)
            _np.add.at(cont_fr, m_idx, o_fr)
            _np.add.at(cont_nt, m_idx, (m_end - m_start))

    out: Dict[int, dict] = {}
    for i, eid in enumerate(eids):
        clean_tot = float(tot_w[i] - cont_tot[i])
        clean_inf = float(fr_w[i][a_e[i]] - cont_fr[i][a_e[i]])
        out[int(eid)] = _elong_evidence(
            float(tot_w[i]),
            int(ncov_w[i]),
            clean_tot,
            clean_inf,
            float(cont_tot[i]),
            list(map(float, cont_fr[i])),
            int(a_e[i]),
            comp_frames.get(i, {}),
            int(cont_nt[i]),
            int(ends[i] - starts[i]),
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
    junction_support: Optional[Dict[int, float]] = None,
    group: str = "aggregate",
    tier: str = "aggregate",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
) -> pl.DataFrame:
    """Reference scorer: score every event in `events` → granular long-form table.

    events columns: event_id, type, start, end, strand, phase.
    overlaps: {elongation event_id → [(o_start, o_end, competitor_phase, competitor_id)]}.
    junction_support: {junction event_id → spanning count}.

    Scalar reference implementation; score_events_vectorised must reproduce this.
    """
    overlaps = overlaps or {}
    junction_support = junction_support or {}
    rows: List[dict] = []
    for r in events.iter_rows(named=True):
        t, eid = r["type"], r["event_id"]
        if t == "init":
            raw = score_initiation_event(r["start"], r["strand"], coverage, thr=thr)
        elif t == "term":
            raw = score_termination_event(r["start"], r["strand"], coverage, thr=thr)
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
        rows.append(event_record(eid, t, group, tier, raw, thr.version))
    return (
        pl.from_dicts(rows, schema=_RECORD_SCHEMA) if rows else pl.DataFrame(schema=_RECORD_SCHEMA)
    )


def score_events_vectorised(
    events: pl.DataFrame,
    cov_df: pl.DataFrame,
    *,
    overlaps_df: Optional[pl.DataFrame] = None,
    junction_support: Optional[Dict[int, float]] = None,
    group: str = "aggregate",
    tier: str = "aggregate",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
) -> pl.DataFrame:
    """Vectorised scorer: elongation via prefix sums, init/term scalar.

    cov_df columns: pos (Int64), count (Float64).
    overlaps_df columns: event_id, other_event_id, overlap_start, overlap_end, comp_phase.

    Reproduces score_events exactly; only elongation summation is vectorised.
    """
    junction_support = junction_support or {}
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
            raw = score_initiation_event(r["start"], r["strand"], coverage, thr=thr)
        elif t == "term":
            raw = score_termination_event(r["start"], r["strand"], coverage, thr=thr)
        elif t == "junction":
            raw = score_junction_event(junction_support.get(eid, {}), thr=thr)
        else:
            continue
        rows.append(event_record(eid, t, group, tier, raw, thr.version))
    return (
        pl.from_dicts(rows, schema=_RECORD_SCHEMA) if rows else pl.DataFrame(schema=_RECORD_SCHEMA)
    )

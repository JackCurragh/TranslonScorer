"""Per-event scoring with frame-resolved identifiability in contended regions.

An elongation event with phase ``φ`` has its in-frame (codon position 0) signal
at genomic positions where ``(φ ± p) % 3 == 0`` — i.e. in a single absolute
genomic frame ``a_E``:

    + strand:  a_E = (-φ) % 3
    - strand:  a_E = ( φ) % 3

A competing event overlapping in a *different* phase occupies a *different*
genomic frame, so in the overlap the coverage partitions cleanly by ``p % 3``
into {this event, each competitor, residual/noise}. Identifiability is then the
fraction of the convoluted signal attributable to this event's frame — no ORF
enumeration needed, just the frame decomposition the matrix already produces.

Clean (un-overlapped) and contended portions are scored separately: a real ORF
should be strongly in-frame where it is *unconfounded*; the contended portion's
identifiability says whether the overlap actually muddies attribution.
"""
from __future__ import annotations

import math
import statistics

from typing import Dict, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Decision layer — kept SEPARATE from measurement. Scorers emit continuous
# evidence + an eligibility (can we judge this at all?); `call` is a thin,
# auditable function of the evidence and these thresholds. Change a threshold,
# re-derive calls, never re-measure.
# ---------------------------------------------------------------------------

from TranslonScorer.model import ScoreThresholds  # noqa: F401 — re-exported for back-compat

DEFAULT_THRESHOLDS = ScoreThresholds()

# Eligibility ∈ {NA, INSUFFICIENT, ELIGIBLE}; call ∈ {SUPPORTED, UNSUPPORTED,
# AMBIGUOUS, None}.  "couldn't test" (INSUFFICIENT/NA) is never "failed".


def _decide_step(ev: dict, *, rise_thr: float, thr: ScoreThresholds) -> Tuple[str, Optional[str]]:
    # Trust the robust median consensus (already discounts flank spikes). The
    # peakiness/stability fields stay in evidence as review flags, but do NOT
    # hard-gate the call — that wrongly overrode strong, clean starts/stops.
    if ev["n_reads"] < thr.min_reads:
        return "INSUFFICIENT", None
    r = ev["consensus_rise"]
    if r >= rise_thr:
        return "ELIGIBLE", "SUPPORTED"
    if r <= 0:
        return "ELIGIBLE", "UNSUPPORTED"
    return "ELIGIBLE", "AMBIGUOUS"     # borderline 0 < rise < threshold


# ---------------------------------------------------------------------------
# Initiation / termination: smoothed, flank-robust relative-change scoring
# ---------------------------------------------------------------------------
# A real start is a *step up* in coverage across the start site; a real stop is
# a *step down* after it. We must not be fooled by a single spurious peak in a
# flank, so each flank is binned into codons (3 nt, periodicity-aware) and we
# take the MEDIAN bin as the robust level. We compute the step at several flank
# lengths and report the consensus + its stability, plus an explicit
# flank-peakiness flag (max bin / median bin) that catches a confounding peak.

def _flank_positions(pos: int, strand: int, length: int, *, body_side: bool) -> range:
    """Genomic positions for a flank of `length` nt on the body or outer side of
    `pos`, in transcript orientation. body_side=True → into the ORF."""
    if strand > 0:
        return range(pos, pos + length) if body_side else range(pos - length, pos)
    return range(pos - length + 1, pos + 1) if body_side else range(pos + 1, pos + length + 1)


def _codon_levels(coverage: Dict[int, float], positions: Sequence[int]) -> Tuple[float, float, float]:
    """(median codon-bin sum, max codon-bin sum, total) over a window. Median is
    robust to a single spike; codon bins fold in 3-nt periodicity."""
    pos = list(positions)
    bins = [sum(coverage.get(p, 0.0) for p in pos[i:i + 3]) for i in range(0, max(len(pos) - 2, 0), 3)]
    total = sum(coverage.get(p, 0.0) for p in pos)
    if not bins:
        return total, total, total
    return statistics.median(bins), max(bins), total


def _step_score(
    pos: int, strand: int, coverage: Dict[int, float],
    flanks: Sequence[int], min_reads: float, alpha: float,
) -> dict:
    """Robust step (body_level - outer_level)/(sum) across `pos`, over several
    flank lengths. Metric = log2 fold-change of body vs outer robust levels with
    a pseudocount α (de-saturates; positive = body higher than the outer flank)."""
    rises: Dict[int, float] = {}
    for L in flanks:
        body_med, _, _ = _codon_levels(coverage, _flank_positions(pos, strand, L, body_side=True))
        out_med, _, _ = _codon_levels(coverage, _flank_positions(pos, strand, L, body_side=False))
        rises[L] = math.log2((body_med + alpha) / (out_med + alpha))
    vals = list(rises.values())
    consensus = statistics.median(vals)
    stability = (max(vals) - min(vals)) if len(vals) > 1 else 0.0  # spread across scales
    refL = sorted(flanks)[len(flanks) // 2]
    out_med, out_max, _ = _codon_levels(coverage, _flank_positions(pos, strand, refL, body_side=False))
    flank_peakiness = (out_max / out_med) if out_med > 0 else (float("inf") if out_max > 0 else 0.0)
    n_reads = _codon_levels(coverage, _flank_positions(pos, strand, min(flanks), body_side=True))[2]
    return {
        "rise_by_flank": rises, "consensus_rise": consensus, "stability": stability,
        "flank_peakiness": flank_peakiness, "n_reads": n_reads,
    }


def score_initiation_event(start_pos: int, strand: int, coverage: Dict[int, float], *,
                           flanks: Sequence[int] = (9, 18, 30, 60),
                           thr: ScoreThresholds = DEFAULT_THRESHOLDS) -> dict:
    """Initiation = smoothed, flank-robust step UP at the start. Continuous
    evidence: log2FC at each flank length, consensus, scale-stability, and a
    flank-peakiness review flag. headline metric = consensus_rise (log2 fold-change)."""
    s = _step_score(start_pos, strand, coverage, flanks, thr.min_reads, thr.step_pseudocount)
    s["metric"] = s["consensus_rise"]
    s["metric_name"] = "init_rise"
    s["eligibility"], s["call"] = _decide_step(s, rise_thr=thr.init_rise, thr=thr)
    return s


def score_termination_event(term_pos: int, strand: int, coverage: Dict[int, float], *,
                            flanks: Sequence[int] = (9, 18, 30, 60),
                            thr: ScoreThresholds = DEFAULT_THRESHOLDS) -> dict:
    """Termination = step DOWN after the stop (body before > UTR after). Same
    robust machinery as initiation. headline metric = consensus log2 fold-change (drop)."""
    alpha = thr.step_pseudocount
    rises: Dict[int, float] = {}
    for L in flanks:
        body_med, _, _ = _codon_levels(coverage, _flank_positions(term_pos, strand, L, body_side=False))
        utr_med, _, _ = _codon_levels(coverage, _flank_positions(term_pos, strand, L, body_side=True))
        rises[L] = math.log2((body_med + alpha) / (utr_med + alpha))
    vals = list(rises.values())
    refL = sorted(flanks)[len(flanks) // 2]
    utr_med, utr_max, _ = _codon_levels(coverage, _flank_positions(term_pos, strand, refL, body_side=True))
    s = {
        "rise_by_flank": rises, "consensus_rise": statistics.median(vals),
        "stability": (max(vals) - min(vals)) if len(vals) > 1 else 0.0,
        "flank_peakiness": (utr_max / utr_med) if utr_med > 0 else (float("inf") if utr_max > 0 else 0.0),
        "n_reads": _codon_levels(coverage, _flank_positions(term_pos, strand, min(flanks), body_side=False))[2],
    }
    s["metric"] = s["consensus_rise"]
    s["metric_name"] = "term_drop"
    s["eligibility"], s["call"] = _decide_step(s, rise_thr=thr.term_drop, thr=thr)
    return s


def abs_frame(phase: int, strand: int) -> int:
    """Genomic frame (p % 3) carrying this event's in-frame signal."""
    return (-phase) % 3 if strand > 0 else phase % 3


def score_elongation_event(
    start: int,
    end: int,
    phase: int,
    strand: int,
    coverage: Dict[int, float],
    overlaps: Optional[List[Tuple[int, int, int, int]]] = None,
    *,
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
) -> dict:
    """Score one elongation event.

    overlaps: list of (overlap_start, overlap_end, competitor_phase,
    competitor_event_id) for the competing events that share genomic span in a
    different frame.

    Returns clean vs contended in-frame fractions, per-overlap identifiability,
    and per-competitor signal shares.
    """
    overlaps = overlaps or []
    a_e = abs_frame(phase, strand)

    contended: Dict[int, None] = {}
    comp_frame: Dict[int, int] = {}        # competitor_id -> its genomic frame
    for os, oe, cph, cid in overlaps:
        for p in range(max(os, start), min(oe, end)):
            contended[p] = None
        comp_frame[cid] = abs_frame(cph, strand)

    n = clean_tot = clean_inf = cont_tot = 0.0
    cont_by_frame = [0.0, 0.0, 0.0]
    covered = 0
    for p in range(start, end):
        c = coverage.get(p, 0.0)
        if c <= 0:
            continue
        covered += 1
        n += c
        if p in contended:
            cont_tot += c
            cont_by_frame[p % 3] += c
        else:
            clean_tot += c
            if p % 3 == a_e:
                clean_inf += c

    return _elong_evidence(n, covered, clean_tot, clean_inf, cont_tot, cont_by_frame,
                           a_e, comp_frame, len(contended), end - start, thr)


def _elong_evidence(n, covered, clean_tot, clean_inf, cont_tot, cont_by_frame,
                    a_e, comp_frame, contended_nt, span_nt, thr) -> dict:
    """Turn per-event frame sums into evidence + eligibility/call. Shared by the
    scalar and vectorised scorers so the decision logic is identical."""
    cont_e = cont_by_frame[a_e]
    clean_in_frame = (clean_inf / clean_tot) if clean_tot else None
    identifiability = (cont_e / cont_tot) if cont_tot else None
    overall_in_frame = ((clean_inf + cont_e) / n) if n else 0.0
    # breadth = fraction of the segment's positions covered — distinguishes a
    # processive ORF (broadly covered) from one in-frame spike (high in-frame,
    # tiny breadth). One-third coverage is the periodicity ceiling.
    breadth = (covered / span_nt) if span_nt else 0.0

    competitor_share: Dict[int, float] = {
        cid: ((cont_by_frame[af] / cont_tot) if cont_tot else 0.0) for cid, af in comp_frame.items()}
    noise_frame = ({0, 1, 2} - {a_e} - set(comp_frame.values()))
    noise_share = (sum(cont_by_frame[f] for f in noise_frame) / cont_tot) if cont_tot else 0.0

    effective_in_frame = clean_in_frame if clean_in_frame is not None else overall_in_frame
    if n < thr.min_reads:
        eligibility, call = "INSUFFICIENT", None
    elif identifiability is not None and identifiability < thr.elong_identifiability:
        eligibility, call = "ELIGIBLE", "AMBIGUOUS"
    elif effective_in_frame >= thr.elong_in_frame and breadth >= thr.elong_breadth:
        eligibility, call = "ELIGIBLE", "SUPPORTED"
    else:
        eligibility, call = "ELIGIBLE", "UNSUPPORTED"   # weak periodicity or single-spike

    return {
        "n_reads": n, "covered_nt": covered, "metric": effective_in_frame,
        "metric_name": "elong_in_frame", "overall_in_frame": overall_in_frame,
        "breadth": breadth, "span_nt": span_nt, "clean_in_frame": clean_in_frame,
        "contended_nt": contended_nt, "identifiability": identifiability,
        "competitor_share": competitor_share, "noise_share": noise_share,
        "eligibility": eligibility, "call": call,
    }


def score_junction_event(support: dict, *, thr: ScoreThresholds = DEFAULT_THRESHOLDS) -> dict:
    """Splice support with an overhang-confidence (identifiability) layer.

    support = {span_conf, span_short, unspliced} weighted read counts.
      span_conf  — intron matches with ≥6 nt aligned both sides (confident map)
      span_short — matches but short overhang (~ambiguous across the junction)
      unspliced  — aligned straight through the donor (intron retention)

    headline metric = confident spanning count; psi = spliced-in ratio.
    """
    conf = float(support.get("span_conf", 0.0))
    short = float(support.get("span_short", 0.0))
    unspl = float(support.get("unspliced", 0.0))
    spanning = conf + short
    psi = (spanning / (spanning + unspl)) if (spanning + unspl) > 0 else None
    if conf >= thr.junc_min_spanning:
        elig, call = "ELIGIBLE", "SUPPORTED"
    elif spanning >= thr.junc_min_spanning:
        elig, call = "ELIGIBLE", "AMBIGUOUS"      # support only via short overhangs
    else:
        elig, call = "INSUFFICIENT", None
    return {"n_reads": spanning, "metric": conf, "metric_name": "junc_confident_spanning",
            "confident_spanning": conf, "short_spanning": short, "unspliced": unspl,
            "psi": psi, "eligibility": elig, "call": call}


# ---------------------------------------------------------------------------
# Granular score records — uniform long-form rows (the persisted reference).
# Every continuous field survives in `evidence`; eligibility/call/thresholds
# are explicit so the decision is re-derivable without re-measuring.
# ---------------------------------------------------------------------------

import json as _json
import math as _math

import polars as pl

_RECORD_SCHEMA = {
    "event_id": pl.UInt64, "aspect": pl.Utf8, "group": pl.Utf8, "tier": pl.Utf8,
    "n_reads": pl.Float64, "metric": pl.Float64, "metric_name": pl.Utf8,
    "eligibility": pl.Utf8, "call": pl.Utf8, "evidence": pl.Utf8, "thresholds_version": pl.Utf8,
}


def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, float) and (_math.isinf(v) or _math.isnan(v)):
        return "inf" if v > 0 else "-inf"
    return v


def event_record(event_id: int, aspect: str, group: str, tier: str, raw: dict,
                 thr_version: str = "v0") -> dict:
    drop = {"eligibility", "call", "metric", "metric_name", "n_reads"}
    evidence = {k: _jsonable(v) for k, v in raw.items() if k not in drop}
    m = raw.get("metric")
    return {
        "event_id": int(event_id), "aspect": aspect, "group": group, "tier": tier,
        "n_reads": float(raw.get("n_reads", 0.0) or 0.0),
        "metric": (None if m is None else float(m)),
        "metric_name": raw.get("metric_name"),
        "eligibility": raw["eligibility"], "call": raw.get("call"),
        "evidence": _json.dumps(evidence), "thresholds_version": thr_version,
    }


import numpy as _np


def _prefix_sums(cov_pos: "_np.ndarray", cov_cnt: "_np.ndarray"):
    order = _np.argsort(cov_pos, kind="stable")
    p = cov_pos[order].astype(_np.int64)
    c = cov_cnt[order].astype(_np.float64)
    fr = _np.mod(p, 3)
    cum_all = _np.concatenate([[0.0], _np.cumsum(c)])
    cum_f = [_np.concatenate([[0.0], _np.cumsum(_np.where(fr == f, c, 0.0))]) for f in range(3)]
    return p, cum_all, cum_f


def _range_sums(p, cum_all, cum_f, starts, ends):
    lo = _np.searchsorted(p, starts, "left")
    hi = _np.searchsorted(p, ends, "left")
    total = cum_all[hi] - cum_all[lo]
    frames = _np.stack([cum_f[f][hi] - cum_f[f][lo] for f in range(3)], axis=-1)  # [n,3]
    return total, frames, (hi - lo)


def score_elongation_batch(elong: pl.DataFrame, overlaps_df: pl.DataFrame,
                           cov_pos: "_np.ndarray", cov_cnt: "_np.ndarray",
                           thr: ScoreThresholds = DEFAULT_THRESHOLDS) -> Dict[int, dict]:
    """Vectorised elongation: per-event per-frame sums via coverage prefix sums,
    then the shared `_elong_evidence`. Identical decision logic to the scalar
    path; only the summation is vectorised."""
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
        od = overlaps_df.with_columns([
            pl.col("event_id").cast(pl.UInt64), pl.col("other_event_id").cast(pl.UInt64),
        ]).filter(pl.col("event_id").is_in(elong["event_id"].cast(pl.UInt64))).sort(
            ["event_id", "overlap_start"])
        if not od.is_empty():
            # Merge each event's overlap intervals (union) so contended coverage is
            # not double-counted where competitors overlap each other. Competitor
            # frames are collected from the raw (unmerged) rows.
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
                        merged.append((cs, ce)); cs, ce = s, e
                merged.append((cs, ce))
                for s, e in merged:
                    m_start.append(s); m_end.append(e); m_idx.append(i)
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
            float(tot_w[i]), int(ncov_w[i]), clean_tot, clean_inf,
            float(cont_tot[i]), list(map(float, cont_fr[i])), int(a_e[i]),
            comp_frames.get(i, {}), int(cont_nt[i]), int(ends[i] - starts[i]), thr)
    return out


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
    """Vectorised scorer: elongation via prefix sums, init/term scalar (small
    flank windows), junction from precomputed support. Reproduces score_events.

    cov_df: (pos, count). overlaps_df: (event_id, other_event_id, overlap_start,
    overlap_end, comp_phase).
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
    return pl.from_dicts(rows, schema=_RECORD_SCHEMA) if rows else pl.DataFrame(schema=_RECORD_SCHEMA)


from TranslonScorer.io.store import persist_scores  # noqa: F401 — re-exported


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

    events: rows with event_id, type, start, end, strand, phase.
    overlaps: {elongation event_id → [(o_start,o_end,competitor_phase,competitor_id)]}.
    junction_support: {junction event_id → spanning count}.
    Scalar reference implementation; the vectorised path must reproduce this.
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
            raw = score_elongation_event(r["start"], r["end"], r["phase"], r["strand"],
                                         coverage, overlaps.get(eid), thr=thr)
        elif t == "junction":
            raw = score_junction_event(junction_support.get(eid, {}), thr=thr)
        else:
            continue
        rows.append(event_record(eid, t, group, tier, raw, thr.version))
    return pl.from_dicts(rows, schema=_RECORD_SCHEMA) if rows else pl.DataFrame(schema=_RECORD_SCHEMA)

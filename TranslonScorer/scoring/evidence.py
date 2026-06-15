"""Evidence builders and eligibility/call decisions.

These are the pure, measurement-free decision functions: they turn pre-computed
frame sums and per-position coverage tallies into structured evidence dicts, then
apply threshold comparisons to derive eligibility and call.  Thresholds live in
ScoreThresholds; changing a threshold lets you re-derive calls without
re-measuring coverage.

Public API
----------
_RECORD_SCHEMA      — Polars schema for the long-form score table
_jsonable           — JSON-safe value coercion (inf → string)
_decide_step        — coverage step evidence → (eligibility, call)
_elong_evidence     — frame-sum tallies → elongation evidence dict
event_record        — raw evidence dict → one long-form record dict
"""

from __future__ import annotations

import json as _json
import math as _math
from typing import Dict, Optional, Tuple

import polars as pl

from TranslonScorer.model import ScoreThresholds

# ---------------------------------------------------------------------------
# Long-form record schema (mirrors ScoreRecord dataclass in model.py)
# ---------------------------------------------------------------------------

_RECORD_SCHEMA = {
    "event_id": pl.UInt64,
    "aspect": pl.Utf8,
    "group": pl.Utf8,
    "tier": pl.Utf8,
    "n_reads": pl.Float64,
    "metric": pl.Float64,
    "metric_name": pl.Utf8,
    "eligibility": pl.Utf8,
    "call": pl.Utf8,
    "evidence": pl.Utf8,
    "thresholds_version": pl.Utf8,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _jsonable(v):
    """Coerce a value to JSON-safe form (replaces inf/nan with strings)."""
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, float) and (_math.isinf(v) or _math.isnan(v)):
        return "inf" if v > 0 else "-inf"
    return v


# ---------------------------------------------------------------------------
# Eligibility / call
# ---------------------------------------------------------------------------


def _decide_step(
    ev: dict,
    *,
    rise_thr: float,
    thr: ScoreThresholds,
) -> Tuple[str, Optional[str]]:
    """Map step-scorer evidence → (eligibility, call).

    Peakiness and stability remain in evidence as review flags but do NOT
    hard-gate the call — that wrongly overrode strong, clean starts/stops.
    """
    if ev["n_reads"] < thr.min_reads:
        return "INSUFFICIENT", None
    r = ev["consensus_rise"]
    if r >= rise_thr:
        return "ELIGIBLE", "SUPPORTED"
    if r <= 0:
        return "ELIGIBLE", "UNSUPPORTED"
    return "ELIGIBLE", "AMBIGUOUS"  # borderline: 0 < rise < threshold


def _elong_evidence(
    n: float,
    covered: int,
    clean_tot: float,
    clean_inf: float,
    cont_tot: float,
    cont_by_frame,  # list/array of length 3
    a_e: int,
    comp_frame: Dict[int, int],
    contended_nt: int,
    span_nt: int,
    thr: ScoreThresholds,
) -> dict:
    """Turn per-event frame sums into evidence + eligibility/call.

    Shared by the scalar (score_elongation_event) and vectorised
    (score_elongation_batch) scorers so decision logic is identical.
    """
    cont_e = cont_by_frame[a_e]
    clean_in_frame = (clean_inf / clean_tot) if clean_tot else None
    identifiability = (cont_e / cont_tot) if cont_tot else None
    overall_in_frame = ((clean_inf + cont_e) / n) if n else 0.0
    breadth = (covered / span_nt) if span_nt else 0.0

    competitor_share: Dict[int, float] = {
        cid: ((cont_by_frame[af] / cont_tot) if cont_tot else 0.0) for cid, af in comp_frame.items()
    }
    noise_frame = {0, 1, 2} - {a_e} - set(comp_frame.values())
    noise_share = (sum(cont_by_frame[f] for f in noise_frame) / cont_tot) if cont_tot else 0.0

    effective_in_frame = clean_in_frame if clean_in_frame is not None else overall_in_frame
    if n < thr.min_reads:
        eligibility, call = "INSUFFICIENT", None
    elif identifiability is not None and identifiability < thr.elong_identifiability:
        eligibility, call = "ELIGIBLE", "AMBIGUOUS"
    elif effective_in_frame >= thr.elong_in_frame and breadth >= thr.elong_breadth:
        eligibility, call = "ELIGIBLE", "SUPPORTED"
    else:
        eligibility, call = "ELIGIBLE", "UNSUPPORTED"

    return {
        "n_reads": n,
        "covered_nt": covered,
        "metric": effective_in_frame,
        "metric_name": "elong_in_frame",
        "overall_in_frame": overall_in_frame,
        "breadth": breadth,
        "span_nt": span_nt,
        "clean_in_frame": clean_in_frame,
        "contended_nt": contended_nt,
        "identifiability": identifiability,
        "competitor_share": competitor_share,
        "noise_share": noise_share,
        "eligibility": eligibility,
        "call": call,
    }


# ---------------------------------------------------------------------------
# Record serialisation
# ---------------------------------------------------------------------------


def event_record(
    event_id: int,
    aspect: str,
    group: str,
    tier: str,
    raw: dict,
    thr_version: str = "v0",
) -> dict:
    """Serialise a raw evidence dict into one long-form record dict."""
    drop = {"eligibility", "call", "metric", "metric_name", "n_reads"}
    evidence = {k: _jsonable(v) for k, v in raw.items() if k not in drop}
    m = raw.get("metric")
    return {
        "event_id": int(event_id),
        "aspect": aspect,
        "group": group,
        "tier": tier,
        "n_reads": float(raw.get("n_reads", 0.0) or 0.0),
        "metric": (None if m is None else float(m)),
        "metric_name": raw.get("metric_name"),
        "eligibility": raw["eligibility"],
        "call": raw.get("call"),
        "evidence": _json.dumps(evidence),
        "thresholds_version": thr_version,
    }

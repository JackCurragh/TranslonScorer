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
# Long-form record schema (the single source of truth for a score row)
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
    # Region-context annotation, NOT the BAM-NH-tag SupportsMappability/
    # mappability_ledger concept (coverage/base.py) — distinct prefix
    # deliberately avoids conflating the two. Never affects eligibility/call
    # (same "review flag, not a gate" precedent as flank_peakiness/stability).
    "map_track_mean": pl.Float64,
    "map_track_low": pl.Boolean,
    # Elongation only; null on init/term/junction rows. First-class (not
    # evidence-JSON-only) so translon-level composition can weight by codon
    # count instead of read count — see report.compose_block_detail and
    # _compose_per_translon's elongation_cif_approx.
    "cif": pl.Float64,
    "n_codons": pl.Int64,
    # docs/significance_testing_plan.md §6 -- fraction of this aspect's lens
    # battery that agreed (scoring/attribution.py's compose_confidence), NOT
    # a count of "how many tests of any kind passed". First-class (not
    # evidence-JSON-only) for the same reason cif/n_codons are: composition
    # needs to weight/aggregate by it (consequential.py's tier_confidence),
    # not just display it. Null where no battery exists yet for that aspect
    # (junction) or no lens could be evaluated for this event.
    "confidence": pl.Float64,
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

    The borderline band (0 < rise < rise_thr) gets one more chance: if the
    frame-resolved periodicity significance test (`_periodicity_significance`
    in aspects.py) is significant at enough flank lengths, the call is
    upgraded to SUPPORTED — distinguishing "genuinely marginal" from
    "borderline depth but statistically real periodicity." This can only
    move AMBIGUOUS -> SUPPORTED: an event already SUPPORTED or UNSUPPORTED
    by depth alone is untouched, so no existing call can flip. Sets
    `ev["periodicity_resolved_ambiguous"]` (only reached, hence only present,
    for borderline events) so the upgrade is auditable in the evidence JSON.
    """
    if ev["n_reads"] < thr.min_reads:
        return "INSUFFICIENT", None
    r = ev["consensus_rise"]
    if r >= rise_thr:
        return "ELIGIBLE", "SUPPORTED"
    if r <= 0:
        return "ELIGIBLE", "UNSUPPORTED"
    # borderline: 0 < rise < threshold
    axes_by_flank = ev.get("boundary_axes_by_flank") or {}
    p_vals = [
        a["periodicity_p"] for a in axes_by_flank.values() if a.get("periodicity_p") is not None
    ]
    resolved = bool(p_vals) and (
        sum(1 for p in p_vals if p < thr.periodicity_significance_alpha) / len(p_vals)
        >= thr.periodicity_min_agree_frac
    )
    ev["periodicity_resolved_ambiguous"] = resolved
    if resolved:
        return "ELIGIBLE", "SUPPORTED"
    return "ELIGIBLE", "AMBIGUOUS"


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
    cif_value: Optional[float] = None,
    n_codons: Optional[int] = None,
    frame_chisq_p: Optional[float] = None,
    cif_significance: Optional[dict] = None,
    cif_contiguity: Optional[dict] = None,
    body_uniformity_p: Optional[float] = None,
) -> dict:
    """Turn per-event frame sums into evidence + eligibility/call.

    Shared by the scalar (score_elongation_event) and vectorised
    (score_elongation_batch) scorers so decision logic is identical.

    Two of the Chothani et al. (2022) signature scores ride along here as
    evidence rather than as headline metrics, because the record schema
    carries exactly one metric per event and `elong_in_frame` already owns it:

    * ``overall_in_frame`` IS their PIF — frame-0 signal over total signal,
      across every position in the span. It is not a separate field because it
      would be a duplicate one; see docs and the CHANGELOG.
    * ``cif`` is their CIF, supplied by the caller (it needs the per-codon
      vector, which neither the prefix-sum kernel nor this function has).
      None when the span is not a whole number of codons or no vector was
      built. ``n_codons`` is the codon count that CIF was actually computed
      over (the caller's trimmed vector length // 3, not ``span_nt // 3`` —
      they differ whenever ``orf_signal_vector`` trims a partial codon at the
      phase-offset boundary), promoted alongside `cif` so translon-level
      composition can weight by codon count rather than read count.

    Note the headline ``metric`` stays ``clean_in_frame`` where available —
    PIF's uncontended-only sibling. The two differ exactly where events
    overlap in different frames, which is the case PIF cannot express.
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
        "overall_in_frame": overall_in_frame,  # == Chothani PIF
        "cif": cif_value,
        "n_codons": n_codons,
        # docs/significance_testing_plan.md §2 "Level" lens -- chi-square
        # goodness-of-fit of the whole-span 3-way frame tally against
        # uniform 1/3. Evidence-only (does not gate eligibility/call here),
        # same precedent as periodicity_p on the boundary axes.
        "frame_chisq_p": frame_chisq_p,
        # docs/significance_testing_plan.md §3 -- per-codon significance
        # (low-depth-aware companion to `cif`) and its contiguity/run-length
        # summary. Both evidence-only dicts, JSON-only (not schema columns):
        # exploratory, not yet something composition weights by.
        "cif_significance": cif_significance,
        "cif_contiguity": cif_contiguity,
        # docs/significance_testing_plan.md §2 "Uniformity" lens -- large p
        # means frame-0 dominance is consistent across the body's two
        # halves, small p flags a local patch. Evidence-only.
        "body_uniformity_p": body_uniformity_p,
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
    drop = {
        "eligibility",
        "call",
        "metric",
        "metric_name",
        "n_reads",
        "map_track_mean",
        "map_track_low",
        "cif",
        "n_codons",
        "confidence",
    }
    evidence = {k: _jsonable(v) for k, v in raw.items() if k not in drop}
    m = raw.get("metric")
    mtm = raw.get("map_track_mean")
    cif_val = raw.get("cif")
    n_codons = raw.get("n_codons")
    confidence = raw.get("confidence")
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
        "map_track_mean": (None if mtm is None else float(mtm)),
        "map_track_low": raw.get("map_track_low"),
        "cif": (None if cif_val is None else float(cif_val)),
        "n_codons": (None if n_codons is None else int(n_codons)),
        "confidence": (None if confidence is None else float(confidence)),
    }

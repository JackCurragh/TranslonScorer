"""Per-aspect event scorers: initiation, elongation, termination, junction.

Each scorer receives pre-extracted coverage and returns a raw evidence dict.
Decision logic (eligibility/call) is applied inside each scorer via helpers
in evidence.py so that the measurement and decision layers are clearly
separated: changing thresholds re-derives calls without re-measuring coverage.

Public API
----------
abs_frame               — phase × strand → genomic frame carrying in-frame signal
score_initiation_event  — P-site coverage step UP at start codon
score_termination_event — A-site coverage step DOWN at stop codon
score_elongation_event  — frame-resolved in-frame fraction + identifiability
score_junction_event    — splice-junction confident-spanning count + psi
"""

from __future__ import annotations

import math
import statistics
from typing import Dict, List, Optional, Sequence, Tuple

from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.evidence import _decide_step, _elong_evidence


# ---------------------------------------------------------------------------
# Frame utility
# ---------------------------------------------------------------------------


def abs_frame(phase: int, strand: int) -> int:
    """Genomic frame (p % 3) carrying this event's in-frame (codon-0) signal."""
    return (-phase) % 3 if strand > 0 else phase % 3


# ---------------------------------------------------------------------------
# Initiation / termination: smoothed, flank-robust step scoring
# ---------------------------------------------------------------------------
# Each flank is binned into codons (3 nt, periodicity-aware) and the MEDIAN
# bin is taken as the robust level.  The step is measured at several flank
# lengths; the consensus (median of log2FCs) is the headline metric.


def _flank_positions(pos: int, strand: int, length: int, *, body_side: bool) -> range:
    """Genomic positions for a flank of `length` nt on the body or outer side of
    `pos`, in transcript orientation. body_side=True → into the ORF."""
    if strand > 0:
        return range(pos, pos + length) if body_side else range(pos - length, pos)
    return range(pos - length + 1, pos + 1) if body_side else range(pos + 1, pos + length + 1)


def _codon_levels(
    coverage: Dict[int, float],
    positions: Sequence[int],
) -> Tuple[float, float, float]:
    """(median codon-bin sum, max codon-bin sum, total) over a window.

    Median is robust to a single spike; codon bins fold in 3-nt periodicity.
    """
    pos = list(positions)
    bins = [
        sum(coverage.get(p, 0.0) for p in pos[i : i + 3]) for i in range(0, max(len(pos) - 2, 0), 3)
    ]
    total = sum(coverage.get(p, 0.0) for p in pos)
    if not bins:
        return total, total, total
    return statistics.median(bins), max(bins), total


def _step_score(
    pos: int,
    strand: int,
    coverage: Dict[int, float],
    flanks: Sequence[int],
    min_reads: float,
    alpha: float,
) -> dict:
    """Robust step (body_level - outer_level) across `pos`, over several flank
    lengths. Metric = log2 fold-change (body vs outer) with pseudocount α."""
    rises: Dict[int, float] = {}
    for L in flanks:
        body_med, _, _ = _codon_levels(coverage, _flank_positions(pos, strand, L, body_side=True))
        out_med, _, _ = _codon_levels(coverage, _flank_positions(pos, strand, L, body_side=False))
        rises[L] = math.log2((body_med + alpha) / (out_med + alpha))
    vals = list(rises.values())
    consensus = statistics.median(vals)
    stability = (max(vals) - min(vals)) if len(vals) > 1 else 0.0
    refL = sorted(flanks)[len(flanks) // 2]
    out_med, out_max, _ = _codon_levels(
        coverage, _flank_positions(pos, strand, refL, body_side=False)
    )
    flank_peakiness = (out_max / out_med) if out_med > 0 else (float("inf") if out_max > 0 else 0.0)
    n_reads = _codon_levels(coverage, _flank_positions(pos, strand, min(flanks), body_side=True))[2]
    return {
        "rise_by_flank": rises,
        "consensus_rise": consensus,
        "stability": stability,
        "flank_peakiness": flank_peakiness,
        "n_reads": n_reads,
    }


_DEFAULT_THR = ScoreThresholds()


def score_initiation_event(
    start_pos: int,
    strand: int,
    coverage: Dict[int, float],
    *,
    flanks: Sequence[int] = (9, 18, 30, 60),
    thr: ScoreThresholds = _DEFAULT_THR,
) -> dict:
    """Initiation = smoothed, flank-robust step UP at the start codon.

    Headline metric: consensus log2 fold-change (body / outer flank).
    """
    s = _step_score(start_pos, strand, coverage, flanks, thr.min_reads, thr.step_pseudocount)
    s["metric"] = s["consensus_rise"]
    s["metric_name"] = "init_rise"
    s["eligibility"], s["call"] = _decide_step(s, rise_thr=thr.init_rise, thr=thr)
    return s


def score_termination_event(
    term_pos: int,
    strand: int,
    coverage: Dict[int, float],
    *,
    flanks: Sequence[int] = (9, 18, 30, 60),
    thr: ScoreThresholds = _DEFAULT_THR,
) -> dict:
    """Termination = step DOWN after the stop codon (body before > UTR after).

    Headline metric: consensus log2 fold-change (body / UTR) — positive = drop.
    """
    alpha = thr.step_pseudocount
    rises: Dict[int, float] = {}
    for L in flanks:
        body_med, _, _ = _codon_levels(
            coverage, _flank_positions(term_pos, strand, L, body_side=False)
        )
        utr_med, _, _ = _codon_levels(
            coverage, _flank_positions(term_pos, strand, L, body_side=True)
        )
        rises[L] = math.log2((body_med + alpha) / (utr_med + alpha))
    vals = list(rises.values())
    refL = sorted(flanks)[len(flanks) // 2]
    utr_med, utr_max, _ = _codon_levels(
        coverage, _flank_positions(term_pos, strand, refL, body_side=True)
    )
    s = {
        "rise_by_flank": rises,
        "consensus_rise": statistics.median(vals),
        "stability": (max(vals) - min(vals)) if len(vals) > 1 else 0.0,
        "flank_peakiness": (
            (utr_max / utr_med) if utr_med > 0 else (float("inf") if utr_max > 0 else 0.0)
        ),
        "n_reads": _codon_levels(
            coverage, _flank_positions(term_pos, strand, min(flanks), body_side=False)
        )[2],
    }
    s["metric"] = s["consensus_rise"]
    s["metric_name"] = "term_drop"
    s["eligibility"], s["call"] = _decide_step(s, rise_thr=thr.term_drop, thr=thr)
    return s


# ---------------------------------------------------------------------------
# Elongation
# ---------------------------------------------------------------------------


def score_elongation_event(
    start: int,
    end: int,
    phase: int,
    strand: int,
    coverage: Dict[int, float],
    overlaps: Optional[List[Tuple[int, int, int, int]]] = None,
    *,
    thr: ScoreThresholds = _DEFAULT_THR,
) -> dict:
    """Score one elongation event.

    overlaps: [(overlap_start, overlap_end, competitor_phase, competitor_event_id)]
    for every competing event sharing genomic span in a different frame.

    Returns clean vs contended in-frame fractions, per-overlap identifiability,
    and per-competitor signal shares via _elong_evidence.
    """
    overlaps = overlaps or []
    a_e = abs_frame(phase, strand)

    contended: Dict[int, None] = {}
    comp_frame: Dict[int, int] = {}
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

    return _elong_evidence(
        n,
        covered,
        clean_tot,
        clean_inf,
        cont_tot,
        cont_by_frame,
        a_e,
        comp_frame,
        len(contended),
        end - start,
        thr,
    )


# ---------------------------------------------------------------------------
# Junction
# ---------------------------------------------------------------------------


def score_junction_event(
    support: dict,
    *,
    thr: ScoreThresholds = _DEFAULT_THR,
) -> dict:
    """Splice support with an overhang-confidence (identifiability) layer.

    support keys: span_conf, span_short, unspliced (weighted read counts).
      span_conf  — intron matches with ≥6 nt aligned both sides (confident)
      span_short — matches with short overhang (~ambiguous)
      unspliced  — aligned straight through the donor (intron retention)

    Headline metric: confident spanning count. psi = spliced-in ratio.
    """
    conf = float(support.get("span_conf", 0.0))
    short = float(support.get("span_short", 0.0))
    unspl = float(support.get("unspliced", 0.0))
    spanning = conf + short
    psi = (spanning / (spanning + unspl)) if (spanning + unspl) > 0 else None
    if conf >= thr.junc_min_spanning:
        elig, call = "ELIGIBLE", "SUPPORTED"
    elif spanning >= thr.junc_min_spanning:
        elig, call = "ELIGIBLE", "AMBIGUOUS"
    else:
        elig, call = "INSUFFICIENT", None
    return {
        "n_reads": spanning,
        "metric": conf,
        "metric_name": "junc_confident_spanning",
        "confident_spanning": conf,
        "short_spanning": short,
        "unspliced": unspl,
        "psi": psi,
        "eligibility": elig,
        "call": call,
    }

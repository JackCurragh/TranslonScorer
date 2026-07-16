"""Periodicity and frame-dominance QC — pure transforms.

These functions operate on pre-computed frame distributions (dicts/DataFrames)
and do not perform I/O, BAM scanning, or offset inference.  The BAM-scanning
and per-partition accumulation lives in TranslonScorer/matrix/qc.py (I/O shell).

Public API
----------
_compute_periodicity_from_frames  — {frame: count} → periodicity score + fracs
_ribometric_frame_scores          — {length: {frame: count}} → weighted score
_nudge_to_frame0                  — nudge per-length offsets ±1 (RiboMetric logic)
_assign_frames_sweep              — sweep-line CDS frame assignment for reads
frame_dominance_matrix            — periodicity_qc output → read_length × sample matrix
_empty_periodicity_schema         — empty DataFrame with the periodicity_qc schema
"""

from __future__ import annotations

import heapq
import json
import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import polars as pl

# ---------------------------------------------------------------------------
# Periodicity score (RiboMetric formula)
# ---------------------------------------------------------------------------


def _compute_periodicity_from_frames(frames: Dict[int, float]) -> Dict:
    """Entropy-reduction periodicity score — exact RiboMetric formula (sqrt of info ratio).

    frames: {0: count, 1: count, 2: count}
    Returns {periodicity_score, f0, f1, f2}.
    """
    total = sum(frames.values())
    if total == 0:
        return {"periodicity_score": 0.0, "f0": 0.0, "f1": 0.0, "f2": 0.0}
    pseudocount = 1e-100
    probs = [(frames.get(f, 0.0) + pseudocount) / (total + pseudocount) for f in range(3)]
    max_ent = math.log2(3)
    entropy = -sum(p * math.log2(p) for p in probs if p > 0)
    score = float(math.sqrt(max(0.0, (max_ent - entropy) / max_ent)))
    return {"periodicity_score": score, "f0": probs[0], "f1": probs[1], "f2": probs[2]}


def _ribometric_frame_scores(
    rfd: Dict[int, Dict[int, float]],
    min_count_frac: float = 0.05,
) -> Dict:
    """Apply RiboMetric's read_frame_information_content + information_metric_cutoff.

    rfd = {length: {0: count, 1: count, 2: count}} — the read_frame_distribution.
    Returns {periodicity_score, f0, f1, f2, n_reads, per_length: {length: score}}.
    """
    pseudocount = 1e-100
    log2_3 = math.log2(3)
    total_reads = sum(sum(fd.values()) for fd in rfd.values())

    per_length: Dict[int, Tuple[float, int]] = {}
    for length, fd in rfd.items():
        n = sum(fd.values())
        probs = [(fd.get(f, 0.0) + pseudocount) / (n + pseudocount) for f in range(3)]
        entropy = -sum(p * math.log2(p) for p in probs if p > 0)
        per_length[length] = (math.sqrt(max(0.0, (log2_3 - entropy) / log2_3)), int(n))

    valid = {l: (s, n) for l, (s, n) in per_length.items() if n > total_reads * min_count_frac}
    if valid:
        w = sum(n for _, n in valid.values())
        global_score = sum(s * n for s, n in valid.values()) / w
    else:
        global_score = 0.0

    agg: Dict[int, float] = {0: 0.0, 1: 0.0, 2: 0.0}
    for fd in rfd.values():
        for f, c in fd.items():
            agg[f] += c
    t = sum(agg.values()) or 1.0

    return {
        "periodicity_score": global_score,
        "f0": agg[0] / t,
        "f1": agg[1] / t,
        "f2": agg[2] / t,
        "n_reads": int(total_reads),
        "per_length": {l: s for l, (s, _) in per_length.items()},
    }


# ---------------------------------------------------------------------------
# Offset nudging
# ---------------------------------------------------------------------------


def _nudge_to_frame0(
    length_frame_dist: Dict[int, Dict[int, float]],
    base_offsets: Dict[int, int],
    *,
    min_reads: int = 50,
    min_fraction: float = 0.60,
) -> Dict[int, int]:
    """Nudge per-length offsets ±1 to put dominant frame back to 0.

    Matches RiboMetric._frame_calibrated_offsets() logic.
    """
    out = dict(base_offsets)
    for length, fd in length_frame_dist.items():
        total = sum(fd.values())
        if total < min_reads:
            continue
        dom = max(fd, key=fd.get)
        if dom == 0 or fd[dom] / total < min_fraction:
            continue
        shift = -1 if dom == 1 else 1
        if int(length) in out:
            out[int(length)] += shift
    return out


# ---------------------------------------------------------------------------
# Frame assignment (sweep-line, pure)
# ---------------------------------------------------------------------------


def _assign_frames_sweep(
    read_ids: List[int],
    asites: "np.ndarray",
    strand: str,
    intervals: List[Tuple[int, int, int]],
) -> Dict[int, int]:
    """Sweep-line O((m+n) log n) CDS frame assignment.

    Iterates reads sorted by A-site position and maintains an active heap of
    intervals sorted by stop. Returns {read_id: frame} for CDS-overlapping reads.

    intervals: [(start, end_exclusive, phase)]
    """
    if not intervals or len(read_ids) == 0:
        return {}

    order = np.argsort(asites)
    sorted_asites = asites[order].tolist()
    sorted_rids = [read_ids[i] for i in order]

    result: Dict[int, int] = {}
    active: list = []  # min-heap: (end_exclusive, phase)
    iv_ptr = 0
    n_ivs = len(intervals)

    for asite, rid in zip(sorted_asites, sorted_rids):
        while iv_ptr < n_ivs and intervals[iv_ptr][0] <= asite:
            es, ee, phase = intervals[iv_ptr]
            heapq.heappush(active, (ee, phase))
            iv_ptr += 1
        while active and active[0][0] <= asite:
            heapq.heappop(active)
        if active:
            _, phase_top = active[0]
            if strand == "+":
                result[rid] = (phase_top + asite) % 3
            else:
                result[rid] = (phase_top - asite) % 3

    return result


# ---------------------------------------------------------------------------
# Frame-dominance matrix (pure pivot over periodicity_qc output)
# ---------------------------------------------------------------------------


def frame_dominance_matrix(
    perio_df: pl.DataFrame,
    *,
    metric: str = "dominance",
    min_reads: int = 50,
) -> pl.DataFrame:
    """Pivot per-sample periodicity output into a read_length × sample matrix.

    perio_df : output of periodicity_qc (must carry the per-length JSON
               columns read_frame_distribution and either per_length_dominance
               or per_length_periodicity).
    metric   : "dominance" → dominant-frame fraction (max frame proportion,
               in [1/3, 1]); "periodicity" → sqrt entropy-reduction score (in [0, 1]).

    Returns a wide DataFrame: one row per read_length (sorted), one column per
    sample, cells holding the chosen metric.  Cells with fewer than min_reads
    reads at that (sample, length) are left null.
    """
    if perio_df.is_empty() or "read_frame_distribution" not in perio_df.columns:
        return pl.DataFrame(schema={"read_length": pl.Int64})

    src_col = "per_length_dominance" if metric == "dominance" else "per_length_periodicity"
    long_rows: list = []
    for row in perio_df.iter_rows(named=True):
        sample = str(row["sample_id"])
        rfd = json.loads(row.get("read_frame_distribution") or "{}")
        vals = json.loads(row.get(src_col) or "{}")
        for length_str, counts in rfd.items():
            n = sum(counts)
            if n < min_reads:
                continue
            long_rows.append(
                {
                    "read_length": int(length_str),
                    "sample_id": sample,
                    "value": float(vals.get(length_str, 0.0)),
                }
            )

    if not long_rows:
        return pl.DataFrame(schema={"read_length": pl.Int64})

    long = pl.from_dicts(long_rows)
    wide = long.pivot(values="value", index="read_length", on="sample_id").sort("read_length")
    return wide


# ---------------------------------------------------------------------------
# Length distribution metrics  (RiboMetric equivalents, pure transforms)
# ---------------------------------------------------------------------------

def _weighted_percentile(lengths: np.ndarray, counts: np.ndarray, q: float) -> float:
    """Weighted percentile via cumulative sum of sorted lengths."""
    order = np.argsort(lengths)
    sl, sc = lengths[order], counts[order]
    cumw = np.cumsum(sc)
    total = cumw[-1]
    idx = np.searchsorted(cumw, q * total)
    return float(sl[min(idx, len(sl) - 1)])


def _length_distribution_metrics(length_counts: Dict[int, float]) -> Dict[str, float]:
    """All RiboMetric read-length-distribution metrics from {length: count}.

    Returns:
        rld_IQR_metric        — 1 − IQR/(P90−P10), higher = tighter peak
        rld_CV_metric         — 1/(1+CV), higher = tighter peak
        rld_normality_metric  — 1 − normaltest p-value (non-normal = good for RPFs)
        rld_max_prop_metric   — fraction of reads at the single peak length
        rld_bimodality        — 1/(1+bimodality_coeff), higher = unimodal
        peak_length           — length with most reads
        mean_length           — weighted mean
        rpf_28_32_prop        — fraction in 28–32 nt window
        total_reads           — sum of all counts
    """
    if not length_counts:
        return {k: 0.0 for k in (
            "rld_IQR_metric", "rld_CV_metric", "rld_normality_metric",
            "rld_max_prop_metric", "rld_bimodality",
            "peak_length", "mean_length", "rpf_28_32_prop", "total_reads",
        )}

    lengths = np.array(sorted(length_counts), dtype=float)
    counts = np.array([length_counts[int(l)] for l in lengths], dtype=float)
    total = float(counts.sum())

    mean_l = float(np.average(lengths, weights=counts))
    std_l = float(np.sqrt(np.average((lengths - mean_l) ** 2, weights=counts)))
    cv = std_l / mean_l if mean_l > 0 else 0.0

    peak_l = int(lengths[int(np.argmax(counts))])
    rpf_mask = (lengths >= 28) & (lengths <= 32)
    rpf_prop = float(counts[rpf_mask].sum() / total) if total > 0 else 0.0

    # IQR metric
    p10 = _weighted_percentile(lengths, counts, 0.10)
    p25 = _weighted_percentile(lengths, counts, 0.25)
    p75 = _weighted_percentile(lengths, counts, 0.75)
    p90 = _weighted_percentile(lengths, counts, 0.90)
    range_10_90 = p90 - p10
    iqr_metric = float(1.0 - (p75 - p25) / range_10_90) if range_10_90 > 0 else 1.0

    # Max-proportion metric
    max_prop = float(counts.max() / total) if total > 0 else 0.0

    # Bimodality coefficient (Sarle's B)
    try:
        from scipy.stats import skew as _skew, kurtosis as _kurt, normaltest as _ntest
        expanded = np.repeat(lengths.astype(int), counts.astype(int).clip(0))
        if len(expanded) >= 20:
            n = len(expanded)
            sk = float(_skew(expanded))
            ku = float(_kurt(expanded))
            denom = ku + 3.0 * ((n - 1) ** 2) / ((n - 2) * (n - 3)) if n > 3 else 1.0
            bm = (sk ** 2 + 1) / denom if denom != 0 else 0.0
            bimodality = float(1.0 / (1.0 + max(bm, 0.0)))
            norm_pval = float(_ntest(expanded).pvalue)
            normality_metric = float(max(0.0, min(1.0, 1.0 - norm_pval)))
        else:
            bimodality = 1.0
            normality_metric = 0.0
    except Exception:
        bimodality = 1.0
        normality_metric = 0.0

    return {
        "rld_IQR_metric": iqr_metric,
        "rld_CV_metric": float(1.0 / (1.0 + cv)),
        "rld_normality_metric": normality_metric,
        "rld_max_prop_metric": max_prop,
        "rld_bimodality": bimodality,
        "peak_length": float(peak_l),
        "mean_length": mean_l,
        "rpf_28_32_prop": rpf_prop,
        "total_reads": total,
    }


# ---------------------------------------------------------------------------
# Ligation bias metrics  (RiboMetric equivalents, pure transforms)
# ---------------------------------------------------------------------------

_ALL_DINUCS = [a + b for a in "ACGT" for b in "ACGT"]
_UNIFORM_DINUC = {d: 1.0 / 16.0 for d in _ALL_DINUCS}


def _kl_divergence(obs: Dict[str, float], exp: Dict[str, float]) -> float:
    """KL divergence D(obs||exp) in bits, clipped to ≥ 0."""
    kl = 0.0
    for d, p in obs.items():
        q = exp.get(d, 0.0)
        if p > 0 and q > 0:
            kl += p * math.log2(p / q)
    return max(0.0, kl)


def _ligation_bias_metrics(
    five_prime: Dict[str, float],
    three_prime: Dict[str, float],
    background_5p: Optional[Dict[str, float]] = None,
    background_3p: Optional[Dict[str, float]] = None,
) -> Dict[str, float]:
    """Per-sample ligation bias metrics.

    Normalises observed dinucleotide frequencies against background (defaults to
    uniform 1/16 when no transcriptome background is provided).

    Returns:
        ligation_bias_KL_5p      — KL divergence 5' vs background (bits)
        ligation_bias_KL_3p      — KL divergence 3' vs background
        ligation_bias_score_5p   — 1/(1+KL), 1.0 = no bias
        ligation_bias_score_3p   — 1/(1+KL), 1.0 = no bias
        ligation_bias_max_abs_5p — max |obs − exp| at 5' end
        ligation_bias_max_abs_3p — max |obs − exp| at 3' end
    """
    bg5 = background_5p or _UNIFORM_DINUC
    bg3 = background_3p or _UNIFORM_DINUC

    def _normalise(counts: Dict[str, float]) -> Dict[str, float]:
        total = sum(counts.values()) or 1.0
        return {k: v / total for k, v in counts.items() if "N" not in k}

    obs5 = _normalise(five_prime)
    obs3 = _normalise(three_prime)

    kl5 = _kl_divergence(obs5, bg5)
    kl3 = _kl_divergence(obs3, bg3)

    max_abs5 = max((abs(obs5.get(d, 0.0) - bg5.get(d, 0.0)) for d in _ALL_DINUCS), default=0.0)
    max_abs3 = max((abs(obs3.get(d, 0.0) - bg3.get(d, 0.0)) for d in _ALL_DINUCS), default=0.0)

    return {
        "ligation_bias_KL_5p": kl5,
        "ligation_bias_KL_3p": kl3,
        "ligation_bias_score_5p": 1.0 / (1.0 + kl5),
        "ligation_bias_score_3p": 1.0 / (1.0 + kl3),
        "ligation_bias_max_abs_5p": max_abs5,
        "ligation_bias_max_abs_3p": max_abs3,
    }


# ---------------------------------------------------------------------------
# Derived / composite metrics
# ---------------------------------------------------------------------------

def _recommend_read_lengths(
    rfd: Dict[int, Dict[int, float]],
    rld: Dict[int, float],
    offsets: Optional[Dict[int, int]] = None,
    min_periodicity: float = 0.5,
    min_read_proportion: float = 0.05,
) -> Dict[str, Any]:
    """RiboMetric recommend_read_lengths — which lengths carry clean periodicity.

    rfd : {length: {frame: count}}
    rld : {length: total_count}
    Returns dict with recommended_lengths, n_recommended, recommended_read_proportion.
    """
    total_reads = sum(rld.values()) or 1
    recommended: List[int] = []
    for length, frames in rfd.items():
        rl = int(length)
        frame_total = sum(frames.values())
        if frame_total == 0:
            continue
        dom_frac = max(frames.values()) / frame_total
        proportion = rld.get(rl, 0) / total_reads
        if dom_frac >= min_periodicity and proportion >= min_read_proportion:
            recommended.append(rl)

    recommended_lengths = sorted(recommended)
    rec_prop = sum(rld.get(rl, 0) for rl in recommended_lengths) / total_reads
    result: Dict[str, Any] = {
        "recommended_lengths": recommended_lengths,
        "n_recommended": len(recommended_lengths),
        "recommended_read_proportion": round(rec_prop, 4),
    }
    if offsets is not None:
        result["recommended_offsets_filtered"] = {
            rl: offsets[rl] for rl in recommended_lengths if rl in offsets
        }
    return result


def _classify_library(
    periodicity: float,
    prop_cds: float,
    min_periodicity: float = 0.4,
    min_cds_proportion: float = 0.5,
) -> str:
    """Simplified RiboMetric classify_library_type (no initiation ratio without metagene)."""
    if periodicity < min_periodicity or prop_cds < min_cds_proportion:
        return "low_quality"
    return "elongation"


# ---------------------------------------------------------------------------
# Empty schema helper
# ---------------------------------------------------------------------------


def _empty_periodicity_schema() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "sample_id": pl.Utf8,
            "periodicity_score": pl.Float64,
            "f0": pl.Float64,
            "f1": pl.Float64,
            "f2": pl.Float64,
            "n_cds_reads": pl.Int64,
            "recommended_offsets": pl.Utf8,
            "read_frame_distribution": pl.Utf8,
            "per_length_periodicity": pl.Utf8,
            "per_length_dominance": pl.Utf8,
        }
    )

"""Periodicity and frame-dominance QC — pure transforms.

These functions operate on pre-computed frame distributions (dicts/DataFrames)
and do not perform I/O, BAM scanning, or offset inference.  The BAM-scanning
and per-partition accumulation lives in TranslonScorer/matrix_qc.py (I/O shell).

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
from typing import Dict, List, Optional, Tuple

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

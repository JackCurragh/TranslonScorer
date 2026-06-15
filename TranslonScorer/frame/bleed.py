from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import polars as pl


def _cds_interior(cds_df: pl.DataFrame, trim_nt: int = 30) -> pl.DataFrame:
    """Return CDS interior intervals per transcript: [start+trim, stop-trim).

    Expects columns: tran_id, tran_start, tran_stop (0-based, half-open).
    Filters out CDS shorter than 2*trim.
    """
    if cds_df.is_empty():
        return cds_df
    # Support both (tran_start,tran_stop) and (start,stop) column names
    start_col = "tran_start" if "tran_start" in cds_df.columns else "start"
    stop_col = "tran_stop" if "tran_stop" in cds_df.columns else "stop"
    df = cds_df.with_columns((pl.col(stop_col) - pl.col(start_col)).alias("cds_len")).filter(
        pl.col("cds_len") >= 2 * trim_nt
    )
    if df.is_empty():
        return pl.DataFrame({"tran_id": [], "start": [], "stop": []})
    return df.select(
        pl.col("tran_id"),
        (pl.col(start_col) + trim_nt).alias("start"),
        (pl.col(stop_col) - trim_nt).alias("stop"),
    )


def learn_confusion(
    profiles: pl.DataFrame,
    cds_df: pl.DataFrame,
    by_length: bool = True,
) -> Dict[Optional[int], np.ndarray]:
    """Learn a 3x3 confusion matrix C per read length (or global) from CDS interiors.

    Profiles schema: tran_id, pos, count[, length]
    Returns dict mapping length (or None) -> matrix C (rows observed frame i, cols true frame j).

    We approximate that true frame is j=0 in CDS interior, then infer the full
    cyclic leakage pattern from observed frames modulo CDS start. A learned
    pattern [0.70, 0.10, 0.20] is therefore kept asymmetric:

        true0 -> observed [0.70, 0.10, 0.20]
        true1 -> observed [0.20, 0.70, 0.10]
        true2 -> observed [0.10, 0.20, 0.70]
    """
    if profiles.is_empty() or cds_df.is_empty():
        return {None: np.eye(3, dtype=float)}

    has_length = ("length" in profiles.columns) and by_length
    # Build joins to annotate CDS-relative frame
    cds_in = _cds_interior(cds_df)
    if cds_in.is_empty():
        return {None: np.eye(3, dtype=float)}

    # Join profiles with CDS interiors per transcript; filter to interior
    p = profiles.join(cds_in, on="tran_id", how="inner")
    p = p.filter((pl.col("pos") >= pl.col("start")) & (pl.col("pos") < pl.col("stop")))
    if p.is_empty():
        return {None: np.eye(3, dtype=float)}

    # Compute frame relative to CDS start
    p = p.with_columns(((pl.col("pos") - pl.col("start")) % 3).alias("frame"))

    if has_length:
        groups = ["length", "frame"]
    else:
        groups = ["frame"]

    agg = p.group_by(groups).agg(pl.col("count").sum().alias("sum")).sort(groups)

    out: Dict[Optional[int], np.ndarray] = {}
    if has_length:
        global_sums = np.array(
            [float(agg.filter(pl.col("frame") == i)["sum"].sum()) for i in range(3)]
        )
        out[None] = _shrink_and_stabilize(_cyclic_confusion_from_offsets(global_sums))
        for L in agg.get_column("length").unique().to_list():
            sub = agg.filter(pl.col("length") == L)
            sums = np.array(
                [float(sub.filter(pl.col("frame") == i)["sum"].sum()) for i in range(3)]
            )
            out[int(L)] = _shrink_and_stabilize(_cyclic_confusion_from_offsets(sums))
    else:
        sums = np.array([float(agg.filter(pl.col("frame") == i)["sum"].sum()) for i in range(3)])
        out[None] = _shrink_and_stabilize(_cyclic_confusion_from_offsets(sums))

    return out


def _cyclic_confusion_from_offsets(offset_counts: np.ndarray) -> np.ndarray:
    """Build C[observed,true] from leakage offsets observed around CDS frame 0."""
    counts = offset_counts.astype(float)
    total = float(counts.sum())
    if total <= 0:
        return np.eye(3, dtype=float)
    probs = counts / total
    return np.array(
        [[probs[(observed - true) % 3] for true in range(3)] for observed in range(3)],
        dtype=float,
    )


def _shrink_and_stabilize(C: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Apply diagonal shrinkage and renormalize columns to sum 1."""
    I = np.eye(3)
    C2 = (1.0 - alpha) * C + alpha * I
    # Normalize columns to 1 (treat columns as true frame j)
    colsum = C2.sum(axis=0, keepdims=True)
    colsum[colsum == 0] = 1.0
    C2 = C2 / colsum
    return C2


def _row_normalize(frame_counts: np.ndarray) -> np.ndarray:
    out = frame_counts.astype(float).copy()
    rowsum = out.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    return out / rowsum


def apply_confusion_counts(
    frame_counts: np.ndarray,
    C: np.ndarray,
    alpha: float = 0.05,
    nnls: bool = False,
    preserve_total: bool = True,
) -> np.ndarray:
    """Invert confusion to estimate adjusted true-frame counts per codon.

    frame_counts: shape (N,3) observed counts in frames 0/1/2
    Returns: shape (N,3) nonnegative adjusted counts. When preserve_total is
    true, each row is rescaled to the original observed count total after
    negative values are clipped.
    """
    if frame_counts.size == 0:
        return frame_counts

    frame_counts = frame_counts.astype(float)

    if nnls:
        try:
            from scipy.optimize import nnls as scipy_nnls

            T = np.vstack([scipy_nnls(C, row)[0] for row in frame_counts])
        except Exception:
            T = _pinv_correct(frame_counts, C, alpha)
    else:
        T = _pinv_correct(frame_counts, C, alpha)

    T[T < 0] = 0.0

    if preserve_total:
        observed_total = frame_counts.sum(axis=1, keepdims=True)
        corrected_total = T.sum(axis=1, keepdims=True)
        nonzero = corrected_total.squeeze() > 0
        if np.any(nonzero):
            T[nonzero] *= observed_total[nonzero] / corrected_total[nonzero]
    return T


def _pinv_correct(frame_counts: np.ndarray, C: np.ndarray, alpha: float) -> np.ndarray:
    # Regularized pseudo-inverse.
    I = np.eye(3)
    Cinvl = np.linalg.pinv(C + alpha * I)
    return frame_counts @ Cinvl.T


def apply_confusion(
    frame_counts: np.ndarray,
    C: np.ndarray,
    alpha: float = 0.05,
    nnls: bool = False,
) -> np.ndarray:
    """Invert confusion to estimate latent frame probabilities per codon.

    frame_counts: shape (N,3) observed counts in frames 0/1/2
    Returns: shape (N,3) nonnegative estimates, row-normalized to sum 1 when nonzero.
    """
    return _row_normalize(apply_confusion_counts(frame_counts, C, alpha=alpha, nnls=nnls))

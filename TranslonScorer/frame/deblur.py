from __future__ import annotations

"""
Positional deblurring for Ribo-seq transcript profiles.

Method 2 from the frame-assignment framework:
  "Observed RPF signal is a blurred version of true ribosome occupancy.
   For each footprint length, learn a blur kernel from aggregate profiles,
   then deconvolve observed signal to recover latent A-site positions."

When length information is absent (BigWig-derived profiles), a single
aggregate kernel is learned across all reads.

Kernel learning
---------------
We build a metagene by aggregating the read count distribution in a
+/-W nucleotide window around every in-frame (frame-0) codon start in
CDS interior regions. The normalised histogram is the empirical blur kernel K.

  K[delta] = P(observed read offset = delta | true ribosome at codon start)

Deconvolution
-------------
Given observed profile O and kernel K, we recover latent A-site profile T
using Richardson-Lucy (RL):

  T^(t+1) = T^(t) * conv(O / conv(T^(t), K), K_flip)

where K_flip = K[::-1].  RL is non-negative and mass-preserving.
"""

from typing import Dict, Optional

import numpy as np
import polars as pl

from .bleed import _cds_interior

_DEFAULT_HALFWIDTH = 9  # +/-9 nt window around each codon start
_RL_ITERS = 20  # Richardson-Lucy iterations


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _delta_kernel(W: int) -> np.ndarray:
    k = np.zeros(2 * W + 1, dtype=np.float64)
    k[W] = 1.0
    return k


def _dense_array(sub: pl.DataFrame, max_pos: int, pad: int = 0) -> np.ndarray:
    arr = np.zeros(max_pos + pad + 1, dtype=np.float64)
    for p, c in zip(sub["pos"].to_list(), sub["count"].to_list()):
        idx = int(p)
        if 0 <= idx < len(arr):
            arr[idx] += float(c)
    return arr


def _accumulate_metagene(
    cds_in: pl.DataFrame,
    prof_by_tran: Dict[str, np.ndarray],
    W: int,
) -> np.ndarray:
    """Sum read windows around all frame-0 codon starts in CDS interior."""
    acc = np.zeros(2 * W + 1, dtype=np.float64)
    for row in cds_in.iter_rows(named=True):
        arr = prof_by_tran.get(str(row["tran_id"]))
        if arr is None:
            continue
        cds_start = int(row["start"])
        cds_stop = int(row["stop"])
        arr_len = len(arr)
        for anchor in range(cds_start, cds_stop, 3):
            lo = anchor - W
            hi = anchor + W + 1
            if lo < 0 or hi > arr_len:
                continue
            acc += arr[lo:hi]
    return acc


def _to_np_kernel(k: object) -> np.ndarray:
    if isinstance(k, np.ndarray):
        return k
    return np.array(k, dtype=np.float64)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def learn_kernel(
    profiles: pl.DataFrame,
    cds_df: pl.DataFrame,
    by_length: bool = True,
    kernel_halfwidth: int = _DEFAULT_HALFWIDTH,
) -> Dict[Optional[int], np.ndarray]:
    """Learn positional blur kernel(s) from the CDS metagene.

    Parameters
    ----------
    profiles : tran_id, pos, count[, length] – transcript-space profiles.
    cds_df   : tran_id, start, stop – CDS interior bounds in *transcript* space.
    by_length: learn per-read-length kernels when a 'length' column exists.
    kernel_halfwidth: window half-width W (nt) for metagene accumulation.

    Returns
    -------
    dict: length (int) or None -> 1-D kernel array (len 2W+1, sums to 1).
    """
    W = kernel_halfwidth
    fallback = _delta_kernel(W)

    cds_in = _cds_interior(cds_df)
    if cds_in.is_empty() or profiles.is_empty():
        return {None: fallback}

    has_length = ("length" in profiles.columns) and by_length

    # Build a CDS lookup once: tran_id -> (start, stop) for fast per-transcript access
    cds_lookup: Dict[str, tuple] = {
        row["tran_id"]: (int(row["start"]), int(row["stop"]))
        for row in cds_in.iter_rows(named=True)
    }
    cds_tran_ids = set(cds_lookup.keys())

    if not has_length:
        acc = np.zeros(2 * W + 1, dtype=np.float64)
        for key, sub in profiles.group_by("tran_id"):
            tid = str(key[0] if isinstance(key, tuple) else key)
            if tid not in cds_tran_ids:
                continue
            cds_start, cds_stop = cds_lookup[tid]
            mxp = int(sub["pos"].max())
            arr = _dense_array(sub, mxp, pad=W + 1)
            for anchor in range(cds_start, cds_stop, 3):
                lo, hi = anchor - W, anchor + W + 1
                if lo < 0 or hi > len(arr):
                    continue
                acc += arr[lo:hi]
        total = acc.sum()
        return {None: acc / total if total > 0 else fallback}

    result: Dict[Optional[int], np.ndarray] = {}
    for L_val in sorted(profiles["length"].unique().to_list()):
        sub_l = profiles.filter(pl.col("length") == L_val)
        acc = np.zeros(2 * W + 1, dtype=np.float64)
        for key, sub in sub_l.group_by("tran_id"):
            tid = str(key[0] if isinstance(key, tuple) else key)
            if tid not in cds_tran_ids:
                continue
            cds_start, cds_stop = cds_lookup[tid]
            mxp = int(sub["pos"].max())
            arr = _dense_array(sub, mxp, pad=W + 1)
            for anchor in range(cds_start, cds_stop, 3):
                lo, hi = anchor - W, anchor + W + 1
                if lo < 0 or hi > len(arr):
                    continue
                acc += arr[lo:hi]
        total = acc.sum()
        result[int(L_val)] = acc / total if total > 0 else fallback.copy()

    return result


def _rl_deconvolve(
    signal: np.ndarray,
    kernel: np.ndarray,
    n_iter: int = _RL_ITERS,
) -> np.ndarray:
    """1-D Richardson-Lucy deconvolution.

    Recovers T from O ~= conv(T, K) via:
        T^(t+1) = T^(t) * conv(O / conv(T^(t), K), K_flip)

    Output has same total mass as ``signal``.
    """
    orig_sum = float(signal.sum())
    if orig_sum == 0:
        return signal.copy()

    W = len(kernel) // 2
    pad = W + 2

    padded = np.pad(signal, pad, mode="edge").astype(np.float64)
    T = padded.copy()
    T[T <= 0] = 1e-9

    K = kernel.astype(np.float64)
    K_flip = K[::-1]

    for _ in range(n_iter):
        O_hat = np.convolve(T, K, mode="same")
        O_hat[O_hat <= 0] = 1e-12
        ratio = padded / O_hat
        correction = np.convolve(ratio, K_flip, mode="same")
        correction[correction < 0] = 0.0
        T = T * correction
        T[T < 0] = 0.0

    T_out = T[pad : pad + len(signal)]
    out_sum = T_out.sum()
    if out_sum > 0:
        T_out *= orig_sum / out_sum
    return T_out


def deconvolve(
    profiles: pl.DataFrame,
    kernels: Dict[Optional[int], object],
    rl_iters: int = _RL_ITERS,
) -> pl.DataFrame:
    """Apply Richardson-Lucy deconvolution to transcript-space profiles.

    Parameters
    ----------
    profiles : tran_id, pos, count[, length]
    kernels  : output of ``learn_kernel`` – length (or None) -> kernel array.
    rl_iters : number of RL iterations.

    Returns
    -------
    DataFrame with same schema but sharpened (deconvolved) counts.
    """
    if profiles.is_empty():
        return profiles

    has_length = "length" in profiles.columns

    kernels_np: Dict[Optional[int], np.ndarray] = {k: _to_np_kernel(v) for k, v in kernels.items()}
    global_kernel = kernels_np.get(None)

    group_cols = ["tran_id"] + (["length"] if has_length else [])

    # Accumulate output as parallel lists then build DataFrame once — avoids
    # the OOM cost of building millions of Python dicts.
    out_tran_ids: list = []
    out_positions: list = []
    out_counts: list = []
    out_lengths: list = []

    for key, sub in profiles.group_by(group_cols):
        if isinstance(key, tuple):
            tran_id = key[0]
            L_val = int(key[1]) if has_length else None
        else:
            tran_id = key
            L_val = None

        K = kernels_np.get(L_val, global_kernel) if has_length else global_kernel

        sub_s = sub.sort("pos")
        positions = sub_s["pos"].to_numpy().astype(int)
        counts = sub_s["count"].to_numpy().astype(np.float64)

        if len(positions) == 0:
            continue

        # Delta or missing kernel: pass-through without deconvolution
        if K is None or (len(K) == 1 and float(K[0]) >= 0.99):
            nz_mask = counts > 0
            nz_pos = positions[nz_mask]
            nz_cnt = counts[nz_mask]
        else:
            min_p = int(positions.min())
            max_p = int(positions.max())
            dense = np.zeros(max_p - min_p + 1, dtype=np.float64)
            np.add.at(dense, positions - min_p, counts)

            deconv = _rl_deconvolve(dense, K, n_iter=rl_iters)
            nz_idx = np.nonzero(deconv > 0)[0]
            nz_pos = nz_idx + min_p
            nz_cnt = deconv[nz_idx]

        n = len(nz_pos)
        if n == 0:
            continue

        out_tran_ids.extend([tran_id] * n)
        out_positions.append(nz_pos)
        out_counts.append(nz_cnt)
        if has_length:
            out_lengths.extend([L_val] * n)

    if not out_tran_ids:
        return profiles

    all_pos = np.concatenate(out_positions)
    all_cnt = np.concatenate(out_counts)

    cols = {
        "tran_id": out_tran_ids,
        "pos": all_pos,
        "count": all_cnt,
    }
    if has_length:
        cols["length"] = out_lengths

    return pl.DataFrame(cols).sort(group_cols)

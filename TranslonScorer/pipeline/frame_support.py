from __future__ import annotations

from typing import Dict, List

import polars as pl
import numpy as np

from ..utils.logging import log_info, log_warning
from .config import Config
from ..frame.bleed import learn_confusion, apply_confusion_counts
from ..frame.hmm import smooth_posteriors
from ..frame.deblur import learn_kernel as learn_deblur_kernel, deconvolve
from ..frame.latent import fit_latent


MAX_FRAME_ENTROPY = float(np.log2(3.0))
QC_LOW_SUPPORT_HIGH_CONF = 1 << 0
QC_MIXED_FRAME_EVIDENCE = 1 << 1


def _profiles_to_codon_frame_counts(profiles: pl.DataFrame) -> pl.DataFrame:
    """Aggregate nucleotide positions into per-codon frame triplets per transcript.

    Input schema: tran_id, pos, count[, length]
    Output: tran_id, codon, f0, f1, f2[, length]
    where codon = pos // 3 and fX sums counts at positions 3k+X.
    """
    cols    = profiles.columns
    has_len = "length" in cols
    df      = profiles.with_columns(
        (pl.col("pos") // 3).alias("codon"),
        (pl.col("pos") %  3).alias("frame"),
    )
    groups = ["tran_id", "codon"] + (["length"] if has_len else [])
    agg = (
        df.group_by(groups + ["frame"])
        .agg(pl.col("count").sum().alias("sum"))
        .pivot(values="sum", index=groups, on="frame")
        .fill_null(0)
        .rename({"0": "f0", "1": "f1", "2": "f2"})
    )
    for col in ("f0", "f1", "f2"):
        if col not in agg.columns:
            agg = agg.with_columns(pl.lit(0.0).alias(col))
    return agg.select(groups + ["f0", "f1", "f2"]).sort(groups)


def _empty_frame_support(has_len: bool = False) -> pl.DataFrame:
    schema = {
        "tran_id": pl.Utf8,
        "codon": pl.Int64,
        "observed_f0": pl.Float64,
        "observed_f1": pl.Float64,
        "observed_f2": pl.Float64,
        "adjusted_f0": pl.Float64,
        "adjusted_f1": pl.Float64,
        "adjusted_f2": pl.Float64,
        "total_count": pl.Float64,
        "p0": pl.Float64,
        "p1": pl.Float64,
        "p2": pl.Float64,
        "entropy": pl.Float64,
        "p_max": pl.Float64,
        "secondary_frame_mass": pl.Float64,
        "frame_periodicity_score": pl.Float64,
        "depth_score": pl.Float64,
        "support_evidence": pl.Float64,
        "support_gated_p0": pl.Float64,
        "support_gated_p1": pl.Float64,
        "support_gated_p2": pl.Float64,
        "method": pl.Utf8,
        "qc_flags": pl.Int32,
    }
    if has_len:
        rest = {k: v for k, v in schema.items() if k not in {"tran_id", "codon"}}
        schema = {"tran_id": pl.Utf8, "codon": pl.Int64, "length": pl.Int64, **rest}
    return pl.DataFrame(schema=schema)


def _posterior_from_counts(adjusted: np.ndarray) -> np.ndarray:
    P = adjusted.astype(float).copy()
    rowsum = P.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    return P / rowsum


def _support_from_arrays(
    grp: pl.DataFrame,
    ids: list[str],
    obs: np.ndarray,
    adjusted: np.ndarray,
    method: str,
) -> pl.DataFrame:
    P = _posterior_from_counts(adjusted)
    ent = _entropy_vec(P)
    total = obs.sum(axis=1)
    return (
        grp.select(ids)
        .with_columns([
            pl.Series("observed_f0", obs[:, 0]),
            pl.Series("observed_f1", obs[:, 1]),
            pl.Series("observed_f2", obs[:, 2]),
            pl.Series("adjusted_f0", adjusted[:, 0]),
            pl.Series("adjusted_f1", adjusted[:, 1]),
            pl.Series("adjusted_f2", adjusted[:, 2]),
            pl.Series("total_count", total),
            pl.Series("p0", P[:, 0]),
            pl.Series("p1", P[:, 1]),
            pl.Series("p2", P[:, 2]),
            pl.Series("entropy", ent),
            pl.lit(method).alias("method"),
            pl.lit(0).cast(pl.Int32).alias("qc_flags"),
        ])
    )


def _hmm_smooth(out: pl.DataFrame, has_len: bool, lam: float) -> pl.DataFrame:
    """Apply HMM forward-backward smoothing per transcript (and per length)."""
    by = ["tran_id"] + (["length"] if has_len else [])
    parts: List[pl.DataFrame] = []
    for _key, df in out.sort(by + ["codon"]).group_by(by):
        arr    = df.select(["p0", "p1", "p2"]).to_numpy()
        P_s, _ = smooth_posteriors(arr, lam)
        exprs = [
            pl.Series("p0", P_s[:, 0]),
            pl.Series("p1", P_s[:, 1]),
            pl.Series("p2", P_s[:, 2]),
            pl.Series("entropy", _entropy_vec(P_s)),
        ]
        if "total_count" in df.columns:
            total = df.get_column("total_count").to_numpy()
            exprs.extend([
                pl.Series("adjusted_f0", P_s[:, 0] * total),
                pl.Series("adjusted_f1", P_s[:, 1] * total),
                pl.Series("adjusted_f2", P_s[:, 2] * total),
            ])
        parts.append(df.with_columns(exprs))
    return pl.concat(parts) if parts else out


def _entropy(p: np.ndarray) -> float:
    pnz = p[p > 0]
    return float(-(pnz * np.log2(pnz)).sum())


def _entropy_vec(P: np.ndarray) -> np.ndarray:
    """Vectorised Shannon entropy for each row of P (shape N×3).  Returns shape (N,)."""
    safe = np.where(P > 0, P, 1.0)
    return -(P * np.log2(safe)).sum(axis=1)


def _support_depth_scale(out: pl.DataFrame) -> float:
    if out.is_empty() or "total_count" not in out.columns:
        return 1.0
    positive = out.filter(pl.col("total_count") > 0)
    if positive.is_empty():
        return 1.0
    try:
        scale = float(positive.select(pl.col("total_count").median()).item())
    except Exception:
        scale = 1.0
    if not np.isfinite(scale) or scale <= 0:
        try:
            scale = float(positive.select(pl.col("total_count").mean()).item())
        except Exception:
            scale = 1.0
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0
    return scale


def _add_support_diagnostics(out: pl.DataFrame) -> pl.DataFrame:
    """Add evidence diagnostics without changing conditional frame posteriors.

    p0/p1/p2 answer "which frame, conditional on translated frame signal".
    support_evidence is an intentionally conservative diagnostic combining
    local depth and periodicity; it is not a calibrated translation posterior.
    """
    if out.is_empty():
        return out

    depth_scale = _support_depth_scale(out)
    pmax = pl.max_horizontal(["p0", "p1", "p2"])
    pmin = pl.min_horizontal(["p0", "p1", "p2"])

    out = (
        out.with_columns([
            pmax.alias("p_max"),
            (pl.col("p0") + pl.col("p1") + pl.col("p2") - pmax - pmin).alias("secondary_frame_mass"),
            (1.0 - (pl.col("entropy") / MAX_FRAME_ENTROPY)).clip(0.0, 1.0).alias("frame_periodicity_score"),
            (1.0 - (-pl.col("total_count") / depth_scale).exp()).clip(0.0, 1.0).alias("depth_score"),
        ])
        .with_columns((pl.col("frame_periodicity_score") * pl.col("depth_score")).alias("support_evidence"))
        .with_columns([
            (pl.col("p0") * pl.col("support_evidence")).alias("support_gated_p0"),
            (pl.col("p1") * pl.col("support_evidence")).alias("support_gated_p1"),
            (pl.col("p2") * pl.col("support_evidence")).alias("support_gated_p2"),
        ])
    )

    low_support_high_conf = (pl.col("total_count") < depth_scale) & (pl.col("p_max") >= 0.80)
    mixed_frame_evidence = (pl.col("total_count") >= depth_scale) & (pl.col("secondary_frame_mass") >= 0.20)
    return out.with_columns(
        (
            low_support_high_conf.cast(pl.Int32) * QC_LOW_SUPPORT_HIGH_CONF
            + mixed_frame_evidence.cast(pl.Int32) * QC_MIXED_FRAME_EVIDENCE
        ).alias("qc_flags")
    )


# ---------------------------------------------------------------------------
# Linear correction path  (methods: linear, linear+hmm, deblur+linear+hmm)
# ---------------------------------------------------------------------------

def _build_linear(
    prof_use:  pl.DataFrame,
    cds_df:    pl.DataFrame,
    has_len:   bool,
    method:    str,
) -> pl.DataFrame:
    """Build frame posteriors via linear confusion-matrix correction."""
    codon_fc  = _profiles_to_codon_frame_counts(prof_use)
    C_by_len  = learn_confusion(prof_use, cds_df, by_length=has_len)

    if has_len:
        # Per-length: batch each unique length together to stay vectorised
        parts: List[pl.DataFrame] = []
        for (length,), grp in codon_fc.group_by(["length"]):
            _c = C_by_len.get(int(length))
            if _c is None:
                _c = C_by_len.get(None)
            C   = _c if _c is not None else np.eye(3)
            obs = grp.select(["f0", "f1", "f2"]).to_numpy().astype(float)
            adj = apply_confusion_counts(obs, C, alpha=0.05)
            parts.append(_support_from_arrays(grp, ["tran_id", "codon", "length"], obs, adj, method))
        return pl.concat(parts) if parts else _empty_frame_support(has_len=True)
    else:
        _c  = C_by_len.get(None)
        C   = _c if _c is not None else np.eye(3)
        obs = codon_fc.select(["f0", "f1", "f2"]).to_numpy().astype(float)
        adj = apply_confusion_counts(obs, C, alpha=0.05)
        return _support_from_arrays(codon_fc, ["tran_id", "codon"], obs, adj, method)


# ---------------------------------------------------------------------------
# Latent EM path  (method: latent)
# ---------------------------------------------------------------------------

def _build_latent(
    prof_use:   pl.DataFrame,
    cds_df:     pl.DataFrame,
    has_len:    bool,
    method:     str,
    background: str,
) -> pl.DataFrame:
    """Build frame posteriors via joint EM over the Poisson mixture model.

    Operates on raw codon frame counts. The confusion matrix learned from
    CDS is held fixed by default, so EM refines per-codon latent frame support
    T without letting sparse or non-canonical regions re-train the leakage
    calibration.
    """
    codon_fc = _profiles_to_codon_frame_counts(prof_use)
    C_by_len = learn_confusion(prof_use, cds_df, by_length=has_len)

    if has_len:
        parts: List[pl.DataFrame] = []
        for (length,), grp in codon_fc.group_by(["length"]):
            _ci    = C_by_len.get(int(length))
            if _ci is None:
                _ci = C_by_len.get(None)
            C_init = _ci if _ci is not None else np.eye(3)
            obs    = grp.select(["f0", "f1", "f2"]).to_numpy().astype(float)
            P, _, _ = fit_latent(obs, C_init, background=background)
            adj = P * obs.sum(axis=1, keepdims=True)
            parts.append(_support_from_arrays(grp, ["tran_id", "codon", "length"], obs, adj, method))
        return pl.concat(parts) if parts else _empty_frame_support(has_len=True)
    else:
        _ci    = C_by_len.get(None)
        C_init = _ci if _ci is not None else np.eye(3)
        obs    = codon_fc.select(["f0", "f1", "f2"]).to_numpy().astype(float)
        P, _, _ = fit_latent(obs, C_init, background=background)
        adj = P * obs.sum(axis=1, keepdims=True)
        return _support_from_arrays(codon_fc, ["tran_id", "codon"], obs, adj, method)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def build_frame_support(profiles: pl.DataFrame, cds_df: pl.DataFrame, cfg: Config) -> pl.DataFrame:
    """Build codon-level frame posteriors and write to Parquet when requested.

    Returns a DataFrame: tran_id, codon, p0, p1, p2, method[, length, entropy, qc_flags].

    Supported methods
    -----------------
    linear            : confusion-matrix linear correction
    linear+hmm        : linear correction + HMM forward-backward smoothing
    deblur+linear+hmm : RL positional deblur → linear correction → HMM smoothing
    latent            : count-weighted EM for Poisson mixture (O = M T + B)
    """
    method = (cfg.frame_method or "none").lower()
    if method == "none":
        log_info("Frame method=none; skipping frame support")
        return _empty_frame_support(has_len=False)

    has_len  = ("length" in profiles.columns) and bool(cfg.frame_by_length)

    # --- Optional RL deblurring pre-step ---
    prof_use = profiles
    if method.startswith("deblur"):
        try:
            kernels  = learn_deblur_kernel(profiles, cds_df, by_length=has_len)
            prof_use = deconvolve(profiles, kernels)
        except Exception as e:
            log_warning(f"Deblur failed ({e}); continuing without deblurring")
            prof_use = profiles

    # --- Route to the appropriate correction method ---
    if method == "latent":
        out = _build_latent(prof_use, cds_df, has_len, method, cfg.frame_background)
    else:
        # linear / linear+hmm / deblur+linear+hmm  (all use linear correction base)
        out = _build_linear(prof_use, cds_df, has_len, method)

    # --- Optional HMM smoothing ---
    if method.endswith("+hmm") and not out.is_empty():
        out = _hmm_smooth(out, has_len, cfg.frame_hmm_lambda)

    out = _add_support_diagnostics(out)

    # --- Write output ---
    if cfg.frame_support_out and not out.is_empty():
        try:
            out.write_parquet(cfg.frame_support_out)
            log_info(f"Frame support written: {cfg.frame_support_out}")
        except Exception as e:
            log_warning(f"Failed to write frame support ({e})")

    return out

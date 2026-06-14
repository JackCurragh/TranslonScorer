"""Re-export shim — implementations live in TranslonScorer.frame_support.

The old public API accepted a Config object; the new pure function takes a
FrameSupportParams.  This shim bridges the two so existing callers continue
to work unchanged until the final cleanup task (T14).
"""
from __future__ import annotations

import polars as pl

from TranslonScorer.model import FrameSupportParams
from TranslonScorer.frame_support import (  # noqa: F401 — re-exported
    MAX_FRAME_ENTROPY,
    QC_LOW_SUPPORT_HIGH_CONF,
    QC_MIXED_FRAME_EVIDENCE,
    _empty_frame_support,
    _profiles_to_codon_frame_counts,
    _posterior_from_counts,
    _entropy_vec,
    _support_from_arrays,
    _hmm_smooth,
    _support_depth_scale,
    _add_support_diagnostics,
    _build_linear,
    _build_latent,
    build_frame_support as _build_frame_support_pure,
)


def build_frame_support(profiles: pl.DataFrame, cds_df: pl.DataFrame, cfg: object) -> pl.DataFrame:
    """Adapter: build_frame_support(profiles, cds_df, cfg: Config) → DataFrame.

    Translates old-style Config into FrameSupportParams, calls the pure
    function, and handles the optional file-write that the pure function omits.
    """
    params = FrameSupportParams(
        frame_method=getattr(cfg, "frame_method", "none"),
        frame_by_length=bool(getattr(cfg, "frame_by_length", False)),
        frame_background=getattr(cfg, "frame_background", "uniform"),
        frame_hmm_lambda=float(getattr(cfg, "frame_hmm_lambda", 1.0)),
    )
    out = _build_frame_support_pure(profiles, cds_df, params)

    out_path = getattr(cfg, "frame_support_out", None)
    if out_path and not out.is_empty():
        try:
            out.write_parquet(out_path)
        except Exception:
            pass

    return out

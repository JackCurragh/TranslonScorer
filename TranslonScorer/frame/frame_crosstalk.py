from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import polars as pl
from ..utils.logging import log_info


def estimate_crosstalk_matrix(profiles: pl.DataFrame, cds_windows: pl.DataFrame, length: int) -> Tuple[np.ndarray, float]:
    """Estimate 3x3 frame confusion matrix C for a given read length.

    Inputs:
      profiles(tran_id,pos,count,length?) corrected to include per-length; if length column missing, caller must subset
      cds_windows(tran_id, start_pos, end_pos): high-confidence CDS regions
    Returns: (C, alpha)
    """
    # Simplified: compute average proportion in frames 0,1,2 relative to frame-0 signal
    # Build counts by frame across confident windows
    if 'length' in profiles.columns:
        p = profiles.filter(pl.col('length') == length)
    else:
        p = profiles

    if p.is_empty() or cds_windows.is_empty():
        return np.eye(3, dtype=float), 0.05

    joined = p.join(cds_windows, on='tran_id', how='inner') \
             .filter((pl.col('pos') >= pl.col('start_pos')) & (pl.col('pos') < pl.col('end_pos'))) \
             .with_columns((pl.col('pos') % 3).alias('frame'))

    agg = joined.group_by('frame').agg(pl.col('count').sum().alias('sum')).sort('frame')
    sums = np.array([float(agg['sum'][i]) if i < agg.height else 0.0 for i in range(3)], dtype=float)
    total = sums.sum() or 1.0
    # Create a diagonal-dominant C with leakage inferred from proportions
    f0 = (sums[0] / total) if total > 0 else 1.0
    leak = max(1e-6, (1 - f0) / 2)
    C = np.array([[f0, leak, leak], [leak, f0, leak], [leak, leak, f0]], dtype=float)
    alpha = 0.05
    return C, alpha


def invert_and_correct(C: np.ndarray, alpha: float, f_counts: np.ndarray) -> np.ndarray:
    """Regularized inversion and nonnegativity clamp.

    f_counts: shape (N,3) counts per codon frame
    """
    I = np.eye(3)
    Cinvl = np.linalg.pinv(C + alpha * I)
    corrected = f_counts @ Cinvl.T
    corrected[corrected < 0] = 0
    return corrected


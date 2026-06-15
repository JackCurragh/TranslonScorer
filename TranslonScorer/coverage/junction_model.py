from __future__ import annotations

from typing import Dict, Tuple

import polars as pl
import numpy as np
from ..utils.logging import log_info


def estimate_junction_expectation(
    genome_bam_df: pl.DataFrame,
    cds_junctions: pl.DataFrame,
    *,
    n_bins: int = 6,
) -> pl.DataFrame:
    """Estimate expected split-read fraction per geometry bin with gene-level shrinkage.

    Inputs:
      genome_bam_df: rows with columns chr,start,stop,strand,length,count and split flag via CIGAR 'N' pre-parsed
      cds_junctions: junction table for confident CDS transcripts (chr,donor,acceptor,gene_id)
    Output schema:
      bin_id, E_split_fraction, var, gene_id (optional), length(optional)
    """
    # Placeholder: compute global baseline by binning junction distance and local coverage
    # In practice, caller should precompute split flags and geometry
    if genome_bam_df.is_empty() or cds_junctions.is_empty():
        return pl.DataFrame({"bin_id": [], "E_split_fraction": [], "var": []})

    # Mock bin: single bin with global rate
    split = genome_bam_df.filter(pl.col("cigar_N") == True)["count"].sum()
    total = genome_bam_df["count"].sum() or 1
    frac = float(split / total)
    return pl.DataFrame(
        {"bin_id": [0], "E_split_fraction": [frac], "var": [frac * (1 - frac) + 1e-6]}
    )


def junction_llr(obs_split: int, exp_frac: float, total: int, var: float) -> float:
    """Simple Gaussian LLR proxy until NB model is added."""
    exp = exp_frac * total
    denom = max(var * total, 1e-6)
    return -0.5 * ((obs_split - exp) ** 2) / denom

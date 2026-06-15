"""Pure P/A-site profile generation from assigned reads + explicit offset table.

Design constraints
------------------
- NO offset inference here.  The caller must supply a pre-computed offset table
  (dict mapping read_length → P-site offset in transcript coordinates) produced
  by `TranslonScorer.offsets.make_offset_table`.
- A-site = P-site + 3 nt by default.  If the offset table carries an explicit
  'A' column the provider may pass those directly; otherwise A = P + 3.
- Aggregation happens *after* per-sample/per-length offsets are applied; the
  caller folds/merges profiles, not this module.
- Prefix-sum helpers for vectorised range queries live here (the scoring layer
  imports them from this module or from `scoring/run.py` which re-exports them).

Public API
----------
apply_offsets(reads_df, offset_table, *, site, keep_length) → tidy profiles
size_factors(profiles, sample_col)                          → {sample: float}
_prefix_sums(cov_pos, cov_cnt)        → (pos_arr, cum_all, cum_frame)
_range_sums(pos_arr, cum_all, cum_f, starts, ends) → (total, frames, covered)
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl


# ---------------------------------------------------------------------------
# Strand-aware genomic site placement (pure)
# ---------------------------------------------------------------------------


def site_position(
    ref_start: int, ref_end: int, is_reverse: bool, p_offset: int, site: str = "A"
) -> Tuple[int, int]:
    """Genomic (strand, position) of the P- or A-site for one aligned read.

    The 5' end of the footprint is `ref_start` on + strand and `ref_end - 1`
    on - strand. The P-site sits `p_offset` nt 3' of the 5' end; the A-site is
    one codon (3 nt) further 3'. Returns (strand∈{1,-1}, genomic_pos).
    """
    a_shift = 3 if site == "A" else 0
    if is_reverse:
        return -1, (ref_end - 1) - p_offset - a_shift
    return 1, ref_start + p_offset + a_shift


# ---------------------------------------------------------------------------
# Core profile builder (pure — no I/O, no offset inference)
# ---------------------------------------------------------------------------


def apply_offsets(
    reads_df: pl.DataFrame,
    offset_table: Dict[int, int],
    *,
    site: str = "A",
    keep_length: bool = False,
    default_offset: int = 12,
) -> pl.DataFrame:
    """Apply an explicit P-site offset table to transcript-mapped reads.

    Parameters
    ----------
    reads_df     : DataFrame with columns tran_id, tran_start_bam, length, count.
                   May also carry sample_id and/or study_id for multi-sample inputs.
    offset_table : {read_length: p_site_offset} produced by offsets.make_offset_table.
                   Lengths absent from the table fall back to default_offset.
    site         : "P" → P-site offset as-is; "A" → P-site offset + 3 nt.
    keep_length  : include the read-length column in the output (for by-length QC).
    default_offset : P-site offset for read lengths not in the table.

    Returns
    -------
    Tidy DataFrame: tran_id, pos, count[, sample_id][, length]
    Positions are grouped and summed so the output has one row per (tran_id, pos)
    (or per (tran_id, pos, sample_id[, length])).
    """
    if site not in {"P", "A"}:
        raise ValueError(f"site must be 'P' or 'A', got {site!r}")
    if reads_df.is_empty():
        cols = {"tran_id": pl.Utf8, "pos": pl.Int64, "count": pl.Float64}
        if keep_length:
            cols["length"] = pl.Int64
        return pl.DataFrame(schema=cols)

    required = {"tran_id", "tran_start_bam", "length", "count"}
    missing = required - set(reads_df.columns)
    if missing:
        raise ValueError(f"reads_df is missing required columns: {sorted(missing)}")

    a_shift = 3 if site == "A" else 0

    def _offset(length: int) -> int:
        return int(offset_table.get(int(length), default_offset)) + a_shift

    sample_cols = [c for c in ("sample_id", "sample_index", "study_id") if c in reads_df.columns]
    group_cols = sample_cols + ["tran_id", "pos"]
    if keep_length and "length" in reads_df.columns:
        group_cols = group_cols + ["length"]

    out = (
        reads_df.with_columns(
            pl.col("length").map_elements(_offset, return_dtype=pl.Int64).alias("_ofs"),
            pl.col("tran_start_bam").cast(pl.Int64),
        )
        .with_columns((pl.col("tran_start_bam") + pl.col("_ofs")).alias("pos"))
        .drop("_ofs", "tran_start_bam")
    )
    if not keep_length:
        out = out.drop("length", strict=False)
    return out.group_by(group_cols).agg(pl.col("count").sum()).sort(group_cols)


# ---------------------------------------------------------------------------
# Size factor estimation (pure — no file I/O)
# ---------------------------------------------------------------------------


def size_factors(
    profiles: pl.DataFrame,
    sample_col: str = "sample_id",
) -> Dict[str, float]:
    """Estimate per-sample depth-normalisation factors (median-ratio method).

    If the profiles DataFrame has no sample_col, returns {"": 1.0}.

    The median-ratio method:
        geometric_mean_per_pos = exp(mean(log(count + ε)))  across samples
        ratio_per_pos          = count / geometric_mean
        size_factor_per_sample = median(ratio_per_pos)
    Falls back to library-size normalisation when too few positions pass
    the geometric mean filter.

    Parameters
    ----------
    profiles  : tidy DataFrame with at minimum (pos, count[, sample_col]).
    sample_col: column name holding the sample identifier.

    Returns
    -------
    {sample_id: float} — factors suitable for dividing raw counts.
    """
    if sample_col not in profiles.columns:
        return {"": 1.0}

    eps = 1e-6
    samples = profiles.get_column(sample_col).unique().sort().to_list()
    if len(samples) == 1:
        return {str(samples[0]): 1.0}

    wide = (
        profiles.select([sample_col, "pos", "count"])
        .pivot(values="count", index="pos", on=sample_col)
        .fill_null(0.0)
    )
    mat = wide.drop("pos").to_numpy().astype(float)
    log_mat = np.log(mat + eps)
    geo_mean = log_mat.mean(axis=1)
    nz = geo_mean > np.log(eps * 10)
    if nz.sum() < 10:
        totals = mat.sum(axis=0)
        med = float(np.median(totals)) or 1.0
        return {str(s): float(t / med) for s, t in zip(samples, totals)}

    ratios = mat[nz] / np.exp(geo_mean[nz, None])
    sf = np.median(ratios, axis=0)
    sf[sf <= 0] = 1.0
    return {str(s): float(f) for s, f in zip(samples, sf)}


# ---------------------------------------------------------------------------
# Prefix-sum helpers for vectorised coverage range queries
# (the scoring layer imports from scoring/run.py which re-exports these)
# ---------------------------------------------------------------------------


def _prefix_sums(
    cov_pos: np.ndarray,
    cov_cnt: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, List]:
    """Build prefix-sum arrays for O(1) frame-aware range queries.

    Parameters
    ----------
    cov_pos : int64 array of positions with coverage (need not be sorted).
    cov_cnt : float64 array of counts at those positions.

    Returns
    -------
    (p, cum_all, cum_f)
    p       : positions sorted ascending (int64).
    cum_all : cumulative sum of all counts; length = n+1, cum_all[0] = 0.
    cum_f   : list of 3 cumulative-sum arrays, one per mod-3 frame.
              cum_f[k] is the cumulative sum of counts at positions where
              pos % 3 == k; same length convention as cum_all.
    """
    order = np.argsort(cov_pos, kind="stable")
    p = cov_pos[order].astype(np.int64)
    c = cov_cnt[order].astype(np.float64)
    fr = np.mod(p, 3)
    cum_all = np.concatenate([[0.0], np.cumsum(c)])
    cum_f = [np.concatenate([[0.0], np.cumsum(np.where(fr == f, c, 0.0))]) for f in range(3)]
    return p, cum_all, cum_f


def _range_sums(
    p: np.ndarray,
    cum_all: np.ndarray,
    cum_f: List,
    starts: np.ndarray,
    ends: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised range-sum query over prefix arrays from _prefix_sums.

    Parameters
    ----------
    p       : sorted positions from _prefix_sums.
    cum_all : cumulative-all from _prefix_sums.
    cum_f   : list of 3 per-frame cumulative arrays from _prefix_sums.
    starts  : int array of query range starts (inclusive).
    ends    : int array of query range ends (exclusive).

    Returns
    -------
    (total, frames, covered)
    total   : float array — total counts in each [start, end) range.
    frames  : float array of shape [n, 3] — per-frame counts in each range.
    covered : int array   — number of covered positions in each range.
    """
    lo = np.searchsorted(p, starts, side="left")
    hi = np.searchsorted(p, ends, side="left")
    total = cum_all[hi] - cum_all[lo]
    frames = np.stack([cum_f[f][hi] - cum_f[f][lo] for f in range(3)], axis=-1)
    covered = hi - lo
    return total, frames, covered

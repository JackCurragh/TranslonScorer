"""P/A-site offset calibration: plausible-range guard, method dispatch, offset tables.

All functions are pure (no I/O side effects) except `file_offsets`, which reads a
single CSV file and is idempotent.  `metagene_offsets` is a stub — the per-BAM
5′-end histogram required for the real implementation is built by BamSetProvider
(T12) and passed in explicitly; it is never inferred inside a locus or profile.

Public API
----------
plausible_offset_range  — (read_length, OffsetParams) → (lo, hi)
usable_read_length      — (read_length, OffsetParams) → bool
psite_to_asite          — P-site offset → A-site offset (P + 3 nt)
metagene_offsets        — 5′-end histogram DataFrame + OffsetParams → {length: offset}  [stub]
file_offsets            — CSV path + OffsetParams → {length: offset}
global_offsets          — OffsetParams [+ Series of lengths] → {length: offset}
make_offset_table       — dispatch by OffsetParams.method → {length: offset}

Contract (all functions)
- Return type is always Dict[int, int]: {read_length: p_site_offset}.
- Keys cover only usable read lengths (via usable_read_length).
- Offsets are clamped to plausible_offset_range.
- Inputs are OffsetParams + a data source; NEVER a locus, transcript, or profile matrix.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import polars as pl

from TranslonScorer.model import OffsetParams


# ---------------------------------------------------------------------------
# Plausibility guard
# ---------------------------------------------------------------------------

def plausible_offset_range(read_length: int, p: OffsetParams) -> Tuple[int, int]:
    """Physically plausible P-site offset window for `read_length`.

    A 25 nt read cannot have an 18 nt offset:
      hi = min(p.offset_max, floor(read_length × p.max_frac))
    """
    hi = min(p.offset_max, int(read_length * p.max_frac))
    return p.offset_min, hi


def usable_read_length(read_length: int, p: OffsetParams) -> bool:
    """True when `read_length` falls within the configured window."""
    return p.min_read_len <= read_length <= p.max_read_len


def psite_to_asite(offset_p: int) -> int:
    """A-site offset = P-site + 3 nt (next codon downstream)."""
    return offset_p + 3


# ---------------------------------------------------------------------------
# Offset methods
# ---------------------------------------------------------------------------

def metagene_offsets(
    five_prime_by_length: pl.DataFrame,
    p: OffsetParams,
) -> Dict[int, int]:
    """P-site offsets from a whole-sample 5′-end metagene (stub — pending T12).

    Parameters
    ----------
    five_prime_by_length:
        DataFrame with columns read_length (int), rel_pos (5′ end position
        relative to annotated start codon), count (float).  Must be built
        from UNIQUE reads across the whole BAM, NOT from a single locus,
        transcript, or profile matrix.
    p:
        Offset parameters controlling the plausible window.

    Returns
    -------
    {read_length: p_site_offset} restricted to usable lengths and the
    plausible window for each length.
    """
    raise NotImplementedError(
        "metagene_offsets: real implementation pending BamSetProvider (T12). "
        "Build the per-BAM per-length 5′-end histogram from the whole BAM "
        "before calling this function."
    )


def file_offsets(path: str, p: OffsetParams) -> Dict[int, int]:
    """Load per-length P-site offsets from a CSV file.

    Required columns: read_length (int), offset (int).
    Optional column:  site (str) — rows where site ∉ {'P', 'p'} are dropped
                      so that A-site files do not silently corrupt the table.

    Results are filtered to usable read lengths and clamped to
    plausible_offset_range for each length.
    """
    df = pl.read_csv(path)
    required = {"read_length", "offset"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Offsets file {path!r} missing columns {missing}; got {list(df.columns)}"
        )
    if "site" in df.columns:
        df = df.filter(pl.col("site").is_in(["P", "p"]))
    result: Dict[int, int] = {}
    for row in df.select(["read_length", "offset"]).iter_rows():
        length, raw_offset = int(row[0]), int(row[1])
        if not usable_read_length(length, p):
            continue
        lo, hi = plausible_offset_range(length, p)
        result[length] = max(lo, min(hi, raw_offset))
    return result


def global_offsets(
    p: OffsetParams,
    read_lengths: Optional[pl.Series] = None,
) -> Dict[int, int]:
    """Return a flat {read_length: p.global_offset} table.

    If `read_lengths` is supplied (a Series of integer lengths observed in the
    data), the table covers exactly the usable subset of those lengths.
    If omitted, the table covers the full usable range
    [p.min_read_len, p.max_read_len].
    """
    if read_lengths is not None:
        lengths = sorted(
            int(l) for l in read_lengths.drop_nulls().unique().to_list()
            if usable_read_length(int(l), p)
        )
    else:
        lengths = list(range(p.min_read_len, p.max_read_len + 1))
    return {length: p.global_offset for length in lengths}


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

def make_offset_table(
    p: OffsetParams,
    read_lengths: Optional[pl.Series] = None,
) -> Dict[int, int]:
    """Dispatch to the offset method named in `p.method`.

    "global"   → global_offsets   (real; uses p.global_offset)
    "file"     → file_offsets     (real; requires p.offsets_file)
    "metagene" → metagene_offsets (stub — caller must build the 5′-end
                                   histogram first and call it directly)

    Returns {read_length: p_site_offset}.
    """
    if p.method == "global":
        return global_offsets(p, read_lengths)
    if p.method == "file":
        if p.offsets_file is None:
            raise ValueError(
                "OffsetParams.method='file' requires offsets_file to be set"
            )
        return file_offsets(p.offsets_file, p)
    if p.method == "metagene":
        raise NotImplementedError(
            "metagene_offsets requires a 5′-end histogram DataFrame; "
            "call metagene_offsets() directly once the whole-BAM histogram "
            "is available from BamSetProvider."
        )
    raise ValueError(
        f"Unknown offset method {p.method!r}; expected 'metagene', 'file', or 'global'"
    )

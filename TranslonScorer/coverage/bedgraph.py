"""BedgraphSetProvider — coverage from one or more bedGraph files.

Exists so the normal scorer can consume the exact input the Chothani et al.
(Molecular Cell 2022) R reference consumes.  Without it, comparing our numbers
against published ones means comparing across two different coverage models as
well as two implementations, and any disagreement is uninterpretable.

Same capability profile as BigwigSetProvider — a bedGraph is per-base depth:
- ``coverage(…, site=…)`` accepts the parameter but IGNORES it (no read lengths,
  so no offset shift is possible).
- ``junction_support()`` raises NotImplementedError (no CIGAR).
- ``size_factors()`` returns 1.0 per file (caller normalises externally).

Coordinates: bedGraph is BED-style, 0-based half-open ``[start, end)``, and
``coverage()`` emits 0-based positions to match BigwigSetProvider.  The R
reference works in 1-based transcript coordinates and tests
``bedgraph[,2] < x & bedgraph[,3] >= x``; for ``x = p + 1`` that is exactly
``start <= p < end``, so the two conventions agree.

Sparse by design: positions with no interval are simply absent from the
returned frame, which is what a bedGraph gap means (zero signal).  Note this is
NOT what the R does — see ``scoring.signature.signal_from_intervals`` for the
``unlist()`` behaviour that drops uncovered positions and shifts frame.

Strand convention matches BigwigSetProvider: a stranded entry is a dict with
'forward'/'reverse' keys naming two files that both hold POSITIVE values.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import polars as pl

from TranslonScorer.model import Region

PathLike = Union[str, Path]

# chrom -> (starts, ends, values), each sorted by start
_ChromIndex = Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]


def _read_bedgraph(path: PathLike) -> _ChromIndex:
    """Load a bedGraph into per-chromosome sorted interval arrays."""
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    chroms: List[str] = []
    starts: List[int] = []
    ends: List[int] = []
    values: List[float] = []
    with opener(path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                raise ValueError(f"malformed bedGraph line in {path}: {line!r}")
            chroms.append(f[0])
            starts.append(int(f[1]))
            ends.append(int(f[2]))
            values.append(float(f[3]))
    if not chroms:
        raise ValueError(f"bedGraph is empty: {path}")

    index: _ChromIndex = {}
    chrom_arr = np.asarray(chroms)
    start_arr = np.asarray(starts, dtype=np.int64)
    end_arr = np.asarray(ends, dtype=np.int64)
    value_arr = np.asarray(values, dtype=float)
    for chrom in np.unique(chrom_arr):
        sel = chrom_arr == chrom
        s, e, v = start_arr[sel], end_arr[sel], value_arr[sel]
        order = np.argsort(s, kind="stable")
        index[str(chrom)] = (s[order], e[order], v[order])
    return index


class BedgraphSetProvider:
    """CoverageProvider backed by one or more bedGraph files.

    Files are sum-merged across the set, which is how the reference combines
    the five atomic runs in a group: the signal vectors are summed
    position-wise *before* scoring, never by averaging already-computed
    PIF/CIF values.  Passing the five runs as five entries with
    ``by_sample=False`` gives exactly that.

    Parameters
    ----------
    bedgraphs   : paths, or dicts with 'forward'/'reverse' keys when stranded.
    sample_names: optional labels (one per entry).
    stranded    : if True every entry must be a {'forward','reverse'} dict and
                  coverage() emits a `strand` column (1 forward, -1 reverse).
    """

    def __init__(
        self,
        bedgraphs: List[Union[PathLike, Dict[str, PathLike]]],
        *,
        sample_names: Optional[List[str]] = None,
        stranded: bool = False,
    ) -> None:
        if not bedgraphs:
            raise ValueError("at least one bedGraph is required")
        self._bedgraphs = bedgraphs
        self._sample_names = sample_names or [
            (Path(b).stem if isinstance(b, (str, Path)) else f"sample_{i}")
            for i, b in enumerate(bedgraphs)
        ]
        self._stranded = stranded
        self._cache: Dict[str, _ChromIndex] = {}

    def _index(self, path: PathLike) -> _ChromIndex:
        key = str(path)
        if key not in self._cache:
            self._cache[key] = _read_bedgraph(path)
        return self._cache[key]

    # ------------------------------------------------------------------
    # CoverageProvider
    # ------------------------------------------------------------------

    def coverage(
        self,
        regions: List[Region],
        *,
        site: str = "A",
        by_sample: bool = False,
    ) -> pl.DataFrame:
        """Per-position coverage summed across bedGraph files.

        ``site`` is accepted and ignored — a bedGraph carries no read-length
        information, so no P/A-site shift can be applied.
        """
        if site not in {"P", "A"}:
            raise ValueError(f"site must be 'P' or 'A', got {site!r}")

        schema: Dict[str, type] = {"pos": pl.Int64, "count": pl.Float64}
        if self._stranded:
            schema["strand"] = pl.Int64
        if by_sample:
            schema["sample_id"] = pl.Utf8

        frames: List[pl.DataFrame] = []
        for entry, sample_id in zip(self._bedgraphs, self._sample_names):
            if self._stranded:
                if not isinstance(entry, dict) or not {"forward", "reverse"} <= entry.keys():
                    raise ValueError(
                        "stranded=True requires every bedGraph entry to be a dict with "
                        f"'forward'/'reverse' keys; got {entry!r} for sample {sample_id!r}"
                    )
                sources: List[Tuple[PathLike, Optional[int]]] = [
                    (entry["forward"], 1),
                    (entry["reverse"], -1),
                ]
            elif isinstance(entry, (str, Path)):
                sources = [(entry, None)]
            else:
                sources = [(v, None) for v in entry.values() if v is not None]

            for path, strand_val in sources:
                index = self._index(path)
                for region in regions:
                    hit = index.get(region.chrom)
                    if hit is None:
                        continue
                    starts, ends, values = hit
                    pos, cnt = _expand_region(starts, ends, values, region.start, region.end)
                    if pos.size == 0:
                        continue
                    frame = pl.DataFrame({"pos": pos, "count": cnt})
                    if self._stranded:
                        frame = frame.with_columns(
                            pl.lit(strand_val, dtype=pl.Int64).alias("strand")
                        )
                    if by_sample:
                        frame = frame.with_columns(pl.lit(sample_id).alias("sample_id"))
                    frames.append(frame)

        if not frames:
            return pl.DataFrame(schema=schema)

        df = pl.concat(frames)
        group_cols = (
            ["pos"] + (["strand"] if self._stranded else []) + (["sample_id"] if by_sample else [])
        )
        return df.group_by(group_cols).agg(pl.col("count").sum()).sort(group_cols)

    def size_factors(self) -> Dict[str, float]:
        """Return 1.0 per file (external normalisation assumed)."""
        return {s: 1.0 for s in self._sample_names}

    def junction_support(self, *args, **kwargs):  # type: ignore[override]
        raise NotImplementedError(
            "BedgraphSetProvider does not support junction_support(): "
            "bedGraph files contain no read alignment structure (CIGAR)."
        )

    def mappability_ledger(self, events: pl.DataFrame) -> pl.DataFrame:
        raise NotImplementedError(
            "BedgraphSetProvider does not support mappability_ledger(): "
            "bedGraph files contain no unique/multimapper information."
        )


def _expand_region(
    starts: np.ndarray,
    ends: np.ndarray,
    values: np.ndarray,
    region_start: int,
    region_end: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Expand the intervals overlapping [region_start, region_end) to positions.

    Only nonzero values are emitted, matching BigwigSetProvider's sparse
    output; absent positions mean zero signal.
    """
    if region_end <= region_start or starts.size == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=float)
    # Intervals with start < region_end and end > region_start.
    lo = 0
    hi = int(np.searchsorted(starts, region_end, side="left"))
    if hi <= lo:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=float)
    s = np.maximum(starts[lo:hi], region_start)
    e = np.minimum(ends[lo:hi], region_end)
    v = values[lo:hi]
    keep = (e > s) & (v != 0)
    if not np.any(keep):
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=float)
    s, e, v = s[keep], e[keep], v[keep]
    lengths = (e - s).astype(np.int64)
    pos = np.repeat(s, lengths) + _within_offsets(lengths)
    cnt = np.repeat(v, lengths)
    return pos, cnt


def _within_offsets(lengths: np.ndarray) -> np.ndarray:
    """0,1,..,n-1 concatenated for each n in lengths, without a Python loop."""
    total = int(lengths.sum())
    if total == 0:
        return np.empty(0, dtype=np.int64)
    offsets = np.repeat(np.cumsum(lengths) - lengths, lengths)
    return np.arange(total, dtype=np.int64) - offsets

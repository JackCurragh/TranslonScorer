from __future__ import annotations
"""
Zarr + unique-read index loader utilities.

This module streams normalized read chunks from a Zarr counts array and a
Parquet alignment index. It emits DataFrames with columns:
  chr, start, stop, length, strand, count

Assumptions:
- Parquet index has columns: read_id (int or str), chr, start, stop, strand, length.
- Zarr root contains an array "counts" shaped (n_samples, n_reads) or (n_reads, n_samples).

These helpers are intentionally simple and functional; they stream in chunks and
avoid materializing the full matrix.
"""

from typing import Iterator, List, Tuple, Optional, Dict

import polars as pl
import importlib

_zarr_mod = None
_zarr_err = None

def _require_zarr():
    global _zarr_mod, _zarr_err
    if _zarr_mod is None and _zarr_err is None:
        try:
            _zarr_mod = importlib.import_module('zarr')
        except Exception as e:
            _zarr_err = e
    if _zarr_mod is None:
        raise ImportError(
            "Zarr is required for Zarr-based profiles. Install with pip: "
            "pip install 'zarr>=2.16' 'numcodecs>=0.12'"
        ) from _zarr_err


def _open_counts(zroot: str):
    """Open and return the counts array from a Zarr root path.

    Handles stores that either:
      - have a root group with a child array named 'counts' (recommended), or
      - point directly to the array path (e.g. ".../global_matrix.zarr/counts"), or
      - have a single top-level array (fallback to the first child).
    """
    _require_zarr()
    # First, try to open as a group and access 'counts'
    try:
        root = _zarr_mod.open(zroot, mode="r")
        # If this is already an Array (path points directly to an array), return it
        try:
            from zarr.core import Array as _ZarrArray  # zarr 2.x
        except Exception:
            _ZarrArray = None  # type: ignore
        if _ZarrArray is not None and isinstance(root, _ZarrArray):
            return root
        # Otherwise, it's a group-like; prefer 'counts' child
        if hasattr(root, "__contains__") and ("counts" in root):
            return root["counts"]
        # If no 'counts', but there are children, pick the first child
        keys = list(root.keys()) if hasattr(root, "keys") else []
        if keys:
            # Prefer a child literally named 'counts' if present
            if 'counts' in keys:
                return root['counts']
            return root[keys[0]]
    except Exception:
        # Could not open as group; fall through to direct array paths
        pass

    # Try an explicit '/counts' child path
    try:
        return _zarr_mod.open(zroot.rstrip('/') + '/counts', mode='r')
    except Exception:
        # Final attempt: open the given path as an array
        arr = _zarr_mod.open(zroot, mode='r')
        return arr


def _counts_is_samples_first(arr) -> bool:
    # Heuristic: assume samples dimension is the smaller one for typical setups
    return arr.shape[0] <= arr.shape[1]


def iter_reads_from_zarr(
    zarr_root: str,
    read_index_parquet: str,
    samples: List[str],
    *,
    sample_to_index: dict | None = None,
    chunk_size: int = 1_000_000,
    include_read_id: bool = False,
) -> Iterator[Tuple[str, pl.DataFrame]]:
    """
    Stream normalized read chunks for each sample from Zarr + Parquet index.

    Yields: (sample_name, DataFrame[chr,start,stop,length,strand,count])
    """
    arr = _open_counts(zarr_root)
    samples_first = _counts_is_samples_first(arr)

    # Map sample names to indices if provided; otherwise assume integer names
    if sample_to_index is None:
        # If no mapping provided, assume samples are addressed by integer indices in order
        sample_to_index = {s: i for i, s in enumerate(samples)}

    # Lazy scan the index to get total rows
    scan = pl.scan_parquet(read_index_parquet).select(pl.len())
    n_rows = scan.collect().item()

    # Iterate in read_id chunks
    for start in range(0, n_rows, chunk_size):
        end = min(start + chunk_size, n_rows)
        # Load index slice with read_id
        idx_df = (
            pl.scan_parquet(read_index_parquet)
            .slice(start, end - start)
            .select(["read_id", "chr", "start", "stop", "strand", "length"])  # keep only needed
            .collect()
        )

        # Skip empty
        if idx_df.is_empty():
            continue

        # For each sample, slice counts and emit
        for s in samples:
            si = sample_to_index[s]
            read_ids = idx_df.get_column('read_id').to_numpy()
            if samples_first:
                counts = arr.get_orthogonal_selection((si, read_ids))
            else:
                counts = arr.get_orthogonal_selection((read_ids, si))

            # Build DF with counts
            if include_read_id:
                df = idx_df.with_columns(pl.Series("count", counts)).select(["read_id","chr","start","stop","strand","length","count"])  # keep read_id
            else:
                df = idx_df.with_columns(pl.Series("count", counts)).select(["chr","start","stop","strand","length","count"])  # drop read_id

            yield s, df

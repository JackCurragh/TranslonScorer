from __future__ import annotations

import importlib
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
import polars as pl

from ..utils import log_info, log_warning

_zarr_mod = None
_zarr_err = None


def _require_zarr():
    global _zarr_mod, _zarr_err
    if _zarr_mod is None and _zarr_err is None:
        try:
            _zarr_mod = importlib.import_module("zarr")
        except Exception as e:
            _zarr_err = e
    if _zarr_mod is None:
        raise ImportError(
            "Zarr is required for locus profile export. Install with pip: "
            "pip install 'zarr>=2.16' 'numcodecs>=0.12'"
        ) from _zarr_err
    return _zarr_mod


def _load_offsets_dict(path: Optional[str], default_offset: int = 15) -> Dict[int, int]:
    if not path:
        return {}
    try:
        df = pl.read_csv(path)
        # Expect columns: length, offset
        return {int(r[0]): int(r[1]) for r in df.select([df.columns[0], df.columns[1]]).iter_rows()}
    except Exception as e:
        log_warning(f"Failed to read offsets file {path}: {e}; using default {default_offset}")
        return {}


def _offset_for(length: int, offsets: Dict[int, int], default_offset: int) -> int:
    try:
        return int(offsets.get(int(length), default_offset))
    except Exception:
        return default_offset


def _a_site_positions(
    starts: np.ndarray,
    stops: np.ndarray,
    strands: List[str],
    offsets: np.ndarray,
) -> np.ndarray:
    """Strand-aware A-site genomic positions.

    The offset is the distance from the read's 5' end, so the genomic
    direction depends on strand: plus-strand reads advance from ``start``,
    minus-strand reads retreat from the genomic-right 5' end (``stop - 1``).
    Using ``start + offset`` for both strands (the pre-fix behaviour) shifts
    every reverse-strand profile by the whole footprint length and destroys
    its frame assignment -- a ~15nt error on a typical RPF, not a rounding
    quirk.
    """
    is_minus = np.asarray([str(s) == "-" for s in strands], dtype=bool)
    return np.where(is_minus, stops - 1 - offsets, starts + offsets)


def _iter_loci_from_bed(bed_path: str) -> Iterator[Tuple[str, str, int, int]]:
    """Yield (locus_id, chr, start, stop) from a 4+ column BED-like file."""
    df = pl.read_csv(bed_path, has_header=False, separator="\t")
    # Columns: chrom, start, end, name, ...
    for row in df.iter_rows():
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        name = str(row[3]) if len(row) > 3 else f"{chrom}:{start}-{end}"
        yield name, chrom, start, end


def build_locus_profiles_zarr(
    zarr_root: str,
    read_index_parquet: str,
    samples: List[str],
    loci_bed: str,
    *,
    sample_to_index: Optional[Dict[str, int]] = None,
    offsets_file: Optional[str] = None,
    default_offset: int = 15,
    out_zarr: str = "locus_profiles.zarr",
) -> str:
    """Build a Zarr with per-locus A-site genomic profiles: profiles[sample, position].

    - X axis (columns) = genomic positions within the locus (contiguous, start..end-1)
    - Y axis (rows) = samples (order matches provided samples list)
    - Values = counts at A-site position

    Zarr layout:
      /samples               (1D utf-8 array of sample names)
      /loci/{locus_id}/pos   (1D int64 genomic positions)
      /loci/{locus_id}/profiles  (2D float32 [n_samples, len(pos)])
    """
    # Open inputs
    zarr = _require_zarr()
    store = zarr.open(zarr_root, mode="r")
    if "counts" not in store:
        raise ValueError("Zarr root does not contain 'counts' array")
    counts = store["counts"]
    samples_first = counts.shape[0] <= counts.shape[1]

    # Map sample names to indices
    if sample_to_index is None:
        sample_to_index = {s: i for i, s in enumerate(samples)}

    # Write output Zarr
    out = zarr.open(out_zarr, mode="a")
    # Save samples once
    if "samples" in out:
        del out["samples"]
    out.create_dataset(
        "samples", data=np.array(samples, dtype=object), dtype=object, overwrite=True
    )

    # Prepare index scan
    scan_idx = pl.scan_parquet(read_index_parquet).select(
        ["read_id", "chr", "start", "stop", "strand", "length"]
    )

    # Load offsets mapping
    offsets = _load_offsets_dict(offsets_file, default_offset=default_offset)

    # Iterate loci
    for locus_id, chrom, locus_start, locus_end in _iter_loci_from_bed(loci_bed):
        log_info(f"Building locus {locus_id} {chrom}:{locus_start}-{locus_end}")
        locus_len = int(locus_end - locus_start)
        if locus_len <= 0:
            continue

        # Create/overwrite datasets for this locus
        grp = out.require_group(f"loci/{locus_id}")
        # Genomic coordinate axis
        pos_vec = np.arange(locus_start, locus_end, dtype=np.int64)
        if "pos" in grp:
            del grp["pos"]
        grp.create_dataset("pos", data=pos_vec, dtype="i8", overwrite=True)

        # Profiles array
        if "profiles" in grp:
            del grp["profiles"]
        prof = grp.create_dataset(
            "profiles",
            shape=(len(samples), locus_len),
            chunks=(min(len(samples), 64), min(locus_len, 16384)),
            dtype="f4",
            overwrite=True,
        )

        # Get read_ids overlapping locus from index (predicate pushdown)
        idx = scan_idx.filter(
            (pl.col("chr") == chrom)
            & (pl.col("start") < locus_end)
            & (pl.col("stop") > locus_start)
        ).collect()
        if idx.is_empty():
            continue

        lengths = idx.get_column("length").to_list()
        ofs = np.array([_offset_for(L, offsets, default_offset) for L in lengths], dtype=np.int64)
        a_site = _a_site_positions(
            idx.get_column("start").to_numpy(),
            idx.get_column("stop").to_numpy(),
            idx.get_column("strand").to_list(),
            ofs,
        )
        # Keep those within locus bounds
        mask = (a_site >= locus_start) & (a_site < locus_end)
        if not mask.any():
            continue
        idx = idx.filter(pl.Series(mask))
        a_site = a_site[mask]

        # Column indices inside locus
        col_idx = (a_site - locus_start).astype(np.int64)

        # Read IDs for slicing counts
        read_ids = idx.get_column("read_id").to_numpy()

        # For each sample, slice counts and scatter-add into profiles
        for s_i, s_name in enumerate(samples):
            si = sample_to_index[s_name]
            if samples_first:
                c = counts.get_orthogonal_selection((si, read_ids))
            else:
                c = counts.get_orthogonal_selection((read_ids, si))
            # Accumulate into locus vector
            # Use np.add.at for scatter add
            vec = np.zeros(locus_len, dtype=np.float32)
            np.add.at(vec, col_idx, c.astype(np.float32))
            prof[s_i, :] = vec

    return out_zarr

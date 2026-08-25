"""BigWig reading: coverage extraction and transcript-space projection.

Reads bigWigs and maps them into transcript coordinates.  Scoring lives on the
spine (``coverage/bigwig.py`` provides coverage to ``scoring/run.py``); this
module deliberately contains no scoring path of its own.
"""

import importlib
from typing import Any, List, Tuple

import polars as pl

_bw_mod = None
_bw_err = None


def _require_pybigwig():
    global _bw_mod, _bw_err
    if _bw_mod is None and _bw_err is None:
        try:
            _bw_mod = importlib.import_module("pyBigWig")
        except Exception as e:
            _bw_err = e
    if _bw_mod is None:
        raise ImportError(
            "pyBigWig failed to import. This often indicates a NumPy ABI mismatch. "
            "For a pip-only setup, install compatible wheels:\n"
            "  pip install --upgrade 'numpy<2' 'pyBigWig>=0.3.22'\n"
            "Then reinstall this package if needed (pip install -e .)."
        ) from _bw_err
    return _bw_mod


import gc
from contextlib import contextmanager
from dataclasses import dataclass

import numpy as np

from ..utils.logging import log_error, log_warning


@dataclass
class ProcessingConfig:
    """Configuration for transcript processing."""

    max_workers: int = 0  # 0 means auto-detect based on CPU count
    batch_size: int = 500  # Number of transcripts to process in each batch
    sru_range: int = 100  # Range for SRU score calculation
    max_region_size: int = 1_000_000  # Maximum region size to process at once
    chunk_size: int = 100_000  # Size of chunks for processing large regions
    use_temp_files: bool = True  # Use temporary files for large intermediate results
    max_regions_per_worker: int = 500  # Maximum regions to process per worker
    prefetch_bigwig: bool = True  # Prefetch BigWig data for common chromosomes
    use_shared_data: bool = True  # Use shared memory for data when possible


@contextmanager
def open_bigwig(path: str):
    """Safely open and close a BigWig file."""
    bw_file = None
    try:
        bw_module = _require_pybigwig()
        bw_file = bw_module.open(path)
        if not bw_file.isBigWig():
            raise ValueError(f"File {path} is not a valid bigWig file")
        yield bw_file
    finally:
        if bw_file:
            bw_file.close()


class BigWigCache:
    """Cache for BigWig data to reduce repeated file access."""

    def __init__(self, bigwig_path, max_cache_size=1000):
        self.bigwig_path = bigwig_path
        self.max_cache_size = max_cache_size
        self.cache = {}
        self.chroms = set()
        self.hits = 0
        self.misses = 0

        # Initialize with chromosome information
        with open_bigwig(bigwig_path) as bwfile:
            self.chroms = set(bwfile.chroms().keys())

    def get_values(self, chrom, start, stop):
        """Get values for a region, using cache if available."""
        if chrom not in self.chroms:
            return None

        # For very large regions, don't cache
        if stop - start > 100000:
            with open_bigwig(self.bigwig_path) as bwfile:
                return bwfile.values(chrom, start, stop)

        # Check cache
        cache_key = (chrom, start, stop)
        if cache_key in self.cache:
            self.hits += 1
            return self.cache[cache_key]

        # Fetch from file
        self.misses += 1
        with open_bigwig(self.bigwig_path) as bwfile:
            values = bwfile.values(chrom, start, stop)

            # Manage cache size with LRU-like behavior
            if len(self.cache) >= self.max_cache_size:
                # Remove oldest item (simple approximation)
                self.cache.pop(next(iter(self.cache)))

            self.cache[cache_key] = values
            return values


def _flatten_list_or_value(value: Any) -> List:
    """
    Helper function to flatten a value that might be a list or a single value.
    Returns a list in either case.
    """
    if isinstance(value, list):
        return value
    else:
        return [value]


def extract_regions(exon_df: pl.DataFrame, orf_tran_ids=None) -> List[Tuple[str, int, int, str]]:
    """
    Extract regions (chr, start, stop, tran_id) from exon DataFrame.
    Handles both list and scalar values.

    Args:
        exon_df: DataFrame with exon data
        orf_tran_ids: Optional set of transcript IDs to filter by

    Returns:
        List of (chr, start, stop, tran_id) tuples
    """
    regions = []

    # Process each row in the exon dataframe
    for row in exon_df.iter_rows(named=True):
        chrom = row["chr"]
        tran_id = row["tran_id"]

        # Skip if we're filtering by transcript ID and this one isn't in the list
        if orf_tran_ids is not None and tran_id not in orf_tran_ids:
            continue

        # Handle both single values and lists for start/stop
        starts = _flatten_list_or_value(row["start"])
        stops = _flatten_list_or_value(row["stop"])

        # Make sure we have matching start/stop pairs
        if len(starts) != len(stops):
            log_warning(
                f"Mismatched start/stop lists for {tran_id} on {chrom}: {len(starts)} starts, {len(stops)} stops"
            )
            continue

        # Process each start/stop pair
        for start, stop in zip(starts, stops):
            # Ensure start and stop are integers
            try:
                start = int(start)
                stop = int(stop)
            except (ValueError, TypeError):
                continue

            if start >= stop:
                continue

            regions.append((chrom, start, stop, tran_id))

    return regions


def process_batch(batch_data):
    """
    Process a batch of regions - module level function for multiprocessing

    Args:
        batch_data: Tuple containing (regions, bigwig_path)

    Returns:
        Dictionary with processing results
    """
    regions, bigwig_path = batch_data

    # Use an optimized cache size based on batch size
    cache_size = min(2000, max(500, len(regions) // 4))
    bw_cache = BigWigCache(bigwig_path, max_cache_size=cache_size)

    # Group regions by transcript AND chromosome for better locality
    regions_by_transcript = {}
    for region in regions:
        chrom, start, stop, tran_id = region
        if tran_id not in regions_by_transcript:
            regions_by_transcript[tran_id] = {}

        if chrom not in regions_by_transcript[tran_id]:
            regions_by_transcript[tran_id][chrom] = []

        regions_by_transcript[tran_id][chrom].append((start, stop))

    batch_results = {}
    processed = 0
    failed = 0

    # Process each transcript's regions
    for tran_id, chroms in regions_by_transcript.items():
        # Store as sparse representation
        positions = []
        values = []
        max_position = 0

        # Process each chromosome separately for better BigWig cache efficiency
        for chrom, regions in chroms.items():
            # Sort regions by position for better sequential access
            regions.sort()

            current_position = len(positions)

            for start, stop in regions:
                try:
                    stop - start
                    raw_values = bw_cache.get_values(chrom, start, stop)

                    if raw_values is None:
                        # Region not in bigwig
                        processed += 1
                        continue

                    # Only store non-zero values to save memory
                    for i, v in enumerate(raw_values):
                        if v is not None and v > 0:
                            positions.append(current_position + i)
                            values.append(v)

                    # Update max position and current position
                    max_position = max(max_position, current_position + len(raw_values))
                    current_position += len(raw_values)
                    processed += 1
                except Exception as e:
                    log_error(f"Error processing region {chrom}:{start}-{stop}: {str(e)}")
                    failed += 1
                    current_position += stop - start  # Still increment position

        # Only store if we have data
        if positions and values:
            batch_results[tran_id] = {
                "positions": positions,
                "values": values,
                "max_pos": max_position,
            }

    # Clear references to help garbage collection
    bw_cache = None
    regions_by_transcript = None
    gc.collect()

    return {"results": batch_results, "processed": processed, "failed": failed}


def transcriptreads(bigwig_file, exon_df, transcript_id=None):
    """
    Extract transcript coverage from a BigWig file, supporting both
    genomic and transcriptomic BigWig inputs.

    Parameters:
    ----------
    bigwig_file : str or pyBigWig.pyBigWig
        Path to a BigWig file or an open BigWig handle
    exon_df : polars.DataFrame
        DataFrame containing exon annotations with columns: chr, start, stop, tran_id
    transcript_id : str, optional
        If provided, only process this specific transcript

    Returns:
    -------
    polars.DataFrame
        DataFrame containing transcript coverage with columns: tran_start, counts
    """
    bw_module = _require_pybigwig()
    import polars as pl

    # Open BigWig file if a path was given
    if isinstance(bigwig_file, str):
        bw_handle = bw_module.open(bigwig_file)
        need_close = True
    else:
        bw_handle = bigwig_file
        need_close = False

    try:
        # Get chromosomes in BigWig
        bw_chroms = set(bw_handle.chroms().keys())

        # Filter exon_df for specific transcript if requested
        if transcript_id:
            exon_df = exon_df.filter(pl.col("tran_id") == transcript_id)

        # Check BigWig type (genomic or transcriptomic)
        exon_trans = set(exon_df["tran_id"].unique())
        exon_chroms = set(exon_df["chr"].unique())

        trans_match = len(bw_chroms.intersection(exon_trans))
        chrom_match = len(bw_chroms.intersection(exon_chroms))

        # Determine if this is a genomic or transcriptomic BigWig
        is_genomic = chrom_match >= trans_match
        if is_genomic:
            return process_genomic_bigwig(bw_handle, exon_df)
        else:
            return process_transcriptomic_bigwig(bw_handle, exon_df)

    finally:
        # Close BigWig handle if we opened it
        if need_close:
            bw_handle.close()


def process_genomic_bigwig(bw_handle, exon_df):
    """
    Process a genomic BigWig file and map to transcript coordinates.

    Parameters:
    ----------
    bw_handle : pyBigWig.pyBigWig
        Open BigWig handle
    exon_df : polars.DataFrame
        DataFrame containing exon annotations

    Returns:
    -------
    polars.DataFrame
        DataFrame with transcript coordinates and coverage values
    """
    import polars as pl

    from ..utils.logging import log_warning

    # Get chromosomes in BigWig
    bw_chroms = set(bw_handle.chroms().keys())

    # Create lists to store transcript data
    all_tran_ids = []
    all_tran_starts = []
    all_counts = []
    # Simple in-memory cache to avoid re-reading identical exon intervals shared across isoforms
    _cache: dict[tuple[str, int, int], np.ndarray] = {}
    _cache_cap = 50000

    # Process each transcript
    for tran_id in exon_df["tran_id"].unique():
        tran_exons = exon_df.filter(pl.col("tran_id") == tran_id)

        # Skip if no exons
        if tran_exons.is_empty():
            continue

        # Get chromosome for this transcript
        chrom = tran_exons["chr"][0]

        # Skip if chromosome not in BigWig
        if chrom not in bw_chroms:
            # Try without 'chr' prefix
            if chrom.startswith("chr") and chrom[3:] in bw_chroms:
                chrom = chrom[3:]
            # Try with 'chr' prefix
            elif f"chr{chrom}" in bw_chroms:
                chrom = f"chr{chrom}"
            else:
                continue

        # Collect sparse positions/values directly (vectorized)
        sparse_pos: list[int] = []
        sparse_vals: list[float] = []

        # Process each exon in this transcript
        for row in tran_exons.iter_rows(named=True):
            # Determine strand for position mapping
            strand = row.get("strand", "+") or "+"

            # Get genomic coordinates
            if isinstance(row["start"], list):
                # Multiple exons case
                for i in range(len(row["start"])):
                    start = int(row["start"][i])
                    stop = int(row["stop"][i])
                    tran_start = int(row["tran_start"][i])
                    # Guard against zero/negative length intervals
                    if stop <= start:
                        log_warning(
                            f"Skipping invalid exon interval {chrom}:{start}-{stop} (zero/negative length)"
                        )
                        continue

                    try:
                        key = (chrom, start, stop)
                        arr = _cache.get(key)
                        if arr is None:
                            vals = bw_handle.values(chrom, start, stop, numpy=True)
                            arr = np.asarray(vals, dtype=np.float32)
                            if len(_cache) >= _cache_cap:
                                _cache.clear()
                            _cache[key] = arr
                        mask = np.isfinite(arr) & (arr != 0.0)
                        if mask.any():
                            idx = np.nonzero(mask)[0]
                            exon_len = stop - start
                            if strand == "-":
                                # For – strand: genomic index k → transcript pos tran_start + (exon_len - 1 - k)
                                tran_idx = (tran_start + (exon_len - 1 - idx)).astype(int)
                            else:
                                tran_idx = (tran_start + idx).astype(int)
                            sparse_pos.extend(tran_idx.tolist())
                            sparse_vals.extend(arr[mask].astype(float).tolist())
                    except Exception as e:
                        log_warning(f"Error reading {chrom}:{start}-{stop}: {str(e)}")
            else:
                # Single exon case
                start = int(row["start"])
                stop = int(row["stop"])
                tran_start = int(row["tran_start"])
                # Guard against zero/negative length intervals
                if stop <= start:
                    log_warning(
                        f"Skipping invalid exon interval {chrom}:{start}-{stop} (zero/negative length)"
                    )
                    continue

                try:
                    key = (chrom, start, stop)
                    arr = _cache.get(key)
                    if arr is None:
                        vals = bw_handle.values(chrom, start, stop, numpy=True)
                        arr = np.asarray(vals, dtype=np.float32)
                        if len(_cache) >= _cache_cap:
                            _cache.clear()
                        _cache[key] = arr
                    mask = np.isfinite(arr) & (arr != 0.0)
                    if mask.any():
                        idx = np.nonzero(mask)[0]
                        exon_len = stop - start
                        if strand == "-":
                            tran_idx = (tran_start + (exon_len - 1 - idx)).astype(int)
                        else:
                            tran_idx = (tran_start + idx).astype(int)
                        sparse_pos.extend(tran_idx.tolist())
                        sparse_vals.extend(arr[mask].astype(float).tolist())
                except Exception as e:
                    log_warning(f"Error reading {chrom}:{start}-{stop}: {str(e)}")

        # Add coverage data for this transcript
        if sparse_pos:
            n = len(sparse_pos)
            all_tran_ids.extend([tran_id] * n)
            all_tran_starts.extend(sparse_pos)
            all_counts.extend(sparse_vals)

    # Create DataFrame from collected data
    if all_tran_ids:
        return pl.DataFrame(
            {"tran_id": all_tran_ids, "tran_start": all_tran_starts, "counts": all_counts}
        )
    else:
        # Return empty DataFrame with correct structure
        return pl.DataFrame({"tran_id": [], "tran_start": [], "counts": []})


def process_transcriptomic_bigwig(bw_handle, exon_df):
    """
    Process a transcriptomic BigWig file.

    Parameters:
    ----------
    bw_handle : pyBigWig.pyBigWig
        Open BigWig handle
    exon_df : polars.DataFrame
        DataFrame containing exon annotations

    Returns:
    -------
    polars.DataFrame
        DataFrame with transcript coordinates and coverage values
    """
    import polars as pl

    from ..utils.logging import log_warning

    # Get transcripts in BigWig
    bw_trans = set(bw_handle.chroms().keys())

    # Create a dictionary to store transcript coverage
    results = []

    # Process each transcript
    for tran_id in exon_df["tran_id"].unique():
        # Skip if transcript not in BigWig
        if tran_id not in bw_trans:
            continue

        try:
            # Get transcript size from BigWig
            tran_size = bw_handle.chroms()[tran_id]
            # Get coverage values for the entire transcript
            arr = np.asarray(bw_handle.values(tran_id, 0, tran_size, numpy=True), dtype=np.float32)
            mask = np.isfinite(arr) & (arr != 0.0)
            if mask.any():
                pos = np.nonzero(mask)[0].astype(int).tolist()
                vals = arr[mask].astype(float).tolist()
                tran_df = pl.DataFrame(
                    {
                        "tran_id": [tran_id] * len(pos),
                        "tran_start": pos,
                        "counts": vals,
                    }
                )
                results.append(tran_df)
        except Exception as e:
            log_warning(f"Error reading transcript {tran_id}: {str(e)}")

    # Combine all transcripts
    if results:
        return pl.concat(results)
    else:
        # Return empty DataFrame with correct structure
        return pl.DataFrame({"tran_id": [], "tran_start": [], "counts": []})


# Add this function outside of process_strand_orfs


# Then modify process_strand_orfs to use this function

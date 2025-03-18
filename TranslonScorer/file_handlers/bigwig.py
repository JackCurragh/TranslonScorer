"""
BigWig file handling functionality for TranslonScorer with high-performance optimization.

This module contains optimized functions for reading and processing BigWig files,
including conversion to other formats and coordinate transformations.
"""

from typing import Dict, List, Optional, Union, Tuple, Any, Set
import polars as pl
import pyBigWig as bw
from ..utils.logging import log_info, log_warning, log_error
from ..core.scoring import oldscoring, newscoring, globalscores, existingscore, assigningscore

import concurrent.futures
from functools import partial
import os
import time
import numpy as np
from dataclasses import dataclass, field
from contextlib import contextmanager
import tempfile
import pickle
import gc
from collections import defaultdict


@dataclass
class ProcessingConfig:
    """Configuration for transcript processing."""
    max_workers: int = 0  # 0 means auto-detect based on CPU count
    batch_size: int = 50  # Number of transcripts to process in each batch
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
        bw_file = bw.open(path)
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
            log_warning(f"Mismatched start/stop lists for {tran_id} on {chrom}: {len(starts)} starts, {len(stops)} stops")
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
                    region_size = stop - start
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
                    current_position += (stop - start)  # Still increment position
        
        # Only store if we have data
        if positions and values:
            batch_results[tran_id] = {
                'positions': positions,
                'values': values,
                'max_pos': max_position
            }
    
    # Clear references to help garbage collection
    bw_cache = None
    regions_by_transcript = None
    gc.collect()
    
    return {"results": batch_results, "processed": processed, "failed": failed}


def score_transcript(args):
    """
    Score a single transcript - module level function for multiprocessing
    
    Args:
        args: Tuple containing (tran_id, transcript_data, orf_data, old_scoring, sru_range)
        
    Returns:
        Scored ORF DataFrame or None
    """
    tran_id, tran_data, orf_data, old_scoring, sru_range = args
    
    try:
        start_time = time.time()  # Start timing the function
        
        # Create DataFrame from transcript data
        if isinstance(tran_data, dict) and 'positions' in tran_data:
            # Convert from sparse representation
            max_pos = tran_data['max_pos'] + 1
            full_array = np.zeros(max_pos, dtype=np.float32)
            
            # Fill in non-zero values
            fill_start_time = time.time()  # Start timing this stage
            for pos, val in zip(tran_data['positions'], tran_data['values']):
                if pos < max_pos:  # Safety check
                    full_array[pos] = val
            fill_duration = time.time() - fill_start_time
            log_info(f"Time taken to fill array for {tran_id}: {fill_duration:.4f} seconds")
            
            # Create DataFrame from full array
            tran_reads = pl.DataFrame({
                "tran_start": pl.arange(0, max_pos, eager=True),
                "counts": full_array
            })
        elif isinstance(tran_data, list):
            # Convert from simple list
            tran_reads = pl.DataFrame({
                "tran_start": pl.arange(0, len(tran_data), eager=True),
                "counts": tran_data
            })
        else:
            # Invalid data format
            return None
        
        # Create ORF DataFrame
        orfs_start_time = time.time()  # Start timing ORF DataFrame creation
        orfs = pl.DataFrame(orf_data)
        orf_duration = time.time() - orfs_start_time
        log_info(f"Time taken to create ORF DataFrame for {tran_id}: {orf_duration:.4f} seconds")
        
        if orfs.is_empty() or tran_reads.is_empty():
            return None
        
        results = []
        
        # Score each ORF type separately
        for typeorf in orfs["type"].unique():
            orfs_filtered = orfs.filter(pl.col("type") == typeorf)
            
            if orfs_filtered.is_empty():
                continue
                
            try:
                scoring_start_time = time.time()  # Start timing scoring
                if old_scoring:
                    # Use classic scoring method
                    scored_orfs = oldscoring(
                        orfs_filtered, tran_reads, sru_range, typeorf
                    )
                    results.append(scored_orfs)
                else:
                    # Use modern scoring method with score caching
                    empty_dict = {"rise_up": defaultdict(float), "step_down": defaultdict(float)}
                    
                    # Find ORFs that need scoring
                    emptyscore_df = existingscore(orfs_filtered, typeorf, empty_dict)
                    
                    if not emptyscore_df.is_empty():
                        # Calculate new scores
                        scoredict = newscoring(
                            emptyscore_df, tran_reads, sru_range, typeorf, empty_dict
                        )
                        
                        # Assign scores to ORFs
                        orfs_filtered = assigningscore(
                            orfs_filtered, scoredict, typeorf
                        )
                        
                        # Calculate global scores
                        orfs_filtered = globalscores(orfs_filtered, tran_reads, typeorf)
                        results.append(orfs_filtered)
                scoring_duration = time.time() - scoring_start_time
                log_info(f"Time taken to score ORF type {typeorf} for {tran_id}: {scoring_duration:.4f} seconds")
            except Exception as e:
                log_error(f"Error scoring ORFs for transcript {tran_id}, type {typeorf}: {str(e)}")
        
        # Combine results for this transcript
        if results:
            total_duration = time.time() - start_time
            log_info(f"Total time taken for transcript {tran_id}: {total_duration:.4f} seconds")
            return pl.concat(results).to_dict(as_series=False)
        
        return None
    except Exception as e:
        log_error(f"Error in score_transcript for {tran_id}: {str(e)}")
        return None


def stream_results_to_disk(results_batch, output_file):
    """
    Stream results to a disk file to save memory.
    
    Args:
        results_batch: List of DataFrames to save
        output_file: Path to save the results
    """
    if not results_batch:
        return
        
    # Check if file already exists
    file_exists = os.path.exists(output_file)
    
    # Concatenate the batch
    combined = pl.concat(results_batch)
    
    if file_exists:
        # Append to existing file
        try:
            combined.write_csv(output_file, mode="a", include_header=False)
        except Exception as e:
            log_error(f"Error appending to CSV: {e}")
    else:
        # Create new file
        try:
            combined.write_csv(output_file)
        except Exception as e:
            log_error(f"Error writing to CSV: {e}")
    
    # Clear the batch
    for i in range(len(results_batch)):
        results_batch[i] = None

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
    import pyBigWig as bw
    import polars as pl
    from ..utils.logging import log_info, log_warning, log_error
    import numpy as np
    
    # Open BigWig file if a path was given
    if isinstance(bigwig_file, str):
        bw_handle = bw.open(bigwig_file)
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
            log_info("Detected genomic BigWig, mapping to transcript coordinates")
            return process_genomic_bigwig(bw_handle, exon_df)
        else:
            log_info("Detected transcriptomic BigWig, using directly")
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
    import numpy as np
    from ..utils.logging import log_info, log_warning
    
    # Get chromosomes in BigWig
    bw_chroms = set(bw_handle.chroms().keys())
    
    # Create a dictionary to store transcript coverage
    transcript_coverage = {}
    
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
            if chrom.startswith('chr') and chrom[3:] in bw_chroms:
                chrom = chrom[3:]
            # Try with 'chr' prefix
            elif f"chr{chrom}" in bw_chroms:
                chrom = f"chr{chrom}"
            else:
                continue
        
        # Initialize coverage array for this transcript
        coverage_dict = {}
        max_tran_pos = 0
        
        # Process each exon in this transcript
        for row in tran_exons.iter_rows(named=True):
            # Get genomic coordinates
            if isinstance(row["start"], list):
                # Multiple exons case
                for i in range(len(row["start"])):
                    start = int(row["start"][i])
                    stop = int(row["stop"][i])
                    tran_start = int(row["tran_start"][i])
                    
                    # Calculate transcript length for this exon
                    exon_length = stop - start
                    max_tran_pos = max(max_tran_pos, tran_start + exon_length)
                    
                    try:
                        # Get coverage values from BigWig
                        values = bw_handle.values(chrom, start, stop)
                        
                        # Map to transcript coordinates
                        for j, value in enumerate(values):
                            tran_pos = tran_start + j
                            if value is not None:
                                coverage_dict[tran_pos] = value
                    except Exception as e:
                        log_warning(f"Error reading {chrom}:{start}-{stop}: {str(e)}")
            else:
                # Single exon case
                start = int(row["start"])
                stop = int(row["stop"])
                tran_start = int(row["tran_start"])
                
                # Calculate transcript length for this exon
                exon_length = stop - start
                max_tran_pos = max(max_tran_pos, tran_start + exon_length)
                
                try:
                    # Get coverage values from BigWig
                    values = bw_handle.values(chrom, start, stop)
                    
                    # Map to transcript coordinates
                    for j, value in enumerate(values):
                        tran_pos = tran_start + j
                        if value is not None:
                            coverage_dict[tran_pos] = value
                except Exception as e:
                    log_warning(f"Error reading {chrom}:{start}-{stop}: {str(e)}")
        
        # Store coverage for this transcript
        if coverage_dict:
            transcript_coverage[tran_id] = (coverage_dict, max_tran_pos)
    
    # Combine all transcripts into a single DataFrame
    results = []
    
    for tran_id, (coverage_dict, max_pos) in transcript_coverage.items():
        # Create array for this transcript
        coverage_array = np.zeros(max_pos + 1)
        
        # Fill in coverage values
        for pos, value in coverage_dict.items():
            if pos < len(coverage_array):
                coverage_array[pos] = value
        
        # Create DataFrame for this transcript
        tran_df = pl.DataFrame({
            "tran_id": [tran_id] * len(coverage_array),
            "tran_start": pl.arange(0, len(coverage_array)),
            "counts": coverage_array
        })
        
        results.append(tran_df)
    
    # Combine all transcripts
    if results:
        return pl.concat(results)
    else:
        # Return empty DataFrame with correct structure
        return pl.DataFrame({
            "tran_id": [],
            "tran_start": [],
            "counts": []
        })


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
    import numpy as np
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
            values = bw_handle.values(tran_id, 0, tran_size)
            
            # Create DataFrame for this transcript
            tran_df = pl.DataFrame({
                "tran_id": [tran_id] * len(values),
                "tran_start": pl.arange(0, len(values)),
                "counts": [v if v is not None else 0.0 for v in values]
            })
            
            results.append(tran_df)
        except Exception as e:
            log_warning(f"Error reading transcript {tran_id}: {str(e)}")
    
    # Combine all transcripts
    if results:
        return pl.concat(results)
    else:
        # Return empty DataFrame with correct structure
        return pl.DataFrame({
            "tran_id": [],
            "tran_start": [],
            "counts": []
        })
def scoring(bigwig, exon, orfs, old_scoring, sru_range, stranded=False, batch_size=50, max_workers=None):
    """
    Score ORFs using bigwig coverage data with memory-efficient streaming.
    Now supports genomic BigWig files and strand-specific analysis.
    
    Args:
        bigwig (str or dict): Path to bigwig file or dict with 'forward' and 'reverse' keys
        exon (str or DataFrame): Path to exon file or DataFrame
        orfs (str or DataFrame): Path to ORFs file or DataFrame
        old_scoring (bool): Whether to use old scoring method
        sru_range (int): Range for SRU score calculation
        stranded (bool): Whether to process strands separately
        batch_size (int): Size of transcript batches for processing
        max_workers (int, optional): Maximum number of worker processes
        
    Returns:
        DataFrame: Scored ORFs
    """
    import os
    import polars as pl
    import pyBigWig as bw
    import time
    import concurrent.futures
    import gc
    from ..utils.logging import log_info, log_warning, log_error
    
    start_time = time.time()
    
    # Handle different bigwig input formats
    if isinstance(bigwig, dict) and 'forward' in bigwig and 'reverse' in bigwig:
        log_info("Using strand-specific BigWig files")
        forward_bw = bigwig['forward']
        reverse_bw = bigwig['reverse']
        stranded = True
    else:
        log_info(f"Using single BigWig file: {bigwig}")
        forward_bw = reverse_bw = bigwig
    
    # Load exon data efficiently
    log_info("Loading exon data...")
    if isinstance(exon, str):
        if os.path.exists(exon) and os.path.getsize(exon) > 1e9:  # 1 GB
            # For large files, stream and select only needed columns
            exon_df = pl.scan_csv(exon).select(["chr", "start", "stop", "tran_id", "tran_start", "tran_stop"]).collect()
        else:
            exon_df = pl.read_csv(exon)
            
            # Ensure exon coordinates are properly formatted (convert strings to lists if needed)
            for col in ["start", "stop", "tran_start", "tran_stop"]:
                if col in exon_df.columns:
                    exon_df = exon_df.with_columns(
                        pl.col(col).map_elements(lambda x: 
                            x.split(",") if isinstance(x, str) else x
                        ).alias(col)
                    )
    else:
        exon_df = exon
    
    # Load ORF data efficiently
    log_info("Loading ORF data...")
    if isinstance(orfs, str):
        if os.path.exists(orfs) and os.path.getsize(orfs) > 1e9:  # 1 GB
            # For large files, stream and select only needed columns
            needed_cols = ["tran_id", "type", "start", "stop", "length"]
            
            # Add strand if available
            orf_scan = pl.scan_csv(orfs)
            if "strand" in orf_scan.columns:
                needed_cols.append("strand")
                
            # Add any scoring-related columns
            for col in orf_scan.columns:
                if "score" in col.lower() or col in ["rise_up", "step_down", "hrf", "avg", "nzc"]:
                    needed_cols.append(col)
                    
            orf_df = orf_scan.select(needed_cols).collect()
        else:
            orf_df = pl.read_csv(orfs)
    else:
        orf_df = orfs
    
    # Add strand info if missing
    if "strand" not in orf_df.columns:
        log_info("Adding strand information to ORFs based on transcript annotation")
        orf_df = add_strand_to_orfs(orf_df, exon_df)
    
    # Split ORFs by strand if using strand-specific data
    if stranded:
        log_info("Processing ORFs by strand")
        
        # Group ORFs by strand
        pos_strand_orfs = orf_df.filter(pl.col("strand") == "+")
        neg_strand_orfs = orf_df.filter(pl.col("strand") == "-")
        unstrand_orfs = orf_df.filter(~pl.col("strand").is_in(["+", "-"]))
        
        log_info(f"Found {len(pos_strand_orfs)} positive strand, {len(neg_strand_orfs)} negative strand, and {len(unstrand_orfs)} unstranded ORFs")
        
        # Process each strand separately
        results = []
        
        # Process positive strand ORFs with forward BigWig
        if not pos_strand_orfs.is_empty():
            log_info("Scoring positive strand ORFs")
            pos_results = process_strand_orfs(forward_bw, exon_df, pos_strand_orfs, old_scoring, sru_range, batch_size, max_workers)
            results.append(pos_results)
        
        # Process negative strand ORFs with reverse BigWig
        if not neg_strand_orfs.is_empty():
            log_info("Scoring negative strand ORFs")
            neg_results = process_strand_orfs(reverse_bw, exon_df, neg_strand_orfs, old_scoring, sru_range, batch_size, max_workers)
            results.append(neg_results)
        
        # Process unstranded ORFs with forward BigWig (default)
        if not unstrand_orfs.is_empty():
            log_info("Scoring unstranded ORFs using forward strand data")
            unstrand_results = process_strand_orfs(forward_bw, exon_df, unstrand_orfs, old_scoring, sru_range, batch_size, max_workers)
            results.append(unstrand_results)
        
        # Combine results
        if results:
            final_df = pl.concat(results)
            log_info(f"Total scored ORFs: {len(final_df)}")
            return final_df
        else:
            return pl.DataFrame()
    else:
        # Process all ORFs with the same BigWig (original behavior)
        return process_strand_orfs(forward_bw, exon_df, orf_df, old_scoring, sru_range, batch_size, max_workers)


def process_strand_orfs(bigwig_path, exon_df, orf_df, old_scoring, sru_range, batch_size=50, max_workers=None):
    """
    Process and score ORFs for a specific strand.
    
    Args:
        bigwig_path (str): Path to bigwig file
        exon_df (DataFrame): DataFrame containing exon information
        orf_df (DataFrame): DataFrame containing ORF information
        old_scoring (bool): Whether to use old scoring method
        sru_range (int): Range for SRU score calculation
        batch_size (int): Size of transcript batches for processing
        max_workers (int, optional): Maximum number of worker processes
        
    Returns:
        DataFrame: Scored ORFs
    """
    import polars as pl
    import pyBigWig as bw
    import os
    import gc
    import time
    import tempfile
    import concurrent.futures
    from ..utils.logging import log_info
    
    start_time = time.time()
    
    # Determine optimal worker count based on system resources
    if max_workers is None:
        try:
            import psutil
            mem_gb = psutil.virtual_memory().available / (1024**3)
            # More conservative worker allocation
            max_workers = max(1, min(os.cpu_count(), int(mem_gb / 2)))
            max_workers = min(max_workers, 6)  # Cap at 6 workers for stability
        except ImportError:
            # Conservative default if psutil not available
            max_workers = min(os.cpu_count() or 4, 4)
    
    log_info(f"Using {max_workers} worker processes")
    
    # Create a temporary file for streaming results
    temp_dir = tempfile.gettempdir()
    temp_results_file = os.path.join(temp_dir, f"translonscorer_results_{int(time.time())}.csv")
    log_info(f"Will stream results to temporary file: {temp_results_file}")
    
    # Get unique transcripts with ORFs
    transcript_ids = set(orf_df["tran_id"].unique())
    log_info(f"Processing {len(transcript_ids)} transcripts with ORFs")
    
    # Process in batches for memory efficiency
    transcript_batches = [list(transcript_ids)[i:i+batch_size] for i in range(0, len(transcript_ids), batch_size)]
    
    processed_count = 0
    for batch_idx, transcript_batch in enumerate(transcript_batches):
        log_info(f"Processing batch {batch_idx+1}/{len(transcript_batches)} ({len(transcript_batch)} transcripts)")
        
        # Process transcripts in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for tran_id in transcript_batch:
                # Get ORFs for this transcript
                tran_orfs = orf_df.filter(pl.col("tran_id") == tran_id)
                
                if tran_orfs.is_empty():
                    continue
                
                # Get exons for this transcript
                tran_exons = exon_df.filter(pl.col("tran_id") == tran_id)
                
                if tran_exons.is_empty():
                    continue
                
                # Submit task
                future = executor.submit(
                    score_single_transcript, 
                    bigwig_path, 
                    tran_exons, 
                    tran_orfs, 
                    old_scoring, 
                    sru_range
                )
                futures.append(future)
            
            # Process results
            results = []
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                        processed_count += 1
                except Exception as e:
                    log_info(f"Error processing transcript: {str(e)}")
            
            # Combine and save batch results
            if results:
                batch_df = pl.concat(results)
                
                # Write to temp file
                if batch_idx == 0:
                    batch_df.write_csv(temp_results_file)
                else:
                    batch_df.write_csv(temp_results_file, mode="a", include_header=False)
                
                log_info(f"Processed {processed_count}/{len(transcript_ids)} transcripts")
        
        # Clear memory between batches
        gc.collect()
    
    # Load and return final results
    if os.path.exists(temp_results_file) and os.path.getsize(temp_results_file) > 0:
        log_info(f"Loading final results from {temp_results_file}")
        
        try:
            # For large files, use lazy loading
            if os.path.getsize(temp_results_file) > 1e9:  # 1 GB
                final_df = pl.scan_csv(temp_results_file).collect()
            else:
                final_df = pl.read_csv(temp_results_file)
                
            log_info(f"Successfully loaded {final_df.height} scored ORFs")
            
            # Clean up temp file
            try:
                os.remove(temp_results_file)
                log_info(f"Temporary file removed: {temp_results_file}")
            except Exception as e:
                log_info(f"Could not remove temporary file: {temp_results_file}")
            
            # Calculate and log performance metrics
            duration = time.time() - start_time
            orfs_per_second = final_df.height / max(0.001, duration)
            
            log_info(f"Scoring completed in {duration:.2f} seconds")
            log_info(f"Performance: {orfs_per_second:.2f} ORFs/second")
            
            return final_df
        except Exception as e:
            log_info(f"Error loading results: {str(e)}")
            return pl.DataFrame()
    else:
        log_info("No results were generated or temporary file is empty")
        return pl.DataFrame()


def score_single_transcript(bigwig_path, exon_df, orf_df, old_scoring, sru_range):
    """
    Score ORFs for a single transcript.
    
    Args:
        bigwig_path (str): Path to bigwig file
        exon_df (DataFrame): DataFrame containing exon information for one transcript
        orf_df (DataFrame): DataFrame containing ORF information for one transcript
        old_scoring (bool): Whether to use old scoring method
        sru_range (int): Range for SRU score calculation
        
    Returns:
        DataFrame: Scored ORFs
    """
    import polars as pl
    import pyBigWig as bw
    from ..utils.logging import log_info
    
    # Get transcript ID (should be the same for all ORFs)
    tran_id = orf_df["tran_id"][0]
    
    # Open BigWig file
    with bw.open(bigwig_path) as bw_file:
        # Get transcript coverage
        coverage_df = transcriptreads(bw_file, exon_df, tran_id)
        
        if coverage_df.is_empty():
            return None
        
        # Process different ORF types separately
        results = []
        
        for orf_type in orf_df["type"].unique():
            type_orfs = orf_df.filter(pl.col("type") == orf_type)
            
            # Score ORFs based on method
            if old_scoring:
                from ..core.scoring import oldscoring
                scored_orfs = oldscoring(type_orfs, coverage_df, sru_range, orf_type)
            else:
                from ..core.scoring import globalscores, existingscore, assigningscore, newscoring
                
                # Create empty score dict
                score_dict = {"rise_up": {}, "step_down": {}}
                
                # Find ORFs that need scoring
                empty_score_df = existingscore(type_orfs, orf_type, score_dict)
                
                if not empty_score_df.is_empty():
                    # Calculate new scores
                    score_dict = newscoring(empty_score_df, coverage_df, sru_range, orf_type, score_dict)
                    
                    # Assign scores to ORFs
                    type_orfs = assigningscore(type_orfs, score_dict, orf_type)
                    
                    # Calculate global scores
                    scored_orfs = globalscores(type_orfs, coverage_df, orf_type)
                else:
                    scored_orfs = type_orfs
            
            results.append(scored_orfs)
        
        # Combine results
        if results:
            return pl.concat(results)
        else:
            return None


def add_strand_to_orfs(orf_df, exon_df):
    """
    Add strand information to ORFs based on transcript annotation.
    
    Args:
        orf_df (DataFrame): DataFrame containing ORF information
        exon_df (DataFrame): DataFrame containing exon information
        
    Returns:
        DataFrame: ORF DataFrame with added strand column
    """
    import polars as pl
    
    # Create a mapping of transcript ID to strand
    tran_strands = {}
    for row in exon_df.iter_rows(named=True):
        tran_id = row["tran_id"]
        
        # Determine strand
        if "strand" in row:
            strand = row["strand"]
            # Handle list case
            if isinstance(strand, list):
                strand = strand[0]
        else:
            # Default to forward strand if not specified
            strand = "+"
            
        tran_strands[tran_id] = strand
    
    # Add strand column to ORF DataFrame
    orf_df = orf_df.with_columns(
        pl.col("tran_id").map_elements(lambda tran_id: 
            tran_strands.get(tran_id, "+")  # Default to + if transcript not found
        ).alias("strand")
    )
    
    return orf_df


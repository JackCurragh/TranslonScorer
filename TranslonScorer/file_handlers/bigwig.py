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
import itertools
import gc
from memory_profiler import profile  # Import the memory profiler


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
            return self.cache[cache_key]
        
        # Fetch from file
        with open_bigwig(self.bigwig_path) as bwfile:
            values = bwfile.values(chrom, start, stop)
            
            # Manage cache size
            if len(self.cache) >= self.max_cache_size:
                # Remove a random item (simple strategy)
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


def extract_regions(exon_df: pl.DataFrame) -> List[Tuple[str, int, int, str]]:
    """
    Extract regions (chr, start, stop, tran_id) from exon DataFrame.
    Handles both list and scalar values.
    
    Returns:
        List of (chr, start, stop, tran_id) tuples
    """
    regions = []
    
    # Process each row in the exon dataframe
    for row in exon_df.iter_rows(named=True):
        chrom = row["chr"]
        tran_id = row["tran_id"]
        
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


def process_df_in_chunks(df, chunk_size=10000, func=None):
    """
    Process a large DataFrame in chunks to reduce memory usage.
    
    Args:
        df: Polars DataFrame
        chunk_size: Number of rows to process at once
        func: Function to apply to each chunk
        
    Returns:
        List of results from processing each chunk
    """
    results = []
    total_rows = df.height
    
    for i in range(0, total_rows, chunk_size):
        end = min(i + chunk_size, total_rows)
        chunk = df.slice(i, end - i)
        
        if func:
            result = func(chunk)
            results.append(result)
        else:
            results.append(chunk)
            
        # Force garbage collection after processing each chunk
        chunk = None
        gc.collect()
        
    return results


def process_batch_memory_efficient(batch_data):
    """
    Process a batch of regions using memory-efficient approach.
    
    Args:
        batch_data: Tuple containing (indices, all_regions, bwfile_path, chunk_size)
        
    Returns:
        Dictionary with processing results
    """
    batch_indices, all_regions, bwfile_path, chunk_size = batch_data
    start_idx, end_idx = batch_indices
    batch_regions = all_regions[start_idx:end_idx]
    
    # Use BigWigCache for faster access
    bw_cache = BigWigCache(bwfile_path, max_cache_size=2000)  # Larger cache
    
    # Group regions by transcript to minimize BigWig access
    regions_by_transcript = {}
    for chrom, start, stop, tran_id in batch_regions:
        if tran_id not in regions_by_transcript:
            regions_by_transcript[tran_id] = []
        regions_by_transcript[tran_id].append((chrom, start, stop))
    
    batch_results = {}
    processed = 0
    failed = 0
    
    # Process each transcript's regions
    for tran_id, regions in regions_by_transcript.items():
        # For small transcripts, use a regular list
        # For larger ones, use a NumPy array directly
        all_reads = []
        total_length = sum(stop - start for _, start, stop in regions)
        
        # If the total length is large, pre-allocate a numpy array
        if total_length > 50000:
            all_reads = np.zeros(total_length, dtype=np.float32)
            current_pos = 0
            
            for chrom, start, stop in regions:
                try:
                    region_size = stop - start
                    values = bw_cache.get_values(chrom, start, stop)
                    
                    if values and any(v is not None for v in values):
                        for i, v in enumerate(values):
                            if v is not None:
                                all_reads[current_pos + i] = v
                    
                    current_pos += region_size
                    processed += 1
                except Exception:
                    failed += 1
                    # Skip filling, already zeros
                    current_pos += stop - start
            
            # Only store if there are non-zero values
            if np.any(all_reads):
                batch_results[tran_id] = all_reads
        else:
            # For small regions, use the original approach (faster for small data)
            for chrom, start, stop in regions:
                try:
                    values = bw_cache.get_values(chrom, start, stop)
                    if values and any(v is not None for v in values):
                        all_reads.extend(v if v is not None else 0.0 for v in values)
                    else:
                        all_reads.extend([0.0] * (stop - start))
                    
                    processed += 1
                except Exception:
                    failed += 1
                    all_reads.extend([0.0] * (stop - start))
            
            # Only keep non-empty results
            if all_reads and any(v > 0 for v in all_reads):
                batch_results[tran_id] = np.array(all_reads, dtype=np.float32)
    
    return {"results": batch_results, "processed": processed, "failed": failed}


def score_transcripts_memory_efficient(batch_data):
    """
    Score a batch of transcripts with memory-efficient approach.
    
    Args:
        batch_data: Tuple containing necessary data
        
    Returns:
        List of scored ORF DataFrames
    """
    batch_indices, available_transcripts, transcript_data, orf_df, old_scoring, sru_range = batch_data
    start_idx, end_idx = batch_indices
    batch_transcripts = available_transcripts[start_idx:end_idx]
    
    batch_results = []
    
    for tran in batch_transcripts:
        if tran not in transcript_data:
            continue
            
        # Convert NumPy array to DataFrame
        tran_array = transcript_data[tran]
        tran_reads = pl.DataFrame({
            "tran_start": pl.arange(0, len(tran_array), eager=True),
            "counts": tran_array
        })
        
        orfs = orf_df.filter(pl.col("tran_id") == tran)
        
        if orfs.is_empty() or tran_reads.is_empty():
            continue
            
        for typeorf in orfs["type"].unique():
            orfs_filtered = orfs.filter(pl.col("type") == typeorf)
            
            if orfs_filtered.is_empty():
                continue
                
            try:
                if old_scoring:
                    orfs_filtered = oldscoring(
                        orfs_filtered, tran_reads, sru_range, typeorf
                    )
                    batch_results.append(orfs_filtered)
                else:
                    emptyscore_df = existingscore(orfs_filtered, typeorf, {"rise_up": {}, "step_down": {}})
                    if not emptyscore_df.is_empty():
                        scoredict = newscoring(
                            emptyscore_df, tran_reads, sru_range, typeorf, {"rise_up": {}, "step_down": {}}
                        )
                        orfs_filtered = assigningscore(
                            orfs_filtered, scoredict, typeorf
                        )
                        orfs_filtered = globalscores(orfs_filtered, tran_reads, typeorf)
                        batch_results.append(orfs_filtered)
            except Exception as e:
                log_error(f"Error scoring ORFs for transcript {tran}, type {typeorf}: {str(e)}")
    
    return batch_results


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
        # Read the existing header first
        try:
            with open(output_file, 'r') as f:
                header = f.readline().strip()
            
            # Write without header
            with open(output_file, 'a') as f:
                combined.write_csv(f, include_header=False)
        except Exception as e:
            log_error(f"Error appending to CSV: {e}")
    else:
        # Create new file
        try:
            with open(output_file, 'w') as f:
                combined.write_csv(f, include_header=True)
        except Exception as e:
            log_error(f"Error writing to CSV: {e}")
    
    # Clear the batch
    for i in range(len(results_batch)):
        results_batch[i] = None

def transcriptreads(bwfile: bw.pyBigWig, exon_df: pl.DataFrame) -> pl.DataFrame:
    """
    Converts a BigWig file to a DataFrame based on provided exon annotation.
    
    Parameters:
    ----------
    bwfile : pyBigWig.pyBigWig
        An open BigWig file handle
    exon_df : polars.DataFrame
        DataFrame containing exon annotations with columns: chr, start, stop
        
    Returns:
    -------
    polars.DataFrame
        DataFrame containing transcript information with columns: tran_start, counts
    """
    config = ProcessingConfig()
    
    if not isinstance(bwfile, bw.pyBigWig):
        raise TypeError("bwfile must be a pyBigWig handle")
    
    if exon_df.is_empty():
        raise ValueError("Empty exon DataFrame provided")
        
    if not all(col in exon_df.columns for col in ["chr", "start", "stop"]):
        raise ValueError("Exon DataFrame missing required columns (chr, start, stop)")
    
    reads = []
    processed_regions = 0
    failed_regions = 0
    
    # Get chromosomes that exist in the bigwig file
    valid_chroms = set(bwfile.chroms().keys())
    
    start_time = time.time()
    
    try:
        # Process each row in the exon dataframe
        for row in exon_df.iter_rows(named=True):
            chrom = row["chr"]
            
            if chrom not in valid_chroms:
                log_warning(f"Chromosome {chrom} not found in bigWig file")
                continue
            
            # Handle both single values and lists for start/stop
            starts = _flatten_list_or_value(row["start"])
            stops = _flatten_list_or_value(row["stop"])
            
            # Make sure we have matching start/stop pairs
            if len(starts) != len(stops):
                log_warning(f"Mismatched start/stop lists: {len(starts)} starts, {len(stops)} stops")
                failed_regions += 1
                continue
            
            # Process each start/stop pair
            for start, stop in zip(starts, stops):
                # Ensure start and stop are integers
                try:
                    start = int(start)
                    stop = int(stop)
                except (ValueError, TypeError) as e:
                    log_warning(f"Invalid start/stop values: {start}, {stop} - {str(e)}")
                    failed_regions += 1
                    continue
                
                if start >= stop:
                    log_warning(f"Invalid region {chrom}:{start}-{stop} (start >= stop)")
                    failed_regions += 1
                    continue
                
                # Check if region is too large
                region_size = stop - start
                if region_size > config.max_region_size:
                    log_warning(f"Very large region detected: {chrom}:{start}-{stop} ({region_size} bp). Chunking.")
                    
                    # Process in chunks to avoid memory issues
                    for chunk_start in range(start, stop, config.chunk_size):
                        chunk_end = min(chunk_start + config.chunk_size, stop)
                        
                        try:
                            values = bwfile.values(chrom, chunk_start, chunk_end)
                            if values and any(v is not None for v in values):
                                reads.extend(v if v is not None else 0.0 for v in values)
                            else:
                                reads.extend([0.0] * (chunk_end - chunk_start))
                        except Exception as e:
                            log_error(f"Error processing chunk {chrom}:{chunk_start}-{chunk_end}: {str(e)}")
                            reads.extend([0.0] * (chunk_end - chunk_start))
                    
                    processed_regions += 1
                else:
                    # Process the whole region at once
                    try:
                        values = bwfile.values(chrom, start, stop)
                        if values and any(v is not None for v in values):
                            reads.extend(v if v is not None else 0.0 for v in values)
                            processed_regions += 1
                        else:
                            reads.extend([0.0] * (stop - start))
                            processed_regions += 1
                    except Exception as e:
                        log_error(f"Could not read values for {chrom}:{start}-{stop}: {str(e)}")
                        failed_regions += 1
                    
    except Exception as e:
        log_error(f"Error processing exon data: {str(e)}")
        raise
    
    if not reads:
        msg = f"No reads extracted. Processed: {processed_regions}, Failed: {failed_regions}"
        log_error(msg, exception_type=ValueError)
        raise ValueError(msg)
        
    duration = time.time() - start_time
    rate = processed_regions / max(0.001, duration)
    log_info(f"Successfully processed {processed_regions} regions ({rate:.2f} regions/sec), {failed_regions} failed")
    
    # Create DataFrame safely
    return pl.DataFrame({
        "tran_start": list(range(len(reads))),
        "counts": reads
    })


@profile
def scoring(bigwig, exon, orfs, old_scoring, sru_range, batch_size=50, max_workers=None):
    """
    Score ORFs using bigwig coverage data with memory-efficient streaming.
    
    Args:
        bigwig (str): Path to bigwig file
        exon (str or DataFrame): Path to exon file or DataFrame
        orfs (str or DataFrame): Path to ORFs file or DataFrame
        old_scoring (bool): Whether to use old scoring method
        sru_range (int): Range for SRU score calculation
        batch_size (int): Size of transcript batches for processing
        max_workers (int, optional): Maximum number of worker processes
        
    Returns:
        DataFrame: Scored ORFs
    """
    import gc
    import os
    import time
    import tempfile
    import concurrent.futures
    import numpy as np
    from functools import partial
    import polars as pl
    from ..utils.logging import log_info, log_warning, log_error
    
    start_time = time.time()
    config = ProcessingConfig(batch_size=batch_size, sru_range=sru_range)
    
    bwfile_path = bigwig
    
    # Load data efficiently
    log_info("Loading exon and ORF data...")
    
    # Load exon data with optimized approach
    if isinstance(exon, str):
        # Use streaming approach for large files
        if os.path.exists(exon) and os.path.getsize(exon) > 1e9:  # 1 GB
            # Stream read with specific columns - much more memory efficient
            exon_df = pl.scan_csv(exon).select(["chr", "start", "stop", "tran_id"]).collect()
        else:
            exon_df = pl.read_csv(exon)
    else:
        exon_df = exon
    
    # For ORFs - similar optimization
    if isinstance(orfs, str):
        if os.path.exists(orfs) and os.path.getsize(orfs) > 1e9:  # 1 GB
            # Create a LazyFrame for large ORF files
            orf_scan = pl.scan_csv(orfs)
            # Only collect essential columns and materialize only once
            needed_cols = ["tran_id", "type", "start", "stop", "length"]
            # Add any scoring-related columns
            all_cols = orf_scan.columns
            for col in all_cols:
                if "score" in col.lower() or col in ["rise_up", "step_down"]:
                    needed_cols.append(col)
            orf_df = orf_scan.select(needed_cols).collect()
        else:
            orf_df = pl.read_csv(orfs)
    else:
        orf_df = orfs

    # Pre-compute and cache ORF transcript IDs to reduce repeated computations
    orf_tran_ids = set(orf_df["tran_id"].unique().to_list())
    
    # Only process exons that have corresponding ORFs to save time
    exon_df = exon_df.filter(pl.col("tran_id").is_in(orf_tran_ids))
    
    # Optimize worker count based on system resources
    if max_workers is None:
        try:
            import psutil
            # Use available memory to determine worker count
            available_memory_gb = psutil.virtual_memory().available / (1024**3)
            # Estimate 1GB per worker as baseline with safety margin
            suggested_workers = max(1, min(os.cpu_count(), int(available_memory_gb * 0.75)))
            max_workers = min(suggested_workers, 8)  # Cap at 8 workers for stability
        except ImportError:
            # If psutil is not available, use a conservative approach
            max_workers = min(os.cpu_count() or 4, 4)  # Lower default cap
    
    log_info(f"Using {max_workers} workers for parallel processing")
    
    # Extract regions more efficiently using Polars' native operations
    log_info("Extracting regions from exon data...")
    
    # Handle different data formats for start/stop
    # This is a critical optimization to avoid exploding memory with Python lists
    flat_regions = []
    
    # Process exon data to extract regions
    for row in exon_df.iter_rows(named=True):
        chrom = row["chr"]
        tran_id = row["tran_id"]
        
        # Skip if transcript has no ORFs (using pre-computed set)
        if tran_id not in orf_tran_ids:
            continue
        
        # Handle both single values and lists for start/stop
        if isinstance(row["start"], list):
            starts = row["start"]
            stops = row["stop"] if isinstance(row["stop"], list) else [row["stop"]] * len(starts)
        else:
            starts = [row["start"]]
            stops = [row["stop"]]
        
        # Ensure equal length
        min_len = min(len(starts), len(stops))
        if min_len == 0:
            continue
            
        starts = starts[:min_len]
        stops = stops[:min_len]
        
        # Process regions
        for start, stop in zip(starts, stops):
            try:
                start = int(start)
                stop = int(stop)
                
                if start >= stop:
                    continue
                    
                flat_regions.append((chrom, start, stop, tran_id))
            except (TypeError, ValueError):
                continue
    
    # Clear exon_df to free memory
    exon_df = None
    gc.collect()
    
    total_regions = len(flat_regions)
    log_info(f"Found {total_regions} regions to process from {len(orf_tran_ids)} transcripts with ORFs")
    
    if total_regions == 0:
        log_warning("No valid regions found to process")
        return pl.DataFrame()
    
    # Group regions by transcript to improve data locality and reduce BigWig access
    regions_by_transcript = {}
    for chrom, start, stop, tran_id in flat_regions:
        if tran_id not in regions_by_transcript:
            regions_by_transcript[tran_id] = []
        regions_by_transcript[tran_id].append((chrom, start, stop))
    
    # Clear flat_regions to free memory
    flat_regions = None
    gc.collect()
    
    # Create a more efficient batching strategy - group by chromosomes for better cache locality
    # This improves BigWig access patterns since the file is indexed by chromosome
    regions_by_chrom = {}
    for tran_id, regions in regions_by_transcript.items():
        for chrom, start, stop in regions:
            if chrom not in regions_by_chrom:
                regions_by_chrom[chrom] = []
            regions_by_chrom[chrom].append((start, stop, tran_id))
    
    # Create balanced batches that keep chromosome locality
    batches = []
    chrom_order = sorted(regions_by_chrom.keys())
    
    # Target a specific number of batches based on worker count
    target_batch_count = max_workers * 3  # Create 3x as many batches as workers for better load balancing
    target_batch_size = max(1, total_regions // target_batch_count)
    
    current_batch = []
    current_size = 0
    
    for chrom in chrom_order:
        chrom_regions = regions_by_chrom[chrom]
        # Sort by position to improve cache locality
        chrom_regions.sort()  # Sort by start position
        
        for region in chrom_regions:
            current_batch.append((chrom,) + region)
            current_size += 1
            
            if current_size >= target_batch_size:
                batches.append(current_batch)
                current_batch = []
                current_size = 0
    
    # Add the last batch if it has any regions
    if current_batch:
        batches.append(current_batch)
    
    # Clear region mappings to free memory
    regions_by_chrom = None
    regions_by_transcript = None
    gc.collect()
    
    log_info(f"Created {len(batches)} balanced processing batches")
    
    # Process batches with improved parallelism
    transcript_data = {}
    
    # Custom worker function with better memory management
    def process_batch(batch):
        # Use BigWigCache only for this batch
        bw_cache = BigWigCache(bwfile_path, max_cache_size=1000)
        
        # Group by transcript again
        local_transcript_regions = {}
        for chrom, start, stop, tran_id in batch:
            if tran_id not in local_transcript_regions:
                local_transcript_regions[tran_id] = []
            local_transcript_regions[tran_id].append((chrom, start, stop))
        
        batch_results = {}
        processed = 0
        failed = 0
        
        # Process each transcript
        for tran_id, regions in local_transcript_regions.items():
            # Skip transcripts that don't have ORFs
            if tran_id not in orf_tran_ids:
                continue
                
            try:
                # Calculate total region size to preallocate array
                total_size = sum(stop - start for _, start, stop in regions)
                
                # Use a NumPy array directly to save memory and improve performance
                position_data = []
                read_data = []
                
                # Process regions
                for chrom, start, stop in regions:
                    try:
                        region_size = stop - start
                        values = bw_cache.get_values(chrom, start, stop)
                        
                        if values is None:
                            # Region not in bigwig, add zeros
                            read_data.extend([0.0] * region_size)
                            position_data.extend(range(len(position_data), len(position_data) + region_size))
                            processed += 1
                            continue
                        
                        # Process valid values
                        current_position = len(position_data)
                        
                        for i, v in enumerate(values):
                            read_data.append(v if v is not None else 0.0)
                            position_data.append(current_position + i)
                        
                        processed += 1
                    except Exception as e:
                        failed += 1
                        # Handle error by filling with zeros
                        read_data.extend([0.0] * region_size)
                        position_data.extend(range(len(position_data), len(position_data) + region_size))
                
                # Only keep non-empty results and convert to NumPy arrays
                if read_data and any(v > 0 for v in read_data):
                    # Store as compressed tuple of position and reads - more memory efficient
                    # Only store non-zero values to save memory
                    non_zero_indices = [i for i, v in enumerate(read_data) if v > 0]
                    if non_zero_indices:
                        # Store sparse representation
                        batch_results[tran_id] = {
                            'positions': [position_data[i] for i in non_zero_indices],
                            'values': [read_data[i] for i in non_zero_indices],
                            'max_pos': max(position_data) if position_data else 0
                        }
            except Exception as e:
                log_error(f"Error processing transcript {tran_id}: {str(e)}")
        
        # Clear local variables
        bw_cache = None
        local_transcript_regions = None
        gc.collect()
        
        return {"results": batch_results, "processed": processed, "failed": failed}
    
    # Process batches with better parallelism
    # Use a context manager for better resource management
    total_processed = 0
    total_failed = 0
    
    # Create a temporary file for streaming results
    temp_results_file = os.path.join(tempfile.gettempdir(), f"translonscorer_results_{int(time.time())}.csv")
    log_info(f"Will stream results to temporary file: {temp_results_file}")
    
    # Process batches in staggered groups to limit memory pressure
    stagger_size = min(len(batches), max_workers * 2)
    
    for batch_start in range(0, len(batches), stagger_size):
        batch_end = min(batch_start + stagger_size, len(batches))
        current_batches = batches[batch_start:batch_end]
        
        log_info(f"Processing batch group {batch_start//stagger_size + 1}/{(len(batches) + stagger_size - 1)//stagger_size}")
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_batch, batch): i for i, batch in enumerate(current_batches)}
            
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    
                    # Update transcript data with results
                    for tran_id, data in result["results"].items():
                        transcript_data[tran_id] = data
                    
                    total_processed += result["processed"]
                    total_failed += result["failed"]
                    
                    # Log progress periodically
                    if len(transcript_data) % 100 == 0:
                        log_info(f"Processed data for {len(transcript_data)} transcripts so far")
                except Exception as e:
                    log_error(f"Error processing batch: {str(e)}")
        
        # Force garbage collection between batch groups
        gc.collect()
    
    log_info(f"Completed BigWig data extraction: {total_processed} regions processed, {total_failed} failed")
    log_info(f"Extracted data for {len(transcript_data)} transcripts")
    
    # Clear batches to free memory
    batches = None
    gc.collect()
    
    # Convert sparse transcript data to DataFrame format for scoring
    def prepare_transcript_reads(tran_id, data):
        # Convert sparse representation to DataFrame
        max_pos = data['max_pos']
        sparse_positions = data['positions']
        sparse_values = data['values']
        
        # Create DataFrame directly with sparse data
        # This avoids creating full arrays with mostly zeros
        tran_reads = pl.DataFrame({
            "tran_start": sparse_positions,
            "counts": sparse_values
        })
        
        return tran_id, tran_reads
    
    # Score ORFs using more efficient batching
    available_transcripts = list(transcript_data.keys())
    log_info(f"Scoring {len(available_transcripts)} transcripts with available data")
    
    # Create a more memory-efficient way to score ORFs
    # Process one transcript at a time to limit memory usage
    def score_transcript(tran_id):
        if tran_id not in transcript_data:
            return None
            
        # Create DataFrame from sparse data
        data = transcript_data[tran_id]
        _, tran_reads = prepare_transcript_reads(tran_id, data)
        
        # Filter ORFs for this transcript
        orfs_for_tran = orf_df.filter(pl.col("tran_id") == tran_id)
        
        if orfs_for_tran.is_empty() or tran_reads.is_empty():
            return None
            
        results = []
        
        # Score each ORF type
        for typeorf in orfs_for_tran["type"].unique():
            orfs_filtered = orfs_for_tran.filter(pl.col("type") == typeorf)
            
            if orfs_filtered.is_empty():
                continue
                
            try:
                if old_scoring:
                    orfs_filtered = oldscoring(
                        orfs_filtered, tran_reads, sru_range, typeorf
                    )
                    results.append(orfs_filtered)
                else:
                    # Optimize dictionary handling - use defaultdict to avoid key checks
                    from collections import defaultdict
                    empty_dict = {"rise_up": defaultdict(float), "step_down": defaultdict(float)}
                    
                    # Check existing scores
                    emptyscore_df = existingscore(orfs_filtered, typeorf, empty_dict)
                    if not emptyscore_df.is_empty():
                        # Calculate new scores
                        scoredict = newscoring(
                            emptyscore_df, tran_reads, sru_range, typeorf, empty_dict
                        )
                        # Assign scores
                        orfs_filtered = assigningscore(
                            orfs_filtered, scoredict, typeorf
                        )
                        # Calculate global scores
                        orfs_filtered = globalscores(orfs_filtered, tran_reads, typeorf)
                        results.append(orfs_filtered)
            except Exception as e:
                log_error(f"Error scoring ORFs for transcript {tran_id}, type {typeorf}: {str(e)}")
        
        # Combine results for this transcript
        if results:
            return pl.concat(results)
        return None
    
    # Process transcripts in batches to control memory usage
    scoring_batch_size = max(5, min(50, (len(available_transcripts) + max_workers - 1) // max_workers))
    log_info(f"Using batch size of {scoring_batch_size} for scoring")
    
    all_results = []
    processed_count = 0
    
    # Use a reduced worker count for scoring to manage memory
    scoring_workers = max(1, min(max_workers, 4))
    log_info(f"Using {scoring_workers} workers for scoring to manage memory usage")
    
    for batch_start in range(0, len(available_transcripts), scoring_batch_size):
        batch_end = min(batch_start + scoring_batch_size, len(available_transcripts))
        batch_transcripts = available_transcripts[batch_start:batch_end]
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=scoring_workers) as executor:
            batch_results = list(filter(None, executor.map(score_transcript, batch_transcripts)))
        
        # Stream results to disk instead of keeping in memory
        if batch_results:
            stream_results_to_disk(batch_results, temp_results_file)
        
        processed_count += len(batch_transcripts)
        completion_pct = (processed_count / len(available_transcripts)) * 100
        log_info(f"Processed {processed_count}/{len(available_transcripts)} transcripts ({completion_pct:.1f}%)")
        
        # Clear batch results to free memory
        batch_results = None
        gc.collect()
    
    # Load the final results
    if os.path.exists(temp_results_file) and os.path.getsize(temp_results_file) > 0:
        log_info("Reading final results from disk...")
        
        # Use streaming for large files
        if os.path.getsize(temp_results_file) > 1e9:  # > 1GB
            final_df = pl.scan_csv(temp_results_file).collect()
        else:
            final_df = pl.read_csv(temp_results_file)
            
        log_info(f"Loaded {final_df.height} scored ORFs successfully")
        
        # Clean up
        try:
            os.remove(temp_results_file)
        except:
            log_warning(f"Could not remove temporary file: {temp_results_file}")
        
        duration = time.time() - start_time
        log_info(f"Scoring completed in {duration:.2f} seconds")
        return final_df
    else:
        log_warning("No ORFs were scored or result file is empty")
        return pl.DataFrame()
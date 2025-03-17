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

import numpy as np
import tempfile
import mmap
import os

def create_memory_mapped_array(data, dtype=np.float32):
    """
    Create a memory-mapped numpy array to store large datasets.
    
    Args:
        data: Data to store in memory-mapped array
        dtype: Data type for the array
        
    Returns:
        tuple: (temp_file, np_array) - Keep the temp_file reference to prevent garbage collection
    """
    # Create a temporary file
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    temp_filename = temp_file.name
    temp_file.close()
    
    # Create a memory-mapped array
    shape = (len(data),)
    np_array = np.memmap(temp_filename, dtype=dtype, mode='w+', shape=shape)
    
    # Copy data to the memory-mapped array
    np_array[:] = data
    np_array.flush()
    
    return temp_filename, np_array

# 2. Chunked dataframe processing for Polars
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

# 3. Optimized BigWig reading with memory mapping
class MemoryEfficientBigWigReader:
    """A memory-efficient BigWig reader that uses memory mapping."""
    
    def __init__(self, bigwig_path):
        self.bigwig_path = bigwig_path
        self.temp_files = []  # Keep track of temp files
        
    def __del__(self):
        # Clean up temp files
        for file_path in self.temp_files:
            try:
                os.unlink(file_path)
            except:
                pass
                
    def get_values(self, chrom, start, stop, chunk_size=100000):
        """
        Get values from BigWig file for a region, using memory mapping for large regions.
        
        Args:
            chrom: Chromosome name
            start: Start position
            stop: End position
            chunk_size: Size of chunks to process
            
        Returns:
            numpy array with values
        """
        region_size = stop - start
        
        # For small regions, just read directly
        if region_size <= chunk_size:
            with open_bigwig(self.bigwig_path) as bwfile:
                values = bwfile.values(chrom, start, stop)
                if values is None:
                    return np.zeros(region_size, dtype=np.float32)
                return np.array(values, dtype=np.float32)
        
        # For large regions, use memory mapping
        result_filename, result_array = create_memory_mapped_array(
            np.zeros(region_size, dtype=np.float32)
        )
        self.temp_files.append(result_filename)
        
        # Process in chunks
        with open_bigwig(self.bigwig_path) as bwfile:
            for chunk_start in range(start, stop, chunk_size):
                chunk_end = min(chunk_start + chunk_size, stop)
                chunk_size_actual = chunk_end - chunk_start
                
                values = bwfile.values(chrom, chunk_start, chunk_end)
                if values is not None and any(v is not None for v in values):
                    # Copy valid values to the memory-mapped array
                    for i, v in enumerate(values):
                        if v is not None:
                            result_array[chunk_start - start + i] = v
                
                # Force flush after each chunk
                result_array.flush()
        
        return result_array

# 4. Improved transcript reads storage with numpy memory mapping
class MemoryEfficientTranscriptStorage:
    """
    Store transcript reads data in a memory-efficient way using memory mapping.
    """
    
    def __init__(self):
        self.transcript_data = {}  # Maps transcript ID to (filename, array) tuple
        self.temp_files = []
        
    def __del__(self):
        # Clean up temp files
        for file_path in self.temp_files:
            try:
                os.unlink(file_path)
            except:
                pass
                
    def add_transcript(self, tran_id, counts):
        """
        Add transcript data using memory mapping.
        
        Args:
            tran_id: Transcript ID
            counts: List of count values
        """
        if not counts:
            return
            
        # Create memory-mapped array for counts
        filename, array = create_memory_mapped_array(counts)
        self.temp_files.append(filename)
        self.transcript_data[tran_id] = (filename, array)
        
    def get_transcript_df(self, tran_id):
        """
        Get transcript data as a Polars DataFrame.
        
        Args:
            tran_id: Transcript ID
            
        Returns:
            Polars DataFrame with transcript data
        """
        if tran_id not in self.transcript_data:
            return pl.DataFrame({"tran_start": [], "counts": []})
            
        _, array = self.transcript_data[tran_id]
        
        # Convert to DataFrame efficiently
        return pl.DataFrame({
            "tran_start": pl.arange(0, len(array), eager=True),
            "counts": array.copy()  # Make a copy to avoid issues with the memmap
        })
        
    def get_transcript_ids(self):
        """Get all transcript IDs in the storage."""
        return list(self.transcript_data.keys())
    

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


def extract_regions(exon_df: pl.DataFrame) -> List[Tuple[str, int, int]]:
    """
    Extract regions (chr, start, stop) from exon DataFrame.
    Handles both list and scalar values.
    
    Returns:
        List of (chr, start, stop) tuples
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


def process_region_batch(region_batch, bigwig_path, config=None):
    """
    Process a batch of regions from the BigWig file.
    
    Args:
        region_batch: List of (chrom, start, stop, tran_id) tuples
        bigwig_path: Path to BigWig file
        config: ProcessingConfig object
        
    Returns:
        Dictionary mapping transcript IDs to reads DataFrames
    """
    if config is None:
        config = ProcessingConfig()
        
    # Create BigWig cache
    bw_cache = BigWigCache(bigwig_path)
    
    # Group regions by transcript
    regions_by_transcript = {}
    for chrom, start, stop, tran_id in region_batch:
        if tran_id not in regions_by_transcript:
            regions_by_transcript[tran_id] = []
        regions_by_transcript[tran_id].append((chrom, start, stop))
    
    results = {}
    processed_regions = 0
    failed_regions = 0
    
    # Process each transcript
    for tran_id, regions in regions_by_transcript.items():
        all_reads = []
        
        for chrom, start, stop in regions:
            # Check if region is too large
            region_size = stop - start
            if region_size > config.max_region_size:
                # Process in chunks
                for chunk_start in range(start, stop, config.chunk_size):
                    chunk_end = min(chunk_start + config.chunk_size, stop)
                    
                    try:
                        values = bw_cache.get_values(chrom, chunk_start, chunk_end)
                        if values and any(v is not None for v in values):
                            all_reads.extend(v if v is not None else 0.0 for v in values)
                        else:
                            all_reads.extend([0.0] * (chunk_end - chunk_start))
                            
                        processed_regions += 1
                    except Exception as e:
                        failed_regions += 1
                        all_reads.extend([0.0] * (chunk_end - chunk_start))
            else:
                # Process the whole region
                try:
                    values = bw_cache.get_values(chrom, start, stop)
                    if values and any(v is not None for v in values):
                        all_reads.extend(v if v is not None else 0.0 for v in values)
                    else:
                        all_reads.extend([0.0] * (stop - start))
                        
                    processed_regions += 1
                except Exception:
                    failed_regions += 1
                    all_reads.extend([0.0] * (stop - start))
        
        # Create DataFrame for this transcript
        if all_reads:
            results[tran_id] = pl.DataFrame({
                "tran_start": list(range(len(all_reads))),
                "counts": all_reads
            })
    
    return {
        "results": results, 
        "processed": processed_regions,
        "failed": failed_regions
    }

def score_transcript_batch(transcript_batch, transcript_reads, orf_df, old_scoring, sru_range):
    """
    Score a batch of transcripts.
    
    Args:
        transcript_batch: List of transcript IDs
        transcript_reads: Dictionary mapping transcript IDs to reads DataFrames
        orf_df: DataFrame containing all ORF data
        old_scoring: Whether to use old scoring method
        sru_range: Range for SRU score calculation
        
    Returns:
        List of scored ORF DataFrames
    """
    batch_results = []
    
    for tran in transcript_batch:
        if tran not in transcript_reads:
            continue
            
        tran_reads = transcript_reads[tran]
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


def process_region_batch_by_indices(indices_and_data):
    """
    Process a batch of regions from the BigWig file using indices.
    This must be defined at the module level for multiprocessing to work.
    
    Args:
        indices_and_data: Tuple containing (indices, all_regions, bigwig_path, config)
    
    Returns:
        Results from process_region_batch
    """
    indices, all_regions, bigwig_path, config = indices_and_data
    start_idx, end_idx = indices
    # Extract the actual data when needed inside the worker process
    batch_data = all_regions[start_idx:end_idx]
    # Call the original function with the extracted data
    return process_region_batch(batch_data, bigwig_path, config)

def score_transcript_batch_by_indices(indices_and_data):
    """
    Score a batch of transcripts using indices.
    This must be defined at the module level for multiprocessing to work.
    
    Args:
        indices_and_data: Tuple containing (indices, transcripts_list, transcript_reads, orf_df, old_scoring, sru_range)
    
    Returns:
        Results from score_transcript_batch
    """
    indices, transcripts_list, transcript_reads, orf_df, old_scoring, sru_range = indices_and_data
    start_idx, end_idx = indices
    batch = transcripts_list[start_idx:end_idx]
    return score_transcript_batch(batch, transcript_reads, orf_df, old_scoring, sru_range)

@profile
def scoring(bigwig, exon, orfs, old_scoring, sru_range, batch_size=50, max_workers=None):
    """
    Score ORFs using bigwig coverage data with high-performance optimization and memory efficiency.
    
    This implementation uses a two-stage approach:
    1. Extract all BigWig data in parallel by region batches using memory mapping
    2. Score all transcripts in parallel using the pre-extracted data
    
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
    import gc  # Add garbage collection
    
    start_time = time.time()
    config = ProcessingConfig(batch_size=batch_size, sru_range=sru_range)
    
    bwfile_path = bigwig  # Store the path instead of opening it here
    
    # Load data with memory-efficient approach
    log_info("Loading exon and ORF data...")
    
    # Check if exon is a file path or a DataFrame - use lazy loading for large files
    if isinstance(exon, str):
        if os.path.getsize(exon) > 1e9:  # 1 GB
            # Use scan_csv with lazy evaluation and filter/select only needed columns
            exon_scan = pl.scan_csv(exon, has_header=True, separator=",")
            # Only select columns we actually need to reduce memory usage
            needed_cols = ["chr", "start", "stop", "tran_id"]
            exon_df = exon_scan.select([col for col in needed_cols if col in exon_scan.columns]).collect()
        else:
            exon_df = pl.read_csv(exon, has_header=True, separator=",")
    else:
        exon_df = exon

    # For ORFs, use lazy loading with streaming for large files
    if isinstance(orfs, str):
        if os.path.getsize(orfs) > 1e9:  # 1 GB
            # Create a scanning DataFrame and filter out any unnecessary columns
            orf_scan = pl.scan_csv(orfs, has_header=True, separator=",")
            orf_df = orf_scan.collect()  
            # Immediately after collecting, clean up any unused columns
            needed_cols = ["tran_id", "type", "start", "stop", "length"]
            all_cols = set(orf_df.columns)
            for col in all_cols:
                if col not in needed_cols and col not in ["rise_up", "step_down"]:
                    orf_df = orf_df.drop(col)
        else:
            orf_df = pl.read_csv(orfs, has_header=True, separator=",")
    else:
        orf_df = orfs

    # Set max workers if not specified
    if max_workers is None:
        max_workers = min(os.cpu_count() or 4, 8)  # Limit to 8 workers to avoid memory pressure
    config.max_workers = max_workers
    
    log_info(f"Using {max_workers} workers for parallel processing")
    
    # STEP 1: Extract all regions from exon data
    log_info("Extracting regions from exon data...")
    all_regions = []
    
    # Process in smaller chunks to avoid memory spikes
    chunk_size = 10000
    total_rows = exon_df.height
    total_chunks = (total_rows + chunk_size - 1) // chunk_size  # Ceiling division
    unique_transcripts = set()
    
    for i in range(total_chunks):
        start_idx = i * chunk_size
        end_idx = min(start_idx + chunk_size, total_rows)
        
        # Extract chunk and process
        chunk = exon_df.slice(start_idx, end_idx - start_idx)
        chunk_regions = extract_regions(chunk)
        all_regions.extend(chunk_regions)
        
        # Track unique transcripts
        for region in chunk_regions:
            unique_transcripts.add(region[3])  # tran_id is the 4th element
        
        # Clear chunk to free memory
        chunk = None
        
        # Force garbage collection
        if i > 0 and i % 5 == 0:  # Every 5 chunks
            gc.collect()
            
    # Get count of regions
    total_regions = len(all_regions)
    total_unique_transcripts = len(unique_transcripts)
    log_info(f"Found {total_regions} regions to process across {total_unique_transcripts} transcripts")
    
    # Free up memory
    unique_transcripts = None
    gc.collect()
    
    # Group regions into batches for better work distribution - using indices only
    region_batches = []
    batch_size = min(config.max_regions_per_worker, 200)  # Limit batch size to avoid memory issues
    for i in range(0, total_regions, batch_size):
        end = min(i + batch_size, total_regions)
        # Store only indices, not data copies
        region_batches.append((i, end))
    
    log_info(f"Distributing regions into {len(region_batches)} balanced batches")
    
    # Create a memory-efficient transcript storage
    log_info("Creating memory-efficient transcript storage...")
    transcript_storage = MemoryEfficientTranscriptStorage()
    
    # Create a memory-efficient BigWig reader
    bw_reader = MemoryEfficientBigWigReader(bwfile_path)
    
    # STEP 2: Process all region batches with memory-efficient approach
    total_processed = 0
    total_failed = 0
    
    # Staggered parallelization constants
    max_concurrent_jobs = min(max_workers, 8)  # Limit concurrent jobs even more
    stagger_size = min(len(region_batches) // 10 + 1, max_concurrent_jobs)  # Process ~10% at a time
    log_interval = max(1, len(region_batches) // 20)  # Log ~20 times during processing
    
    # Process batches in parallel with memory-efficient worker function
    def process_batch_memory_efficient(batch_indices):
        start_idx, end_idx = batch_indices
        batch_regions = all_regions[start_idx:end_idx]
        
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
            all_reads = []
            
            for chrom, start, stop in regions:
                try:
                    # Use the memory-efficient reader
                    values = bw_reader.get_values(chrom, start, stop)
                    if values is not None and len(values) > 0:
                        all_reads.extend(values)
                    else:
                        all_reads.extend([0.0] * (stop - start))
                    
                    processed += 1
                except Exception as e:
                    failed += 1
                    all_reads.extend([0.0] * (stop - start))
            
            # Only keep non-empty results
            if all_reads:
                batch_results[tran_id] = all_reads
        
        return {"results": batch_results, "processed": processed, "failed": failed}
    
    log_info(f"Processing {len(region_batches)} region batches...")
    
    # Process batches in sequential staggered groups to limit memory usage
    for start_batch in range(0, len(region_batches), stagger_size):
        end_batch = min(start_batch + stagger_size, len(region_batches))
        batch_group = region_batches[start_batch:end_batch]
        
        log_info(f"Processing batch group {start_batch//stagger_size + 1}: "
                 f"batches {start_batch} to {end_batch-1}")
        
        # Process this group in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(process_batch_memory_efficient, indices) for indices in batch_group]
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    
                    # Add transcript data to storage
                    for tran_id, reads in result["results"].items():
                        transcript_storage.add_transcript(tran_id, reads)
                    
                    # Track statistics
                    total_processed += result["processed"]
                    total_failed += result["failed"]
                    
                except Exception as exc:
                    log_error(f"Batch processing error: {exc}")
        
        # Force garbage collection after each stagger group
        gc.collect()
    
    # Clear regions data as it's no longer needed
    all_regions = None
    region_batches = None
    gc.collect()
    
    log_info(f"Completed BigWig data extraction: {total_processed} regions processed, {total_failed} failed")
    
    # Get list of transcript IDs with data
    available_transcripts = transcript_storage.get_transcript_ids()
    log_info(f"Extracted data for {len(available_transcripts)} transcripts")
    
    # STEP 3: Score transcripts using the extracted reads
    log_info(f"Scoring {len(available_transcripts)} transcripts with available data")
    
    # Create batches for scoring - indices only
    transcript_indices = []
    batch_size = min(config.batch_size, 25)  # Smaller batches to reduce memory pressure
    for i in range(0, len(available_transcripts), batch_size):
        end = min(i + batch_size, len(available_transcripts))
        transcript_indices.append((i, end))
    
    # Define a memory-efficient function to score transcripts
    def score_transcripts_memory_efficient(batch_indices):
        start_idx, end_idx = batch_indices
        batch_transcripts = available_transcripts[start_idx:end_idx]
        
        # Score each transcript
        batch_results = []
        for tran in batch_transcripts:
            # Get transcript data as DataFrame
            tran_reads = transcript_storage.get_transcript_df(tran)
            if tran_reads.is_empty():
                continue
                
            # Filter ORFs for this transcript
            orfs_for_tran = orf_df.filter(pl.col("tran_id") == tran)
            if orfs_for_tran.is_empty():
                continue
                
            for typeorf in orfs_for_tran["type"].unique():
                orfs_filtered = orfs_for_tran.filter(pl.col("type") == typeorf)
                
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
                
                # Free memory after processing each ORF type
                orfs_filtered = None
                gc.collect()
        
        return batch_results
    
    # Score in parallel with staggered execution to limit memory
    all_results = []
    scoring_stagger_size = min(len(transcript_indices) // 4 + 1, max_workers)  # Process ~25% at a time
    
    for start_batch in range(0, len(transcript_indices), scoring_stagger_size):
        end_batch = min(start_batch + scoring_stagger_size, len(transcript_indices))
        batch_group = transcript_indices[start_batch:end_batch]
        
        log_info(f"Processing scoring batch group {start_batch//scoring_stagger_size + 1}: "
                 f"batches {start_batch} to {end_batch-1}")
        
        # Process this group in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(score_transcripts_memory_efficient, indices) for indices in batch_group]
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(futures):
                try:
                    batch_results = future.result()
                    if batch_results:
                        all_results.extend(batch_results)
                except Exception as exc:
                    log_error(f"Transcript batch scoring error: {exc}")
        
        # Force garbage collection after each stagger group
        gc.collect()
    
    # Combine all results
    if not all_results:
        log_warning("No ORFs were scored")
        return pl.DataFrame()
    
    try:
        # Combine results in chunks to avoid memory issues
        chunk_size = min(100, len(all_results))
        final_dfs = []
        
        for i in range(0, len(all_results), chunk_size):
            end = min(i + chunk_size, len(all_results))
            chunk_df = pl.concat(all_results[i:end])
            final_dfs.append(chunk_df)
            
            # Clear processed results to free memory
            for j in range(i, end):
                all_results[j] = None
            
            gc.collect()
        
        # Final concatenation
        final_df = pl.concat(final_dfs)
        
        duration = time.time() - start_time
        transcripts_per_second = len(available_transcripts) / duration
        log_info(f"Scoring completed in {duration:.2f} seconds ({transcripts_per_second:.2f} transcripts/sec)")
        return final_df
    except Exception as exc:
        log_error(f"Error combining all results: {exc}")
        return pl.DataFrame()
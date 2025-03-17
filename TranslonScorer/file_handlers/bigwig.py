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

@profile  # Add the memory profiler decorator
def scoring(bigwig, exon, orfs, old_scoring, sru_range, batch_size=50, max_workers=None):
    """
    Score ORFs using bigwig coverage data with high-performance optimization.
    
    This implementation uses a two-stage approach:
    1. Extract all BigWig data in parallel by region batches
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
    
    # Check if exon is a file path or a DataFrame
    if isinstance(exon, str):
        if os.path.getsize(exon) > 1e9:  # 1 GB
            exon_df = pl.scan_csv(exon, has_header=True, separator=",").collect()
        else:
            exon_df = pl.read_csv(exon, has_header=True, separator=",")
    else:
        exon_df = exon

    # Similar check for ORFs
    if isinstance(orfs, str):
        if os.path.getsize(orfs) > 1e9:  # 1 GB
            orf_df = pl.scan_csv(orfs, has_header=True, separator=",").collect()
        else:
            orf_df = pl.read_csv(orfs, has_header=True, separator=",")
    else:
        orf_df = orfs

    # Set max workers if not specified
    if max_workers is None:
        max_workers = os.cpu_count() or 4
    config.max_workers = max_workers
    
    log_info(f"Using {max_workers} workers for parallel processing")
    
    # STEP 1: Extract all regions from exon data
    all_regions = extract_regions(exon_df)
    total_regions = len(all_regions)
    log_info(f"Found {total_regions} regions to process across all transcripts")
    
    # Calculate unique transcripts without creating a large intermediate list
    unique_transcripts = set()
    # Process in smaller chunks to avoid memory spikes
    chunk_size = 10000
    total_chunks = (len(all_regions) + chunk_size - 1) // chunk_size  # Ceiling division
    
    for i in range(total_chunks):
        start_idx = i * chunk_size
        end_idx = min(start_idx + chunk_size, len(all_regions))
        
        # Process directly from all_regions without creating a new list
        unique_transcripts.update(all_regions[j][3] for j in range(start_idx, end_idx))
        
        # Optional: periodically force garbage collection for very large datasets
        if i > 0 and i % 10 == 0:  # Every 10 chunks
            gc.collect()

    total_unique_transcripts = len(unique_transcripts)
    unique_transcripts = None  # Free memory
    gc.collect()
    
    # Group regions into batches for better work distribution - using indices only
    region_batches = []
    batch_size = config.max_regions_per_worker
    for i in range(0, total_regions, batch_size):
        end = min(i + batch_size, total_regions)
        # Store only indices, not data copies
        region_batches.append((i, end))
    
    log_info(f"Distributing regions into {len(region_batches)} balanced batches")
    
    # STEP 2: Process all region batches in parallel with staggered execution
    transcript_reads = {}  # Will hold reads for all transcripts
    total_processed = 0
    total_failed = 0
    
    # Free up memory before parallel processing
    gc.collect()
    
    # Modified process_region_batch function to work with indices instead of data slices
    def process_region_batch_by_indices(indices, all_regions, bwfile_path, config):
        start_idx, end_idx = indices
        # Extract the actual data when needed inside the worker process
        batch_data = all_regions[start_idx:end_idx]
        # Call the original function with the extracted data
        return process_region_batch(batch_data, bwfile_path, config)
    
    # Staggered parallelization constants
    max_concurrent_jobs = min(max_workers * 2, 16)  # Limit concurrent jobs
    stagger_size = min(len(region_batches) // 10 + 1, max_concurrent_jobs)  # Process ~10% at a time
    log_interval = max(1, len(region_batches) // 20)  # Log ~20 times during processing
    
    completed_batches = 0
    total_batches = len(region_batches)
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Process batches in a staggered fashion
        for start_batch in range(0, total_batches, stagger_size):
            end_batch = min(start_batch + stagger_size, total_batches)
            current_batch_count = end_batch - start_batch
            
            log_info(f"Processing batch group {start_batch//stagger_size + 1}: batches {start_batch} to {end_batch-1}")
            
            # Submit a group of batches
            futures = []
            for i in range(start_batch, end_batch):
                indices = region_batches[i]
                futures.append(executor.submit(
                    process_region_batch_by_indices, 
                    indices, 
                    all_regions, 
                    bwfile_path, 
                    config
                ))
            
            # Process results as they complete for this group
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    completed_batches += 1
                    
                    # Merge results
                    transcript_reads.update(result["results"])
                    total_processed += result["processed"]
                    total_failed += result["failed"]
                    
                    # Log progress less frequently
                    if completed_batches % log_interval == 0 or completed_batches == total_batches:
                        progress = completed_batches / total_batches * 100
                        log_info(f"Processed {completed_batches}/{total_batches} region batches ({progress:.1f}%) "
                                 f"- {len(transcript_reads)}/{total_unique_transcripts} transcripts with data")
                        
                except Exception as exc:
                    log_error(f"Region batch processing error: {exc}")
            
            # Optional: force garbage collection after each stagger group
            gc.collect()
    
    # Clear regions data as it's no longer needed
    all_regions = None
    region_batches = None
    gc.collect()
    
    log_info(f"Completed BigWig data extraction: {total_processed} regions processed, {total_failed} failed")
    log_info(f"Extracted data for {len(transcript_reads)} transcripts")
    
    # STEP 3: Score transcripts using the extracted reads - Also with staggered parallelization
    # Get transcripts with available reads
    available_transcripts = list(transcript_reads.keys())
    log_info(f"Scoring {len(available_transcripts)} transcripts with available data")
    
    # Create batches for scoring - indices only
    transcript_indices = []
    batch_size = config.batch_size
    for i in range(0, len(available_transcripts), batch_size):
        end = min(i + batch_size, len(available_transcripts))
        transcript_indices.append((i, end))
    
    # Limit to test batch count if needed
    test_batch_count = 3  # Ensure at least one batch
    transcript_indices = transcript_indices[:test_batch_count]
    
    log_info(f"Processing {len(transcript_indices)} transcript batches for testing.")
    
    # Score in parallel - with staggered execution
    all_results = []
    
    # Free memory before second parallel processing
    gc.collect()
    
    # Define function that works with indices
    def score_transcript_batch_by_indices(indices, transcripts_list, transcript_reads, orf_df, old_scoring, sru_range):
        start_idx, end_idx = indices
        batch = transcripts_list[start_idx:end_idx]
        return score_transcript_batch(batch, transcript_reads, orf_df, old_scoring, sru_range)
    
    # Staggered parallelization for scoring
    scoring_stagger_size = min(len(transcript_indices) // 4 + 1, max_workers)  # Process ~25% at a time
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Process transcript batches in a staggered fashion
        for start_batch in range(0, len(transcript_indices), scoring_stagger_size):
            end_batch = min(start_batch + scoring_stagger_size, len(transcript_indices))
            
            log_info(f"Processing scoring batch group {start_batch//scoring_stagger_size + 1}: "
                     f"batches {start_batch} to {end_batch-1}")
            
            # Submit a group of scoring batches
            futures = []
            for i in range(start_batch, end_batch):
                indices = transcript_indices[i]
                futures.append(executor.submit(
                    score_transcript_batch_by_indices,
                    indices,
                    available_transcripts,
                    transcript_reads,
                    orf_df,
                    old_scoring,
                    sru_range
                ))
            
            # Process results as they complete for this group
            completed = 0
            total_futures = len(futures)
            
            for future in concurrent.futures.as_completed(futures):
                try:
                    batch_results = future.result()
                    if batch_results:
                        all_results.extend(batch_results)
                    
                    completed += 1
                    progress = (completed / total_futures) * 100
                    log_info(f"Scored batch {completed}/{total_futures} ({progress:.1f}%)")
                        
                except Exception as exc:
                    log_error(f"Transcript batch scoring error: {exc}")
            
            # Force garbage collection after each stagger group
            gc.collect()
    
    # Free memory after processing
    transcript_reads = None
    gc.collect()
    
    # Combine all results
    if not all_results:
        log_warning("No ORFs were scored")
        return pl.DataFrame()
    
    try:
        final_df = pl.concat(all_results)
        duration = time.time() - start_time
        transcripts_per_second = len(available_transcripts) / duration
        log_info(f"Scoring completed in {duration:.2f} seconds ({transcripts_per_second:.2f} transcripts/sec)")
        return final_df
    except Exception as exc:
        log_error(f"Error combining all results: {exc}")
        return pl.DataFrame()
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
                empty_dict = {"rise_up": defaultdict(float), "step_down": defaultdict(float)}
                
                # Find ORFs that need scoring
                find_orfs_start_time = time.time()
                emptyscore_df = existingscore(orfs_filtered, typeorf, empty_dict)
                find_orfs_duration = time.time() - find_orfs_start_time
                log_info(f"Time taken to find ORFs for type {typeorf} in {tran_id}: {find_orfs_duration:.4f} seconds")
                
                if not emptyscore_df.is_empty():
                    # Calculate new scores
                    calculate_scores_start_time = time.time()
                    scoredict = newscoring(
                        emptyscore_df, tran_reads, sru_range, typeorf, empty_dict
                    )
                    calculate_scores_duration = time.time() - calculate_scores_start_time
                    log_info(f"Time taken to calculate scores for type {typeorf} in {tran_id}: {calculate_scores_duration:.4f} seconds")
                    
                    # Assign scores to ORFs
                    assign_scores_start_time = time.time()
                    orfs_filtered = assigningscore(
                        orfs_filtered, scoredict, typeorf
                    )
                    assign_scores_duration = time.time() - assign_scores_start_time
                    log_info(f"Time taken to assign scores for type {typeorf} in {tran_id}: {assign_scores_duration:.4f} seconds")
                    
                    # Calculate global scores
                    global_scores_start_time = time.time()
                    orfs_filtered = globalscores(orfs_filtered, tran_reads, typeorf)
                    global_scores_duration = time.time() - global_scores_start_time
                    log_info(f"Time taken to calculate global scores for type {typeorf} in {tran_id}: {global_scores_duration:.4f} seconds")
                    
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
    start_time = time.time()
    log_info(f"Starting scoring with BigWig file: {bigwig}")
    
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
    
    # Load exon data efficiently
    log_info("Loading exon data...")
    if isinstance(exon, str):
        if os.path.exists(exon) and os.path.getsize(exon) > 1e9:  # 1 GB
            # For large files, stream and select only needed columns
            exon_df = pl.scan_csv(exon).select(["chr", "start", "stop", "tran_id"]).collect()
        else:
            exon_df = pl.read_csv(exon)
    else:
        exon_df = exon
    
    # Load ORF data efficiently
    log_info("Loading ORF data...")
    if isinstance(orfs, str):
        if os.path.exists(orfs) and os.path.getsize(orfs) > 1e9:  # 1 GB
            # For large files, stream and select only needed columns
            needed_cols = ["tran_id", "type", "start", "stop", "length"]
            orf_scan = pl.scan_csv(orfs)
            
            # Add any scoring-related columns
            for col in orf_scan.columns:
                if "score" in col.lower() or col in ["rise_up", "step_down", "hrf", "avg", "nzc"]:
                    needed_cols.append(col)
                    
            orf_df = orf_scan.select(needed_cols).collect()
        else:
            orf_df = pl.read_csv(orfs)
    else:
        orf_df = orfs
    
    # Get transcript IDs from ORFs for filtering
    log_info("Identifying transcripts with ORFs...")
    orf_tran_ids = set(orf_df["tran_id"].unique().to_list())
    log_info(f"Found {len(orf_tran_ids)} transcripts with ORFs")
    
    # Filter exons to only include transcripts with ORFs
    exon_df = exon_df.filter(pl.col("tran_id").is_in(orf_tran_ids))
    log_info(f"Filtered exon data to {exon_df.height} rows with relevant transcripts")
    
    # Extract regions from exon data
    log_info("Extracting genomic regions...")
    regions = extract_regions(exon_df, orf_tran_ids)
    log_info(f"Extracted {len(regions)} regions to process")
    
    # Free memory from exon data as it's no longer needed
    exon_df = None
    gc.collect()
    
    if not regions:
        log_warning("No valid regions found to process")
        return pl.DataFrame()
    
    # Create balanced batches for parallel processing
    log_info("Creating balanced processing batches...")
    
    # Group regions by chromosome for better BigWig performance
    regions_by_chrom = {}
    for region in regions:
        chrom = region[0]
        if chrom not in regions_by_chrom:
            regions_by_chrom[chrom] = []
        regions_by_chrom[chrom].append(region)
    
    # Create batches that preserve chromosome locality
    batches = []
    target_size = max(1, len(regions) // (max_workers * 3))  # Aim for 3x batches per worker
    
    for chrom, chrom_regions in regions_by_chrom.items():
        # Sort regions by start position for better locality
        chrom_regions.sort(key=lambda r: r[1])
        
        # Create batches of appropriate size
        for i in range(0, len(chrom_regions), target_size):
            batch = chrom_regions[i:min(i+target_size, len(chrom_regions))]
            if batch:
                batches.append((batch, bigwig))
    
    log_info(f"Created {len(batches)} processing batches")
    
    # Free memory from regions data
    regions = None
    regions_by_chrom = None
    gc.collect()
    
    # Process batches in parallel
    log_info("Starting BigWig data extraction...")
    transcript_data = {}
    total_processed = 0
    total_failed = 0
    
    # Process in smaller groups to manage memory
    group_size = min(len(batches), max_workers * 2)
    batch_groups = [batches[i:i+group_size] for i in range(0, len(batches), group_size)]
    
    for group_idx, batch_group in enumerate(batch_groups):
        log_info(f"Processing batch group {group_idx+1}/{len(batch_groups)}")
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_batch = {executor.submit(process_batch, batch_data): i 
                              for i, batch_data in enumerate(batch_group)}
            
            for future in concurrent.futures.as_completed(future_to_batch):
                try:
                    result = future.result()
                    
                    # Update transcript data
                    for tran_id, data in result["results"].items():
                        transcript_data[tran_id] = data
                    
                    total_processed += result["processed"]
                    total_failed += result["failed"]
                except Exception as e:
                    log_error(f"Error processing batch: {str(e)}")
        
        # Force garbage collection between groups
        gc.collect()
        
        # Log progress
        current_transcripts = len(transcript_data)
        log_info(f"Processed {current_transcripts} transcripts so far "
                f"({total_processed} regions, {total_failed} failed)")
    
    log_info(f"Completed BigWig data extraction for {len(transcript_data)} transcripts")
    
    # Free memory
    batches = None
    batch_groups = None
    gc.collect()
    
        # Add after creating transcript_data but before scoring:
    log_info("Filtering transcripts by coverage...")
    min_nonzero_reads = 10  # Minimum number of positions with non-zero coverage
    min_mean_coverage = 0.5  # Minimum mean coverage across transcript

    filtered_transcripts = {}
    for tran_id, tran_data in transcript_data.items():
        # For sparse representation
        if isinstance(tran_data, dict) and 'positions' in tran_data:
            # Check if we have enough non-zero positions
            if len(tran_data['positions']) < min_nonzero_reads:
                continue
                
            # Check mean coverage (only consider positions with data)
            if np.mean(tran_data['values']) < min_mean_coverage:
                continue
                
            filtered_transcripts[tran_id] = tran_data
        # For list representation
        elif isinstance(tran_data, list):
            nonzero_count = sum(1 for v in tran_data if v > 0)
            if nonzero_count < min_nonzero_reads:
                continue
                
            mean_coverage = sum(tran_data) / max(1, len(tran_data))
            if mean_coverage < min_mean_coverage:
                continue
                
            filtered_transcripts[tran_id] = tran_data

    log_info(f"Filtered out {len(transcript_data) - len(filtered_transcripts)} transcripts with insufficient coverage")
    log_info(f"Proceeding with {len(filtered_transcripts)} well-covered transcripts")

    # Then replace transcript_data with filtered_transcripts
    transcript_data = filtered_transcripts

    # Create a temporary file for streaming results
    temp_dir = tempfile.gettempdir()
    temp_results_file = os.path.join(temp_dir, f"translonscorer_results_{int(time.time())}.csv")
    log_info(f"Will stream results to temporary file: {temp_results_file}")
    
    # Prepare scoring tasks
    log_info("Preparing scoring tasks...")
    
    # Group ORFs by transcript and type
    log_info("Organizing ORF data by transcript...")
    orfs_by_tran = {}
    
    for row in orf_df.iter_rows(named=True):
        tran_id = row["tran_id"]
        
        # Skip transcripts with no BigWig data
        if tran_id not in transcript_data:
            continue
            
        if tran_id not in orfs_by_tran:
            orfs_by_tran[tran_id] = {}
        
        # Group by ORF type
        typeorf = row["type"]
        if typeorf not in orfs_by_tran[tran_id]:
            orfs_by_tran[tran_id][typeorf] = []
            
        orfs_by_tran[tran_id][typeorf].append(row)
    
    # Create scoring tasks
    scoring_tasks = []
    
    for tran_id in transcript_data.keys():
        if tran_id not in orfs_by_tran:
            continue
            
        # Combine ORFs for this transcript
        orf_rows = []
        for typeorf, type_rows in orfs_by_tran[tran_id].items():
            orf_rows.extend(type_rows)
        
        if not orf_rows:
            continue
            
        # Convert to dictionary format for serialization
        orf_dict = {}
        for key in orf_rows[0].keys():
            orf_dict[key] = [row[key] for row in orf_rows]
            
        # Add task
        scoring_tasks.append((tran_id, transcript_data[tran_id], orf_dict, old_scoring, sru_range))
    
    # Free memory
    orfs_by_tran = None
    orf_df = None
    gc.collect()
    
    log_info(f"Created {len(scoring_tasks)} scoring tasks")
    
    # Score in batches to manage memory
    log_info("Starting ORF scoring...")
    
    # Use fewer workers for scoring to manage memory
    scoring_workers = max(1, min(max_workers, 4))
    log_info(f"Using {scoring_workers} workers for scoring phase")
    
    # Process in smaller batches
    score_batch_size = min(100, max(10, len(scoring_tasks) // (scoring_workers * 4)))
    task_batches = [scoring_tasks[i:i+score_batch_size] 
                   for i in range(0, len(scoring_tasks), score_batch_size)]
    
    total_orfs = 0
    
    for batch_idx, task_batch in enumerate(task_batches):
        log_info(f"Processing scoring batch {batch_idx+1}/{len(task_batches)}")
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=scoring_workers) as executor:
            results = list(filter(None, executor.map(score_transcript, task_batch)))
            
        # Convert results to DataFrames
        if results:
            result_dfs = [pl.DataFrame(res) for res in results if res]
            
            if result_dfs:
                # Stream to disk
                if batch_idx == 0:
                    # First batch - create file
                    pl.concat(result_dfs).write_csv(temp_results_file)
                else:
                    # Append to existing file - open in append mode and write without header
                    with open(temp_results_file, 'a') as f:
                        pl.concat(result_dfs).write_csv(f, include_header=False)
                # Update counts
                total_orfs += sum(df.height for df in result_dfs)
                
                # Log progress
                log_info(f"Processed {total_orfs} ORFs so far")
                
        # Force garbage collection between batches
        result_dfs = None
        results = None
        gc.collect()
    
    # Free memory
    scoring_tasks = None
    task_batches = None
    transcript_data = None
    gc.collect()
    
    # Load final results
    if os.path.exists(temp_results_file) and os.path.getsize(temp_results_file) > 0:
        log_info(f"Loading final results from {temp_results_file}")
        
        try:
            # For large files, use lazy loading
            if os.path.getsize(temp_results_file) > 1e9:  # 1 GB
                final_df = pl.scan_csv(temp_results_file).collect()
            else:
                final_df = pl.read_csv(temp_results_file)
                
            log_info(f"Successfully loaded {final_df.height} scored ORFs")
            
            # Clean up temp file# Clean up temp file
            try:
                os.remove(temp_results_file)
                log_info(f"Temporary file removed: {temp_results_file}")
            except Exception as e:
                log_warning(f"Could not remove temporary file: {temp_results_file}")
            
            # Calculate and log performance metrics
            duration = time.time() - start_time
            orfs_per_second = final_df.height / max(0.001, duration)
            
            log_info(f"Scoring completed in {duration:.2f} seconds")
            log_info(f"Performance: {orfs_per_second:.2f} ORFs/second")
            
            return final_df
        except Exception as e:
            log_error(f"Error loading results: {str(e)}")
            return pl.DataFrame()
    else:
        log_warning("No results were generated or temporary file is empty")
        return pl.DataFrame()
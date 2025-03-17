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
    
    start_time = time.time()
    config = ProcessingConfig(batch_size=batch_size, sru_range=sru_range)
    
    bwfile_path = bigwig
    
    # Load data efficiently
    log_info("Loading exon and ORF data...")
    
    # Check if exon is a file path or a DataFrame
    if isinstance(exon, str):
        if os.path.exists(exon) and os.path.getsize(exon) > 1e9:  # 1 GB
            # For large files, only load needed columns
            exon_df = pl.read_csv(exon, columns=["chr", "start", "stop", "tran_id"])
        else:
            exon_df = pl.read_csv(exon)
    else:
        exon_df = exon

    # For ORFs
    if isinstance(orfs, str):
        if os.path.exists(orfs) and os.path.getsize(orfs) > 1e9:  # 1 GB
            # Create a LazyFrame for large ORF files
            orf_scan = pl.scan_csv(orfs)
            # Only collect essential columns
            needed_cols = ["tran_id", "type", "start", "stop", "length"]
            # Add any scoring-related columns
            for col in orf_scan.columns:
                if "score" in col.lower() or col in ["rise_up", "step_down"]:
                    needed_cols.append(col)
            orf_df = orf_scan.select(needed_cols).collect()
        else:
            orf_df = pl.read_csv(orfs)
    else:
        orf_df = orfs

    # Optimize worker count based on system resources
    if max_workers is None:
        try:
            import psutil
            # Use available memory to determine worker count
            available_memory_gb = psutil.virtual_memory().available / (1024**3)
            # Estimate 1GB per worker as baseline
            suggested_workers = max(1, min(os.cpu_count(), int(available_memory_gb)))
            max_workers = min(suggested_workers, 12)  # Cap at 12 workers
        except ImportError:
            # If psutil is not available, use a conservative approach
            max_workers = min(os.cpu_count() or 4, 6)  # Use at most 6 workers
    
    log_info(f"Using {max_workers} workers for parallel processing")
    
    # STEP 1: Extract regions in chunks to reduce memory pressure
    log_info("Extracting regions from exon data...")
    all_regions = []
    
    # Process in manageable chunks
    chunk_size = min(10000, exon_df.height // 10 + 1)  # Adaptive chunk size
    for chunk_df in process_df_in_chunks(exon_df, chunk_size):
        chunk_regions = extract_regions(chunk_df)
        all_regions.extend(chunk_regions)
        gc.collect()  # Clean up after each chunk
    
    total_regions = len(all_regions)
    log_info(f"Found {total_regions} regions to process")
    
    # STEP 2: Process region batches more efficiently
    # Create balanced batches for better throughput
    region_batch_size = min(config.max_regions_per_worker * 2, 500)  # Larger batch size for better performance
    region_batches = []
    for i in range(0, total_regions, region_batch_size):
        end = min(i + region_batch_size, total_regions)
        region_batches.append((i, end))
    
    log_info(f"Distributing regions into {len(region_batches)} batches")
    
    # Use a dictionary to store transcript data
    transcript_data = {}
    total_processed = 0
    total_failed = 0
    
    # Process in staggered groups
    stagger_size = min(len(region_batches) // 5 + 1, max_workers * 2)  # Process ~20% at a time
    
    for start_batch in range(0, len(region_batches), stagger_size):
        end_batch = min(start_batch + stagger_size, len(region_batches))
        batch_group = region_batches[start_batch:end_batch]
        
        log_info(f"Processing BigWig batch group {start_batch//stagger_size + 1}/{(len(region_batches) + stagger_size - 1)//stagger_size}: "
                 f"batches {start_batch} to {end_batch-1}")
        
        # Process this group with better parallelism
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Include chunk size parameter for performance tuning
            chunk_size = config.chunk_size
            batch_data_list = [(indices, all_regions, bwfile_path, chunk_size) for indices in batch_group]
            
            # Use a map instead of submitting individual futures for better scheduling
            for result in executor.map(process_batch_memory_efficient, batch_data_list):
                # Add transcript data directly to the dictionary
                for tran_id, reads in result["results"].items():
                    transcript_data[tran_id] = reads
                
                total_processed += result["processed"]
                total_failed += result["failed"]
        
        # Only do garbage collection between batch groups
        gc.collect()
    
    # Clear regions data
    all_regions = None
    region_batches = None
    gc.collect()
    
    log_info(f"Completed BigWig data extraction: {total_processed} regions processed, {total_failed} failed")
    log_info(f"Extracted data for {len(transcript_data)} transcripts")
    
    # STEP 3: Score transcripts with streaming to disk
    available_transcripts = list(transcript_data.keys())
    log_info(f"Scoring {len(available_transcripts)} transcripts with available data")
    
    # Create a temporary file for results
    temp_results_file = os.path.join(tempfile.gettempdir(), f"translonscorer_results_{int(time.time())}.csv")
    log_info(f"Streaming results to temporary file: {temp_results_file}")
    
    # Use adaptive batch size
    scoring_batch_size = max(5, min(50, len(available_transcripts) // (max_workers * 3) + 1))
    transcript_indices = []
    for i in range(0, len(available_transcripts), scoring_batch_size):
        end = min(i + scoring_batch_size, len(available_transcripts))
        transcript_indices.append((i, end))
    
    # Reduce worker count for scoring to minimize memory pressure
    scoring_workers = max(2, min(max_workers, 4))  # Use at most 4 workers for scoring
    log_info(f"Using {scoring_workers} workers for scoring to reduce memory pressure")
    
    # Track progress
    total_batches = len(transcript_indices)
    processed_batches = 0
    
    try:
        # Process in smaller staggered groups
        scoring_stagger_size = min(len(transcript_indices) // 5 + 1, scoring_workers)
        
        for start_batch in range(0, len(transcript_indices), scoring_stagger_size):
            end_batch = min(start_batch + scoring_stagger_size, len(transcript_indices))
            batch_group = transcript_indices[start_batch:end_batch]
            
            progress = (processed_batches / total_batches) * 100
            log_info(f"Processing scoring batch group {start_batch//scoring_stagger_size + 1}/{(len(transcript_indices) + scoring_stagger_size - 1)//scoring_stagger_size}: "
                     f"batches {start_batch} to {end_batch-1} ({progress:.1f}% complete)")
            
            # Process this group with fewer workers
            with concurrent.futures.ProcessPoolExecutor(max_workers=scoring_workers) as executor:
                batch_data_list = [
                    (indices, available_transcripts, transcript_data, orf_df, old_scoring, sru_range) 
                    for indices in batch_group
                ]
                
                # Process batches and stream results as they complete
                futures = [executor.submit(score_transcripts_memory_efficient, data) for data in batch_data_list]
                
                for future in concurrent.futures.as_completed(futures):
                    try:
                        batch_results = future.result()
                        if batch_results:
                            # Stream results to disk instead of keeping in memory
                            stream_results_to_disk(batch_results, temp_results_file)
                        processed_batches += 1
                    except Exception as exc:
                        log_error(f"Transcript batch scoring error: {exc}")
            
            # Force garbage collection between groups
            gc.collect()
        
        log_info(f"Finished scoring. Reading results from temporary file.")
        
        # Read the final results from disk
        if os.path.exists(temp_results_file) and os.path.getsize(temp_results_file) > 0:
            # Load in chunks if file is large
            if os.path.getsize(temp_results_file) > 1e9:  # > 1GB
                log_info("Result file is large, reading in chunks...")
                final_df = pl.scan_csv(temp_results_file).collect()
            else:
                final_df = pl.read_csv(temp_results_file)
                
            log_info(f"Successfully loaded {final_df.height} scored ORFs")
            
            # Clean up the temporary file
            try:
                os.remove(temp_results_file)
                log_info("Temporary file removed")
            except:
                log_warning(f"Could not remove temporary file: {temp_results_file}")
                
            duration = time.time() - start_time
            transcripts_per_second = len(available_transcripts) / duration
            log_info(f"Scoring completed in {duration:.2f} seconds ({transcripts_per_second:.2f} transcripts/sec)")
            return final_df
        else:
            log_warning("No ORFs were scored or temporary file is empty")
            return pl.DataFrame()
            
    except Exception as exc:
        log_error(f"Error during scoring process: {exc}")
        return pl.DataFrame()
    finally:
        # Ensure temp file is removed even on error
        if os.path.exists(temp_results_file):
            try:
                os.remove(temp_results_file)
            except:
                pass
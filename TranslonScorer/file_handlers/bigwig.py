"""
BigWig file handling functionality for TranslonScorer with optimized performance.

This module contains optimized functions for reading and processing BigWig files,
including conversion to other formats and coordinate transformations.
"""

from typing import Dict, List, Optional, Union, Tuple
import polars as pl
import pyBigWig as bw
from ..utils.logging import log_info, log_warning, log_error
from ..core.scoring import oldscoring, newscoring, globalscores, existingscore, assigningscore

import concurrent.futures
from functools import partial
import os
import time
import numpy as np
from dataclasses import dataclass
from contextlib import contextmanager


@dataclass
class ProcessingConfig:
    """Configuration for transcript processing."""
    max_workers: int = 0  # 0 means auto-detect based on CPU count
    batch_size: int = 50  # Number of transcripts to process in each batch
    sru_range: int = 100  # Range for SRU score calculation
    region_batch_size: int = 50  # Number of regions to process in each batch
    use_stats_when_possible: bool = True  # Use BigWig stats() for efficiency when appropriate


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


def transcriptreads_optimized(bwfile: bw.pyBigWig, exon_df: pl.DataFrame, 
                              config: ProcessingConfig = None) -> pl.DataFrame:
    """
    Optimized version that converts a BigWig file to a DataFrame based on provided exon annotation.

    Parameters:
    ----------
    bwfile : pyBigWig.pyBigWig
        An open BigWig file handle
    exon_df : polars.DataFrame
        DataFrame containing exon annotations with columns: chr, start, stop
    config : ProcessingConfig, optional
        Configuration for processing

    Returns:
    -------
    polars.DataFrame
        DataFrame containing transcript information with columns: tran_start, counts
    """
    if config is None:
        config = ProcessingConfig()
        
    if not isinstance(bwfile, bw.pyBigWig):
        raise TypeError("bwfile must be a pyBigWig handle")
    
    if exon_df.is_empty():
        raise ValueError("Empty exon DataFrame provided")
        
    if not all(col in exon_df.columns for col in ["chr", "start", "stop"]):
        raise ValueError("Exon DataFrame missing required columns (chr, start, stop)")

    # Calculate total bases to process for better memory pre-allocation
    approx_total_bases = exon_df.select(
        pl.sum((pl.col("stop") - pl.col("start")).alias("bases"))
    ).item()
    
    # Pre-allocate a numpy array for better memory efficiency
    reads = np.zeros(approx_total_bases, dtype=np.float32)
    current_pos = 0
    
    processed_regions = 0
    failed_regions = 0
    
    # Get chromosomes that exist in the bigwig file
    valid_chroms = set(bwfile.chroms().keys())
    
    start_time = time.time()
    
    try:
        # Explode the lists into rows and sort by chromosome and start position
        exon_exploded = exon_df.with_columns([
            pl.col("start").alias("start_list"),
            pl.col("stop").alias("stop_list")
        ]).explode(["start_list", "stop_list"])
        
        # Sort by chromosome and start position
        exon_exploded = exon_exploded.sort(["chr", "start_list"])
        
        # Process each chromosome separately to maintain order
        for chrom in exon_exploded["chr"].unique():
            if chrom not in valid_chroms:
                log_warning(f"Chromosome {chrom} not found in bigWig file")
                continue
                
            chrom_data = exon_exploded.filter(pl.col("chr") == chrom)
            
            # Get sorted positions for this chromosome
            starts = chrom_data["start_list"].to_list()
            stops = chrom_data["stop_list"].to_list()
            
            # Process regions in batches
            for i in range(0, len(starts), config.region_batch_size):
                batch_starts = starts[i:i+config.region_batch_size]
                batch_stops = stops[i:i+config.region_batch_size]
                
                # Filter out invalid regions
                valid_regions = [(s, e) for s, e in zip(batch_starts, batch_stops) if s < e]
                if not valid_regions:
                    continue
                
                for start, stop in valid_regions:
                    region_size = stop - start
                    
                    try:
                        # Use stats method for small regions when appropriate
                        if region_size < 100 and config.use_stats_when_possible:
                            stats = bwfile.stats(chrom, start, stop, type="mean")
                            if stats and stats[0] is not None:
                                # Ensure we have enough space in our array
                                if current_pos + region_size > len(reads):
                                    reads = np.resize(reads, max(len(reads)*2, current_pos + region_size))
                                
                                # Fill with the mean value
                                reads[current_pos:current_pos+region_size] = stats[0]
                                current_pos += region_size
                                processed_regions += 1
                            else:
                                # Fill with zeros
                                if current_pos + region_size > len(reads):
                                    reads = np.resize(reads, max(len(reads)*2, current_pos + region_size))
                                # Already zeros by initialization
                                current_pos += region_size
                                processed_regions += 1
                        else:
                            # For larger regions, use values()
                            values = bwfile.values(chrom, start, stop)
                            if values and any(v is not None for v in values):
                                # Convert to numpy array for faster processing
                                values_array = np.array([v if v is not None else 0.0 for v in values], dtype=np.float32)
                                
                                # Ensure we have enough space
                                if current_pos + len(values_array) > len(reads):
                                    reads = np.resize(reads, max(len(reads)*2, current_pos + len(values_array)))
                                
                                reads[current_pos:current_pos+len(values_array)] = values_array
                                current_pos += len(values_array)
                                processed_regions += 1
                            else:
                                # Fill with zeros
                                if current_pos + region_size > len(reads):
                                    reads = np.resize(reads, max(len(reads)*2, current_pos + region_size))
                                # Already zeros by initialization
                                current_pos += region_size
                                processed_regions += 1
                    except RuntimeError as e:
                        log_error(f"Could not read values for {chrom}:{start}-{stop}: {str(e)}")
                        failed_regions += 1
                        
    except Exception as e:
        log_error(f"Error processing exon data: {str(e)}")
        raise
    
    # Truncate the array to the actual data size
    reads = reads[:current_pos]
    
    if len(reads) == 0:
        msg = f"No reads extracted. Processed: {processed_regions}, Failed: {failed_regions}"
        log_error(msg, exception_type=ValueError)
        raise ValueError(msg)
    
    duration = time.time() - start_time
    rate = processed_regions / max(0.001, duration)
    log_info(f"Successfully processed {processed_regions} regions ({rate:.2f} regions/sec), {failed_regions} failed")
        
    return pl.DataFrame({
        "tran_start": np.arange(len(reads)),
        "counts": reads
    })


def process_transcript_batch(transcript_batch, exon_partitions, orf_partitions, bwfile_path, 
                             old_scoring, sru_range, config=None):
    """
    Process a batch of transcripts at once to reduce overhead.
    
    Args:
        transcript_batch: List of transcript IDs to process in this batch
        exon_partitions: List of DataFrames for exons
        orf_partitions: List of DataFrames for ORFs
        bwfile_path: Path to the BigWig file
        old_scoring: Whether to use old scoring method
        sru_range: Range for SRU score calculation
        config: ProcessingConfig object or None
        
    Returns:
        List of scored ORF DataFrames for all transcripts in the batch
    """
    if config is None:
        config = ProcessingConfig(sru_range=sru_range)
        
    batch_results = []
    
    # Open the BigWig file once for the entire batch
    with open_bigwig(bwfile_path) as bwfile:
        for tran in transcript_batch:
            exons = next((df for df in exon_partitions if tran in df["tran_id"].unique()), pl.DataFrame())
            orfs = next((df for df in orf_partitions if tran in df["tran_id"].unique()), pl.DataFrame())
            
            if exons.is_empty():
                log_warning(f"No exon data found for transcript {tran}")
                continue
            
            try:
                # Use optimized transcript reads function
                tran_reads = transcriptreads_optimized(bwfile, exons, config)
                
                if tran_reads.is_empty():
                    log_warning(f"No transcript reads found for {tran}")
                    continue
                
                for typeorf in orfs["type"].unique():
                    orfs_filtered = orfs.filter(pl.col("type") == typeorf)
                    
                    if orfs_filtered.is_empty():
                        continue
                    
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
                log_error(f"Error processing transcript {tran}: {str(e)}")
    
    return batch_results


def scoring_optimized(bigwig, exon, orfs, old_scoring, sru_range, config=None):
    """
    Optimized version of the scoring function with better performance.
    
    Args:
        bigwig (str): Path to bigwig file
        exon (str or DataFrame): Path to exon file or DataFrame
        orfs (str or DataFrame): Path to ORFs file or DataFrame
        old_scoring (bool): Whether to use old scoring method
        sru_range (int): Range for SRU score calculation
        config (ProcessingConfig, optional): Configuration for processing
        
    Returns:
        DataFrame: Scored ORFs
    """
    start_time = time.time()
    
    if config is None:
        config = ProcessingConfig(sru_range=sru_range)
    
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

    # Get unique transcripts
    unique_transcripts = orf_df["tran_id"].unique()
    total_transcripts = len(unique_transcripts)
    
    log_info(f"Processing {total_transcripts} unique transcripts")
    
    # Pre-group data for faster access
    exon_partitions = exon_df.partition_by("tran_id")
    orf_partitions = orf_df.partition_by("tran_id")
    
    # Set max workers if not specified
    if config.max_workers <= 0:
        config.max_workers = os.cpu_count() or 4
        
    log_info(f"Using {config.max_workers} workers with batch size {config.batch_size}")
    
    # Create batches for worker processes
    transcript_batches = []
    for i in range(0, total_transcripts, config.batch_size):
        end = min(i + config.batch_size, total_transcripts)
        transcript_batches.append(unique_transcripts[i:end])
    
    # Process batches in parallel
    all_results = []
    processed_batches = 0
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=config.max_workers) as executor:
        # Submit all batches
        future_to_batch = {
            executor.submit(
                process_transcript_batch, 
                batch, 
                exon_partitions, 
                orf_partitions, 
                bwfile_path, 
                old_scoring, 
                sru_range,
                config
            ): i for i, batch in enumerate(transcript_batches)
        }
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_batch):
            processed_batches += 1
            batch_index = future_to_batch[future]
            
            try:
                batch_results = future.result()
                all_results.extend(batch_results)
                
                # Log progress periodically
                if processed_batches % 5 == 0 or processed_batches == len(transcript_batches):
                    progress = processed_batches / len(transcript_batches) * 100
                    elapsed = time.time() - start_time
                    estimated_total = elapsed / (processed_batches / len(transcript_batches))
                    remaining = max(0, estimated_total - elapsed)
                    
                    log_info(f"Progress: {progress:.1f}% ({processed_batches}/{len(transcript_batches)} batches) "
                             f"- Elapsed: {elapsed:.1f}s, Remaining: {remaining:.1f}s")
                
            except Exception as exc:
                log_error(f"Batch {batch_index} processing generated an exception: {exc}")
    
    # Combine all results
    if not all_results:
        log_warning("No ORFs were scored")
        return pl.DataFrame()
    
    try:
        final_df = pl.concat(all_results)
        duration = time.time() - start_time
        log_info(f"Scoring completed in {duration:.2f} seconds "
                f"({total_transcripts/duration:.2f} transcripts/sec)")
        return final_df
    except Exception as exc:
        log_error(f"Error combining all results: {exc}")
        return pl.DataFrame()


# Replace the original functions with optimized versions
# You can comment these out if you want to keep both versions separately
transcriptreads = transcriptreads_optimized
process_transcript = process_transcript_batch  
scoring = scoring_optimized
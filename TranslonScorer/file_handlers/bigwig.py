"""
BigWig file handling functionality for TranslonScorer.

This module contains functions for reading and processing BigWig files,
including conversion to other formats and coordinate transformations.
"""

from typing import Dict, List, Optional, Union, Tuple, Any
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
    max_region_size: int = 1_000_000  # Maximum region size to process at once
    chunk_size: int = 50_000  # Size of chunks for processing large regions


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


def _flatten_list_or_value(value: Any) -> List:
    """
    Helper function to flatten a value that might be a list or a single value.
    Returns a list in either case.
    """
    if isinstance(value, list):
        return value
    else:
        return [value]


def transcriptreads(bwfile: bw.pyBigWig, exon_df: pl.DataFrame) -> pl.DataFrame:
    """
    Converts a BigWig file to a DataFrame based on provided exon annotation.
    
    Parameters:
    ----------
    bwfile : pyBigWig.pyBigWig
        An open BigWig file handle
    exon_df : polars.DataFrame
        DataFrame containing exon annotations with columns: chr, start, stop
        The start and stop columns may contain lists of positions
        
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
                log_warning(f"Mismatched start/stop lists for {chrom}: {len(starts)} starts, {len(stops)} stops")
                failed_regions += 1
                continue
            
            # Process each start/stop pair
            for start, stop in zip(starts, stops):
                # Ensure start and stop are integers
                try:
                    start = int(start)
                    stop = int(stop)
                except (ValueError, TypeError) as e:
                    log_warning(f"Invalid start/stop values for {chrom}: {start}, {stop} - {str(e)}")
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


def process_transcript(tran, exon_df, orf_df, bwfile_path, old_scoring, sru_range):
    """
    Process a single transcript.
    
    Args:
        tran: Transcript ID
        exon_df: DataFrame containing all exon data
        orf_df: DataFrame containing all ORF data
        bwfile_path: Path to the BigWig file
        old_scoring: Whether to use old scoring method
        sru_range: Range for SRU score calculation
        
    Returns:
        List of scored ORF DataFrames for this transcript
    """
    # Filter dataframes for this transcript
    exons = exon_df.filter(pl.col("tran_id") == tran)
    orfs = orf_df.filter(pl.col("tran_id") == tran)
    
    transcript_results = []
    
    if exons.is_empty():
        log_warning(f"No exon data found for transcript {tran}")
        return transcript_results
    
    try:
        with open_bigwig(bwfile_path) as bwfile:
            tran_reads = transcriptreads(bwfile, exons)

        if tran_reads.is_empty():
            log_warning(f"No transcript reads found for {tran}")
            return transcript_results
        
        for typeorf in orfs["type"].unique():
            orfs_filtered = orfs.filter(pl.col("type") == typeorf)

            if orfs_filtered.is_empty():
                log_warning(f"No ORFs found for type {typeorf} in transcript {tran}")
                continue

            if old_scoring:
                orfs_filtered = oldscoring(
                    orfs_filtered, tran_reads, sru_range, typeorf
                )
                transcript_results.append(orfs_filtered)
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
                    transcript_results.append(orfs_filtered)
    except Exception as e:
        log_error(f"Error processing transcript {tran}: {str(e)}")
    
    return transcript_results


def scoring(bigwig, exon, orfs, old_scoring, sru_range, batch_size=50, max_workers=None):
    """
    Score ORFs using bigwig coverage data with optimized performance.
    
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
    unique_transcripts = orf_df["tran_id"].unique().to_list()
    total_transcripts = len(unique_transcripts)
    
    log_info(f"Processing {total_transcripts} unique transcripts")
    
    # Set max workers if not specified
    if max_workers is None:
        max_workers = os.cpu_count() or 4
    
    log_info(f"Using {max_workers} workers with batch size {batch_size}")
    
    # Process transcripts in batches
    all_results = []
    processed_count = 0
    failed_count = 0
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        
        for batch_start in range(0, total_transcripts, batch_size):
            batch_end = min(batch_start + batch_size, total_transcripts)
            batch_transcripts = unique_transcripts[batch_start:batch_end]
            
            for tran in batch_transcripts:
                # Submit each transcript as a separate task
                futures.append(executor.submit(
                    process_transcript, 
                    tran, 
                    exon_df, 
                    orf_df, 
                    bwfile_path, 
                    old_scoring, 
                    sru_range
                ))

        # Process results as they complete
        completed = 0
        for future in concurrent.futures.as_completed(futures):
            try:
                transcript_results = future.result()
                if transcript_results:
                    all_results.extend(transcript_results)
                    processed_count += 1
                
                completed += 1
                if completed % 100 == 0 or completed == len(futures):
                    progress = completed / len(futures) * 100
                    elapsed = time.time() - start_time
                    rate = completed / elapsed
                    
                    log_info(f"Processed {completed}/{len(futures)} transcripts ({progress:.1f}%) - "
                             f"Rate: {rate:.2f} transcripts/sec")
                    
            except Exception as exc:
                log_error(f"Transcript processing generated an exception: {exc}")
                failed_count += 1

    # Log final processing statistics
    log_info(f"Completed processing: {processed_count} successful, {failed_count} failed")
    
    # Combine all results
    if not all_results:
        log_warning("No ORFs were scored")
        return pl.DataFrame()
    
    try:
        final_df = pl.concat(all_results)
        duration = time.time() - start_time
        log_info(f"Scoring completed in {duration:.2f} seconds ({total_transcripts/duration:.2f} transcripts/sec)")
        return final_df
    except Exception as exc:
        log_error(f"Error combining all results: {exc}")
        return pl.DataFrame()
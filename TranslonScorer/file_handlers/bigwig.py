"""
BigWig file handling functionality for TranslonScorer.

This module contains functions for reading and processing BigWig files,
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
from dataclasses import dataclass
from contextlib import contextmanager

@dataclass
class ProcessingConfig:
    """Configuration for transcript processing."""
    max_workers: int
    batch_size: int
    sru_range: int

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

    Raises:
    ------
    ValueError
        If exon_df is empty or no reads could be extracted
    RuntimeError
        If there are issues reading values from the bigWig file
    """
    if not isinstance(bwfile, bw.pyBigWig):
        raise TypeError("bwfile must be a pyBigWig handle")
    
    if exon_df.is_empty():
        raise ValueError("Empty exon DataFrame provided")
        
    if not all(col in exon_df.columns for col in ["chr", "start", "stop"]):
        raise ValueError("Exon DataFrame missing required columns (chr, start, stop)")

    reads: List[float] = []
    processed_regions = 0
    failed_regions = 0
    
    # Explode the lists into rows and sort by chromosome and start position
    try:
        exon_exploded = exon_df.with_columns([
            pl.col("start").alias("start_list"),
            pl.col("stop").alias("stop_list")
        ]).explode(["start_list", "stop_list"])
        
        # Sort by chromosome and start position
        exon_exploded = exon_exploded.sort(["chr", "start_list"])
        
        # Process each chromosome separately to maintain order
        for chrom in exon_exploded["chr"].unique():
            chrom_data = exon_exploded.filter(pl.col("chr") == chrom)
            
            # Get sorted positions for this chromosome
            starts = chrom_data["start_list"].to_list()
            stops = chrom_data["stop_list"].to_list()
            
            # Validate chromosome exists in bigwig file
            if not bwfile.chroms().get(chrom):
                log_warning(f"Chromosome {chrom} not found in bigWig file")
                continue
                
            # Process regions for this chromosome
            for start, stop in zip(starts, stops):
                try:
                    if start >= stop:
                        log_warning(f"Invalid region {chrom}:{start}-{stop} (start >= stop)")
                        failed_regions += 1
                        continue
                        
                    values = bwfile.values(chrom, start, stop)
                    if values and any(v is not None for v in values):
                        reads.extend(v if v is not None else 0.0 for v in values)
                        processed_regions += 1
                    else:
                        failed_regions += 1
                except RuntimeError as e:
                    log_error(f"Could not read values for {chrom}:{start}-{stop}: {str(e)}")
                    failed_regions += 1
                    
    except Exception as e:
        log_error(f"Error processing exon data: {str(e)}")
        raise
    
    if not reads:
        msg = f"No reads extracted. Processed: {processed_regions}, Failed: {failed_regions}"
        log_error(msg, exception_type=ValueError)
        raise ValueError(msg)
        
    log_info(f"Successfully processed {processed_regions} regions, {failed_regions} failed")
        
    return pl.DataFrame({
        "tran_start": range(len(reads)),
        "counts": reads
    })

def process_transcript(tran, exon_partitions, orf_partitions, bwfile_path, old_scoring, sru_range):
    """
    Process a single transcript.
    
    Args:
        tran: Transcript ID
        exon_partitions: List of DataFrames for exons
        orf_partitions: List of DataFrames for ORFs
        bwfile_path: Path to the BigWig file
        old_scoring: Whether to use old scoring method
        sru_range: Range for SRU score calculation
        
    Returns:
        List of scored ORF DataFrames for this transcript
    """
    exons = next((df for df in exon_partitions if tran in df["tran_id"].unique()), pl.DataFrame())
    orfs = next((df for df in orf_partitions if tran in df["tran_id"].unique()), pl.DataFrame())
    
    transcript_results = []
    
    if exons.is_empty():
        log_warning(f"No exon data found for transcript {tran}")
        return transcript_results
    
    with bw.open(bwfile_path) as bwfile:
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
    
    return transcript_results

def scoring(bigwig, exon, orfs, old_scoring, sru_range, batch_size=1000, max_workers=None):
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
    
    # Pre-group data for faster access
    exon_partitions = exon_df.partition_by("tran_id")
    orf_partitions = orf_df.partition_by("tran_id")
    
    # Pre-read all necessary data from the BigWig file into memory
    with bw.open(bwfile_path) as bwfile:
        all_reads = {}
        for chrom in exon_df["chr"].unique():
            # Read all regions for this chromosome at once
            chrom_exons = exon_df.filter(pl.col("chr") == chrom)
            starts = chrom_exons["start"].to_list()
            stops = chrom_exons["stop"].to_list()

            # Ensure starts and stops are numeric
            starts = [int(start) for start in starts]
            stops = [int(stop) for stop in stops]

            all_reads[chrom] = bwfile.values(chrom, starts, stops)

    # Process transcripts in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for tran in unique_transcripts:
            futures.append(executor.submit(process_transcript, tran, exon_partitions, orf_partitions, all_reads, old_scoring, sru_range))

        # Print progress updates
        for i, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            try:
                transcript_results = future.result()
                if transcript_results:
                    all_results.extend(transcript_results)
            except Exception as exc:
                log_error(f"Transcript processing generated an exception: {exc}")
            print(f"Processed {i}/{len(futures)} transcripts...")  # Progress update

    # Combine all results
    if not all_results:
        log_warning("No ORFs were scored")
        return pl.DataFrame()
    
    try:
        final_df = pl.concat(all_results)
        return final_df
    except Exception as exc:
        log_error(f"Error combining all results: {exc}")
        return pl.DataFrame()
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




def process_transcript(tran, exon_dict, orf_dict, bwfile, old_scoring, sru_range):
    """
    Process a single transcript.
    
    Args:
        tran: Transcript ID
        exon_dict: Dictionary mapping transcript IDs to exon DataFrames
        orf_dict: Dictionary mapping transcript IDs to ORF DataFrames
        bwfile: BigWig file handle
        old_scoring: Whether to use old scoring method
        sru_range: Range for SRU score calculation
        
    Returns:
        List of scored ORF DataFrames for this transcript
    """
    exons = exon_dict.get(tran, pl.DataFrame())
    orfs = orf_dict.get(tran, pl.DataFrame())
    
    transcript_results = []
    
    if exons.is_empty():
        log_warning(f"No exon data found for transcript {tran}")
        return transcript_results
    
    tran_reads = transcriptreads(bwfile, exons)
    if tran_reads.is_empty():
        return transcript_results
        
    for typeorf in orfs["type"].unique():
        orfs_filtered = orfs.filter(pl.col("type") == typeorf)
        
        if orfs_filtered.is_empty():
            continue
            
        if old_scoring:
            orfs_filtered = oldscoring(
                orfs_filtered, tran_reads, sru_range, typeorf
            )
            # Ensure orfs_filtered is a DataFrame
            if isinstance(orfs_filtered, dict):
                orfs_filtered = pl.DataFrame(orfs_filtered)
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
                # Ensure orfs_filtered is a DataFrame
                if isinstance(orfs_filtered, dict):
                    orfs_filtered = pl.DataFrame(orfs_filtered)
                transcript_results.append(orfs_filtered)
    
    return transcript_results

def scoring(bigwig, exon, orfs, old_scoring, sru_range, max_workers=None, batch_size=1000):
    """
    Score ORFs using bigwig coverage data with optimized performance.
    
    Args:
        bigwig (str): Path to bigwig file
        exon (str or DataFrame): Path to exon file or DataFrame
        orfs (str or DataFrame): Path to ORFs file or DataFrame
        old_scoring (bool): Whether to use old scoring method
        sru_range (int): Range for SRU score calculation
        max_workers (int, optional): Maximum number of worker processes
        batch_size (int): Size of transcript batches for processing
        
    Returns:
        DataFrame: Scored ORFs
    """
    log_info("Opening bigWig file")
    bwfile = bw.open(bigwig)
    
    if not bwfile.isBigWig():
        error_msg = f"File {bigwig} is not a valid bigWig file"
        log_error(error_msg)
        raise ValueError(error_msg)

    log_info("Loading exon and ORF data")
    
    # Check if exon is a file path or a DataFrame
    if isinstance(exon, str):
        log_info(f"Reading exon data from file: {exon}")
        # Use scan_csv for lazy evaluation if file is large
        if os.path.getsize(exon) > 1e9:  # 1 GB
            exon_df = pl.scan_csv(exon, has_header=True, separator=",").collect()
        else:
            exon_df = pl.read_csv(exon, has_header=True, separator=",")
    else:
        log_info("Using provided exon DataFrame")
        exon_df = exon

    # Optimize string conversions - more efficient batch approach
    string_columns = ["start", "stop", "tran_start", "tran_stop"]
    conversions = []
    for col in string_columns:
        if col in exon_df.columns:
            # Get the column's dtype from the DataFrame schema
            col_dtype = exon_df.schema[col]
            if col_dtype == pl.Utf8:  # Check if the column is string type
                # Convert string to list of integers
                conversions.append(
                    pl.col(col)
                    .str.split(",")
                    .map_elements(lambda x: [int(i) for i in x])
                    .alias(col)
                )
    
    # Apply all conversions in one operation if any exist
    if conversions:
        exon_df = exon_df.with_columns(conversions)
    
    # Similar check for ORFs
    if isinstance(orfs, str):
        log_info(f"Reading ORFs data from file: {orfs}")
        if os.path.getsize(orfs) > 1e9:  # 1 GB
            orf_df = pl.scan_csv(orfs, has_header=True, separator=",").collect()
        else:
            orf_df = pl.read_csv(orfs, has_header=True, separator=",")
    else:
        log_info("Using provided ORFs DataFrame")
        orf_df = orfs

    # Get unique transcripts
    unique_transcripts = orf_df["tran_id"].unique()
    total_transcripts = len(unique_transcripts)
    log_info(f"Scoring {total_transcripts} transcripts")
    
    # Pre-group data for faster access
    log_info("Pre-grouping data for faster access")

    exon_dict = {group[0]: group[1] for group in enumerate(exon_df.partition_by("tran_id"))}
    orf_dict = {group[0]: group[1] for group in enumerate(orf_df.partition_by("tran_id"))}
    
    # Process in batches to manage memory
    all_results = []
    
    # Determine max_workers based on CPU count if not specified
    if max_workers is None:
        max_workers = min(os.cpu_count() or 4, 8)  # Reasonable default
    
    log_info(f"Using {max_workers} parallel workers")
    
    for batch_start in range(0, total_transcripts, batch_size):
        batch_end = min(batch_start + batch_size, total_transcripts)
        batch_transcripts = unique_transcripts[batch_start:batch_end]
        
        log_info(f"Processing batch {batch_start//batch_size + 1}/{(total_transcripts + batch_size - 1)//batch_size}: "
                 f"transcripts {batch_start+1} to {batch_end}")
        
        # Create partial function with fixed arguments
        process_func = partial(
            process_transcript, 
            exon_dict=exon_dict, 
            orf_dict=orf_dict, 
            bwfile=bwfile, 
            old_scoring=old_scoring, 
            sru_range=sru_range
        )
        
        batch_results = []
        # Process batch in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_tran = {executor.submit(process_func, tran): tran for tran in batch_transcripts}
            
            for future in concurrent.futures.as_completed(future_to_tran):
                tran = future_to_tran[future]
                try:
                    transcript_results = future.result()
                    if transcript_results:
                        batch_results.extend(transcript_results)
                except Exception as exc:
                    log_error(f"Transcript {tran} generated an exception: {exc}")
        
        # Combine batch results and free memory periodically
        if batch_results:
            try:
                combined_batch = pl.concat(batch_results)
                all_results.append(combined_batch)
                # Free memory
                del batch_results
            except Exception as exc:
                log_error(f"Error combining batch results: {exc}")
    
    # Close BigWig file when done with all processing
    bwfile.close()
    
    # Combine all results
    if not all_results:
        log_warning("No ORFs were scored")
        return pl.DataFrame()
    
    try:
        log_info("Combining all scored ORFs")
        final_df = pl.concat(all_results)
        log_info(f"Scored {len(final_df)} ORFs in total")
        return final_df
    except Exception as exc:
        log_error(f"Error combining all results: {exc}")
        return pl.DataFrame()
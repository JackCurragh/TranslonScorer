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




def process_transcript(
    tran: str,
    exon_dict: Dict[str, pl.DataFrame],
    orf_dict: Dict[str, pl.DataFrame],
    bwfile: bw.pyBigWig,
    old_scoring: bool,
    sru_range: int
) -> List[pl.DataFrame]:
    """
    Process a single transcript.
    
    Parameters
    ----------
    tran : str
        Transcript ID
    exon_dict : Dict[str, pl.DataFrame]
        Dictionary mapping transcript IDs to exon DataFrames
    orf_dict : Dict[str, pl.DataFrame]
        Dictionary mapping transcript IDs to ORF DataFrames
    bwfile : pyBigWig.pyBigWig
        BigWig file handle
    old_scoring : bool
        Whether to use old scoring method
    sru_range : int
        Range for SRU score calculation
        
    Returns
    -------
    List[pl.DataFrame]
        List of scored ORF DataFrames for this transcript
        
    Raises
    ------
    ValueError
        If sru_range is invalid or required data is missing
    TypeError
        If input types are incorrect
    """
    # Input validation
    if not isinstance(tran, str):
        raise TypeError("Transcript ID must be a string")
    if not isinstance(bwfile, bw.pyBigWig):
        raise TypeError("bwfile must be a pyBigWig handle")
    if not isinstance(sru_range, int) or sru_range <= 0:
        raise ValueError("sru_range must be a positive integer")
    if not isinstance(old_scoring, bool):
        raise TypeError("old_scoring must be a boolean")
        
    transcript_results: List[pl.DataFrame] = []
    
    # Get and validate exon data
    exons = exon_dict.get(tran, pl.DataFrame())
    if exons.is_empty():
        log_warning(f"No exon data found for transcript {tran}")
        return transcript_results
        
    # Get and validate ORF data
    orfs = orf_dict.get(tran, pl.DataFrame())
    if orfs.is_empty():
        log_warning(f"No ORF data found for transcript {tran}")
        return transcript_results
        
    # Initialize score cache
    score_cache = {"rise_up": {}, "step_down": {}}
    
    try:
        # Get transcript reads
        tran_reads = transcriptreads(bwfile, exons)
        if tran_reads.is_empty():
            return transcript_results
            
        # Process each ORF type
        for typeorf in orfs["type"].unique():
            try:
                orfs_filtered = orfs.filter(pl.col("type") == typeorf)
                
                if orfs_filtered.is_empty():
                    continue
                    
                if old_scoring:
                    scored_orfs = oldscoring(orfs_filtered, tran_reads, sru_range, typeorf)
                    if isinstance(scored_orfs, dict):
                        scored_orfs = pl.DataFrame(scored_orfs)
                    if not scored_orfs.is_empty():
                        transcript_results.append(scored_orfs)
                else:
                    # Use cached scoring approach
                    emptyscore_df = existingscore(orfs_filtered, typeorf, score_cache)
                    if not emptyscore_df.is_empty():
                        score_cache = newscoring(
                            emptyscore_df, tran_reads, sru_range, typeorf, score_cache
                        )
                        scored_orfs = assigningscore(orfs_filtered, score_cache, typeorf)
                        scored_orfs = globalscores(scored_orfs, tran_reads, typeorf)
                        if isinstance(scored_orfs, dict):
                            scored_orfs = pl.DataFrame(scored_orfs)
                        if not scored_orfs.is_empty():
                            transcript_results.append(scored_orfs)
                            
            except Exception as e:
                log_error(f"Error processing ORF type {typeorf} for transcript {tran}: {str(e)}")
                continue
                
    except Exception as e:
        log_error(f"Error processing transcript {tran}: {str(e)}")
        
    return transcript_results

def scoring(
    bigwig: str,
    exon: Union[str, pl.DataFrame],
    orfs: Union[str, pl.DataFrame],
    old_scoring: bool,
    sru_range: int,
    max_workers: Optional[int] = None,
    batch_size: int = 1000,
    timeout: int = 3600  # 1 hour timeout per batch
) -> pl.DataFrame:
    """
    Score ORFs using bigwig coverage data with optimized performance.
    
    Parameters
    ----------
    bigwig : str
        Path to bigwig file
    exon : Union[str, pl.DataFrame]
        Path to exon file or DataFrame
    orfs : Union[str, pl.DataFrame]
        Path to ORFs file or DataFrame
    old_scoring : bool
        Whether to use old scoring method
    sru_range : int
        Range for SRU score calculation
    max_workers : Optional[int]
        Maximum number of worker processes
    batch_size : int
        Size of transcript batches for processing
    timeout : int
        Maximum time in seconds to wait for a batch to complete
        
    Returns
    -------
    pl.DataFrame
        Scored ORFs DataFrame
        
    Raises
    ------
    ValueError
        If input parameters are invalid or processing fails
    TypeError
        If input types are incorrect
    TimeoutError
        If batch processing exceeds timeout
    """
    # Input validation
    if not isinstance(bigwig, str):
        raise TypeError("bigwig must be a string path")
    if not isinstance(old_scoring, bool):
        raise TypeError("old_scoring must be a boolean")
    if not isinstance(sru_range, int) or sru_range <= 0:
        raise ValueError("sru_range must be a positive integer")
    if not isinstance(batch_size, int) or batch_size <= 0:
        raise ValueError("batch_size must be a positive integer")
    if max_workers is not None and (not isinstance(max_workers, int) or max_workers <= 0):
        raise ValueError("max_workers must be a positive integer")

    # Use context manager for safe file handling
    with open_bigwig(bigwig) as bwfile:
        log_info("Loading exon and ORF data")
        
        try:
            # Load and validate exon data
            exon_df = _load_data(exon, "exon")
            if exon_df.is_empty():
                raise ValueError("Empty exon DataFrame")
                
            # Load and validate ORF data
            orf_df = _load_data(orfs, "ORF")
            if orf_df.is_empty():
                raise ValueError("Empty ORF DataFrame")
                
            # Convert string columns to proper types
            exon_df = _convert_string_columns(exon_df)
            
            # Get unique transcripts and create processing config
            unique_transcripts = orf_df["tran_id"].unique()
            total_transcripts = len(unique_transcripts)
            log_info(f"Scoring {total_transcripts} transcripts")
            
            if total_transcripts == 0:
                raise ValueError("No transcripts found in ORF data")
            
            # Configure processing parameters
            if max_workers is None:
                max_workers = min(os.cpu_count() or 4, 8)
            config = ProcessingConfig(max_workers=max_workers, batch_size=batch_size, sru_range=sru_range)
            
            # Pre-group data for faster access
            log_info("Pre-grouping data for faster access")
            exon_dict = {group[0]: group[1] for group in exon_df.partition_by("tran_id")}
            orf_dict = {group[0]: group[1] for group in orf_df.partition_by("tran_id")}
            
            # Process batches with generator to manage memory
            return _process_batches(
                unique_transcripts=unique_transcripts,
                config=config,
                exon_dict=exon_dict,
                orf_dict=orf_dict,
                bwfile=bwfile,
                old_scoring=old_scoring,
                timeout=timeout
            )
            
        except Exception as e:
            log_error(f"Error in scoring process: {str(e)}")
            raise

def _load_data(data: Union[str, pl.DataFrame], data_type: str) -> pl.DataFrame:
    """Helper function to load data from file or DataFrame."""
    if isinstance(data, str):
        log_info(f"Reading {data_type} data from file: {data}")
        try:
            if os.path.getsize(data) > 1e9:  # 1 GB
                return pl.scan_csv(data, has_header=True, separator=",").collect()
            return pl.read_csv(data, has_header=True, separator=",")
        except Exception as e:
            raise ValueError(f"Error reading {data_type} file: {str(e)}")
    else:
        log_info(f"Using provided {data_type} DataFrame")
        if not isinstance(data, pl.DataFrame):
            raise TypeError(f"{data_type} must be a string path or Polars DataFrame")
        return data

def _convert_string_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Helper function to convert string columns to proper types."""
    string_columns = ["start", "stop", "tran_start", "tran_stop"]
    conversions = []
    
    for col in string_columns:
        if col in df.columns and df.schema[col] == pl.Utf8:
            conversions.append(
                pl.col(col)
                .str.split(",")
                .map_elements(lambda x: [int(i) for i in x])
                .alias(col)
            )
    
    return df.with_columns(conversions) if conversions else df

def _process_batches(
    unique_transcripts: pl.Series,
    config: ProcessingConfig,
    exon_dict: Dict[str, pl.DataFrame],
    orf_dict: Dict[str, pl.DataFrame],
    bwfile: bw.pyBigWig,
    old_scoring: bool,
    timeout: int
) -> pl.DataFrame:
    """Process transcripts in batches with memory management."""
    all_results = []
    total_transcripts = len(unique_transcripts)
    
    log_info(f"Using {config.max_workers} parallel workers")
    
    for batch_start in range(0, total_transcripts, config.batch_size):
        batch_end = min(batch_start + config.batch_size, total_transcripts)
        batch_transcripts = unique_transcripts[batch_start:batch_end]
        
        log_info(f"Processing batch {batch_start//config.batch_size + 1}/"
                f"{(total_transcripts + config.batch_size - 1)//config.batch_size}")
        
        try:
            batch_results = _process_batch(
                batch_transcripts=batch_transcripts,
                config=config,
                exon_dict=exon_dict,
                orf_dict=orf_dict,
                bwfile=bwfile,
                old_scoring=old_scoring,
                timeout=timeout
            )
            
            if batch_results is not None and not batch_results.is_empty():
                all_results.append(batch_results)
                
        except TimeoutError as e:
            log_error(f"Batch processing timeout: {str(e)}")
            continue
        except Exception as e:
            log_error(f"Error processing batch: {str(e)}")
            continue
    
    # Combine all results
    if not all_results:
        log_warning("No ORFs were scored")
        return pl.DataFrame()
    
    try:
        log_info("Combining all scored ORFs")
        final_df = pl.concat(all_results)
        log_info(f"Scored {len(final_df)} ORFs in total")
        return final_df
    except Exception as e:
        log_error(f"Error combining all results: {str(e)}")
        return pl.DataFrame()

def _process_batch(
    batch_transcripts: pl.Series,
    config: ProcessingConfig,
    exon_dict: Dict[str, pl.DataFrame],
    orf_dict: Dict[str, pl.DataFrame],
    bwfile: bw.pyBigWig,
    old_scoring: bool,
    timeout: int
) -> Optional[pl.DataFrame]:
    """Process a single batch of transcripts."""
    process_func = partial(
        process_transcript,
        exon_dict=exon_dict,
        orf_dict=orf_dict,
        bwfile=bwfile,
        old_scoring=old_scoring,
        sru_range=config.sru_range
    )
    
    batch_results = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=config.max_workers) as executor:
        future_to_tran = {
            executor.submit(process_func, tran): tran
            for tran in batch_transcripts
        }
        
        try:
            for future in concurrent.futures.as_completed(future_to_tran, timeout=timeout):
                tran = future_to_tran[future]
                try:
                    transcript_results = future.result()
                    if transcript_results:
                        batch_results.extend(transcript_results)
                except Exception as e:
                    log_error(f"Transcript {tran} generated an exception: {str(e)}")
                    continue
                    
        except concurrent.futures.TimeoutError:
            log_error(f"Batch processing exceeded timeout of {timeout} seconds")
            executor.shutdown(wait=False)
            raise TimeoutError(f"Batch processing exceeded timeout of {timeout} seconds")
    
    if not batch_results:
        return None
        
    try:
        return pl.concat(batch_results)
    except Exception as e:
        log_error(f"Error combining batch results: {str(e)}")
        return None
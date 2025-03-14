"""
Core scoring functionality for TranslonScorer - Optimized Version.

This module contains all functions related to scoring ORFs, including:
- Memory-efficient scoring methods
- SRU (Start/Stop Ramp Up/Down) calculations
- Region scoring calculations
- ORF position classification

All implementations use memory-efficient techniques:
- Reduced DataFrame copies
- Chunked processing for large datasets
- Efficient memory management with garbage collection
"""

import polars as pl
import gc
from ..utils.logging import log_info, log_warning, log_error
from .orffinder import classify_orf
from .orffinder import getexons_and_cds


def sru_score(position, tran_reads, sru_range, direction):
    """
    Calculate the SRU (Start/Stop Ramp Up/Down) score for a given position.
    Memory-efficient implementation.
    
    Args:
        position (int): Position to calculate score for
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for calculation
        direction (int): 0 for start (up), 1 for stop (down)
        
    Returns:
        float: Calculated SRU score
    """
    try:
        # Define filter expressions based on direction
        if direction == 0:
            # For start (rise up)
            before_expr = (
                (pl.col("tran_start") >= position - sru_range) &
                (pl.col("tran_start") < position)
            )
            after_expr = (
                (pl.col("tran_start") >= position) &
                (pl.col("tran_start") < position + sru_range)
            )
        else:
            # For stop (step down)
            before_expr = (
                (pl.col("tran_start") >= position) &
                (pl.col("tran_start") < position + sru_range)
            )
            after_expr = (
                (pl.col("tran_start") >= position + sru_range) &
                (pl.col("tran_start") < position + (2 * sru_range))
            )

        # Calculate means directly on filtered DataFrames without creating full copies
        before_mean = tran_reads.filter(before_expr)["counts"].mean() or 0
        after_mean = tran_reads.filter(after_expr)["counts"].mean() or 0
        
        # Return score based on direction
        return (after_mean - before_mean) if direction == 0 else (before_mean - after_mean)
            
    except Exception as e:
        log_error(f"Error calculating SRU score: {str(e)}")
        return 0.0


def calculate_scores(start, stop, tran_reads):
    """
    Calculate HRF, average, and non-zero coverage scores for a region.
    Memory-efficient implementation.
    
    Args:
        start (int): Start position
        stop (int): Stop position
        tran_reads (DataFrame): Transcript reads data
        
    Returns:
        tuple: (HRF score, average score, non-zero coverage score)
    """
    try:
        # Filter data with a single expression
        region_expr = (pl.col("tran_start") >= start) & (pl.col("tran_start") <= stop)
        
        # Get counts directly without creating a full region DataFrame
        filtered_counts = tran_reads.filter(region_expr)["counts"]
        
        if len(filtered_counts) == 0:
            return (0.0, 0.0, 0.0)
        
        # Calculate scores directly from filtered counts
        hrf = filtered_counts.max() or 0
        avg = filtered_counts.mean() or 0
        nzc = (filtered_counts > 0).sum() / len(filtered_counts) if len(filtered_counts) > 0 else 0
        
        return (hrf, avg, nzc)
        
    except Exception as e:
        log_error(f"Error calculating region scores: {str(e)}")
        return (0.0, 0.0, 0.0)


def process_orfs(df, tran_reads, sru_range, typeorf, batch_size=50, chunk_size=500):
    """
    Process ORFs in chunks, with each chunk processed in batches.
    This approach reduces memory usage significantly for large datasets.
    
    Args:
        df (DataFrame): ORF DataFrame to process
        tran_reads (DataFrame): Transcript read data
        sru_range (int): Range for SRU calculation
        typeorf (str): Type of ORF
        batch_size (int): Number of ORFs to process at once within a chunk
        chunk_size (int): Number of ORFs to process in a single chunk before returning
        
    Returns:
        DataFrame: Processed ORF DataFrame with scores
    """
    log_info(f"Starting memory-efficient scoring for {len(df)} ORFs")
    
    # Process df in chunks to avoid holding all results in memory at once
    result_chunks = []
    total_rows = len(df)
    
    for chunk_start in range(0, total_rows, chunk_size):
        chunk_end = min(chunk_start + chunk_size, total_rows)
        chunk_df = df.slice(chunk_start, chunk_end - chunk_start)
        
        # Initialize result columns for this chunk only
        rise_up_values = [0.0] * (chunk_end - chunk_start)
        step_down_values = [0.0] * (chunk_end - chunk_start)
        hrf_values = [0.0] * (chunk_end - chunk_start)
        avg_values = [0.0] * (chunk_end - chunk_start)
        nzc_values = [0.0] * (chunk_end - chunk_start)
        
        # Process chunk in smaller batches
        for i in range(0, len(chunk_df), batch_size):
            end = min(i + batch_size, len(chunk_df))
            batch = chunk_df.slice(i, end - i)
            
            # Process each ORF in the batch
            for j, row in enumerate(batch.iter_rows(named=True)):
                idx = i + j
                start = row["start"]
                stop = row["stop"]
                
                # Calculate SRU scores based on ORF type
                if typeorf == "uoORF" or typeorf not in ("uoORF", "doORF"):
                    rise_up_values[idx] = sru_score(start, tran_reads, sru_range, 0)
                    
                if typeorf == "doORF" or typeorf not in ("uoORF", "doORF"):
                    step_down_values[idx] = sru_score(stop, tran_reads, sru_range, 1)
                
                # Calculate region scores
                hrf, avg, nzc = calculate_scores(start, stop, tran_reads)
                hrf_values[idx] = hrf
                avg_values[idx] = avg
                nzc_values[idx] = nzc
            
            # Free memory after each batch
            del batch
            gc.collect()
        
        # Add score columns to the chunk DataFrame
        chunk_result = chunk_df.with_columns([
            pl.Series("rise_up", rise_up_values),
            pl.Series("step_down", step_down_values),
            pl.Series("hrf", hrf_values),
            pl.Series("avg", avg_values),
            pl.Series("nzc", nzc_values)
        ])
        
        # Calculate total score
        chunk_result = chunk_result.with_columns(
            score=pl.sum_horizontal("rise_up", "step_down", "hrf", "avg", "nzc")
        )
        
        # Add this chunk to results and free memory
        result_chunks.append(chunk_result)
        del chunk_df, rise_up_values, step_down_values, hrf_values, avg_values, nzc_values
        gc.collect()
        
        log_info(f"Processed chunk {chunk_start}-{chunk_end} of {total_rows} rows")
    
    # Combine all chunks into final result
    if result_chunks:
        result_df = pl.concat(result_chunks)
        # Clean up
        del result_chunks
        gc.collect()
        return result_df
    else:
        return pl.DataFrame()


def capped_dict_update(dict_obj, key, value, max_size=100000):
    """
    Update a dictionary with a size limit to prevent unbounded memory growth.
    If dictionary exceeds max_size, oldest items are removed.
    
    Args:
        dict_obj (dict): Dictionary to update
        key: Key to add/update
        value: Value to set
        max_size (int): Maximum dictionary size
        
    Returns:
        dict: Updated dictionary
    """
    # Add new key-value pair
    dict_obj[key] = value
    
    # Check size and trim if needed
    if len(dict_obj) > max_size:
        # Remove oldest items (first 10% of keys)
        keys_to_remove = list(dict_obj.keys())[:max(1, int(max_size * 0.1))]
        for old_key in keys_to_remove:
            del dict_obj[old_key]
            
    return dict_obj


def cache_scores(df, tran_reads, sru_range, typeorf, scoredict, chunk_size=200):
    """
    Cache SRU scores in a dictionary for faster reuse.
    Memory-efficient implementation using chunking.

    Args:
        df (DataFrame/Series): Input data with ORF information
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for SRU score calculation
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)
        scoredict (dict): Dictionary to cache scores
        chunk_size (int): Size of chunks to process at once

    Returns:
        dict: Updated score dictionary
    """
    try:
        log_info(f"Starting chunked score caching for ORF type: {typeorf}")
        
        # Initialize dictionaries if not present
        if "rise_up" not in scoredict:
            scoredict["rise_up"] = {}
        if "step_down" not in scoredict:
            scoredict["step_down"] = {}
            
        # Process start positions if needed
        if not typeorf == "doORF" and not isinstance(df, pl.Series):
            # Get unique start positions and process in chunks
            startvalues = df.get_column("start").unique()
            
            for chunk_start in range(0, len(startvalues), chunk_size):
                chunk_end = min(chunk_start + chunk_size, len(startvalues))
                chunk = startvalues[chunk_start:chunk_end]
                
                # Process only positions not in cache
                for pos in chunk:
                    if pos not in scoredict["rise_up"]:
                        score = sru_score(pos, tran_reads, sru_range, 0)
                        scoredict["rise_up"] = capped_dict_update(scoredict["rise_up"], pos, score)
                
                # Clean up
                log_info(f"Processed start positions chunk {chunk_start}-{chunk_end} of {len(startvalues)}")
                del chunk
                gc.collect()
                
        # Process stop positions if needed
        if not typeorf == "uoORF" and not isinstance(df, pl.Series):
            # Get unique stop positions and process in chunks
            stopvalues = df.get_column("stop").unique()
            
            for chunk_start in range(0, len(stopvalues), chunk_size):
                chunk_end = min(chunk_start + chunk_size, len(stopvalues))
                chunk = stopvalues[chunk_start:chunk_end]
                
                # Process only positions not in cache
                for pos in chunk:
                    if pos not in scoredict["step_down"]:
                        score = sru_score(pos, tran_reads, sru_range, 1)
                        scoredict["step_down"] = capped_dict_update(scoredict["step_down"], pos, score)
                
                # Clean up
                log_info(f"Processed stop positions chunk {chunk_start}-{chunk_end} of {len(stopvalues)}")
                del chunk
                gc.collect()
                
        return scoredict
    except Exception as e:
        log_error(f"Error in memory-efficient score caching: {str(e)}")
        return {"rise_up": {}, "step_down": {}}


def calculate_global_scores(df, tran_reads, typeorf, chunk_size=200, batch_size=20):
    """
    Calculate global scores (HRF, average, NZC) for ORFs using chunking.
    Memory-efficient implementation.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)
        chunk_size (int): Size of chunks to process at once
        batch_size (int): Size of batches within each chunk

    Returns:
        DataFrame: DataFrame with calculated global scores
    """
    try:
        log_info(f"Starting chunked global scores for {len(df)} ORFs")
        
        total_rows = len(df)
        result_chunks = []
        
        for chunk_start in range(0, total_rows, chunk_size):
            chunk_end = min(chunk_start + chunk_size, total_rows)
            chunk = df.slice(chunk_start, chunk_end - chunk_start)
            
            # Initialize result arrays for this chunk only
            hrf_values = [0.0] * (chunk_end - chunk_start)
            avg_values = [0.0] * (chunk_end - chunk_start)
            nzc_values = [0.0] * (chunk_end - chunk_start)
            
            # Process in batches
            for i in range(0, len(chunk), batch_size):
                end = min(i + batch_size, len(chunk))
                batch = chunk.slice(i, end - i)
                
                # Process each ORF in the batch
                for j, row in enumerate(batch.iter_rows(named=True)):
                    idx = i + j
                    start = row["start"]
                    stop = row["stop"]
                    
                    # Calculate region scores
                    hrf, avg, nzc = calculate_scores(start, stop, tran_reads)
                    hrf_values[idx] = hrf
                    avg_values[idx] = avg
                    nzc_values[idx] = nzc
                
                # Clean up
                del batch
                gc.collect()
            
            # Add score columns to the chunk
            chunk_result = chunk.with_columns([
                pl.Series("hrf", hrf_values),
                pl.Series("avg", avg_values),
                pl.Series("nzc", nzc_values)
            ])
            
            # Calculate total score based on ORF type
            if typeorf == "doORF":
                chunk_result = chunk_result.with_columns(
                    score=pl.col("step_down") + pl.col("hrf") + pl.col("avg") + pl.col("nzc")
                )
            elif typeorf == "uoORF":
                chunk_result = chunk_result.with_columns(
                    score=pl.col("rise_up") + pl.col("hrf") + pl.col("avg") + pl.col("nzc")
                )
            else:
                chunk_result = chunk_result.with_columns(
                    score=pl.col("rise_up") + pl.col("step_down") + pl.col("hrf") + pl.col("avg") + pl.col("nzc")
                )
            
            # Add to results and free memory
            result_chunks.append(chunk_result)
            log_info(f"Processed global scores chunk {chunk_start}-{chunk_end} of {total_rows}")
            
            # Clean up
            del chunk, hrf_values, avg_values, nzc_values
            gc.collect()
        
        # Combine results
        if result_chunks:
            result_df = pl.concat(result_chunks)
            del result_chunks
            gc.collect()
            return result_df
        else:
            return pl.DataFrame()
    except Exception as e:
        log_error(f"Error in memory-efficient global scores: {str(e)}")
        return pl.DataFrame()


def filter_existing_scores(df, typeorf, scoredict):
    """
    Filter out ORFs that already have scores in the cache.
    Memory-efficient implementation.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)
        scoredict (dict): Dictionary of cached scores

    Returns:
        DataFrame/Series: Filtered data containing only ORFs needing scoring
    """
    try:
        # Initialize dictionaries if not present
        if "rise_up" not in scoredict:
            scoredict["rise_up"] = {}
        if "step_down" not in scoredict:
            scoredict["step_down"] = {}
            
        if typeorf == "uoORF":
            # Use expressions to avoid extra copies
            return df.filter(
                ~pl.col("start").is_in(list(scoredict["rise_up"].keys()))
            )["start"]
        elif typeorf == "doORF":
            return df.filter(
                ~pl.col("stop").is_in(list(scoredict["step_down"].keys()))
            )["stop"]
        else:
            # More efficiently filter both start and stop
            return df.filter(
                (~pl.col("start").is_in(list(scoredict["rise_up"].keys()))) |
                (~pl.col("stop").is_in(list(scoredict["step_down"].keys())))
            ).select(["start", "stop"])
    except Exception as e:
        log_error(f"Error checking existing scores: {str(e)}")
        return pl.DataFrame()


def assign_cached_scores(df, scoredict, typeorf, chunk_size=1000):
    """
    Assign cached scores to ORFs using chunking.
    Memory-efficient implementation.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        scoredict (dict): Dictionary of cached scores
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)
        chunk_size (int): Size of chunks to process at once

    Returns:
        DataFrame: DataFrame with scores assigned from cache
    """
    try:
        log_info(f"Starting chunked score assignment for {len(df)} ORFs")
        
        total_rows = len(df)
        result_chunks = []
        
        for chunk_start in range(0, total_rows, chunk_size):
            chunk_end = min(chunk_start + chunk_size, total_rows)
            chunk = df.slice(chunk_start, chunk_end - chunk_start)
            
            # Initialize result lists for this chunk only
            rise_up_values = [0.0] * (chunk_end - chunk_start)
            step_down_values = [0.0] * (chunk_end - chunk_start)
            
            # Process ORFs in chunk
            for j, row in enumerate(chunk.iter_rows(named=True)):
                if typeorf == "uoORF":
                    rise_up_values[j] = scoredict["rise_up"].get(row["start"], 0.0)
                    step_down_values[j] = 0.0
                elif typeorf == "doORF":
                    rise_up_values[j] = 0.0
                    step_down_values[j] = scoredict["step_down"].get(row["stop"], 0.0)
                else:
                    rise_up_values[j] = scoredict["rise_up"].get(row["start"], 0.0)
                    step_down_values[j] = scoredict["step_down"].get(row["stop"], 0.0)
            
            # Add score columns to chunk
            chunk_result = chunk.with_columns([
                pl.Series("rise_up", rise_up_values),
                pl.Series("step_down", step_down_values)
            ])
            
            # Add to results and free memory
            result_chunks.append(chunk_result)
            log_info(f"Assigned scores to chunk {chunk_start}-{chunk_end} of {total_rows}")
            
            # Clean up
            del chunk, rise_up_values, step_down_values
            gc.collect()
        
        # Combine results
        if result_chunks:
            result_df = pl.concat(result_chunks)
            del result_chunks
            gc.collect()
            return result_df
        else:
            return pl.DataFrame()
    except Exception as e:
        log_error(f"Error assigning scores: {str(e)}")
        return pl.DataFrame()


def classify_orf_positions(annotation, df, cds_df, chunk_size=1000):
    """
    Determines the relative position of ORFs to coding sequences (CDS) using chunking.

    Parameters:
        annotation (str): Path to the genome annotation file in BED/GFF/GTF format.
        df (polars.DataFrame): DataFrame containing ORF coordinates.
        cds_df (polars.DataFrame): DataFrame containing CDS coordinates.
        chunk_size (int): Size of chunks to process at once

    Returns:
        tuple: A tuple containing two polars DataFrames:
               - The first DataFrame contains ORF coordinates with an additional column
                 'type' indicating the relative position of each ORF to CDS.
               - The second DataFrame contains exon coordinates.
    """
    log_info(f"Starting chunked ORF relative position for {len(df)} ORFs")
    
    if "cdsdf" not in globals():
        cds_df, exon_coords = getexons_and_cds(annotation, list(df["tran_id"].unique()))

    tranids = cds_df["tran_id"].unique().to_list()

    # Join df with cds_df to include cds_start and cds_stop
    df = df.join(cds_df.select(["tran_id", "tran_start", "tran_stop"]), on="tran_id", how="left")

    # Process in chunks for memory efficiency
    total_rows = len(df)
    result_chunks = []
    
    for chunk_start in range(0, total_rows, chunk_size):
        chunk_end = min(chunk_start + chunk_size, total_rows)
        chunk = df.slice(chunk_start, chunk_end - chunk_start)
        
        # Initialize type values for this chunk only
        type_values = [""] * (chunk_end - chunk_start)
        
        for j, row in enumerate(chunk.iter_rows(named=True)):
            if row["tran_id"] in tranids:
                try:
                    type_values[j] = classify_orf({
                        "start": row["start"],
                        "stop": row["stop"],
                        "tran_start": row["tran_start"],
                        "tran_stop": row["tran_stop"]
                    })
                except:
                    type_values[j] = "Non Coding"
            else:
                type_values[j] = "Non Coding"
        
        # Add type column to chunk
        chunk_result = chunk.with_columns(pl.Series("type", type_values))
        
        # Add to results and free memory
        result_chunks.append(chunk_result)
        log_info(f"Processed ORF position chunk {chunk_start}-{chunk_end} of {total_rows}")
        
        # Clean up
        del chunk, type_values
        gc.collect()
    
    # Combine results
    if result_chunks:
        result_df = pl.concat(result_chunks)
        
        # Filter CDS
        cdslist = result_df.filter(pl.col("type") == "CDS")["tran_id"].unique().to_list()
        
        del result_chunks
        gc.collect()
        
        return result_df, exon_coords
    else:
        return pl.DataFrame(), exon_coords
"""
Core scoring functionality for TranslonScorer.

This module contains all functions related to scoring ORFs, including:
- Classic (old) scoring method
- Modern (new) scoring method
- Global scoring calculations
- Score assignment and management
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
        if direction == 0:
            # Calculate rise up score - more efficient filtering
            before = tran_reads.filter(
                (pl.col("tran_start") >= position - sru_range)
                & (pl.col("tran_start") < position)
            )
            after = tran_reads.filter(
                (pl.col("tran_start") >= position)
                & (pl.col("tran_start") < position + sru_range)
            )
        else:
            # Calculate step down score - more efficient filtering
            before = tran_reads.filter(
                (pl.col("tran_start") >= position)
                & (pl.col("tran_start") < position + sru_range)
            )
            after = tran_reads.filter(
                (pl.col("tran_start") >= position + sru_range)
                & (pl.col("tran_start") < position + (2 * sru_range))
            )

        # Use native Polars mean() for better performance
        before_mean = before["counts"].mean() or 0
        after_mean = after["counts"].mean() or 0
        
        # Clean up to reduce memory usage
        del before, after
        
        if direction == 0:
            return after_mean - before_mean
        else:
            return before_mean - after_mean
            
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
        # Single filter operation for better performance
        region = tran_reads.filter(
            (pl.col("tran_start") >= start) & (pl.col("tran_start") <= stop)
        )
        
        if region.is_empty():
            return (0.0, 0.0, 0.0)
        
        # Use native Polars operations where possible
        counts = region["counts"]
        
        # Highest Reading Frame - use max() directly
        hrf = counts.max() or 0
        
        # Average - use mean() directly
        avg = counts.mean() or 0
        
        # Non-zero Coverage - filter once and calculate
        nzc = (counts > 0).sum() / len(counts) if len(counts) > 0 else 0
        
        # Clean up to reduce memory usage
        del region, counts
        
        return (hrf, avg, nzc)
        
    except Exception as e:
        log_error(f"Error calculating region scores: {str(e)}")
        return (0.0, 0.0, 0.0)


def batch_process_orfs(df, tran_reads, sru_range, typeorf, batch_size=100):
    """
    Process ORFs in small batches to reduce memory usage.
    
    Args:
        df (DataFrame): ORF DataFrame to process
        tran_reads (DataFrame): Transcript read data
        sru_range (int): Range for SRU calculation
        typeorf (str): Type of ORF
        batch_size (int): Number of ORFs to process at once
        
    Returns:
        DataFrame: Processed ORF DataFrame with scores
    """
    # Initialize result columns
    rise_up_values = [0.0] * len(df)
    step_down_values = [0.0] * len(df)
    hrf_values = [0.0] * len(df)
    avg_values = [0.0] * len(df)
    nzc_values = [0.0] * len(df)
    
    # Process in batches
    for i in range(0, len(df), batch_size):
        end = min(i + batch_size, len(df))
        batch = df.slice(i, end - i)
        
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
        
        # Force garbage collection after each batch
        gc.collect()
    
    # Add score columns to the DataFrame
    result_df = df.with_columns([
        pl.Series("rise_up", rise_up_values),
        pl.Series("step_down", step_down_values),
        pl.Series("hrf", hrf_values),
        pl.Series("avg", avg_values),
        pl.Series("nzc", nzc_values)
    ])
    
    # Calculate total score
    result_df = result_df.with_columns(
        score=pl.sum_horizontal("rise_up", "step_down", "hrf", "avg", "nzc")
    )
    
    return result_df


def oldscoring(df, tran_reads, sru_range, typeorf):
    """
    Classic scoring method with memory-efficient implementation.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for SRU score calculation
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with calculated scores
    """
    try:
        # Process in small batches to avoid memory issues
        return batch_process_orfs(df, tran_reads, sru_range, typeorf)
    except Exception as e:
        log_error(f"Error in memory-efficient scoring: {str(e)}")
        return pl.DataFrame()


def newscoring(df, tran_reads, sru_range, typeorf, scoredict):
    """
    Modern scoring method that caches scores in a dictionary.
    Memory-efficient implementation.

    Args:
        df (DataFrame/Series): Input data with ORF information
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for SRU score calculation
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)
        scoredict (dict): Dictionary to cache scores

    Returns:
        dict: Updated score dictionary
    """
    try:
        batch_size = 50  # Process in small batches
        
        if not typeorf == "doORF":
            if not isinstance(df, pl.Series):
                # Get unique start positions
                startvalues = df.get_column("start").unique()
                
                # Process in batches
                for i in range(0, len(startvalues), batch_size):
                    batch_end = min(i + batch_size, len(startvalues))
                    batch = startvalues[i:batch_end]
                    
                    # Process only positions not in cache
                    new_positions = [pos for pos in batch if pos not in scoredict["rise_up"]]
                    
                    if new_positions:
                        for pos in new_positions:
                            scoredict["rise_up"][pos] = sru_score(pos, tran_reads, sru_range, 0)
                    
                    # Clean up
                    gc.collect()
            else:
                # Handle Series case
                startvalues = df.unique()
                
                # Process in batches
                for i in range(0, len(startvalues), batch_size):
                    batch_end = min(i + batch_size, len(startvalues))
                    batch = startvalues[i:batch_end]
                    
                    # Process only positions not in cache
                    new_positions = [pos for pos in batch if pos not in scoredict["rise_up"]]
                    
                    if new_positions:
                        for pos in new_positions:
                            scoredict["rise_up"][pos] = sru_score(pos, tran_reads, sru_range, 0)
                    
                    # Clean up
                    gc.collect()

        if not typeorf == "uoORF":
            if not isinstance(df, pl.Series):
                # Get unique stop positions
                stopvalues = df.get_column("stop").unique()
                
                # Process in batches
                for i in range(0, len(stopvalues), batch_size):
                    batch_end = min(i + batch_size, len(stopvalues))
                    batch = stopvalues[i:batch_end]
                    
                    # Process only positions not in cache
                    new_positions = [pos for pos in batch if pos not in scoredict["step_down"]]
                    
                    if new_positions:
                        for pos in new_positions:
                            scoredict["step_down"][pos] = sru_score(pos, tran_reads, sru_range, 1)
                    
                    # Clean up
                    gc.collect()
            else:
                # Handle Series case
                stopvalues = df.unique()
                
                # Process in batches
                for i in range(0, len(stopvalues), batch_size):
                    batch_end = min(i + batch_size, len(stopvalues))
                    batch = stopvalues[i:batch_end]
                    
                    # Process only positions not in cache
                    new_positions = [pos for pos in batch if pos not in scoredict["step_down"]]
                    
                    if new_positions:
                        for pos in new_positions:
                            scoredict["step_down"][pos] = sru_score(pos, tran_reads, sru_range, 1)
                    
                    # Clean up
                    gc.collect()
                
        return scoredict
    except Exception as e:
        log_error(f"Error in memory-efficient new scoring: {str(e)}")
        return {"rise_up": {}, "step_down": {}}


def globalscores(df, tran_reads, typeorf):
    """
    Calculate global scores (HRF, average, NZC) for ORFs.
    Memory-efficient implementation.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with calculated global scores
    """
    try:
        batch_size = 100  # Process in small batches
        total_rows = len(df)
        
        # Initialize result lists
        hrf_values = [0.0] * total_rows
        avg_values = [0.0] * total_rows
        nzc_values = [0.0] * total_rows
        
        # Process in batches
        for i in range(0, total_rows, batch_size):
            end = min(i + batch_size, total_rows)
            batch = df.slice(i, end - i)
            
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
            gc.collect()
        
        # Add score columns to the DataFrame
        result_df = df.with_columns([
            pl.Series("hrf", hrf_values),
            pl.Series("avg", avg_values),
            pl.Series("nzc", nzc_values)
        ])
        
        # Calculate total score based on ORF type
        if typeorf == "doORF":
            result_df = result_df.with_columns(
                score=pl.sum_horizontal("step_down", "hrf", "avg", "nzc")
            )
        elif typeorf == "uoORF":
            result_df = result_df.with_columns(
                score=pl.sum_horizontal("rise_up", "hrf", "avg", "nzc")
            )
        else:
            result_df = result_df.with_columns(
                score=pl.sum_horizontal("rise_up", "step_down", "hrf", "avg", "nzc")
            )
        
        return result_df
    except Exception as e:
        log_error(f"Error in memory-efficient global scores: {str(e)}")
        return pl.DataFrame()


def existingscore(df, typeorf, scoredict):
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
        if typeorf == "uoORF":
            df = (
                df["start"]
                .map_elements(lambda x: x if x not in scoredict["rise_up"] else None)
                .drop_nulls()
            )
            return df
        elif typeorf == "doORF":
            df = (
                df["stop"]
                .map_elements(lambda x: x if x not in scoredict["step_down"] else None)
                .drop_nulls()
            )
            return df
        else:
            df = df.select(["start", "stop"]).with_columns([
                pl.col("start")
                .apply(lambda x: x if x not in scoredict["rise_up"] else None)
                .alias("in_ru"),
                
                pl.col("stop")
                .apply(lambda x: x if x not in scoredict["step_down"] else None)
                .alias("in_sd")
            ])
            
            df = df.filter(
                (pl.col("in_ru").is_not_null()) | (pl.col("in_sd").is_not_null())
            ).select(["start", "stop"])
            
            return df
    except Exception as e:
        log_error(f"Error checking existing scores: {str(e)}")
        return pl.DataFrame()


def assigningscore(df, scoredict, typeorf):
    """
    Assign cached scores to ORFs.
    Memory-efficient implementation.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        scoredict (dict): Dictionary of cached scores
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with scores assigned from cache
    """
    try:
        batch_size = 100  # Process in small batches
        total_rows = len(df)
        
        # Initialize result lists
        rise_up_values = [0.0] * total_rows
        step_down_values = [0.0] * total_rows
        
        # Process in batches
        for i in range(0, total_rows, batch_size):
            end = min(i + batch_size, total_rows)
            batch = df.slice(i, end - i)
            
            # Process each ORF in the batch
            for j, row in enumerate(batch.iter_rows(named=True)):
                idx = i + j
                
                if typeorf == "uoORF":
                    rise_up_values[idx] = scoredict["rise_up"].get(row["start"], 0.0)
                    step_down_values[idx] = 0.0
                elif typeorf == "doORF":
                    rise_up_values[idx] = 0.0
                    step_down_values[idx] = scoredict["step_down"].get(row["stop"], 0.0)
                else:
                    rise_up_values[idx] = scoredict["rise_up"].get(row["start"], 0.0)
                    step_down_values[idx] = scoredict["step_down"].get(row["stop"], 0.0)
            
            # Clean up
            gc.collect()
        
        # Add score columns to the DataFrame
        result_df = df.with_columns([
            pl.Series("rise_up", rise_up_values),
            pl.Series("step_down", step_down_values)
        ])
        
        return result_df
    except Exception as e:
        log_error(f"Error assigning scores: {str(e)}")
        return pl.DataFrame()


def orfrelativeposition(annotation, df, cds_df):
    """
    Determines the relative position of ORFs to coding sequences (CDS).

    Parameters:
        annotation (str): Path to the genome annotation file in BED/GFF/GTF format.
        df (polars.DataFrame): DataFrame containing ORF coordinates.
        cds_df (polars.DataFrame): DataFrame containing CDS coordinates.

    Returns:
        tuple: A tuple containing two polars DataFrames:
               - The first DataFrame contains ORF coordinates with an additional column
                 'type' indicating the relative position of each ORF to CDS.
               - The second DataFrame contains exon coordinates.
    """
    if not "cdsdf" in globals():
        cds_df, exon_coords = getexons_and_cds(annotation, list(df["tran_id"].unique()))

    tranids = cds_df["tran_id"].unique().to_list()

    print(cds_df.head())
    print(df.head())

    # Join df with cds_df to include cds_start and cds_stop
    df = df.join(cds_df.select(["tran_id", "tran_start", "tran_stop"]), on="tran_id", how="left")

    # Process in batches for memory efficiency
    batch_size = 1000
    total_rows = len(df)
    type_values = [""] * total_rows
    
    for i in range(0, total_rows, batch_size):
        end = min(i + batch_size, total_rows)
        batch = df.slice(i, end - i)
        
        for j, row in enumerate(batch.iter_rows(named=True)):
            idx = i + j
            if row["tran_id"] in tranids:
                try:
                    type_values[idx] = classify_orf({
                        "start": row["start"],
                        "stop": row["stop"],
                        "tran_start": row["tran_start"],
                        "tran_stop": row["tran_stop"]
                    })
                except:
                    type_values[idx] = "Non Coding"
            else:
                type_values[idx] = "Non Coding"
        
        # Clean up
        gc.collect()
    
    # Add type column
    df = df.with_columns(pl.Series("type", type_values))

    # Filter CDS
    cdslist = df.filter(pl.col("type") == "CDS")["tran_id"].unique().to_list()

    return df, exon_coords
"""
Core scoring functionality for TranslonScorer (Optimized Version).

This module contains all functions related to scoring ORFs, including:
- Classic (old) scoring method
- Modern (new) scoring method
- Global scoring calculations
- Score assignment and management
"""

import polars as pl
import numpy as np
from ..utils.logging import log_info, log_warning, log_error
from .orffinder import classify_orf
from .orffinder import getexons_and_cds


def _precompute_tran_reads(tran_reads):
    """
    Precompute arrays for faster access during scoring.
    
    Args:
        tran_reads (DataFrame): Transcript reads data
        
    Returns:
        tuple: (positions array, counts array, total positions)
    """
    positions = tran_reads["tran_start"].to_numpy()
    counts = tran_reads["counts"].to_numpy()
    total_positions = len(positions)
    return positions, counts, total_positions


def sru_score(positions, tran_reads_data, sru_range, direction=0):
    """
    Vectorized version of SRU score calculation for multiple positions at once.
    
    Args:
        positions (list/array): Positions to calculate scores for
        tran_reads_data (tuple): Precomputed transcript reads data
        sru_range (int): Range for calculation
        direction (int): 0 for start (up), 1 for stop (down)
        
    Returns:
        list: Calculated SRU scores for each position
    """
    try:
        tran_positions, tran_counts, _ = tran_reads_data
        
        # Initialize results array
        results = np.zeros(len(positions))
        
        for i, position in enumerate(positions):
            if direction == 0:
                # Calculate rise up score (start)
                before_mask = (tran_positions >= position - sru_range) & (tran_positions < position)
                after_mask = (tran_positions >= position) & (tran_positions < position + sru_range)
            else:
                # Calculate step down score (stop)
                before_mask = (tran_positions >= position) & (tran_positions < position + sru_range)
                after_mask = (tran_positions >= position + sru_range) & (tran_positions < position + (2 * sru_range))
            
            before_vals = tran_counts[before_mask]
            after_vals = tran_counts[after_mask]
            
            before_mean = np.mean(before_vals) if len(before_vals) > 0 else 0
            after_mean = np.mean(after_vals) if len(after_vals) > 0 else 0
            
            if direction == 0:
                results[i] = after_mean - before_mean
            else:
                results[i] = before_mean - after_mean
        
        return results.tolist()
    except Exception as e:
        log_error(f"Error in vectorized SRU score calculation: {str(e)}")
        return [0.0] * len(positions)


def calculate_scores(starts, stops, tran_reads_data):
    """
    Vectorized version to calculate HRF, average, and non-zero coverage scores for multiple regions.
    
    Args:
        starts (list/array): Start positions
        stops (list/array): Stop positions
        tran_reads_data (tuple): Precomputed transcript reads data
        
    Returns:
        tuple: (HRF scores, average scores, non-zero coverage scores)
    """
    try:
        tran_positions, tran_counts, _ = tran_reads_data
        
        # Initialize result arrays
        hrf_scores = np.zeros(len(starts))
        avg_scores = np.zeros(len(starts))
        nzc_scores = np.zeros(len(starts))
        
        for i, (start, stop) in enumerate(zip(starts, stops)):
            # Get counts for this region
            region_mask = (tran_positions >= start) & (tran_positions <= stop)
            region_counts = tran_counts[region_mask]
            
            if len(region_counts) > 0:
                # Highest Reading Frame
                hrf_scores[i] = np.max(region_counts)
                
                # Average
                avg_scores[i] = np.mean(region_counts)
                
                # Non-zero Coverage
                nzc_scores[i] = np.sum(region_counts > 0) / len(region_counts)
        
        return hrf_scores.tolist(), avg_scores.tolist(), nzc_scores.tolist()
    except Exception as e:
        log_error(f"Error in vectorized region score calculation: {str(e)}")
        return ([0.0] * len(starts), [0.0] * len(starts), [0.0] * len(starts))


def oldscoring(df, tran_reads, sru_range, typeorf):
    """
    Optimized classic scoring method using vectorized operations.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for SRU score calculation
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with calculated scores
    """
    try:
        # Precompute transcript reads data for faster access
        tran_reads_data = _precompute_tran_reads(tran_reads)
        
        # Extract start and stop positions
        starts = df["start"].to_list()
        stops = df["stop"].to_list()
        
        # Initialize scores
        rise_up_scores = [0.0] * len(df)
        step_down_scores = [0.0] * len(df)
        
        # Calculate SRU scores based on ORF type
        if typeorf == "uoORF" or typeorf not in ("uoORF", "doORF"):
            rise_up_scores = sru_score(starts, tran_reads_data, sru_range, 0)
            
        if typeorf == "doORF" or typeorf not in ("uoORF", "doORF"):
            step_down_scores = sru_score(stops, tran_reads_data, sru_range, 1)
        
        # Calculate region scores
        hrf_scores, avg_scores, nzc_scores = calculate_scores_vectorized(
            starts, stops, tran_reads_data
        )
        
        # Create result DataFrame
        result_df = df.with_columns([
            pl.Series("rise_up", rise_up_scores),
            pl.Series("step_down", step_down_scores),
            pl.Series("hrf", hrf_scores),
            pl.Series("avg", avg_scores),
            pl.Series("nzc", nzc_scores)
        ])
        
        # Calculate total score
        result_df = result_df.with_columns(
            score=pl.sum_horizontal("rise_up", "step_down", "hrf", "avg", "nzc")
        )
        
        return result_df
    except Exception as e:
        log_error(f"Error in optimized old scoring method: {str(e)}")
        return pl.DataFrame()


def compute_score_dict(unique_positions, tran_reads_data, sru_range, direction):
    """
    Efficiently compute scores for unique positions and create a dictionary.
    
    Args:
        unique_positions (list): Unique positions to calculate scores for
        tran_reads_data (tuple): Precomputed transcript reads data
        sru_range (int): Range for calculation
        direction (int): 0 for start (up), 1 for stop (down)
        
    Returns:
        dict: Dictionary mapping positions to their scores
    """
    scores = sru_score(unique_positions, tran_reads_data, sru_range, direction)
    # Create dictionary in one go instead of repeated updates
    return dict(zip(unique_positions, scores))


def newscoring(df, tran_reads, sru_range, typeorf, scoredict):
    """
    Optimized modern scoring method that efficiently caches scores in a dictionary.

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
        # Precompute transcript reads data
        tran_reads_data = _precompute_tran_reads(tran_reads)
        
        # Update rise_up scores if needed
        if typeorf != "doORF" and not isinstance(df, pl.Series):
            # Get unique start positions not already in the cache
            start_values = df["start"].unique().to_list()
            new_starts = [s for s in start_values if s not in scoredict["rise_up"]]
            
            if new_starts:
                # Calculate scores for new positions only
                rise_up_dict = compute_score_dict(new_starts, tran_reads_data, sru_range, 0)
                # Update the cache in one batch operation
                scoredict["rise_up"].update(rise_up_dict)
                
        elif typeorf != "doORF" and isinstance(df, pl.Series):
            # Handle Series case
            start_values = df.to_list()
            new_starts = [s for s in start_values if s not in scoredict["rise_up"]]
            
            if new_starts:
                rise_up_dict = compute_score_dict(new_starts, tran_reads_data, sru_range, 0)
                scoredict["rise_up"].update(rise_up_dict)
        
        # Update step_down scores if needed
        if typeorf != "uoORF" and not isinstance(df, pl.Series):
            # Get unique stop positions not already in the cache
            stop_values = df["stop"].unique().to_list()
            new_stops = [s for s in stop_values if s not in scoredict["step_down"]]
            
            if new_stops:
                # Calculate scores for new positions only
                step_down_dict = compute_score_dict(new_stops, tran_reads_data, sru_range, 1)
                # Update the cache in one batch operation
                scoredict["step_down"].update(step_down_dict)
                
        elif typeorf != "uoORF" and isinstance(df, pl.Series):
            # Handle Series case
            stop_values = df.to_list()
            new_stops = [s for s in stop_values if s not in scoredict["step_down"]]
            
            if new_stops:
                step_down_dict = compute_score_dict(new_stops, tran_reads_data, sru_range, 1)
                scoredict["step_down"].update(step_down_dict)
        
        return scoredict
    except Exception as e:
        log_error(f"Error in optimized new scoring method: {str(e)}")
        return {"rise_up": {}, "step_down": {}}


def globalscores(df, tran_reads, typeorf):
    """
    Calculate global scores (HRF, average, NZC) for ORFs using optimized vectorized operations.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with calculated global scores
    """
    try:
        # Precompute transcript reads data
        tran_reads_data = _precompute_tran_reads(tran_reads)
        
        # Extract start and stop positions
        starts = df["start"].to_list()
        stops = df["stop"].to_list()
        
        # Calculate region scores in one vectorized operation
        hrf_scores, avg_scores, nzc_scores = calculate_scores(
            starts, stops, tran_reads_data
        )
        
        # Add columns to DataFrame
        result_df = df.with_columns([
            pl.Series("hrf", hrf_scores),
            pl.Series("avg", avg_scores),
            pl.Series("nzc", nzc_scores)
        ])
        
        # Calculate total score
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
        log_error(f"Error in optimized global scores calculation: {str(e)}")
        return pl.DataFrame()


def existingscore(df, typeorf, scoredict):
    """
    Optimized version to filter out ORFs that already have scores in the cache.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)
        scoredict (dict): Dictionary of cached scores

    Returns:
        DataFrame/Series: Filtered data containing only ORFs needing scoring
    """
    try:
        # Convert dict keys to sets for faster membership testing
        rise_up_keys = set(scoredict["rise_up"].keys())
        step_down_keys = set(scoredict["step_down"].keys())
        
        if typeorf == "uoORF":
            # Filter in one vectorized operation
            return df["start"].filter(~pl.col("start").is_in(rise_up_keys))
            
        elif typeorf == "doORF":
            # Filter in one vectorized operation
            return df["stop"].filter(~pl.col("stop").is_in(step_down_keys))
            
        else:
            # For mixed types, check both start and stop
            filtered_df = df.select(["start", "stop"]).filter(
                (~pl.col("start").is_in(rise_up_keys)) | 
                (~pl.col("stop").is_in(step_down_keys))
            )
            return filtered_df
            
    except Exception as e:
        log_error(f"Error in optimized existing score check: {str(e)}")
        return pl.DataFrame()


def assigningscore(df, scoredict, typeorf):
    """
    Optimized version to assign cached scores to ORFs using vectorized operations.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        scoredict (dict): Dictionary of cached scores
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with scores assigned from cache
    """
    try:
        # Convert dictionaries to Series for faster lookups
        rise_up_series = pl.Series(list(scoredict["rise_up"].keys()), 
                                  list(scoredict["rise_up"].values()))
        step_down_series = pl.Series(list(scoredict["step_down"].keys()), 
                                    list(scoredict["step_down"].values()))
        
        # Create lookup expressions
        if typeorf == "uoORF":
            result_df = df.with_columns([
                # Use map_dict for faster lookups
                pl.col("start").map_dict(scoredict["rise_up"], default=0.0).alias("rise_up"),
                pl.lit(0.0).alias("step_down")
            ])
            
        elif typeorf == "doORF":
            result_df = df.with_columns([
                pl.lit(0.0).alias("rise_up"),
                pl.col("stop").map_dict(scoredict["step_down"], default=0.0).alias("step_down")
            ])
            
        else:
            result_df = df.with_columns([
                pl.col("start").map_dict(scoredict["rise_up"], default=0.0).alias("rise_up"),
                pl.col("stop").map_dict(scoredict["step_down"], default=0.0).alias("step_down")
            ])
            
        return result_df
        
    except Exception as e:
        log_error(f"Error in optimized score assignment: {str(e)}")
        return pl.DataFrame()


def orfrelativeposition(annotation, df, cds_df):
    """
    Determines the relative position of ORFs to coding sequences (CDS).
    
    This function takes a genome annotation file and a DataFrame containing ORF coordinates,
    and determines the relative position of each ORF with respect to coding sequences (CDS).
    It classifies each ORF into different categories based on its relationship with CDS.

    Parameters:
        annotation (str): Path to the genome annotation file in BED/GFF/GTF format.
        df (polars.DataFrame): DataFrame containing ORF coordinates. It must have columns
                               'tran_id', 'pos', and 'end' representing transcript ID, start
                               position, and end position of each ORF respectively.

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

    # Vectorized operation to classify ORFs
    df = df.with_columns(
        pl.when(pl.col("tran_id").is_in(tranids))
        .then(
            pl.struct(["start", "stop", "tran_start", "tran_stop"])
            .apply(lambda row: classify_orf(row))
        )
        .otherwise(pl.lit("Non Coding"))
        .alias("type")
    )

    # Filter CDS
    cdslist = df.filter(pl.col("type") == "CDS")["tran_id"].unique().to_list()

    return df, exon_coords
'''Script to score the ORF's'''
import polars as pl
from numpy import log as ln
import numpy as np

def sru_score(start, bigwig_df, rng, invert):
    '''
    docstring
    '''
    positions = pl.Series("positions", np.arange(start - rng, start + 3 + rng))
    # Filter positions to match the correct frame
    positions = positions.filter(positions % 3 == start % 3)
    # Split positions into before and after the start
    before_positions = positions.filter(positions < start)
    after_positions = positions.filter(positions > start)
    # Filter the bigwig_df based on these positions
    before_counts = bigwig_df.filter(pl.col("tran_start").is_in(before_positions))["counts"].sum()
    after_counts = bigwig_df.filter(pl.col("tran_start").is_in(after_positions))["counts"].sum()
    numerator = 1 + after_counts
    denominator = 1 + before_counts

    # Check for Start rise up or step down score
    sru = ln(numerator / denominator) if invert == 0 else ln(denominator / numerator)

    return float(sru)


def calculate_scores(start, stop, bigwig_df):
    '''
    Calculate HRF, average, and NZC scores.
    
    Parameters:
        start (int): Start position.
        stop (int): Stop position.
        bigwig_df (polars DataFrame): DataFrame containing relevant data.
        
    Returns:
        dict: A dictionary containing HRF, average, and NZC scores.
    '''
    # Filter the dataframe for the range of interest
    relevant_df = bigwig_df.filter((pl.col("tran_start") >= start) & (pl.col("tran_start") <= stop))
    # Create a new column for frame based on tran_start
    relevant_df = relevant_df.with_columns((pl.col("tran_start") % 3).alias("frame"))
    codons_df = relevant_df.with_columns(
        pl.when(pl.col("frame") == 0).then(pl.col("counts")).otherwise(0).alias("frame0"),
        pl.when(pl.col("frame") == 1).then(pl.col("counts")).otherwise(0).alias("frame1"),
        pl.when(pl.col("frame") == 2).then(pl.col("counts")).otherwise(0).alias("frame2")
    )
    relevant_df = relevant_df.group_by("frame").agg(pl.sum("counts").alias("frame_counts")).to_dict(as_series=False)
    # Extract frame counts
    framecounts = {frame: relevant_df["frame_counts"][i] if i < len(relevant_df["frame_counts"]) else 0 for i, frame in enumerate([0, 1, 2])}
    # Calculate codons_f0 and total_codons
    codons_f0 = codons_df.filter((pl.col("frame0") > pl.col("frame1")) & (pl.col("frame0") > pl.col("frame2"))).height
    total_codons = codons_df.filter((pl.col("frame0") > 0) | (pl.col("frame1") > 0) | (pl.col("frame2") > 0)).height
    # Calculate hrf, avg, and nzc
    hrf = framecounts[0] / max(framecounts[1], framecounts[2]) if framecounts[1] and framecounts[2] != 0 else 0
    avg = framecounts[0] / ((stop-start) / 3) if framecounts[0] else 0
    nzc = codons_f0 / total_codons if total_codons > 0 else 0

    return float(hrf), float(avg), float(nzc)


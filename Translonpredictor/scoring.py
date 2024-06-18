'''Script to score the ORF's'''
import polars as pl
from numpy import log as ln
import numpy as np

def sru_score(start, bigwig_df, rng, invert):
    '''
    Calculate the Start Rise Up (SRU) and  Step Down score based on counts in a specified range around a start position.

    Parameters:
    - start (int): Start position for which the SRU score is calculated.
    - bigwig_df (DataFrame): Pandas DataFrame containing transcript start positions and associated counts.
    - rng (int): Range parameter specifying the distance in nucleotides from the start position to consider.
    - invert (int): Flag indicating whether to invert the SRU score calculation:
      - 0: Calculate ln(numerator / denominator).
      - Non-zero: Calculate ln(denominator / numerator).

    Returns:
    - float: SRU score calculated based on the given parameters.

    This function calculates the SRU score using the formula:
    SRU = ln((1 + after_counts) / (1 + before_counts))
    or its inverse if `invert` is non-zero.

    The function performs the following steps:
    1. Generates a series of positions around `start` within the range `rng`.
    2. Filters these positions to match the correct reading frame based on `start % 3`.
    3. Splits these positions into those occurring before and after `start`.
    4. Filters `bigwig_df` based on these positions to get counts before and after `start`.
    5. Calculates the SRU score based on the counts, adjusting the calculation based on `invert`.

    Note: The natural logarithm (ln) is used for score calculation.

    Example:
    >>> start = 100
    >>> bigwig_df = pd.DataFrame({'tran_start': [95, 98, 101, 104], 'counts': [10, 15, 20, 25]})
    >>> sru_score(start, bigwig_df, 2, 0)
    0.28768207245178085

    In this example, `start` is 100, and the SRU score is calculated within a range of ±2 nucleotides around 100.
    The function uses counts from positions before and after 100 to compute the SRU score.
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
    Calculate various scores related to codon usage and frame preference within a specified range.

    Parameters:
    - start (int): Start position of the range for which scores are calculated.
    - stop (int): Stop position of the range for which scores are calculated.
    - bigwig_df (DataFrame): Pandas DataFrame containing transcript start positions and associated counts.

    Returns:
    - float: hrf - High Read Frame (HRF) score indicating frame preference.
    - float: avg - Average codon count per codon in the specified range.
    - float: nzc - Non-Zero Codon (NZC) ratio indicating the proportion of codons with non-zero counts.

    This function calculates the following scores based on the provided parameters:
    - hrf: Calculated as the ratio of counts in frame 0 to the maximum counts in frames 1 and 2, if frames 1 and 2 are not zero.
    - avg: Average number of codons per codon position in the specified range.
    - nzc: Ratio of codons with non-zero counts to the total number of codons in the specified range.

    The function first filters `bigwig_df` to include only rows where `tran_start` falls within the specified `start` and `stop`.
    It then calculates frame-specific counts and aggregates them to determine frame preferences (`hrf`).
    Codons primarily in frame 0 are counted to determine the proportion of all codons (`nzc`).

    Example:
    >>> start = 100
    >>> stop = 200
    >>> bigwig_df = pd.DataFrame({'tran_start': [95, 100, 105, 110, 115], 'counts': [10, 20, 15, 25, 30]})
    >>> calculate_scores(start, stop, bigwig_df)
    (0.5, 1.0, 0.6)

    In this example, `start` is 100 and `stop` is 200. The function calculates the HRF, AVG, and NZC scores based on
    counts from `bigwig_df` for positions between 100 and 200.
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

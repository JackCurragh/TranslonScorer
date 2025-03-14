"""This script contains functions to read BigWig files and calculate transcriptomic coordinates."""

import polars as pl
import pyBigWig as bw
from .scoring import sru_score, calculate_scores
from .logging_config import log_info, log_warning, log_error


def transcriptreads(bwfile, exon_df):
    """
    Converts a BigWig file to a DataFrame based on provided exon annotation.

    Parameters:
    - bwfile (pyBigWig): An open BigWig file handle.
    - exon_df (DataFrame): DataFrame containing exon annotations.

    Returns:
    - df_tran (DataFrame): DataFrame containing transcript information derived from the BigWig file.

    Raises:
    - RuntimeError: If there are issues reading values from the BigWig file.
    - ValueError: If no reads could be extracted from any chromosome.
    """
    reads = []
    # Explode the lists into rows and sort by chromosome and start position
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
        
        try:
            # Get values for all positions in this chromosome
            for start, stop in zip(starts, stops):
                try:
                    values = bwfile.values(chrom, start, stop)
                    if values:
                        reads.extend(values)
                except RuntimeError as e:
                    log_error(f"Could not read values for {chrom}:{start}-{stop}: {str(e)}")
        except Exception as e:
            log_error(f"Error processing chromosome {chrom}: {str(e)}")
    
    if not reads:
        log_error("No reads could be extracted from any chromosome", exception_type=ValueError)
        
    return pl.DataFrame({
        "tran_start": range(len(reads)),
        "counts": reads
    })


def oldscoring(df, tran_reads, sru_range, typeorf):
    """
    Optimized classic scoring method using vectorized operations.

    Args:
        df (DataFrame): Input DataFrame with ORF information.
        tran_reads (DataFrame): Transcript reads data.
        sru_range (int): Range for SRU score calculation.
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other).

    Returns:
        DataFrame: DataFrame with calculated scores.
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
        
        # Create result DataFrame with calculated scores
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
        unique_positions (list): Unique positions to calculate scores for.
        tran_reads_data (tuple): Precomputed transcript reads data.
        sru_range (int): Range for calculation.
        direction (int): 0 for start (up), 1 for stop (down).

    Returns:
        dict: Dictionary mapping positions to their scores.
    """
    scores = sru_score(unique_positions, tran_reads_data, sru_range, direction)
    # Create dictionary in one go instead of repeated updates
    return dict(zip(unique_positions, scores))


def newscoring(df, tran_reads, sru_range, typeorf, scoredict):
    """
    Optimized modern scoring method that efficiently caches scores in a dictionary.

    Args:
        df (DataFrame/Series): Input data with ORF information.
        tran_reads (DataFrame): Transcript reads data.
        sru_range (int): Range for SRU score calculation.
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other).
        scoredict (dict): Dictionary to cache scores.

    Returns:
        dict: Updated score dictionary.
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
        df (DataFrame): Input DataFrame with ORF information.
        tran_reads (DataFrame): Transcript reads data.
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other).

    Returns:
        DataFrame: DataFrame with calculated global scores.
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
        log_error(f"Error in optimized global scores calculation: {str(e)}")
        return pl.DataFrame()


def existingscore(df, typeorf, scoredict):
    """
    Optimized version to filter out ORFs that already have scores in the cache.

    Args:
        df (DataFrame): Input DataFrame with ORF information.
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other).
        scoredict (dict): Dictionary of cached scores.

    Returns:
        DataFrame/Series: Filtered data containing only ORFs needing scoring.
    """
    try:
        # Convert dict keys to sets for faster membership testing
        rise_up_keys = set(scoredict["rise_up"].keys())
        step_down_keys = set(scoredict["step_down"].keys())
        
        if typeorf == "uoORF":
            # Filter in one vectorized operation
            return df.filter(~pl.col("start").is_in(rise_up_keys))
            
        elif typeorf == "doORF":
            # Filter in one vectorized operation
            return df.filter(~pl.col("stop").is_in(step_down_keys))
            
        else:
            # For mixed types, check both start and stop
            filtered_df = df.filter(
                (~pl.col("start").is_in(rise_up_keys)) | 
                (~pl.col("stop").is_in(step_down_keys))
            )
            return filtered_df
            
    except Exception as e:
        log_error(f"Error in optimized existing score check: {str(e)}")
        return pl.DataFrame()


def assigningscore(df, scoredict, typeorf):
    """
    Assigns scores from a dictionary to a DataFrame based on the type of ORF.

    This function updates the DataFrame `df` by assigning scores from the `scoredict` to the
    'rise_up' and 'step_down' columns based on the 'start' and 'stop' values. The type of ORF (`typeorf`)
    determines which scores are assigned.

    Parameters:
    df (pl.DataFrame): The input DataFrame containing 'start' and 'stop' columns.
    scoredict (dict): Dictionary containing the scores for 'rise_up' and 'step_down'.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.

    Returns:
    pl.DataFrame: The modified DataFrame with assigned scores.

    Notes:
    - For 'uoORF', assigns 'rise_up' scores from `scoredict` based on 'start' values and sets 'step_down' to 0.0.
    - For 'doORF', assigns 'step_down' scores from `scoredict` based on 'stop' values and sets 'rise_up' to 0.0.
    - For other types, assigns both 'rise_up' and 'step_down' scores from `scoredict` based on 'start' and 'stop' values.
    """
    try:
        if typeorf == "uoORF":
            df = df.with_columns(
                (pl.col("start").apply(lambda x: scoredict["rise_up"].get(x, 0.0)).alias("rise_up")),
                (pl.lit(0.0).alias("step_down")),
            )

        elif typeorf == "doORF":
            df = df.with_columns(
                (
                    pl.col("stop")
                    .apply(lambda x: scoredict["step_down"].get(x, 0.0))
                    .alias("step_down")
                ),
                (pl.lit(0.0).alias("rise_up")),
            )
        else:
            df = df.with_columns(
                (pl.col("start").apply(lambda x: scoredict["rise_up"].get(x, 0.0)).alias("rise_up")),
                (
                    pl.col("stop")
                    .apply(lambda x: scoredict["step_down"].get(x, 0.0))
                    .alias("step_down")
                ),
            )
        return df
    except Exception as e:
        log_error(f"Error assigning scores: {str(e)}")
        return pl.DataFrame()


def scoring(bigwig, exon, orfs, old_scoring, sru_range):
    """
    Score ORFs using bigwig coverage data.
    
    Args:
        bigwig (str): Path to bigwig file
        exon (str): Path to exon file
        orfs (str): Path to ORFs file
        old_scoring (bool): Whether to use old scoring method
        sru_range (int): Range for SRU score calculation
        
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
    exon_df = pl.read_csv(exon, has_header=True, separator=",")
    # Split and convert to integers
    exon_df = exon_df.with_columns([
        pl.col("start").apply(lambda x: [int(i) for i in x.split(",")]),
        pl.col("stop").apply(lambda x: [int(i) for i in x.split(",")]),
        pl.col("tran_start").apply(lambda x: [int(i) for i in x.split(",")]),
        pl.col("tran_stop").apply(lambda x: [int(i) for i in x.split(",")])
    ])
    orf_df = pl.read_csv(orfs, has_header=True, separator=",")

    counter = 0
    orfscores = []
    total_transcripts = len(orf_df["tran_id"].unique())
    log_info(f"Processing {total_transcripts} transcripts")

    for tran in orf_df["tran_id"].unique():
        if counter % 1000 == 0:
            log_info(f"Processed {counter}/{total_transcripts} transcripts")

        exons = exon_df.filter(pl.col("tran_id") == tran)
        orfs = orf_df.filter(pl.col("tran_id") == tran)

        if exons.is_empty():
            log_warning(f"No exon data found for transcript {tran}")
            continue

        tran_reads = transcriptreads(bwfile, exons)
        if not tran_reads.is_empty():
            for typeorf in orfs["type"].unique():
                orfs_filtered = orfs.filter(pl.col("type") == typeorf)

                if old_scoring:
                    orfs_filtered = oldscoring(
                        orfs_filtered, tran_reads, sru_range, typeorf
                    )
                    orfscores.append(orfs_filtered)
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
                        orfscores.append(orfs_filtered)
        counter += 1

    if not orfscores:
        log_warning("No ORFs were scored")
        return pl.DataFrame()

    log_info("Combining scored ORFs")
    final_df = pl.concat(orfscores)
    log_info(f"Scored {len(final_df)} ORFs in total")
    
    return final_df

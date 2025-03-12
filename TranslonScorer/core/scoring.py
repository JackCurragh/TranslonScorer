"""
Core scoring functionality for TranslonScorer.

This module contains all functions related to scoring ORFs, including:
- Classic (old) scoring method
- Modern (new) scoring method
- Global scoring calculations
- Score assignment and management
"""

import polars as pl
from ..utils.logging import log_info, log_warning, log_error
from .orf_classification import classify_orf
from .orf_classification import getexons_and_cds



def sru_score(position, tran_reads, sru_range, direction):
    """
    Calculate the SRU (Start/Stop Ramp Up/Down) score for a given position.
    
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
            # Calculate rise up score
            before = tran_reads.filter(
                (pl.col("tran_start") >= position - sru_range)
                & (pl.col("tran_start") < position)
            )
            after = tran_reads.filter(
                (pl.col("tran_start") >= position)
                & (pl.col("tran_start") < position + sru_range)
            )
        else:
            # Calculate step down score
            before = tran_reads.filter(
                (pl.col("tran_start") >= position)
                & (pl.col("tran_start") < position + sru_range)
            )
            after = tran_reads.filter(
                (pl.col("tran_start") >= position + sru_range)
                & (pl.col("tran_start") < position + (2 * sru_range))
            )

        before_mean = before["counts"].mean() or 0
        after_mean = after["counts"].mean() or 0
        
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
    
    Args:
        start (int): Start position
        stop (int): Stop position
        tran_reads (DataFrame): Transcript reads data
        
    Returns:
        tuple: (HRF score, average score, non-zero coverage score)
    """
    try:
        region = tran_reads.filter(
            (pl.col("tran_start") >= start) & (pl.col("tran_start") <= stop)
        )
        
        if region.is_empty():
            return (0.0, 0.0, 0.0)
            
        counts = region["counts"].to_list()
        
        # Highest Reading Frame
        hrf = max(counts) if counts else 0
        
        # Average
        avg = sum(counts) / len(counts) if counts else 0
        
        # Non-zero Coverage
        nzc = len([x for x in counts if x > 0]) / len(counts) if counts else 0
        
        return (hrf, avg, nzc)
        
    except Exception as e:
        log_error(f"Error calculating region scores: {str(e)}")
        return (0.0, 0.0, 0.0)


def oldscoring(df, tran_reads, sru_range, typeorf):
    """
    Classic scoring method that calculates all scores at once.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for SRU score calculation
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        dict: DataFrame as dictionary with calculated scores
    """
    try:
        # Initialize rise_up and step_down columns with default values
        df = df.with_columns([
            pl.lit(0.0).alias("rise_up"),
            pl.lit(0.0).alias("step_down")
        ])

        if typeorf == "uoORF":
            df = df.with_columns(
                pl.struct(["start"])
                .apply(lambda x: sru_score(x["start"], tran_reads, sru_range, 0))
                .alias("rise_up")
            )
        elif typeorf == "doORF":
            df = df.with_columns(
                pl.struct(["stop"])
                .apply(lambda x: sru_score(x["stop"], tran_reads, sru_range, 1))
                .alias("step_down")
            )
        else:
            # For other types, calculate both rise_up and step_down
            df = df.with_columns([
                pl.struct(["start"])
                .apply(lambda x: sru_score(x["start"], tran_reads, sru_range, 0))
                .alias("rise_up"),
                pl.struct(["stop"])
                .apply(lambda x: sru_score(x["stop"], tran_reads, sru_range, 1))
                .alias("step_down")
            ])

        # Calculate other scores
        df = df.with_columns(
            (
                pl.struct(["start", "stop"])
                .apply(lambda x: calculate_scores(x["start"], x["stop"], tran_reads))
                .alias("list_scores")
            )
        )
        df = df.with_columns(
            (pl.col("list_scores").apply(lambda x: x[0]).alias("hrf")),
            (pl.col("list_scores").apply(lambda x: x[1]).alias("avg")),
            (pl.col("list_scores").apply(lambda x: x[2]).alias("nzc")),
        )

        df = (
            df.with_columns(
                score=pl.sum_horizontal("rise_up", "step_down", "hrf", "avg", "nzc")
            )
            .select(pl.all().exclude("list_scores"))
            .to_dict(as_series=False)
        )
        return df
    except Exception as e:
        log_error(f"Error in old scoring method: {str(e)}")
        return pl.DataFrame()


def newscoring(df, tran_reads, sru_range, typeorf, scoredict):
    """
    Modern scoring method that caches scores in a dictionary.

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
        if not typeorf == "doORF":
            if type(df) != type(pl.Series()):
                startvalues = df.get_column("start").unique()
            else:
                startvalues = df
            startscore = startvalues.map_elements(
                lambda x: sru_score(x, tran_reads, sru_range, 0)
            ).to_list()
            for i, value in enumerate(startscore):
                scoredict["rise_up"].update({startvalues[i]: startscore[i]})

        if not typeorf == "uoORF":
            if type(df) != type(pl.Series()):
                stopvalues = df.get_column("stop").unique()
            else:
                stopvalues = df
            stopscore = stopvalues.map_elements(
                lambda x: sru_score(x, tran_reads, sru_range, 1)
            ).to_list()
            for i, value in enumerate(stopscore):
                scoredict["step_down"].update({stopvalues[i]: stopscore[i]})
        return scoredict
    except Exception as e:
        log_error(f"Error in new scoring method: {str(e)}")
        return {"rise_up": {}, "step_down": {}}


def globalscores(df, tran_reads, typeorf):
    """
    Calculate global scores (HRF, average, NZC) for ORFs.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        dict: DataFrame as dictionary with calculated global scores
    """
    try:
        df = df.with_columns(
            (
                pl.struct(["start", "stop"])
                .apply(lambda x: calculate_scores(x["start"], x["stop"], tran_reads))
                .alias("list_scores")
            )
        )
        df = df.with_columns(
            (pl.col("list_scores").apply(lambda x: x[0]).alias("hrf")),
            (pl.col("list_scores").apply(lambda x: x[1]).alias("avg")),
            (pl.col("list_scores").apply(lambda x: x[2]).alias("nzc")),
        ).select(pl.all().exclude("list_scores"))

        if typeorf == "doORF":
            df = df.with_columns(
                score=pl.sum_horizontal("step_down", "hrf", "avg", "nzc")
            ).to_dict(as_series=False)
        elif typeorf == "uoORF":
            df = df.with_columns(
                score=pl.sum_horizontal("rise_up", "hrf", "avg", "nzc")
            ).to_dict(as_series=False)
        else:
            df = df.with_columns(
                score=pl.sum_horizontal("rise_up", "step_down", "hrf", "avg", "nzc")
            ).to_dict(as_series=False)
        return df
    except Exception as e:
        log_error(f"Error calculating global scores: {str(e)}")
        return pl.DataFrame()


def existingscore(df, typeorf, scoredict):
    """
    Filter out ORFs that already have scores in the cache.

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
                .map_elements(lambda x: x if not x in scoredict["rise_up"] else None)
                .drop_nulls()
            )
            return df
        elif typeorf == "doORF":
            df = (
                df["stop"]
                .map_elements(lambda x: x if not x in scoredict["step_down"] else None)
                .drop_nulls()
            )
            return df
        else:
            df = df.select(["start", "stop"]).with_columns(
                (
                    pl.col("start")
                    .apply(lambda x: x if not x in scoredict["rise_up"] else None)
                    .alias("in_ru")
                ),
                (
                    pl.col("stop")
                    .apply(lambda x: x if not x in scoredict["step_down"] else None)
                    .alias("in_sd")
                ),
            )
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

    Args:
        df (DataFrame): Input DataFrame with ORF information
        scoredict (dict): Dictionary of cached scores
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with scores assigned from cache
    """
    try:
        if typeorf == "uoORF":
            df = df.with_columns(
                (pl.col("start").apply(lambda x: scoredict["rise_up"][x]).alias("rise_up")),
                (pl.lit(0.0).alias("step_down")),
            )

        elif typeorf == "doORF":
            df = df.with_columns(
                (
                    pl.col("stop")
                    .apply(lambda x: scoredict["step_down"][x])
                    .alias("step_down")
                ),
                (pl.lit(0.0).alias("rise_up")),
            )
        else:
            df = df.with_columns(
                (pl.col("start").apply(lambda x: scoredict["rise_up"][x]).alias("rise_up")),
                (
                    pl.col("stop")
                    .apply(lambda x: scoredict["step_down"][x])
                    .alias("step_down")
                ),
            )
        return df
    except Exception as e:
        log_error(f"Error assigning scores: {str(e)}")
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

    Example:
        orf_df, exon_coords = orfrelativeposition("annotation.gff", orf_df)
    """
    if not "cdsdf" in globals():
        cds_df, exon_coords = getexons_and_cds(annotation, list(df["tran_id"].unique()))

    print("Typing ORFS")
    tranids = cds_df["tran_id"].unique().to_list()

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
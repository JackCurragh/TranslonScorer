"""This script contains functions to  and calculate the transcriptomic coordinates"""

import polars as pl
import pyBigWig as bw
from .scoring import sru_score, calculate_scores
from .logging_config import log_info, log_warning, log_error


def transcriptreads(bwfile, exon_df):
    """
    Converts a BigWig file to a DataFrame based on provided exon annotation.

    Parameters:
    - bigwig (str): Path to the BigWig file.
    - exon (str): Path to the exon annotation file.

    Returns:
    - df_tran (DataFrame): DataFrame containing transcript information derived from the BigWig file.

    Raises:
    - RuntimeError: If there are issues reading values from the bigWig file
    - ValueError: If no reads could be extracted from any chromosome

    This function reads a BigWig file and an exon annotation file. It performs various operations to extract transcript information
    from the BigWig file based on the exon coordinates. The resulting transcript information is stored in a DataFrame named `df_tran`.
    The DataFrame includes columns for transcript ID, transcript start and stop coordinates, and counts.

    The function first checks if the given file is a BigWig file. It then reads the exon annotation file and extracts necessary
    information, such as the chromosome notation. It ensures that the chromosome notation in the BigWig file matches the exon annotation.
    Then, it iterates over each exon in the annotation file and retrieves intervals from the BigWig file that correspond to the exon
    coordinates. It calculates the transcript start and stop coordinates for each interval and stores the information in corresponding lists.
    Finally, it constructs the `df_tran` DataFrame using the extracted transcript information and returns it.
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
    Applies specific scoring calculations to a DataFrame based on the type of ORF.

    This function modifies the given DataFrame `df` by adding new columns based on the
    specified type of ORF (`uoORF`, `doORF`, or other). It calculates scores using
    provided `tran_reads` and `sru_range`, and then computes a final score.

    Parameters:
    df (pl.DataFrame): The input Data frame containing 'start' and 'stop' columns.
    tran_reads (df): Data frame containing reads on transcriptomic level required for scoring functions.
    sru_range (int): Range parameter required for SRU scoring functions.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.

    Returns:
    dict: A dictionary representation of the modified DataFrame.
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
    Updates a scoring dictionary with 'rise_up' and 'step_down' scores based on the DataFrame and ORF type.

    This function calculates and updates the `scoredict` with 'rise_up' and 'step_down' scores
    based on unique 'start' and 'stop' values from the DataFrame `df`. The scores are determined
    using the `sru_score` function, and the type of ORF (`typeorf`) dictates which scores are computed.

    Parameters:
    df (pl.DataFrame or pl.Series): The input DataFrame or Series containing 'start' and 'stop' columns.
    tran_reads (df): Data frame containing reads on transcriptomic level required for scoring functions.
    sru_range (int): Range parameter required for SRU scoring functions.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.
    scoredict (dict): Dictionary to store the computed 'rise_up' and 'step_down' scores.

    Returns:
    dict: The updated scoring dictionary with 'rise_up' and 'step_down' scores.

    Notes:
    - For types other than 'doORF', it calculates 'rise_up' scores based on unique 'start' values.
    - For types other than 'uoORF', it calculates 'step_down' scores based on unique 'stop' values.
    - The `scoredict` is updated with these scores, where keys are unique 'start' or 'stop' values
      and values are the corresponding scores from the `sru_score` function.
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
    Computes global scores for a DataFrame based on 'start' and 'stop' values and the type of ORF.

    This function modifies the given DataFrame `df` by calculating a list of scores from 'start'
    and 'stop' values and adds columns for specific scores (Highes reading frame/hrf, Average/avg, Non-zero coverage/nzc).
    It then computes a final score by summing up relevant columns based on the type of ORF (`typeorf`).

    Parameters:
    df (pl.DataFrame): The input DataFrame containing 'start' and 'stop' columns.
    tran_reads (df): Data frame containing reads on transcriptomic level required for scoring functions.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.

    Returns:
    dict: A dictionary representation of the modified DataFrame with computed scores.

    Notes:
    - Computes 'hrf', 'avg', and 'nzc' scores from 'start' and 'stop' columns using `calculate_scores`.
    - For 'doORF', the final score is the sum of 'step_down', 'hrf', 'avg', and 'nzc' columns.
    - For 'uoORF', the final score is the sum of 'rise_up', 'hrf', 'avg', and 'nzc' columns.
    - For other types, the final score is the sum of 'rise_up', 'step_down', 'hrf', 'avg', and 'nzc' columns.
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
    Filters the DataFrame to exclude rows with existing scores in the scoring dictionary.

    This function filters out 'start' and 'stop' values from the DataFrame `df` that already have
    corresponding scores in the `scoredict`. The filtering behavior depends on the type of ORF (`typeorf`).

    Parameters:
    df (pl.DataFrame): The input DataFrame containing 'start' and 'stop' columns.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.
    scoredict (dict): Dictionary containing the existing scores for 'rise_up' and 'step_down'.

    Returns:
    pl.DataFrame or pl.Series: A filtered DataFrame or Series excluding rows with existing scores.

    Notes:
    - For 'uoORF', filters out 'start' values that exist in `scoredict['rise_up']`.
    - For 'doORF', filters out 'stop' values that exist in `scoredict['step_down']`.
    - For other types, filters out rows where either 'start' is in `scoredict['rise_up']` or
      'stop' is in `scoredict['step_down']`, and retains rows where at least one of these conditions is met.
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

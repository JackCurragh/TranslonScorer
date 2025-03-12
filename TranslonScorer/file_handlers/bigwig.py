"""
BigWig file handling functionality for TranslonScorer.

This module contains functions for reading and processing BigWig files,
including conversion to other formats and coordinate transformations.
"""

import polars as pl
import pyBigWig as bw
from ..core.orffinder import orfrelativeposition
from ..utils.logging import log_info, log_warning, log_error
from ..core.scoring import oldscoring, newscoring, globalscores, existingscore, assigningscore


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


def scoring(bigwig, exon, orfs, old_scoring, sru_range):
    """
    Score ORFs using bigwig coverage data.
    
    Args:
        bigwig (str): Path to bigwig file
        exon (str or DataFrame): Path to exon file or DataFrame
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
    # Check if exon is a file path or a DataFrame
    if isinstance(exon, str):
        log_info(f"Reading exon data from file: {exon}")
        exon_df = pl.read_csv(exon, has_header=True, separator=",")
    else:
        log_info("Using provided exon DataFrame")
        exon_df = exon

    # Split and convert to integers
    exon_df = exon_df.with_columns([
        pl.col("start").apply(lambda x: [int(i) for i in x.split(",")] if isinstance(x, str) else x),
        pl.col("stop").apply(lambda x: [int(i) for i in x.split(",")] if isinstance(x, str) else x),
        pl.col("tran_start").apply(lambda x: [int(i) for i in x.split(",")] if isinstance(x, str) else x),
        pl.col("tran_stop").apply(lambda x: [int(i) for i in x.split(",")] if isinstance(x, str) else x)
    ])

    print(orfs)
    orf_df =orfs 
    # = pl.read_csv(orfs, has_header=True, separator=",")

    # Determine the relative position of ORFs to CDS
    orf_df, exon_coords = orfrelativeposition(annotation, orf_df, cds_df)

    # Proceed with scoring
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
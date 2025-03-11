"""
Visualization functionality for TranslonScorer.

This module contains functions for generating plots and visualizations,
including metagene profiles and transcript-specific plots.
"""

import polars as pl
from ..utils.logging import log_info, log_warning, log_error


def plottop10(df, bigwig, exon, range_param, filename, parameters):
    """
    Generate plots and tables summarizing top 10 ORFs per type and metagene profiles.

    Args:
        df (str): Path to CSV file containing ORF information
        bigwig (str): Path to BigWig file for transcript read counts
        exon (str): Path to CSV file containing exon information
        range_param (int): Range for relative coordinates around exon boundaries
        filename (str): Output filename for report
        parameters (dict): Additional parameters for report generation

    This function generates:
    1. Top 10 ORFs plots per type
    2. Metagene profiles based on transcript read counts
    3. Individual transcript plots
    4. Summary tables
    """
    range_list = list(range(-range_param, range_param + 1))
    
    # Read input files
    df = pl.read_csv(df)
    bwfile = bw.open(bigwig)
    
    # Process exon data
    exon_df = pl.read_csv(exon)
    exon_df = exon_df.with_columns([
        pl.col("start").apply(lambda x: [int(i) for i in x.split(",")]),
        pl.col("stop").apply(lambda x: [int(i) for i in x.split(",")]),
        pl.col("tran_start").apply(lambda x: [int(i) for i in x.split(",")]),
        pl.col("tran_stop").apply(lambda x: [int(i) for i in x.split(",")])
    ])
    
    # Generate plots
    plotlist = metageneplot(df, bwfile, exon_df, range_list)
    tranplot, table, pertranscript = pertranscriptplot(df, bwfile, exon_df)
    
    # Generate report
    generate_report(plotlist, tranplot, table, pertranscript, parameters, filename) 
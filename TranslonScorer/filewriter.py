"""Script for writing files that need to be stored during analysis"""

import os
import polars as pl


def saveorfsandexons(orf_df, exon_df, filename):
    """
    Saves annotated ORFs and exons dataframes to CSV files.

    This function takes two polars DataFrames containing information about annotated ORFs
    and exons respectively and saves them to CSV files.

    Parameters:
        orf_df (polars.DataFrame): DataFrame containing annotated ORFs data.
        exon_df (polars.DataFrame): DataFrame containing exon data.
        filename (str): Base name for the output files.

    Returns:
        tuple: Paths to the created ORF and exon CSV files.
    """
    orf_df.write_csv(f"{filename}_annotated_orfs.csv")

    # Keep original chromosome names and concatenate coordinate columns
    exon_df = exon_df.with_columns(
        pl.col("start", "stop", "tran_start", "tran_stop").apply(
            lambda x: ",".join(map(str, x))
        )
    )
    exon_df.write_csv(f"{filename}_exons.csv")
    return f"{filename}_annotated_orfs.csv", f"{filename}_exons.csv"

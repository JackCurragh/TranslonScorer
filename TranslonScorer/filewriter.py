"""Script for writing files that need to be stored during analysis"""

import os
import polars as pl


def saveorfsandexons(orf_df, exon_df, filename):
    """
    Saves annotated ORFs and exons dataframes to CSV files.

    This function takes two polars DataFrames containing information about annotated ORFs
    and exons respectively and saves them to CSV files in the 'data/files' directory.

    Parameters:
        orf_df (polars.DataFrame): DataFrame containing annotated ORFs data.
        exon_df (polars.DataFrame): DataFrame containing exon data.

    Returns:
        None

    Notes:
        - If the 'data/files' directory does not exist, it will be created.
        - The ORFs DataFrame will be saved as 'data/files/annotated_orfs.csv'.
        - The exon DataFrame will be saved as 'data/files/exons.csv'.
        - The 'chr' column in the exon DataFrame will be fixed by removing any extra characters,
          and columns 'start', 'stop', 'tran_start', 'tran_stop' will be concatenated into a single column.

    Example:
        saveorfsandexons(orf_df, exon_df)
    """

    orf_df.write_csv(f"{filename}_annotated_orfs.csv")

    fixchromcol = exon_df.with_columns(pl.col("chr").apply(lambda x: x[0]))

    exon_df = fixchromcol.with_columns(
        pl.col("start", "stop", "tran_start", "tran_stop").apply(
            lambda x: ",".join(map(str, x))
        )
    )
    exon_df.write_csv(f"{filename}_exons.csv")
    return f"{filename}_annotated_orfs.csv", f"{filename}_exons.csv"

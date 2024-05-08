"""This script contains functions to transform data frames into different file types"""

import polars as pl
import os


def asitecalc(df, offsets):
    """
    Calculates A-site positions and aggregates counts based on offset values.

    Parameters:
    - df (DataFrame): Input DataFrame containing 'length' and 'pos' columns for A-site calculation.
    - offsets (dict): Dictionary containing offset values for each 'length' value.

    Returns:
    - df_bed (DataFrame): DataFrame containing aggregated information of A-site positions, their counts, and chromosome information.

    This function calculates the A-site positions based on the provided offsets for different 'length' values. 
    It iterates through the input DataFrame 'df' and calculates the A-site positions using the 'pos' column and the corresponding offset value.
    The calculated A-site positions are stored in a new DataFrame 'df_asite' by concatenating the individual A-site DataFrames.

    The function then groups the 'df_asite' DataFrame by chromosome and A-site position, and aggregates the counts of A-site occurrences using the 'agg' function.
    Subsequently, it creates a new DataFrame 'df_tobed' by adding 1 to the A-site position and sorts the DataFrame based on the chromosome and A-site position.
    Finally, the function constructs the 'df_bed' DataFrame by selecting the 'chr', 'A-site', 'end', and 'count' columns, and returns it.

    """
    # GROUP TO CALCULATE A-SITE
    y = []
    for value, data in df.group_by("length"):
        x = df.filter(pl.col("length") == value).with_columns(
            (pl.col("pos") + offsets[value]).alias("A-site")
        )
        y.append(x)
    df_asite = pl.concat(y)

    # GROUP ON A-SITE
    df_final = df_asite.group_by("chr", "A-site").agg(pl.col("count").sum())
    df_tobed = df_final.with_columns((pl.col("A-site") + 1).alias("end"))
    df_bed = df_tobed.sort(["chr", "A-site"]).select(["chr", "A-site", "end", "count"])

    return df_bed


def dftobed(df):
    """
    Converts a DataFrame to BED format with A-site calculation and offset values.

    Parameters:
    - df (DataFrame): Input DataFrame to be converted to BED format.

    Returns:
    - bed (DataFrame): DataFrame representing the BED format data with A-site positions and aggregated counts.

    This function takes an input DataFrame 'df' and performs the following operations:
    1. Creates a new DataFrame 'df_filtered' by calculating the length using the 'end' and 'pos' columns and excluding the 'end' column from the selection.
    2. Splits the 'qname' column values specifically on "_x", applies a split function to extract the integer value, renames the columns to 'count' and 'chr', and stores the results in a new DataFrame 'df_namesplit'.
    3. Obtains the length values from the 'df_namesplit' DataFrame and creates offsets dictionary with a hardcoded offset value of 15 for each length value.
    4. Calls the 'asitecalc' function to calculate A-site positions using the 'df_namesplit' DataFrame and the generated offsets, and stores the resulting BED format data in the 'bed' DataFrame.

    """
    df_filtered = df.with_columns(
        (pl.col("end") - pl.col("pos")).alias("length")
    ).select(pl.all().exclude("end"))

    # split specifically on _x -> implement more!!!
    split_func = lambda s: int(s.split("_x")[1])
    df_namesplit = df_filtered.with_columns(pl.col("qname").apply(split_func)).rename(
        {"qname": "count", "rname": "chr"}
    )

    # OFFSETS -> DIFFERENT FUNCTION!!
    length_values = df_namesplit.get_column("length")
    offsets = {i: 15 for i in length_values}

    bed = asitecalc(df_namesplit, offsets)

    return bed

def bedtobigwig(file, chromsize):
    """
    Converts a bedGraph file to a bigWig file using the bedGraphToBigWig utility.

    Parameters:
        file (str): The path to the input bedGraph file.
        chromsize (str): The path to the chromosome sizes file.

    Returns:
        None

    Notes:
        - The output bigWig file will be named 'file.bw' and will be created in the same directory as the input file.
        - Requires the bedGraphToBigWig utility to be installed and accessible in the system path.

    Example:
        bedtobigwig("input.bedGraph", "chromsizes.txt")
    """
    os.system(f"bedGraphToBigWig {file} {chromsize} file.bw")

import polars as pl
import numpy as np

def extract_transcript_id(attr_str):
    """
    Extracts transcript ID from a GTF/GFF attribute string.

    This function takes a string representing attributes in GTF/GFF format
    and extracts the transcript ID from it. The transcript ID is typically
    found in attributes such as 'Parent=transcript:', 'ID=transcript:', 
    'transcript_id=', or ' transcript_id '.

    Parameters:
        attr_str (str): A string containing attributes in GTF/GFF format.

    Returns:
        str: The extracted transcript ID, or an empty string if not found.

    Example:
        transcript_id = extract_transcript_id('gene_id="ENSG00000223972"; transcript_id="ENST00000456328"; ')
    """
    for attr in attr_str.split(";"):
        if attr.startswith("Parent=transcript:") or attr.startswith("ID=transcript:"):
            return attr.split(":")[1]
        elif attr.startswith("transcript_id="):
            return attr.split("=")[1]
        elif attr.startswith(" transcript_id "):
            return attr.split(" ")[2].replace('"', "")
    return ""


def getexons_and_cds(annotation_file, exondf, tran=[]):
    """
    Extracts CDS and exon coordinates from an annotation file.

    This function reads an annotation file in GTF/GFF format and extracts the
    coordinates of coding sequences (CDS) and exons. It then processes these
    coordinates to obtain transcript-level coordinates for exons and returns
    the results.

    Parameters:
        annotation_file (str): The path to the annotation file in GTF/GFF format.
        tran (list): A list of transcript IDs to filter. Only coordinates corresponding
                     to these transcripts will be extracted if provided. Default is [].

    Returns:
        tuple: A tuple containing two polars DataFrames:
               - The first DataFrame contains CDS coordinates.
               - The second DataFrame contains exon coordinates.

    Notes:
        - This function assumes the annotation file has columns separated by tabs ('\t').
        - The annotation file is expected to have no header, with comment lines starting with '#'.
        - The following columns are expected in the annotation file: 'chr', 'type', 'start', 'stop',
          'strand', 'attributes'.
        - The 'attributes' column is expected to contain transcript IDs.
        - The function 'extract_transcript_id' is used to extract transcript IDs from the 'attributes' column.
        - If 'tran' is provided, only coordinates corresponding to the specified transcripts will be extracted.

    Example:
        cds_coords, exon_coords = getexons_and_cds("annotation.gff", tran=['ENST00000223972', 'ENST00000456328'])
    """
    df = (
        pl.read_csv(
            annotation_file,
            separator="\t",
            ignore_errors=True,
            has_header=False,
            truncate_ragged_lines=True,
            comment_prefix="#",
        )
        .select(
            ["column_1", "column_3", "column_4", "column_5", "column_7", "column_9"]
        )
        .rename(
            {
                "column_1": "chr",
                "column_3": "type",
                "column_4": "start",
                "column_5": "stop",
                "column_7": "strand",
                "column_9": "attributes",
            }
        )
    )
    df = df.with_columns(
        pl.col("attributes")
        .apply(lambda attributes: extract_transcript_id(attributes))
        .alias("tran_id")
    ).select(pl.all().exclude("attributes"))
    if tran:
        df = df.filter((pl.col("tran_id").is_in(tran)))

    # Getting CDS
    coding_regions = df.filter((pl.col("type") == "CDS"))

    groupedcds = (
        coding_regions.group_by("tran_id")
        .agg(
            pl.col("start"), pl.col("stop"), pl.col("chr")
        )
        .select(["chr", "tran_id", "start", "stop"])
    )

    # Getting exons
    if exondf.is_empty():
        exon_regions = df.filter((pl.col("type") == "exon"))

        pos_exons, neg_exons = procesexons(exon_regions)

        exon_coords_plus = exontranscriptcoords(pos_exons, posstrand=True)
        # column names switched to calculate inverse of positions for negative strands
        exon_coords_neg = exontranscriptcoords(neg_exons, posstrand=False)

        exondf = pl.concat([exon_coords_plus, exon_coords_neg]).select(
            pl.all().exclude("strand")
        )
    
    cds_coords = gettranscriptcoords(groupedcds, exondf)
    return cds_coords, exondf


def procesexons(df):
    """
    Processes exon data by separating them based on strand orientation.

    This function takes a DataFrame containing exon data and separates the exons
    based on their strand orientation (positive or negative). It groups the exons
    by transcript ID and aggregates the start, stop, strand, and chromosome information
    for each group.

    Parameters:
        df (polars.DataFrame): DataFrame containing exon data.

    Returns:
        tuple: A tuple containing two polars DataFrames:
               - The first DataFrame contains exons on the positive strand.
               - The second DataFrame contains exons on the negative strand.

    Example:
        pos_exons, neg_exons = procesexons(exon_df)
    """    
    exonplus = df.filter((pl.col("strand") == "+"))
    exonneg = df.filter((pl.col("strand") == "-"))

    groupedexonspos = (
        exonplus.group_by("tran_id")
        .agg(pl.col("start"), pl.col("stop"), pl.col("strand"), pl.col("chr"))
        .select(["chr", "tran_id", "start", "stop", "strand"])
    )
    groupedexonsneg = (
        exonneg.group_by("tran_id")
        .agg(pl.col("start"), pl.col("stop"), pl.col("strand"), pl.col("chr"))
        .select(["chr", "tran_id", "start", "stop", "strand"])
    )

    return groupedexonspos, groupedexonsneg


def exontranscriptcoords(df: pl.DataFrame, posstrand=True) -> pl.DataFrame:
    """
    Calculates transcript-level coordinates for exons.

    This function takes a DataFrame containing exon coordinates and calculates
    the transcript-level coordinates for each exon. For exons on the positive strand,
    the transcript coordinates start from 0 and increase. For exons on the negative
    strand, the transcript coordinates are calculated in reverse order.

    Parameters:
        df (polars.DataFrame): DataFrame containing exon coordinates.
        posstrand (bool): Flag indicating whether the exons are on the positive strand.
                          If True, transcript coordinates are calculated assuming exons
                          are on the positive strand. If False, transcript coordinates
                          are calculated assuming exons are on the negative strand.
                          Default is True.

    Returns:
        polars.DataFrame: DataFrame containing transcript-level exon coordinates.

    Example:
        exon_transcript_coords = exontranscriptcoords(exon_df, posstrand=True)
    """
    # Initialize new columns
    new_column_1 = []
    new_column_2 = []
    # Iterate over rows
    start_column = "start"
    end_column = "stop"

    for i in range(len(df)):
        start_values = (
            df[start_column][i]
            if posstrand
            else sorted(df[start_column][i], reverse=True)
        )
        end_values = (
            df[end_column][i] if posstrand else sorted(df[end_column][i], reverse=True)
        )
        new_start_values = []  # Starting value is 0
        new_stop_values = []
        for j in range(len(start_values)):
            if j == 0:
                new_start = 0
            else:
                new_start = new_stop_values[j - 1] + 1
            # Calculate stop coordinate
            stop_coordinate = end_values[j] - start_values[j]
            new_start_values.append(new_start)
            new_stop_values.append(new_start + stop_coordinate)
        new_column_1.append(new_start_values)
        new_column_2.append(new_stop_values)

    # Add new columns to the dataframe
    df = df.with_columns((pl.Series(new_column_1)).alias("tran_start"))
    df = df.with_columns((pl.Series(new_column_2)).alias("tran_stop"))
    return df


def gettranscriptcoords(cds_df, exon_df):
    """
    Calculates transcript-level coordinates for coding sequences (CDS).

    This function takes a DataFrame containing CDS coordinates and a DataFrame
    containing exon coordinates. It then calculates the transcript-level coordinates
    for the CDS based on the exon coordinates.

    Parameters:
        cds_df (polars.DataFrame): DataFrame containing CDS coordinates.
        exon_df (polars.DataFrame): DataFrame containing exon coordinates.

    Returns:
        polars.DataFrame: DataFrame containing transcript-level CDS coordinates.

    Example:
        transcript_cds_coords = gettranscriptcoords(cds_df, exon_df)
    """
    # Explode exon start and stop lists to individual rows
    exploded_exon_df = exon_df.explode(["start", "stop", "tran_start", "tran_stop"])
    exploded_cds_df = cds_df.explode(['start', 'stop'])
    # Join the DataFrames on the transcript ID
    combined_df = exploded_cds_df.join(exploded_exon_df, on="tran_id", how="inner")
    # Calculate transcript-level start and stop coordinates
    combined_df = combined_df.with_columns([
        (pl.when((pl.col("start") >= pl.col("start_right")) & (pl.col("start") <= pl.col("stop_right")))
         .then(pl.col("tran_start") + (pl.col("start") - pl.col("start_right")))
         .otherwise(None)).alias("tran_start_cd"),
        
        (pl.when((pl.col("stop") >= pl.col("start_right")) & (pl.col("stop") <= pl.col("stop_right")))
         .then(pl.col("tran_stop") - (pl.col("stop_right") - pl.col("stop")))
         .otherwise(None)).alias("tran_stop_cd")
    ])
    # Drop rows with None values in calculated columns
    combined_df = combined_df.drop_nulls(["tran_start_cd", "tran_stop_cd"])
    
    # Select and rename relevant columns
    result_df = combined_df.select([
        pl.col("tran_id"),
        pl.col("start"),
        pl.col("stop"),
        pl.col("tran_start_cd").alias("tran_start"),
        pl.col("tran_stop_cd").alias("tran_stop")
    ])
    return result_df

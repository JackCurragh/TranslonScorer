"""This script contains functions to transform data frames into different file types"""

import polars as pl
import os

from findexonscds import extract_transcript_id, exontranscriptcoords, procesexons
def getexons(annotation):
    '''
    docstring
    '''
    df = (
        pl.read_csv(
            annotation,
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
    
    exon_regions = df.filter((pl.col("type") == "exon"))

    pos_exons, neg_exons = procesexons(exon_regions)

    exon_coords_plus = exontranscriptcoords(pos_exons, posstrand=True)
    # column names switched to calculate inverse of positions for negative strands
    exon_coords_neg = exontranscriptcoords(neg_exons, posstrand=False)

    exon_coords = pl.concat([exon_coords_plus, exon_coords_neg]).select(
        pl.all().exclude("strand")
    )
    return exon_coords
    
def get_bam_tran(bam_df, exon_df):
    exon_df = exon_df.rename({'start': 'start_left', 'stop': 'stop_left'})

    df_joined = bam_df.join(exon_df, on='chr')

    df_filtered = df_joined.filter(
        (pl.col('start') >= pl.col('start_left')) &
        (pl.col('stop') <= pl.col('stop_left'))
    )
    df_filtered = df_filtered.with_columns(
        (pl.col('tran_start')+(pl.col('start')-pl.col('start_left'))).alias('tran_start_bam'),
        (pl.col('tran_stop')-(pl.col('stop_left')-pl.col('stop'))).alias('tran_stop_bam'))\
        .select(pl.all().exclude('stop_left', 'start_left', 'tran_start','tran_stop'))

    bam_tran_dict = df_filtered.to_dict(as_series=False)

    return bam_tran_dict


def bamtranscript(bam_df, exon_df):
    '''
    docstring
    '''
    exon_flattened = exon_df.with_columns(pl.col('chr').apply(lambda x: x[0].split('r')[1]))
    
    uniquechr_bam = list(bam_df['chr'].unique())
    uniquechr_exon = list(exon_flattened['chr'].unique())

    bam_df = bam_df.with_columns(shared=pl.col('chr').is_in(uniquechr_exon))\
            .filter(pl.col('shared') == True)\
            .select(pl.all().exclude('shared'))
    exon_flattened = exon_flattened.with_columns(shared=pl.col('chr').is_in(uniquechr_bam))\
                .filter(pl.col('shared') == True)\
                .select(pl.all().exclude('shared'))\
                .explode(['start', 'stop', 'tran_start', 'tran_stop'])
    
    results = []
    for chr in bam_df['chr'].unique():
        # Filter exons that overlap with the current bam range
        bam_chr = bam_df.filter(pl.col('chr') == chr)
        exon_flattened = exon_flattened.filter(pl.col('chr') == chr)
        min_val = min(exon_flattened['start'])
        max_val = max(exon_flattened['stop'])
        print(min_val, max_val)
        for i in range(min_val, max_val, 20000000):
            rng = i + 20000000
            exon_chr = exon_flattened.filter((pl.col('start') >= i) & (pl.col('stop') <= rng))

            result_dict = get_bam_tran(bam_chr, exon_chr)
            results.append(result_dict)
    
    bam_df = pl.from_dicts(results)
    print(bam_df)
    return bam_df



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
            (pl.col("start") + offsets[value]).alias("A-site")
        )
        y.append(x)
    df_asite = pl.concat(y)

    # GROUP ON A-SITE
    df_final = df_asite.group_by("chr", "A-site").agg(pl.col("count").sum())
    df_tobed = df_final.with_columns((pl.col("A-site") + 1).alias("stop"))
    df_bed = df_tobed.sort(["chr", "A-site"]).select(["chr", "A-site", "stop", "count"])

    return df_bed

def dftobed(df, annotation):
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
    ).select(pl.all().exclude('flag', 'qual', 'tags','mapq','tlen','seq','pnext','rnext','cigar'))

    # split specifically on _x -> implement more!!!
    split_func = lambda s: int(s.split("_x")[1])
    df_namesplit = df_filtered.with_columns(pl.col("qname").apply(split_func)).rename(
        {"qname": "count", "rname": "chr", 'pos':'start', 'end':'stop'}
    ).cast({'chr':pl.String})

    exon_df = getexons(annotation)
    #calculate transcriptomic coordinates
    bam_tran = bamtranscript(df_namesplit, exon_df)
    print(bam_tran)

    # OFFSETS -> DIFFERENT FUNCTION!!
    length_values = df_namesplit.get_column("length")
    offsets = {i: 15 for i in length_values}

    bed = asitecalc(df_namesplit, offsets)

    return bed, exon_df

def bedtobigwig(bedfile, chromsize):
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
    os.system(f"bedGraphToBigWig {bedfile} {chromsize} file.bw")

    return 'data/file.bw'

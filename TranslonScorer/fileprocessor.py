"""This script contains functions to transform data frames into different file types"""

import polars as pl
import pandas as pd
import os
from typing import Optional
from .logging_config import log_info, log_warning, log_error

from .findexonscds import getexons_and_cds


def get_bam_tran(bam_df, exon_df):
    """
    Merge BAM and exon DataFrames based on chromosome alignment and calculate transcript coordinates in BAM.

    Parameters:
    - bam_df (DataFrame): DataFrame containing BAM file data with chromosome positions.
    - exon_df (DataFrame): DataFrame containing exon annotations with chromosome positions.

    Returns:
    - dict: A dictionary where keys are column names and values are lists of corresponding values,
            representing transcript coordinates aligned with BAM data.

    This function merges `bam_df` and `exon_df` on the 'chr' column and filters rows where the BAM
    positions fall within exon boundaries (`start_right` to `stop_right`). It then calculates
    transcript coordinates (`tran_start_bam` and `tran_stop_bam`) in the BAM file by adjusting
    coordinates relative to exon boundaries.

    The resulting DataFrame (`df_filtered`) excludes redundant columns (`stop_right`, `start_right`,
    `tran_start`, `tran_stop`) and is converted to a dictionary (`bam_tran_dict`) where each key
    corresponds to a column name and the associated value is a list of values from `df_filtered`.

    Note: This function assumes the use of a library like `pandas` (abbreviated here as `pl`) for
    DataFrame operations.
    """
    df_joined = bam_df.join(exon_df, on="chr")
    df_filtered = df_joined.filter(
        (pl.col("start") >= pl.col("start_right"))
        & (pl.col("stop") <= pl.col("stop_right"))
    )
    df_filtered = df_filtered.with_columns(
        (pl.col("tran_start") + (pl.col("start") - pl.col("start_right"))).alias(
            "tran_start_bam"
        )
    )

    df_filtered = df_filtered.with_columns(
        (pl.col("tran_stop") - (pl.col("stop_right") - pl.col("stop"))).alias(
            "tran_stop_bam"
        )
    ).select(pl.all().exclude("stop_right", "start_right", "tran_start", "tran_stop"))
    bam_tran_dict = df_filtered.to_dict(as_series=False)

    return bam_tran_dict


def bamtranscript(bam_df, exon_df):
    """
    Filter BAM and exon DataFrames based on shared chromosome information and flatten exon annotations.

    Parameters:
    - bam_df (DataFrame): DataFrame containing BAM file data with chromosome positions.
    - exon_df (DataFrame): DataFrame containing exon annotations with chromosome positions.

    Returns:
    - DataFrame: A modified version of `bam_df` and `exon_df` after filtering based on shared
                 chromosome information and exploding exon coordinates.

    Raises:
    - ValueError: If there's no overlap between BAM and exon chromosomes
    - ValueError: If no overlapping regions are found between BAM and exon coordinates
    """
    log_info("Starting BAM to transcript conversion")
    exon_flattened = exon_df.with_columns(pl.col("chr"))

    # Get unique chromosomes and check overlap
    uniquechr_bam = set(bam_df["chr"].unique())
    uniquechr_exon = set(exon_flattened["chr"].unique())
    
    log_info(f"Found {len(uniquechr_bam)} unique chromosomes in BAM and {len(uniquechr_exon)} in annotation")
    
    # First try direct matching
    common_chr = uniquechr_bam.intersection(uniquechr_exon)
    
    if not common_chr:
        log_info("No direct chromosome matches found, attempting chromosome name normalization")
        # Check if this might be due to chromosome naming
        bam_has_chr = any(c.startswith('chr') for c in uniquechr_bam)
        exon_has_chr = any(c.startswith('chr') for c in uniquechr_exon)
        
        if bam_has_chr != exon_has_chr:
            if exon_has_chr:
                log_info("Adding 'chr' prefix to BAM chromosomes")
                bam_df = bam_df.with_columns(
                    pl.when(pl.col("chr").str.starts_with("chr"))
                    .then(pl.col("chr"))
                    .otherwise(pl.concat_str([pl.lit("chr"), pl.col("chr")]))
                    .alias("chr")
                )
                uniquechr_bam = set(bam_df["chr"].unique())
            else:
                log_info("Removing 'chr' prefix from BAM chromosomes")
                bam_df = bam_df.with_columns(
                    pl.col("chr").str.replace("^chr", "").alias("chr")
                )
                uniquechr_bam = set(bam_df["chr"].unique())
            
            # Try matching again after normalization
            common_chr = uniquechr_bam.intersection(uniquechr_exon)
            
            if common_chr:
                log_info(f"After normalization, found {len(common_chr)} matching chromosomes")
            else:
                # Create stripped versions for error message
                bam_chr_stripped = {c.replace("chr", "") for c in uniquechr_bam}
                exon_chr_stripped = {c.replace("chr", "") for c in uniquechr_exon}
                error_msg = (
                    "No overlapping chromosomes found between BAM and annotation files after normalization.\n"
                    f"BAM chromosomes: {sorted(uniquechr_bam)[:5]}\n"
                    f"Annotation chromosomes: {sorted(uniquechr_exon)[:5]}\n"
                    f"BAM stripped: {sorted(bam_chr_stripped)[:5]}\n"
                    f"Annotation stripped: {sorted(exon_chr_stripped)[:5]}"
                )
                log_error(error_msg)
                raise ValueError(error_msg)

    # Filter to matching chromosomes
    bam_df = (
        bam_df.with_columns(shared=pl.col("chr").is_in(uniquechr_exon))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )
    exon_flattened = (
        exon_flattened.with_columns(shared=pl.col("chr").is_in(uniquechr_bam))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
        .explode(["start", "stop", "tran_start", "tran_stop"])
    )
    
    if bam_df.is_empty():
        error_msg = (
            "No overlapping chromosomes found between BAM and annotation after name normalization.\n"
            f"BAM chromosomes: {sorted(uniquechr_bam)}\n"
            f"Annotation chromosomes: {sorted(uniquechr_exon)}"
        )
        log_error(error_msg)
        raise ValueError(error_msg)
    
    # Process matching chromosomes
    results = []
    total_chr = len(set(bam_df["chr"].unique()))
    log_info(f"\nProcessing {total_chr} chromosomes")
    
    for idx, chr in enumerate(set(bam_df["chr"].unique()), 1):
        if idx % 5 == 0 or idx == total_chr:
            log_info(f"Processing chromosome {chr} ({idx}/{total_chr})")
            
        bam_with_chr = bam_df.filter(pl.col("chr") == chr)
        exon_with_chr = exon_flattened.filter(pl.col("chr") == chr)
        
        if exon_with_chr.is_empty() or bam_with_chr.is_empty():
            log_warning(f"No data found for chromosome {chr}")
            continue
            
        min_val = min(exon_with_chr["start"])
        max_val = max(exon_with_chr["stop"])
        total_chunks = (max_val - min_val) // 225000 + 1
        
        if total_chunks > 100:  # Only show chunk progress for larger chromosomes
            log_info(f"  Processing {total_chunks} chunks for chromosome {chr}")
        
        for chunk_idx, i in enumerate(range(min_val, max_val, 225000), 1):
            if total_chunks > 100 and (chunk_idx % 50 == 0 or chunk_idx == total_chunks):
                log_info(f"    Chunk {chunk_idx}/{total_chunks}")
                
            rng = i + 225000
            exon_chr = exon_with_chr.filter(
                (pl.col("start") >= i) & (pl.col("stop") <= rng)
            )
            bam_chr = bam_with_chr.filter(
                (pl.col("start") >= i) & (pl.col("stop") <= rng)
            )
            if not exon_chr.is_empty() and not bam_chr.is_empty():
                result_dict = get_bam_tran(bam_chr, exon_chr)
                results.append(result_dict)

    if not results:
        error_msg = (
            "No overlapping regions found between BAM and annotation files.\n"
            "This could be due to:\n"
            "1. Mismatched coordinates\n"
            "2. BAM file aligned to different genome version than annotation\n"
            "3. No reads mapping to annotated regions"
        )
        log_error(error_msg)
        raise ValueError(error_msg)

    log_info("\nCreating final BAM transcript DataFrame")
    bam_df = pl.from_dicts(results)
    bam_df = bam_df.explode(
        [
            "count",
            "chr",
            "start",
            "stop",
            "length",
            "tran_id",
            "tran_start_bam",
            "tran_stop_bam",
        ]
    )
    log_info("BAM transcript conversion complete")
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
        x = (
            df.filter(pl.col("length") == value)
            .with_columns((pl.col("start") + offsets[value]).alias("A-site"))
            .select(
                pl.all().exclude("start", "stop", "length", "tran_id", "bamcds_start")
            )
            .to_dict(as_series=False)
        )
        y.append(x)
    df_asite = pl.from_dicts(y)
    df_asite = df_asite.explode(["chr", "count", "A-site"])
    # GROUP ON A-SITE
    df_asite = df_asite.group_by("chr", "A-site").agg(pl.col("count").sum())
    df_asite = df_asite.with_columns((pl.col("A-site") + 1).alias("stop"))
    df_bed = df_asite.sort(["chr", "A-site"]).select(["chr", "A-site", "stop", "count"])
    return df_bed


def calculate_differences(start, start_dict):
    """
    Calculate the difference between a given start value and each value in a dictionary of start values.

    Parameters:
    - start (int or float): The starting value for which differences are calculated.
    - start_dict (int): A values extracted from a dictionary containing start values (int or float).

    Returns:
    - dict: A dictionary where keys correspond to identifiers from `start_dict` and values represent
            the differences between `start` and each corresponding value in `start_dict`.

    This function computes the difference (`start - start_dict[key]`) for each key-value pair in `start_dict`.
    The result is returned as a dictionary where keys are from `start_dict` and values are the calculated differences.

    Note: This function assumes `start_dict` contains numerical values (integers or floats).
    """
    # Handle polars Series objects
    if isinstance(start, pl.Series):
        start = start[0]  # Take the first value if it's a Series
    start_diff = start - start_dict

    return start_diff


def bamrelativetocds(bamdf, cdsdf):
    """
    Filter BAM and CDS DataFrames based on shared transcript IDs and calculate relative start positions.

    Parameters:
    - bamdf (DataFrame): DataFrame containing BAM file data with transcript information.
    - cdsdf (DataFrame): DataFrame containing CDS annotations with transcript information.

    Returns:
    - DataFrame: A modified version of `bamdf` after calculating relative start positions relative to CDS.

    This function processes `bamdf` and `cdsdf` to filter rows where transcript IDs are shared between
    the two datasets. First, it identifies unique transcript IDs present in both `bamdf` and `cdsdf`.
    Then, it filters `bamdf` to include only rows where the transcript ID matches those found in `cdsdf`,
    and vice versa for `cdsdf`.

    Next, it creates a dictionary (`start_dict`) mapping each transcript ID to its corresponding start
    position in `cdsdf`. Using this dictionary, it calculates relative start positions (`bamcds_start`)
    in `bamdf` by subtracting the start position in `cdsdf` from the start position in `bamdf`.

    The modified `bamdf` DataFrame (`bam_df2`) excludes redundant columns (`tran_start_bam`, `tran_stop_bam`)
    and retains only the calculated `bamcds_start` values for each row.

    Note: This function assumes the use of a library like `pandas` (abbreviated here as `pl`) for
    DataFrame operations, and it depends on a helper function `calculate_differences` to perform
    numerical calculations.
    """
    uniquetran_bam = list(bamdf["tran_id"].unique())
    uniquetran_cds = list(cdsdf["tran_id"].unique())

    bam_df = (
        bamdf.with_columns(shared=pl.col("tran_id").is_in(uniquetran_cds))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )
    cds_df = (
        cdsdf.with_columns(shared=pl.col("tran_id").is_in(uniquetran_bam))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )

    tran_dict = cds_df.to_dict(as_series=False)
    start_dict = dict(zip(tran_dict["tran_id"], tran_dict["tran_start"]))

    bam_df2 = bam_df.with_columns(
        [
            pl.struct(["tran_id", "tran_start_bam"])
            .apply(
                lambda x: calculate_differences(
                    x["tran_start_bam"], start_dict[x["tran_id"]]
                )
            )
            .alias("bamcds_start")
        ]
    )

    bam_df2 = bam_df2.select(pl.all().exclude("tran_start_bam", "tran_stop_bam"))
    return bam_df2


def change_point_analysis(offset_df):
    """
    Calculate the change point for the metagene profile using vectorized operations.
    """
    log_info(f"Processing offset analysis for {len(offset_df['length'].unique())} different read lengths")
    
    offset_dict = {}
    total_lengths = len(offset_df["length"].unique())
    
    for idx, length in enumerate(offset_df["length"].unique(), 1):
        if idx % 10 == 0 or idx == total_lengths:
            log_info(f"Processing length {length} ({idx}/{total_lengths})")
            
        # Get data for this length and sort
        offset_df_len = offset_df.filter(pl.col("length") == length).sort("bamcds_start")
        
        if offset_df_len.is_empty():
            log_warning(f"No data found for length {length}, using default offset of 15")
            offset_dict[length] = 15  # default offset
            continue
            
        # Pre-calculate all positions we need
        positions = list(range(-30, 11))
        max_shift = 0
        max_shift_position = None
        
        # Create a lookup dictionary for counts at each position
        count_dict = dict(zip(
            offset_df_len["bamcds_start"].to_list(),
            offset_df_len["count"].to_list()
        ))
        
        # Vectorized calculation of shifts
        for i in positions:
            left_positions = range(i - 3, i + 1)
            right_positions = range(i + 1, i + 5)
            
            # Get counts, defaulting to 0 for missing positions
            left_counts = [count_dict.get(pos, 0) for pos in left_positions]
            right_counts = [count_dict.get(pos, 0) for pos in right_positions]
            
            mean_left = sum(left_counts) / 4
            mean_right = sum(right_counts) / 4
            shift = abs(mean_right - mean_left)
            
            if shift > max_shift:
                max_shift = shift
                max_shift_position = i
        
        offset_dict[length] = max_shift_position if max_shift_position is not None else 15
    
    log_info("Offset analysis complete")
    return offset_dict


def detect_bam_type(df, exon_df):
    """
    Detect whether a BAM file is genomic or transcriptomic by checking chromosome/transcript ID patterns.
    
    Parameters:
    - df (DataFrame): BAM DataFrame with chromosome/transcript information
    - exon_df (DataFrame): Exon DataFrame with both chromosome and transcript IDs
    
    Returns:
    - str: 'genomic' or 'transcriptomic'
    
    Raises:
    - ValueError: If BAM type cannot be determined or if no matching IDs found
    """
    bam_ids = set(df["chr"].unique())
    exon_chroms = set(exon_df["chr"].unique())
    exon_trans = set(exon_df["tran_id"].unique())
    
    # Check for genomic BAM (chromosome matches)
    chrom_match = len(bam_ids.intersection(exon_chroms))

    if chrom_match == 0:
        exon_chroms = {chr.replace("chr", "") for chr in exon_chroms}
        bam_ids = {chr.replace("chr", "") for chr in bam_ids}
        chrom_match = len(bam_ids.intersection(exon_chroms))
    
    # Check for transcriptomic BAM (transcript matches)
    trans_match = len(bam_ids.intersection(exon_trans))
    
    if chrom_match > trans_match:
        return 'genomic', {}
    elif trans_match > chrom_match:
        return 'transcriptomic', {}
    else:
        raise ValueError(
            "Unable to determine BAM type. No significant matches found with either:\n"
            f"Chromosomes in annotation: {sorted(exon_chroms)[:5]}\n"
            f"Transcript IDs in annotation: {sorted(exon_trans)[:5]}\n"
            f"IDs in BAM: {sorted(bam_ids)[:5]}"
        )


def process_transcriptomic_bam(df_namesplit, cds_df):
    """
    Process a transcriptomic BAM file where reads are already aligned to transcripts.
    
    Parameters:
    - df_namesplit (DataFrame): BAM DataFrame with transcript alignments
    - cds_df (DataFrame): CDS DataFrame with transcript information
    
    Returns:
    - DataFrame: Processed BAM data ready for A-site calculation
    """
    # Rename chr column to tran_id since it contains transcript IDs
    df_with_tran = df_namesplit.rename({"chr": "tran_id"})
    
    # Filter for transcripts present in CDS annotation
    bam_to_cds = (
        df_with_tran.with_columns(shared=pl.col("tran_id").is_in(cds_df["tran_id"]))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )
    
    # Calculate position relative to CDS start
    tran_dict = cds_df.to_dict(as_series=False)
    start_dict = dict(zip(tran_dict["tran_id"], tran_dict["tran_start"]))
    
    bam_to_cds = bam_to_cds.with_columns(
        [
            pl.struct(["tran_id", "start"])
            .apply(lambda x: calculate_differences(x["start"], start_dict[x["tran_id"]]))
            .alias("bamcds_start")
        ]
    )
    
    return bam_to_cds


def dftobed(df, annotation, offsets):
    """
    Converts a DataFrame to BED format with A-site calculation and offset values.
    Automatically detects and handles both genomic and transcriptomic BAM files.

    Parameters:
    - df (DataFrame): Input DataFrame to be converted to BED format
    - annotation (str): Path to annotation file
    - offsets (dict, optional): Pre-calculated offsets for A-site calculation

    Returns:
    - tuple: (bed DataFrame, exon DataFrame, CDS DataFrame)
    """
    log_info("Starting BAM to BED conversion")
    log_info("Filtering and preparing BAM data")
    df_filtered = df.with_columns(
        (pl.col("end") - pl.col("pos")).alias("length")
    ).select(
        pl.all().exclude(
            "flag", "qual", "tags", "mapq", "tlen", "seq", "pnext", "rnext", "cigar"
        )
    )

    # Split read names to get counts
    log_info("Processing read names")
    split_func = lambda s: int(s.split("_x")[1])
    df_namesplit = (
        df_filtered.with_columns(pl.col("qname").apply(split_func))
        .rename({"qname": "count", "rname": "chr", "pos": "start", "end": "stop"})
        .cast({"chr": pl.String})
    )

    # Get annotations
    log_info("Loading annotations")
    cds_df, exon_df = getexons_and_cds(annotation)
    
    # Detect BAM type and get matching info
    log_info("Detecting BAM type")
    bam_type, info = detect_bam_type(df_namesplit, exon_df)
    log_info(f"Detected BAM type: {bam_type}")
    
    # Process based on BAM type
    if bam_type == 'genomic':
        log_info("Processing genomic BAM")
        bam_tran = bamtranscript(df_namesplit, exon_df)
        log_info("Converting BAM coordinates to CDS coordinates")
        bam_to_cds = bamrelativetocds(bam_tran, cds_df)
    else:  # transcriptomic
        log_info("Processing transcriptomic BAM")
        bam_to_cds = process_transcriptomic_bam(df_namesplit, cds_df)
    
    # Calculate offsets if not provided
    if not offsets:
        log_info("Calculating read length offsets")
        bam_offsets = bam_to_cds.group_by("bamcds_start", "length").agg(
            pl.col("count").sum()
        )
        log_info("Performing change point analysis")
        offsets = change_point_analysis(bam_offsets)
    else:
        log_info("Using provided offsets")
    
    # A-site calculation
    log_info("Calculating A-sites")
    bed = asitecalc(bam_to_cds, offsets)
    log_info("BAM to BED conversion complete")

    return bed, exon_df, cds_df


def bedtobigwig(bedfile, chromsize, filename):
    """
    Converts a bedGraph file to a bigWig file using pyBigWig.

    Parameters:
        bedfile (str): The path to the input bedGraph file.
        chromsize (str): The path to the chromosome sizes file.
        filename (str): The name for the generated file.

    Returns:
        str: Path to the created bigWig file

    Notes:
        - The output bigWig file will be named 'filename.bw'
        - Uses pyBigWig library instead of external kent utils
        - Handles chromosome sizes file reading internally
        - Ensures data is sorted by chromosome and position before writing
        - Uses span=1 for single-nucleotide resolution (appropriate for ribo-seq data)
    """
    import pyBigWig as bw
    import polars as pl
    
    log_info("Reading chromosome sizes")
    # Read chromosome sizes file
    chrom_sizes = {}
    with open(chromsize, 'r') as f:
        for line in f:
            chrom, size = line.strip().split('\t')
            chrom_sizes[chrom] = int(size)
    
    log_info("Reading bedGraph data")
    # Read bedGraph data with proper column types
    bed_data = pl.read_csv(
        bedfile, 
        separator='\t', 
        has_header=False,
        dtypes={
            "column_1": pl.String,
            "column_2": pl.Int64,
            "column_3": pl.Int64,
            "column_4": pl.Float64
        }
    )
    bed_data = bed_data.rename({
        "column_1": "chrom",
        "column_2": "start",
        "column_3": "end",
        "column_4": "value"
    })
    
    # Sort the data by chromosome and start position
    bed_data = bed_data.sort(["chrom", "start"])
    
    # Create bigWig file
    log_info("Creating bigWig file")
    bw_file = bw.open(f"{filename}.bw", "w")
    
    # Add header with chromosome sizes
    bw_file.addHeader(list(chrom_sizes.items()))
    
    # Write data chromosome by chromosome
    for chrom in bed_data["chrom"].unique():
        chrom_data = bed_data.filter(pl.col("chrom") == chrom)
        if not chrom_data.is_empty():
            try:
                # Convert polars Series to lists (no need to cast since types are already correct)
                starts = chrom_data["start"].to_list()
                ends = chrom_data["end"].to_list()
                values = chrom_data["value"].to_list()
                
                # Ensure all lists have the same length
                if len(starts) != len(ends) or len(starts) != len(values):
                    log_warning(f"Skipping chromosome {chrom} due to mismatched data lengths")
                    continue
                    
                # Ensure all positions are valid
                if any(end <= start for start, end in zip(starts, ends)):
                    log_warning(f"Skipping chromosome {chrom} due to invalid positions (end <= start)")
                    continue
                
                # Ensure all values are valid numbers using polars
                if chrom_data["value"].is_null().any():
                    log_warning(f"Skipping chromosome {chrom} due to null values")
                    continue
                
                # Create a list of chromosome names matching the length of other lists
                chromosomes = [chrom] * len(starts)
                
                # Add entries to bigWig file with span=1 for single-nucleotide resolution
                bw_file.addEntries(
                    chromosomes,
                    starts,
                    ends=ends,
                    values=values,
                    span=1  # Use single-nucleotide resolution for ribo-seq data
                )
            except Exception as e:
                log_warning(f"Error processing chromosome {chrom}: {str(e)}")
                continue
    
    bw_file.close()
    log_info(f"Successfully created {filename}.bw")
    return f"{filename}.bw"

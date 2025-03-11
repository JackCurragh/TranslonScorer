"""
BAM file handling functionality for TranslonScorer.

This module contains functions for processing BAM files and converting
them to other formats, including coordinate transformations.
"""

import pysam
import polars as pl
import oxbow as ox
from ..utils.logging import log_info, log_error, log_warning


def readbam(bampath):
    """
    Reads a given BAM file, extracts relevant information, and returns it as a DataFrame.

    Parameters:
    - bampath (str): Path to the BAM file to be processed.

    Returns:
    - df (DataFrame): Polars DataFrame containing the extracted information from the BAM file.
                     Contains columns: chr, start, stop, length, strand, count
    """
    log_info("Indexing BAM file")
    pysam.index(bampath)
    log_info("BAM file indexed successfully")
    
    # Read BAM file using pysam
    bam = pysam.AlignmentFile(bampath, "rb")
    
    # Extract relevant information
    records = []
    for read in bam.fetch():
        if read.is_unmapped:
            continue
            
        records.append({
            'chr': bam.get_reference_name(read.reference_id),
            'start': read.reference_start,
            'stop': read.reference_end,
            'length': read.query_length,
            'strand': '-' if read.is_reverse else '+',
            'count': 1
        })
    
    # Convert to DataFrame
    df = pl.DataFrame(records)
    
    # Group by position to get counts
    df = df.group_by(['chr', 'start', 'stop', 'length', 'strand']).agg(
        pl.col('count').sum()
    ).sort(['chr', 'start'])
    
    log_info(f"Processed {len(df)} unique read positions")
    return df


def getexons_and_cds(annotation_file, tran=[]):
    """
    Extract CDS and exon coordinates from an annotation file.

    Args:
        annotation_file (str): Path to annotation file in GTF/GFF format
        tran (list): Optional list of transcript IDs to filter

    Returns:
        tuple: (CDS DataFrame, exon DataFrame)
    """
    log_info("Reading annotation file")
    
    # Read annotation file
    df = pl.read_csv(
        annotation_file,
        has_header=False,
        separator="\t",
        comment_char="#",
        columns=["column_1", "column_3", "column_4", "column_5", "column_7", "column_9"],
    ).rename({
        "column_1": "chr",
        "column_3": "type",
        "column_4": "start",
        "column_5": "stop",
        "column_7": "strand",
        "column_9": "attributes",
    })
    
    # Extract transcript IDs
    df = df.with_columns(
        pl.col("attributes")
        .str.extract(r'transcript_id "([^"]*)"')
        .alias("tran_id")
    )
    
    # Filter by transcript IDs if provided
    if tran:
        df = df.filter(pl.col("tran_id").is_in(tran))
    
    # Split into CDS and exon DataFrames
    cds_df = df.filter(pl.col("type") == "CDS")
    exon_df = df.filter(pl.col("type") == "exon")
    
    if cds_df.is_empty() or exon_df.is_empty():
        log_error("No CDS or exon features found in annotation file")
    
    # Group and sort coordinates
    cds_df = (
        cds_df.group_by("tran_id")
        .agg([
            pl.col("chr").first(),
            pl.col("start").sort(),
            pl.col("stop").sort(),
            pl.col("strand").first(),
        ])
    )
    
    exon_df = (
        exon_df.group_by("tran_id")
        .agg([
            pl.col("chr").first(),
            pl.col("start").sort(),
            pl.col("stop").sort(),
            pl.col("strand").first(),
        ])
    )
    
    # Calculate transcript coordinates
    exon_df = exon_df.with_columns([
        pl.col("start").apply(lambda x: list(range(len(x)))).alias("tran_start"),
        pl.col("stop").apply(lambda x: list(range(len(x)))).alias("tran_stop"),
    ])
    
    log_info(f"Found {len(cds_df)} CDS and {len(exon_df)} exon features")
    return cds_df, exon_df


def get_bam_tran(bam_df, exon_df):
    """
    Merge BAM and exon DataFrames based on chromosome alignment and calculate transcript coordinates in BAM.

    Parameters:
    - bam_df (DataFrame): DataFrame containing BAM file data with chromosome positions.
    - exon_df (DataFrame): DataFrame containing exon annotations with chromosome positions.

    Returns:
    - dict: A dictionary where keys are column names and values are lists of corresponding values,
            representing transcript coordinates aligned with BAM data.
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


def calculate_differences(start, start_dict):
    """
    Calculate the difference between a given start value and each value in a dictionary of start values.

    Parameters:
    - start (int or float): The starting value for which differences are calculated.
    - start_dict (int): A values extracted from a dictionary containing start values (int or float).

    Returns:
    - dict: A dictionary where keys correspond to identifiers from `start_dict` and values represent
            the differences between `start` and each corresponding value in `start_dict`.
    """
    # Handle polars Series objects
    if isinstance(start, pl.Series):
        start = start[0]  # Take the first value if it's a Series
    start_diff = start - start_dict

    return start_diff


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
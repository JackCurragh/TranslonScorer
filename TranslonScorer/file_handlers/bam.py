"""
BAM file handling functionality for TranslonScorer.

This module contains functions for processing BAM files and converting
them to other formats, including coordinate transformations.
"""

import pysam
import polars as pl
import oxbow as ox
from ..utils.logging import log_info, log_error, log_warning
import os


def readbam(bampath):
    """
    Reads a given BAM file, extracts relevant information, and returns it as a DataFrame.

    Parameters:
    - bampath (str): Path to the BAM file to be processed.

    Returns:
    - df (DataFrame): Polars DataFrame containing the extracted information from the BAM file.
                     Contains columns: chr, start, stop, length, strand, count
    """
    # Check if index exists
    if not (os.path.exists(f"{bampath}.bai") or os.path.exists(bampath.replace(".bam", ".bai"))):
        log_info("BAM index not found, creating index...")
        pysam.index(bampath)
        log_info("BAM file indexed successfully")
    else:
        log_info("Using existing BAM index")
    
    # Read BAM file using oxbow for speed
    bamfile = ox.read_bam(bampath)
    log_info("BAM file read successfully")
    
    # Convert to DataFrame
    df = pl.read_ipc(bamfile)
    log_info(f"Available columns: {df.columns}")
    
    # Map oxbow column names to our expected names
    column_mapping = {
        'rname': 'chr',
        'pos': 'start',
        'end': 'stop'
    }
    
    # Rename columns that exist
    for old_name, new_name in column_mapping.items():
        if old_name in df.columns:
            df = df.rename({old_name: new_name})
    
    # Calculate length from sequence
    df = df.with_columns(
        pl.col('seq').str.lengths().alias('length')
    )
    
    # Add strand based on SAM flag (0x10 is the reverse strand bit)
    df = df.with_columns(
        pl.when(pl.col('flag') & 0x10 > 0)
        .then(pl.lit('-'))
        .otherwise(pl.lit('+'))
        .alias('strand')
    ).drop('flag')
    
    # Keep only needed columns and add count
    df = df.select(['chr', 'start', 'stop', 'length', 'strand']).with_columns(count=pl.lit(1))
        
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
            # Take first start position (5' most for + strand, 3' most for - strand)
            pl.col("start").sort().first().alias("start"),
            pl.col("stop").sort().first().alias("stop"),
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
    Map BAM reads to transcript coordinates using exon information.
    Uses interval-based filtering to efficiently find overlaps.

    Parameters:
    - bam_df (DataFrame): DataFrame containing BAM file data with chromosome positions
    - exon_df (DataFrame): DataFrame containing exon annotations with chromosome positions

    Returns:
    - DataFrame: BAM data with added transcript coordinates
    """
    # Ensure chromosome types match (both categorical)
    if exon_df["chr"].dtype != pl.Categorical:
        exon_df = exon_df.with_columns(pl.col("chr").cast(pl.Categorical))
    
    # Sort both DataFrames by start position for efficient overlap checking
    bam_df = bam_df.sort("start")
    exon_df = exon_df.sort("start")
    
    # Create overlapping windows for efficient filtering
    window_size = 10000  # 10kb windows
    results = []
    
    # Process in windows to avoid memory issues
    for chr in bam_df["chr"].unique():
        bam_with_chr = bam_df.filter(pl.col("chr") == chr)
        exon_with_chr = exon_df.filter(pl.col("chr") == chr)
        
        if exon_with_chr.is_empty() or bam_with_chr.is_empty():
            continue
            
        chr_min = bam_with_chr["start"].min()
        chr_max = bam_with_chr["stop"].max()
        
        for start in range(chr_min, chr_max, window_size):
            end = start + window_size
            
            # Get reads and exons in this window
            window_reads = bam_with_chr.filter(
                (pl.col("start") >= start) & (pl.col("start") < end)
            )
            if window_reads.is_empty():
                continue
                
            window_exons = exon_with_chr.filter(
                (pl.col("start") < end) & (pl.col("stop") > start)
            )
            if window_exons.is_empty():
                continue
            
            # Find overlaps within the window using a more efficient join
            overlaps = window_reads.join(
                window_exons,
                on="chr",
                how="inner"
            ).filter(
                (pl.col("start") >= pl.col("start_right")) &
                (pl.col("stop") <= pl.col("stop_right"))
            )
            
            if not overlaps.is_empty():
                # Calculate transcript coordinates
                mapped = overlaps.with_columns([
                    (pl.col("tran_start") + (pl.col("start") - pl.col("start_right")))
                    .alias("tran_start_bam"),
                    (pl.col("tran_stop") - (pl.col("stop_right") - pl.col("stop")))
                    .alias("tran_stop_bam")
                ]).select(
                    pl.all().exclude("start_right", "stop_right", "tran_start", "tran_stop")
                )
                results.append(mapped)
                
            # Clear memory
            del window_reads
            del window_exons
            if 'overlaps' in locals():
                del overlaps
            if 'mapped' in locals():
                del mapped
    
    if not results:
        return pl.DataFrame()
    
    # Concatenate results in chunks to avoid memory issues
    chunk_size = 100
    final_results = []
    for i in range(0, len(results), chunk_size):
        chunk = pl.concat(results[i:i + chunk_size])
        chunk = chunk.unique(
            subset=['chr', 'start', 'stop', 'length', 'strand', 'count', 'tran_id'],
            maintain_order=True
        )
        final_results.append(chunk)
        del chunk
    
    return pl.concat(final_results)


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
    
    # Get unique chromosomes and check overlap
    uniquechr_bam = set(bam_df["chr"].unique().cast(pl.Utf8))
    uniquechr_exon = set(exon_df["chr"].unique().cast(pl.Utf8))
    
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
                # Convert to strings temporarily for the operation
                bam_df = bam_df.with_columns([
                    pl.col("chr").cast(pl.Utf8).map_elements(
                        lambda x: f"chr{x}" if not x.startswith("chr") else x
                    ).cast(pl.Categorical).alias("chr")
                ])
                uniquechr_bam = set(bam_df["chr"].unique().cast(pl.Utf8))
            else:
                log_info("Removing 'chr' prefix from BAM chromosomes")
                bam_df = bam_df.with_columns([
                    pl.col("chr").cast(pl.Utf8).map_elements(
                        lambda x: x.replace("chr", "")
                    ).cast(pl.Categorical).alias("chr")
                ])
                uniquechr_bam = set(bam_df["chr"].unique().cast(pl.Utf8))
            
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

    # Process matching chromosomes
    results = []
    total_chr = len(common_chr)
    log_info(f"\nProcessing {total_chr} chromosomes")
    
    for idx, chr in enumerate(sorted(common_chr), 1):
        if idx % 5 == 0 or idx == total_chr:
            log_info(f"Processing chromosome {chr} ({idx}/{total_chr})")
        
        # Filter BAM data for this chromosome
        bam_with_chr = bam_df.filter(pl.col("chr") == chr)
        if bam_with_chr.is_empty():
            continue
            
        # Get and flatten exon data for this chromosome
        exon_with_chr = (
            exon_df.filter(pl.col("chr") == chr)
            .explode(["start", "stop", "tran_start", "tran_stop"])
        )
        if exon_with_chr.is_empty():
            continue
        
        # Process this chromosome's data
        result_df = get_bam_tran(bam_with_chr, exon_with_chr)
        if not result_df.is_empty():
            results.append(result_df)
            
        # Clear memory
        del bam_with_chr
        del exon_with_chr
        if 'result_df' in locals():
            del result_df

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

    # Concatenate results in chunks to avoid memory issues
    chunk_size = 100
    final_results = []
    for i in range(0, len(results), chunk_size):
        chunk = pl.concat(results[i:i + chunk_size])
        final_results.append(chunk)
        del chunk
    
    log_info("\nCreating final BAM transcript DataFrame")
    final_df = pl.concat(final_results)
    log_info("BAM transcript conversion complete")
    return final_df


def process_transcriptomic_bam(df_namesplit, cds_df):
    """
    Process a transcriptomic BAM file or output from bamtranscript where reads are mapped to transcripts.
    
    Parameters:
    - df_namesplit (DataFrame): BAM DataFrame with transcript alignments, either:
        - From transcriptomic BAM: 'chr' column contains transcript IDs
        - From bamtranscript: has both 'chr' and 'tran_id' columns
    - cds_df (DataFrame): CDS DataFrame with transcript information
    
    Returns:
    - DataFrame: Processed BAM data ready for A-site calculation
    """
    # Check if this is output from bamtranscript (has both chr and tran_id)
    if "tran_id" in df_namesplit.columns:
        df_with_tran = df_namesplit
    else:
        # For direct transcriptomic BAM, rename chr to tran_id
        df_with_tran = df_namesplit.rename({"chr": "tran_id"})
    
    # Filter for transcripts present in CDS annotation
    bam_to_cds = (
        df_with_tran.with_columns(shared=pl.col("tran_id").is_in(cds_df["tran_id"]))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )
    
    # Join with CDS DataFrame to get start positions
    bam_to_cds = bam_to_cds.join(
        cds_df.select(["tran_id", "start"]).rename({"start": "cds_start"}),
        on="tran_id",
        how="left"
    )
    
    # Calculate position relative to CDS start
    bam_to_cds = bam_to_cds.with_columns(
        bamcds_start=pl.col("start") - pl.col("cds_start")
    ).drop("cds_start")
    
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
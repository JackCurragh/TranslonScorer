"""This script contains functions to transform data frames into different file types"""

import polars as pl
import os
import logging
from typing import Optional

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create console handler with formatting
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# Global flag to control logging
ENABLE_LOGGING = False

def enable_logging(enable: bool = True, log_file: Optional[str] = None) -> None:
    """
    Enable or disable logging for the fileprocessor module.
    
    Parameters:
    - enable (bool): Whether to enable logging
    - log_file (str, optional): Path to log file. If provided, logs will be written to this file
                               in addition to console output.
    """
    global ENABLE_LOGGING
    ENABLE_LOGGING = enable
    
    if enable:
        logger.setLevel(logging.INFO)
        if log_file:
            # Add file handler if log file is specified
            fh = logging.FileHandler(log_file)
            fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(fh)
    else:
        logger.setLevel(logging.WARNING)

def log_info(message: str) -> None:
    """Helper function to log messages only if logging is enabled."""
    if ENABLE_LOGGING:
        logger.info(message)

def log_warning(message: str) -> None:
    """Helper function to log warnings only if logging is enabled."""
    if ENABLE_LOGGING:
        logger.warning(message)

def log_error(message: str) -> None:
    """Helper function to log errors only if logging is enabled."""
    if ENABLE_LOGGING:
        logger.error(message)

def normalize_chrom_name(chrom):
    """
    Normalize chromosome names by removing 'chr' prefix and any alternative/patch suffixes.
    
    Parameters:
    - chrom (str): Chromosome name to normalize
    
    Returns:
    - str: Normalized chromosome name
    """
    # Remove 'chr' prefix if present
    norm_chrom = chrom.lower().replace('chr', '')
    
    # Remove any alternative/patch suffixes (e.g., _alt, _fix, _random)
    norm_chrom = norm_chrom.split('_')[0]
    
    return norm_chrom

def get_matching_chromosomes(bam_ids, exon_chroms):
    """
    Find matching chromosomes between BAM and annotation, handling naming variations.
    
    Parameters:
    - bam_ids (set): Set of chromosome IDs from BAM
    - exon_chroms (set): Set of chromosome IDs from annotation
    
    Returns:
    - tuple: (bool for match found, bool for normalization needed, set of common chromosomes)
    """
    # Try direct matching first
    common_chr = bam_ids.intersection(exon_chroms)
    if common_chr:
        return True, False, common_chr
    
    # Try normalized matching
    norm_bam_chroms = {normalize_chrom_name(c) for c in bam_ids}
    norm_exon_chroms = {normalize_chrom_name(c) for c in exon_chroms}
    
    common_chr = norm_bam_chroms.intersection(norm_exon_chroms)
    if common_chr:
        # Map normalized chromosomes back to original exon names
        chrom_map = {normalize_chrom_name(c): c for c in exon_chroms}
        original_common = {chrom_map[c] for c in common_chr}
        return True, True, original_common
        
    return False, False, set()

def detect_bam_type(df, exon_df):
    """
    Detect whether a BAM file is genomic or transcriptomic by checking chromosome/transcript ID patterns.
    Handles 'chr' prefix variations in chromosome names.
    
    Parameters:
    - df (DataFrame): BAM DataFrame with chromosome/transcript information
    - exon_df (DataFrame): Exon DataFrame with both chromosome and transcript IDs
    
    Returns:
    - tuple: (str for BAM type, dict with match info)
    """
    bam_ids = set(df["chr"].unique())
    exon_chroms = set(exon_df["chr"].unique())
    exon_trans = set(exon_df["tran_id"].unique())
    
    # First check for transcript matches
    trans_matches = bam_ids.intersection(exon_trans)
    if trans_matches:
        return 'transcriptomic', {'needs_normalization': False}
    
    # Check for direct chromosome matches
    chrom_matches = bam_ids.intersection(exon_chroms)
    if chrom_matches:
        return 'genomic', {'needs_normalization': False}
    
    # Try matching after removing 'chr' prefix
    bam_no_chr = {c.replace('chr', '') for c in bam_ids}
    exon_no_chr = {c.replace('chr', '') for c in exon_chroms}
    
    if bam_no_chr.intersection(exon_no_chr):
        return 'genomic', {'needs_normalization': True}
    
    # If still no matches, show helpful error
    raise ValueError(
        "Unable to determine BAM type. No matches found between:\n"
        f"BAM chromosomes: {sorted(bam_ids)}\n"
        f"Annotation chromosomes: {sorted(exon_chroms)}\n"
        f"Annotation transcripts: {sorted(exon_trans)}\n"
        "Note: Tried matching with and without 'chr' prefix"
    )

def process_genomic_bam(df_namesplit, exon_df, needs_normalization=False):
    """
    Process a genomic BAM file, handling chromosome name normalization if needed.
    
    Parameters:
    - df_namesplit (DataFrame): BAM DataFrame with genomic alignments
    - exon_df (DataFrame): Exon DataFrame with chromosome information
    - needs_normalization (bool): Whether chromosome names need to be normalized
    
    Returns:
    - DataFrame: Processed BAM data
    """
    if needs_normalization:
        # Normalize chromosome names in both DataFrames
        df_namesplit = df_namesplit.with_columns(
            pl.col("chr").apply(normalize_chrom_name).alias("chr")
        )
        exon_df = exon_df.with_columns(
            pl.col("chr").apply(normalize_chrom_name).alias("chr")
        )
    
    # Filter to matching chromosomes
    bam_ids = set(df_namesplit["chr"].unique())
    exon_chroms = set(exon_df["chr"].unique())
    _, _, common_chr = get_matching_chromosomes(bam_ids, exon_chroms)
    
    df_namesplit = (
        df_namesplit.with_columns(shared=pl.col("chr").is_in(common_chr))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )
    
    exon_df = (
        exon_df.with_columns(shared=pl.col("chr").is_in(common_chr))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )
    
    return bamtranscript(df_namesplit, exon_df)

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

def bamtranscript(bam_df, exon_df):
    """
    Filter BAM and exon DataFrames based on shared chromosome information and flatten exon annotations.
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
                    .otherwise(pl.concat_str(["chr", pl.col("chr")]))
                    .alias("chr")
                )
                uniquechr_bam = set(bam_df["chr"].unique())
            else:
                log_info("Removing 'chr' prefix from BAM chromosomes")
                bam_df = bam_df.with_columns(
                    pl.col("chr").str.replace("chr", "").alias("chr")
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
        bam_tran = process_genomic_bam(df_namesplit, exon_df, info['needs_normalization'])
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
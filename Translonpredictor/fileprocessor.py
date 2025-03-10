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

def detect_bam_type(df, exon_df):
    """
    Detect whether a BAM file is genomic or transcriptomic by checking chromosome/transcript ID patterns.
    Handles variations in chromosome naming (with/without 'chr' prefix).
    
    Parameters:
    - df (DataFrame): BAM DataFrame with chromosome/transcript information
    - exon_df (DataFrame): Exon DataFrame with both chromosome and transcript IDs
    
    Returns:
    - str: 'genomic' or 'transcriptomic'
    - dict: Additional information about the match (e.g., whether chromosome names need normalization)
    
    Raises:
    - ValueError: If BAM type cannot be determined or if no matching IDs found
    """
    bam_ids = set(df["chr"].unique())
    exon_chroms = set(exon_df["chr"].unique())
    exon_trans = set(exon_df["tran_id"].unique())
    
    # Try direct chromosome matching first
    chrom_match = len(bam_ids.intersection(exon_chroms))
    trans_match = len(bam_ids.intersection(exon_trans))
    
    if chrom_match > trans_match:
        return 'genomic', {'needs_normalization': False}
    elif trans_match > chrom_match:
        return 'transcriptomic', {'needs_normalization': False}
    
    # If no direct matches, try normalized chromosome names
    norm_bam_chroms = {normalize_chrom_name(c) for c in bam_ids}
    norm_exon_chroms = {normalize_chrom_name(c) for c in exon_chroms}
    
    norm_chrom_match = len(norm_bam_chroms.intersection(norm_exon_chroms))
    
    if norm_chrom_match > trans_match:
        return 'genomic', {'needs_normalization': True}
    elif trans_match > norm_chrom_match:
        return 'transcriptomic', {'needs_normalization': False}
    else:
        raise ValueError(
            "Unable to determine BAM type. No significant matches found with either:\n"
            f"Chromosomes in annotation (original): {sorted(exon_chroms)}\n"
            f"Chromosomes in annotation (normalized): {sorted(norm_exon_chroms)}\n"
            f"Transcript IDs in annotation: {sorted(exon_trans)}\n"
            f"IDs in BAM (original): {sorted(bam_ids)}\n"
            f"IDs in BAM (normalized): {sorted(norm_bam_chroms)}"
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
    
    return bamtranscript(df_namesplit, exon_df)

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
    df_filtered = df.with_columns(
        (pl.col("end") - pl.col("pos")).alias("length")
    ).select(
        pl.all().exclude(
            "flag", "qual", "tags", "mapq", "tlen", "seq", "pnext", "rnext", "cigar"
        )
    )

    # Split read names to get counts
    split_func = lambda s: int(s.split("_x")[1])
    df_namesplit = (
        df_filtered.with_columns(pl.col("qname").apply(split_func))
        .rename({"qname": "count", "rname": "chr", "pos": "start", "end": "stop"})
        .cast({"chr": pl.String})
    )

    # Get annotations
    cds_df, exon_df = getexons_and_cds(annotation)
    
    # Detect BAM type
    bam_type, info = detect_bam_type(df_namesplit, exon_df)
    
    # Process based on BAM type
    if bam_type == 'genomic':
        # Process genomic BAM with chromosome name handling
        bam_tran = process_genomic_bam(df_namesplit, exon_df, info['needs_normalization'])
        bam_to_cds = bamrelativetocds(bam_tran, cds_df)
    else:  # transcriptomic
        bam_to_cds = process_transcriptomic_bam(df_namesplit, cds_df)
    
    # Calculate offsets if not provided
    if not offsets:
        bam_offsets = bam_to_cds.group_by("bamcds_start", "length").agg(
            pl.col("count").sum()
        )
        offsets = change_point_analysis(bam_offsets)
    
    # A-site calculation
    bed = asitecalc(bam_to_cds, offsets)

    return bed, exon_df, cds_df 
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
    
    # Detect BAM type and get matching info
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

def bamtranscript(bam_df, exon_df):
    """
    Filter BAM and exon DataFrames based on shared chromosome information and flatten exon annotations.

    Parameters:
    - bam_df (DataFrame): DataFrame containing BAM file data with chromosome positions
    - exon_df (DataFrame): DataFrame containing exon annotations with chromosome positions

    Returns:
    - DataFrame: A modified version of bam_df after filtering and coordinate conversion
    """
    # Get unique chromosomes from both datasets
    bam_chroms = set(bam_df["chr"].unique())
    exon_chroms = set(exon_df["chr"].unique())
    
    # Check if we need to add/remove 'chr' prefix
    bam_has_chr = any(c.startswith('chr') for c in bam_chroms)
    exon_has_chr = any(c.startswith('chr') for c in exon_chroms)
    
    if bam_has_chr != exon_has_chr:
        if exon_has_chr:
            # Add 'chr' prefix to BAM chromosomes
            bam_df = bam_df.with_columns(
                pl.when(pl.col("chr").str.starts_with("chr"))
                .then(pl.col("chr"))
                .otherwise(pl.concat_str(["chr", pl.col("chr")]))
                .alias("chr")
            )
        else:
            # Remove 'chr' prefix from BAM chromosomes
            bam_df = bam_df.with_columns(
                pl.col("chr").str.replace("chr", "").alias("chr")
            )
    
    # Now filter for matching chromosomes
    bam_df = (
        bam_df.with_columns(shared=pl.col("chr").is_in(exon_chroms))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )
    
    if bam_df.is_empty():
        raise ValueError(
            "No overlapping chromosomes found between BAM and annotation after name normalization.\n"
            f"BAM chromosomes: {sorted(bam_chroms)}\n"
            f"Annotation chromosomes: {sorted(exon_chroms)}"
        )
    
    # Process matching chromosomes
    results = []
    for chr in set(bam_df["chr"].unique()):
        bam_with_chr = bam_df.filter(pl.col("chr") == chr)
        exon_with_chr = exon_df.filter(pl.col("chr") == chr)
        
        if exon_with_chr.is_empty() or bam_with_chr.is_empty():
            continue
            
        min_val = min(exon_with_chr["start"])
        max_val = max(exon_with_chr["stop"])
        
        for i in range(min_val, max_val, 225000):
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
        raise ValueError(
            "No overlapping regions found between BAM and annotation files.\n"
            "This could be due to:\n"
            "1. Mismatched coordinates\n"
            "2. BAM file aligned to different genome version than annotation\n"
            "3. No reads mapping to annotated regions"
        )

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
    return bam_df 
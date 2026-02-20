"""
BAM file handling functionality for TranslonScorer.

This module contains functions for processing BAM files and converting
them to other formats, including coordinate transformations.
"""

import pysam
import polars as pl
try:
    import oxbow as ox  # optional fast BAM reader
except Exception:
    ox = None
from ..utils.logging import log_info, log_error, log_warning
from typing import Optional
import os
import sys

# Only import memory_profiler if PROFILE environment variable is set
PROFILE = os.environ.get('PROFILE', '0') == '1'
if PROFILE:
    from memory_profiler import profile
else:
    def profile(precision=None):  # Dummy decorator when not profiling
        def wrapper(func):
            return func
        return wrapper


def readbam(
    bampath,
    *,
    collapsed: bool = False,
    count_from: Optional[str] = None,
    count_pattern: Optional[str] = None,
    count_tag: Optional[str] = None,
    unique: bool = False,
    include_qname: bool = False,
):
    """
    Read a BAM file into a normalized table.

    Args:
        bampath: Path to BAM file.
        collapsed: Treat BAM as collapsed (counts per read encoded in name or tag).
        count_from: One of {None, 'name', 'tag'}. If collapsed, where to parse count.
        count_pattern: Regex with named group 'count' when count_from='name', e.g. r'.*_x(?P<count>\\d+)$'.
        count_tag: SAM tag name when count_from='tag', e.g. 'RC'.
        unique: If True, do not group/aggregate; return one row per alignment.
        include_qname: If True, include 'qname' column (useful for Zarr linking).

    Returns:
        Polars DataFrame with columns: chr, start, stop, length, strand, count[, qname]
    """
    # Check if index exists
    if not (os.path.exists(f"{bampath}.bai") or os.path.exists(bampath.replace(".bam", ".bai"))):
        log_info("BAM index not found, creating index...")
        pysam.index(bampath)
        log_info("BAM file indexed successfully")
    else:
        log_info("Using existing BAM index")
    
    # Read BAM
    if ox is not None:
        bamfile = ox.read_bam(bampath)
        log_info("BAM file read via oxbow")
        df = pl.read_ipc(bamfile)
        log_info(f"Available columns: {df.columns}")
    else:
        log_info("oxbow not available; falling back to pysam (slower)")
        rows = []
        with pysam.AlignmentFile(bampath, 'rb') as bam:
            for aln in bam.fetch(until_eof=True):
                if aln.is_unmapped:
                    continue
                rows.append({
                    'rname': bam.get_reference_name(aln.reference_id),
                    'pos': int(aln.reference_start),
                    'end': int(aln.reference_end),
                    'seq': aln.query_sequence or '',
                    'flag': int(aln.flag),
                    'qname': aln.query_name,
                })
        if not rows:
            return pl.DataFrame({
                'chr': [], 'start': [], 'stop': [], 'length': [], 'strand': [], 'count': []
            })
        df = pl.from_dicts(rows)
    
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
    
    # Keep only needed columns; optionally keep qname
    keep_cols = ['chr', 'start', 'stop', 'length', 'strand']
    if include_qname and 'qname' in df.columns:
        keep_cols.append('qname')
    df = df.select(keep_cols)

    # Add count column
    if collapsed:
        if count_from == 'name':
            import re
            pattern = re.compile(count_pattern or r'.*_x(?P<count>\d+)$')
            if 'qname' not in df.columns:
                log_warning("Collapsed BAM count_from=name requested but qname not available; defaulting counts to 1")
                df = df.with_columns(count=pl.lit(1))
            else:
                df = df.with_columns(
                    pl.col('qname')
                    .map_elements(lambda s: int(pattern.match(s).group('count')) if pattern.match(s) else 1)
                    .alias('count')
                )
        elif count_from == 'tag':
            # Fallback to pysam one-pass to fetch tag efficiently if oxbow didn't include it
            if not count_tag:
                log_warning("Collapsed BAM count_from=tag requested but count_tag not provided; defaulting counts to 1")
                df = df.with_columns(count=pl.lit(1))
            else:
                # Build qname -> count map using pysam
                try:
                    import pysam as _pysam
                    q2c = {}
                    with _pysam.AlignmentFile(bampath, 'rb') as bam:
                        for aln in bam.fetch(until_eof=True):
                            try:
                                q2c[aln.query_name] = int(aln.get_tag(count_tag))
                            except Exception:
                                q2c[aln.query_name] = 1
                    if 'qname' in df.columns:
                        df = df.with_columns(
                            pl.col('qname').map_elements(lambda s: q2c.get(s, 1)).alias('count')
                        )
                    else:
                        log_warning("qname not present; cannot apply tag-derived counts. Using count=1")
                        df = df.with_columns(count=pl.lit(1))
                except Exception as e:
                    log_warning(f"Failed to parse tag-derived counts: {e}; using count=1")
                    df = df.with_columns(count=pl.lit(1))
        else:
            df = df.with_columns(count=pl.lit(1))
    else:
        df = df.with_columns(count=pl.lit(1))
        
    # Aggregate unless unique requested
    if not unique:
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
        # Polars 1.36 uses comment_prefix instead of comment_char
        comment_prefix="#",
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

    # Calculate transcript-space exon offsets (strand-aware cumulative coordinates)
    # For each transcript, produce per-exon transcript starts and stops where
    # tran_start[i] = sum_{k < i} (stop[k]-start[k])
    # tran_stop[i]  = tran_start[i] + (stop[i]-start[i])
    def _tran_coords(starts: list[int], stops: list[int]) -> tuple[list[int], list[int]]:
        lens = [int(b) - int(a) for a, b in zip(starts, stops)]
        tran_starts = []
        acc = 0
        for L in lens:
            tran_starts.append(acc)
            acc += L
        tran_stops = [ts + L for ts, L in zip(tran_starts, lens)]
        return tran_starts, tran_stops

    exon_df = exon_df.with_columns([
        pl.struct(["start", "stop"]).map_elements(lambda s: _tran_coords(s["start"], s["stop"]))
        .alias("_tc")
    ])
    exon_df = exon_df.with_columns([
        pl.col("_tc").map_elements(lambda t: t[0]).alias("tran_start"),
        pl.col("_tc").map_elements(lambda t: t[1]).alias("tran_stop"),
    ]).drop("_tc")
    
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
    if bam_df["chr"].dtype != pl.Categorical:
        bam_df = bam_df.with_columns(pl.col("chr").cast(pl.Categorical))
    
    # Create overlapping windows for efficient filtering
    window_size = 225000  # 225kb windows
    results = []
    
    # Process in windows to avoid memory issues
    chr_min = bam_df["start"].min()
    chr_max = bam_df["stop"].max()
    
    for start in range(chr_min, chr_max, window_size):
        end = start + window_size
        
        # Get reads and exons in this window
        window_reads = bam_df.filter(
            (pl.col("start") >= start) & (pl.col("start") < end)
        )
        if window_reads.is_empty():
            continue
            
        window_exons = exon_df.filter(
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
    
    # Concatenate results
    return pl.concat(results).unique(
        subset=['chr', 'start', 'stop', 'length', 'strand', 'count', 'tran_id'],
        maintain_order=True
    )


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
                # Use Polars' native concatenation and ensure categorical type
                bam_df = bam_df.with_columns([
                    pl.when(pl.col("chr").cast(pl.Utf8).str.contains("^chr"))
                    .then(pl.col("chr"))
                    .otherwise(pl.concat_list([
                        pl.Series("chr", ["chr"]),
                        pl.col("chr")
                    ]).list.join(""))
                    .cast(pl.Categorical)
                    .alias("chr")
                ])
            else:
                log_info("Removing 'chr' prefix from BAM chromosomes")
                bam_df = bam_df.with_columns([
                    pl.col("chr").cast(pl.Utf8).str.slice(3).cast(pl.Categorical).alias("chr")
                ])
            
            # Update unique chromosomes after normalization
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

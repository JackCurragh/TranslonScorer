"""
BED file handling functionality for TranslonScorer.

This module contains functions for processing BED files and converting
them to other formats, including coordinate transformations.
"""

import polars as pl
import pyBigWig as bw
from ..utils.logging import log_info, log_warning, log_error


def saveorfsandexons(orf_df, exon_df, filename):
    """
    Save annotated ORFs and exons DataFrames to CSV files.

    Args:
        orf_df (DataFrame): DataFrame containing annotated ORFs data
        exon_df (DataFrame): DataFrame containing exon data
        filename (str): Base name for output files

    Returns:
        tuple: Paths to created ORF and exon CSV files
    """
    orf_df.write_csv(f"{filename}_annotated_orfs.csv")

    # Keep original chromosome names and concatenate coordinate columns
    exon_df = exon_df.with_columns(
        pl.col("start", "stop", "tran_start", "tran_stop").apply(
            lambda x: ",".join(map(str, x))
        )
    )
    exon_df.write_csv(f"{filename}_exons.csv")
    return f"{filename}_annotated_orfs.csv", f"{filename}_exons.csv"


def asitecalc(df, offsets):
    """
    Calculates A-site positions and aggregates counts based on offset values.

    Parameters:
    - df (DataFrame): Input DataFrame containing 'length' and 'pos' columns for A-site calculation.
    - offsets (dict): Dictionary containing offset values for each 'length' value.

    Returns:
    - df_bed (DataFrame): DataFrame containing aggregated information of A-site positions, their counts, and chromosome information.
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


def bedtobigwig(bedfile, chromsize, filename):
    """
    Converts a bedGraph file to a bigWig file using pyBigWig.

    Parameters:
        bedfile (str): The path to the input bedGraph file.
        chromsize (str): The path to the chromosome sizes file.
        filename (str): The name for the generated file.

    Returns:
        str: Path to the created bigWig file
    """
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
    
    # Get unique chromosomes from bedGraph
    bed_chroms = set(bed_data["chrom"].unique())
    bed_has_chr = any(c.startswith('chr') for c in bed_chroms)
    
    log_info("Reading chromosome sizes")
    # Read chromosome sizes file and normalize to match bedGraph format
    chrom_sizes = {}
    with open(chromsize, 'r') as f:
        for line in f:
            chrom, size = line.strip().split('\t')
            # Normalize chromosome name to match bedGraph format
            if bed_has_chr and not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            elif not bed_has_chr and chrom.startswith('chr'):
                chrom = chrom[3:]  # Remove 'chr' prefix
            chrom_sizes[chrom] = int(size)
    
    log_info(f"Found {len(chrom_sizes)} chromosomes in sizes file")
    
    # Filter out chromosomes that aren't in the chromosome sizes file
    bed_data = bed_data.filter(pl.col("chrom").is_in(chrom_sizes.keys()))
    
    if bed_data.is_empty():
        error_msg = (
            "No matching chromosomes found between bedGraph and chromosome sizes after normalization.\n"
            f"BedGraph chromosomes: {sorted(bed_chroms)[:5]}\n"
            f"Chromosome sizes chromosomes: {sorted(chrom_sizes.keys())[:5]}"
        )
        log_error(error_msg)
        raise ValueError(error_msg)
    
    # Ensure all positions are valid
    bed_data = bed_data.filter(pl.col("end") > pl.col("start"))
    
    # Remove any rows with NaN values
    bed_data = bed_data.drop_nulls()
    
    # Sort the data by chromosome and start position
    bed_data = bed_data.sort(["chrom", "start"])
    
    log_info(f"Processed bedGraph data: {len(bed_data)} entries")
    
    # Create bigWig file
    log_info("Creating bigWig file")
    bw_file = bw.open(f"{filename}.bw", "w")
    
    # Add header with chromosome sizes
    bw_file.addHeader(list(chrom_sizes.items()))
    
    # Process chromosomes in a specific order (match chromosome sizes order)
    chroms_to_process = [chrom for chrom in chrom_sizes.keys() if chrom in bed_data["chrom"].unique()]
    if not chroms_to_process:
        log_error("No matched chromosomes to process - check that chromsizes matches BAM/bedGraph")
    
    for chrom in chroms_to_process:
        chrom_data = bed_data.filter(pl.col("chrom") == chrom)
        if len(chrom_data) == 0:
            continue
            
        log_info(f"Processing chromosome {chrom} with {len(chrom_data)} entries")
        
        try:
            # Convert polars Series to lists
            starts = chrom_data["start"].to_list()
            ends = chrom_data["end"].to_list()
            values = chrom_data["value"].to_list()
            
            # Validate positions against chromosome size
            max_pos = max(ends)
            if max_pos > chrom_sizes[chrom]:
                log_warning(f"Trimming entries exceeding chromosome {chrom} size ({max_pos} > {chrom_sizes[chrom]})")
                
                # Filter entries to be within chromosome size
                valid_entries = [(s, e, v) for s, e, v in zip(starts, ends, values) if e <= chrom_sizes[chrom]]
                if not valid_entries:
                    log_warning(f"No valid entries for chromosome {chrom} after trimming")
                    continue
                    
                starts, ends, values = zip(*valid_entries)
            
            # Add entries chromosome by chromosome
            bw_file.addEntries(
                [chrom] * len(starts),
                starts,
                ends=ends,
                values=values
            )
            
            log_info(f"Successfully added {len(starts)} entries for chromosome {chrom}")
            
        except Exception as e:
            log_error(f"Error processing chromosome {chrom}: {str(e)}")
            # Continue with next chromosome instead of terminating the whole process
            continue
    
    bw_file.close()
    log_info(f"Successfully created {filename}.bw")
    return f"{filename}.bw" 
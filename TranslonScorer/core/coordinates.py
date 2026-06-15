"""
Coordinate transformation functionality for TranslonScorer.

This module provides functions for coordinate transformations and ORF classification
relative to transcript features.
"""

import polars as pl
from pyfaidx import Fasta

from ..file_handlers.bam import getexons_and_cds
from ..utils.logging import log_info, log_warning


def change_point_analysis(offset_df):
    """
    Calculate the change point for the metagene profile using vectorized operations.

    Args:
        offset_df (DataFrame): DataFrame containing length and position information

    Returns:
        dict: Dictionary mapping read lengths to their optimal offsets
    """
    # Get unique lengths and process in chunks
    unique_lengths = offset_df["length"].unique().sort()
    total_lengths = len(unique_lengths)
    log_info(f"Processing offset analysis for {total_lengths} different read lengths")

    offset_dict = {}
    chunk_size = 10  # Process 10 lengths at a time

    for chunk_start in range(0, total_lengths, chunk_size):
        chunk_end = min(chunk_start + chunk_size, total_lengths)
        chunk_lengths = unique_lengths[chunk_start:chunk_end]

        # Filter DataFrame for current chunk of lengths
        chunk_df = offset_df.filter(pl.col("length").is_in(chunk_lengths))

        for length in chunk_lengths:
            # Get data for this length and sort
            offset_df_len = chunk_df.filter(pl.col("length") == length).sort("bamcds_start")

            if offset_df_len.is_empty():
                log_warning(f"No data found for length {length}, using default offset of 15")
                offset_dict[length] = 15  # default offset
                continue

            # Pre-calculate all positions we need
            positions = list(range(-30, 11))
            max_shift = 0
            max_shift_position = None

            # Create a lookup dictionary for counts at each position
            count_dict = dict(
                zip(offset_df_len["bamcds_start"].to_list(), offset_df_len["count"].to_list())
            )

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

            # Clean up memory
            del offset_df_len
            del count_dict

        # Log progress
        if chunk_end % 10 == 0 or chunk_end == total_lengths:
            log_info(f"Processed lengths up to {chunk_end}/{total_lengths}")

        # Clean up chunk memory
        del chunk_df

    log_info("Offset analysis complete")
    return offset_dict


def classify_orf(row):
    """
    Classify an ORF based on its relative position to transcript start and stop sites.

    Args:
        row (Series or dict-like): Row containing ORF information with 'start', 'stop',
                                  'tran_start', and 'tran_stop' values

    Returns:
        str: ORF classification:
            - "uORF": Upstream ORF
            - "CDS": Coding Sequence
            - "dORF": Downstream ORF
            - "uoORF": Upstream Overlapping ORF
            - "doORF": Downstream Overlapping ORF
            - "iORF": Internal ORF
            - "eoORF": Encapsulated Overlapping ORF
            - "extORF": Extended ORF
    """
    start = row["start"]
    stop = row["stop"]
    tran_start = row["tran_start"]
    tran_stop = row["tran_stop"]

    if stop < tran_start:
        return "uORF"
    elif start == tran_start and stop == tran_stop:
        return "CDS"
    elif start > tran_stop:
        return "dORF"
    elif start < tran_start and stop == tran_stop:
        return "extORF"
    elif start < tran_start and stop >= tran_start and stop <= tran_stop:
        return "uoORF"
    elif start >= tran_start and start <= tran_stop and stop > tran_stop:
        return "doORF"
    elif start >= tran_start and stop <= tran_stop:
        return "iORF"
    elif start < tran_start and stop > tran_stop:
        return "eoORF"
    else:
        log_warning(
            f"Unexpected ORF coordinates: start={start}, stop={stop}, "
            f"tran_start={tran_start}, tran_stop={tran_stop}"
        )
        return "Unexpected"


def orfrelativeposition(annotation, df, cds_df=None):
    """
    Determine relative positions of ORFs to coding sequences.

    Args:
        annotation (str): Path to genome annotation file
        df (DataFrame): DataFrame containing ORF coordinates
        cds_df (DataFrame, optional): Pre-loaded CDS DataFrame

    Returns:
        tuple: (orf_df, exon_df) where:
            - orf_df: DataFrame with ORFs and their classifications
            - exon_df: DataFrame with exon coordinates
    """
    log_info("Determining ORF positions relative to CDS...")

    # Get CDS and exon data if not provided
    if cds_df is None:
        cds_df, exon_df = getexons_and_cds(annotation, list(df["tran_id"].unique()))

    # Process coding transcripts
    tranids = list(cds_df["tran_id"].unique())

    # Classify coding ORFs
    codingorfs = (
        df.with_columns(shared=pl.col("tran_id").is_in(tranids))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )

    codingorfs = codingorfs.join(cds_df, on="tran_id")
    codingorfs = codingorfs.with_columns(
        pl.struct(["start", "stop", "tran_start", "tran_stop"])
        .apply(lambda row: classify_orf(row))
        .alias("type")
    ).select(pl.all().exclude("tran_start", "tran_stop"))

    # Process non-coding ORFs
    noncodingorfs = (
        df.with_columns(shared=pl.col("tran_id").is_in(tranids))
        .filter(pl.col("shared") == False)
        .select(pl.all().exclude("shared"))
        .with_columns(type=pl.lit("Non Coding"))
    )

    # Combine results
    orf_df = pl.concat([codingorfs, noncodingorfs])

    log_info(f"Classified {len(orf_df)} ORFs")
    return orf_df, exon_df


def gettranscripts(seq, annotation, outfilename):
    """
    Extracts transcript sequences from a genome annotation file.

    This function takes a genome sequence file and a genome annotation file,
    extracts transcript sequences from the annotation file, and saves them
    to a FASTA file.

    Parameters:
        seq (str): Path to the genome sequence file in FASTA format.
        annotation (str): Path to the genome annotation file in BED/GFF/GTF format.
        outfile (str): Path to save the output transcript sequences in FASTA format.
                       Default is 'transcripts.fa'.

    Returns:
        str: Path to the output FASTA file containing transcript sequences.

    Example:
        output_file = gettranscripts("genome.fa", "annotation.gff", outfile="transcripts.fa")
    """
    # Import pyranges lazily to avoid hard dependency at module import time
    import pyranges as pr

    # First, get available chromosomes from the genome.fa file
    log_info("Reading genome file to get available chromosomes")
    genome = Fasta(seq)
    available_chroms = set(genome.keys())
    log_info(f"Found {len(available_chroms)} chromosomes in genome file")

    # Read and filter annotation data
    log_info("Reading annotation file")
    ann = pr.read_gff(annotation, ignore_bad=True)
    transcripts = ann[ann.Feature == "exon"]

    # Get unique chromosomes from annotation
    ann_chroms = set(transcripts.Chromosome.unique())
    log_info(f"Found {len(ann_chroms)} chromosomes in annotation")

    # Find chromosomes that exist in both genome and annotation
    common_chroms = available_chroms.intersection(ann_chroms)
    if not common_chroms:
        # Try without 'chr' prefix
        ann_chroms_no_chr = {c.replace("chr", "") for c in ann_chroms}
        available_chroms_no_chr = {c.replace("chr", "") for c in available_chroms}
        common_chroms = available_chroms_no_chr.intersection(ann_chroms_no_chr)

        if not common_chroms:
            raise ValueError(
                "No matching chromosomes found between genome and annotation files.\n"
                f"Genome chromosomes: {sorted(list(available_chroms))[:5]}\n"
                f"Annotation chromosomes: {sorted(list(ann_chroms))[:5]}"
            )

    log_info(f"Found {len(common_chroms)} common chromosomes")

    # Filter transcripts to only include those on common chromosomes
    transcripts = transcripts[transcripts.Chromosome.isin(common_chroms)]

    # Get transcript sequences
    log_info("Extracting transcript sequences")
    tran_seq = transcripts.get_transcript_sequence(transcript_id="transcript_id", path=seq)

    # Write output
    log_info("Writing transcript sequences to file")
    with open(f"{outfilename}_transcripts.fa", "w") as fw:
        for index, id, seq in tran_seq.itertuples():
            fw.write(f">{id}\n{seq}\n")

    return f"{outfilename}_transcripts.fa"

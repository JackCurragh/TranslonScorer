"""
Coordinate transformation functionality for TranslonScorer.

This module contains functions for converting between different coordinate systems
(genomic, transcriptomic, CDS-relative) and handling coordinate-related calculations.

This module provides functions for coordinate transformations and ORF classification
relative to transcript features.
"""

import polars as pl
import ahocorasick
from Bio import SeqIO
from ..utils.logging import log_info, log_warning, log_error
from ..file_handlers.bam import getexons_and_cds


def change_point_analysis(offset_df):
    """
    Calculate the change point for the metagene profile using vectorized operations.
    
    Args:
        offset_df (DataFrame): DataFrame containing length and position information
        
    Returns:
        dict: Dictionary mapping read lengths to their optimal offsets
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
    elif start < tran_start and stop >= tran_start and stop <= tran_stop:
        return "uoORF"
    elif start >= tran_start and start <= tran_stop and stop > tran_stop:
        return "doORF"
    elif start >= tran_start and stop <= tran_stop:
        return "iORF"
    elif start < tran_start and stop > tran_stop:
        return "eoORF"
    elif start < tran_start and stop == tran_stop:
        return "extORF"
    else:
        log_warning(
            f"Unexpected ORF coordinates: start={start}, stop={stop}, "
            f"tran_start={tran_start}, tran_stop={tran_stop}"
        )
        return "Unexpected"


def find_all_positions(sequence, automaton):
    """
    Find positions of all occurrences of patterns from an Aho-Corasick automaton.

    Args:
        sequence (str): Input sequence to search for patterns
        automaton (ahocorasick.Automaton): Aho-Corasick automaton with patterns

    Returns:
        tuple: (frames dict with positions, codons dict with patterns)
    """
    frames = {0: [], 1: [], 2: []}
    codons = {}
    for i, (order, codon) in automaton.iter(sequence):
        # position of last nucleotide of codon returned
        frames[(i - 2) % 3].append(i)
        codons[i] = codon
    return frames, codons


def find_orfs(sequence, tran_id, startautomaton, stopautomaton, minlength=0, maxlength=1000000):
    """
    Predict Open Reading Frames in a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence
        tran_id (str): Transcript identifier
        startautomaton (ahocorasick.Automaton): Automaton for start codons
        stopautomaton (ahocorasick.Automaton): Automaton for stop codons
        minlength (int): Minimum ORF length
        maxlength (int): Maximum ORF length

    Returns:
        list: List of dictionaries containing ORF information
    """
    orf_list = []
    startpositions, start_codons = find_all_positions(sequence, startautomaton)
    stoppositions, stop_codons = find_all_positions(sequence, stopautomaton)

    for frame, startpositions in startpositions.items():
        for position in startpositions:
            valid_stops = [i for i in stoppositions[frame] if i > position]
            if valid_stops:
                stopposition = min(valid_stops)
                stopcodon = stop_codons[stopposition]
            else:
                stopposition = len(sequence)
                stopcodon = sequence[-3:]
            if stopcodon != "TAA" or stopcodon != "TAG" or stopcodon != "TGA":
                orf_data = {
                    "tran_id": tran_id,
                    "start": position - 2,
                    "stop": stopposition,
                    "length": stopposition - position,
                    "startorf": start_codons[position],
                    "stoporf": stopcodon,
                }
            else:
                orf_data = {
                    "tran_id": tran_id,
                    "start": position - 2,
                    "stop": stopposition - 3,
                    "length": stopposition - position,
                    "startorf": start_codons[position],
                    "stoporf": stopcodon,
                }
            if orf_data["length"] < maxlength and orf_data["length"] > minlength:
                orf_list.append(orf_data)
    return orf_list


def preporfs(transcript, starts, stops, minlength, maxlength):
    """
    Predict ORFs from transcript sequences using start and stop codon patterns.

    Args:
        transcript (str): Path to FASTA file with transcript sequences
        starts (list): List of start codon patterns
        stops (list): List of stop codon patterns
        minlength (int): Minimum ORF length
        maxlength (int): Maximum ORF length

    Returns:
        DataFrame: Predicted ORFs for each transcript sequence
    """
    dict_list = []
    counter = 0
    
    # Create Aho-Corasick automata
    startautomaton = ahocorasick.Automaton()
    stopautomaton = ahocorasick.Automaton()
    
    for start in starts:
        startautomaton.add_word(start, start)
    startautomaton.make_automaton()
    
    for stop in stops:
        stopautomaton.add_word(stop, stop)
    stopautomaton.make_automaton()
    
    # Process each transcript
    for record in SeqIO.parse(transcript, "fasta"):
        if counter % 1000 == 0:
            log_info(f"Processed {counter} transcripts")
            
        sequence = str(record.seq).upper()
        tran_id = record.id
        
        orfs = find_orfs(
            sequence,
            tran_id,
            startautomaton,
            stopautomaton,
            minlength,
            maxlength
        )
        dict_list.extend(orfs)
        counter += 1
        
    log_info(f"Found {len(dict_list)} ORFs in {counter} transcripts")
    return pl.DataFrame(dict_list)


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
    codingorfs = (
        codingorfs.with_columns(
            pl.struct(["start", "stop", "tran_start", "tran_stop"])
            .apply(lambda row: classify_orf(row))
            .alias("type")
        )
        .select(pl.all().exclude("tran_start", "tran_stop"))
    )

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
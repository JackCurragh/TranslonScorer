"""
ORF prediction functionality for TranslonScorer.

This module provides functions for finding Open Reading Frames (ORFs) in nucleotide sequences
using Aho-Corasick pattern matching for efficient codon identification.
"""

import ahocorasick
from ..utils.logging import log_info, log_warning
from pyfaidx import Fasta
import polars as pl


def find_all_positions(sequence, automaton):
    """
    Find positions of all occurrences of patterns from an Aho-Corasick automaton.

    Args:
        sequence (str): Input sequence to search for patterns
        automaton (ahocorasick.Automaton): Aho-Corasick automaton with patterns

    Returns:
        tuple: (frames, codons) where:
            - frames (dict): Maps frame indices (0,1,2) to lists of positions
            - codons (dict): Maps positions to identified codons
    """
    frames = {0: [], 1: [], 2: []}
    codons = {}
    
    for end_pos, pattern in automaton.iter(sequence):
        # position of last nucleotide of codon returned
        frames[(end_pos - 2) % 3].append(end_pos)
        codons[end_pos] = pattern
        
    return frames, codons


def find_orfs(sequence, tran_id, startautomaton, stopautomaton, minlength=0, maxlength=1000000):
    """
    Predict Open Reading Frames (ORFs) in a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence to analyze
        tran_id (str): Transcript identifier
        startautomaton (ahocorasick.Automaton): Automaton for start codons
        stopautomaton (ahocorasick.Automaton): Automaton for stop codons
        minlength (int, optional): Minimum ORF length. Defaults to 0
        maxlength (int, optional): Maximum ORF length. Defaults to 1000000

    Returns:
        list: List of dictionaries, each containing:
            - tran_id: Transcript identifier
            - start: Start position (0-based)
            - stop: Stop position
            - length: ORF length
            - startorf: Start codon sequence
            - stoporf: Stop codon sequence
    """
    log_info(f"Finding ORFs in transcript {tran_id}")
    
    orf_list = []
    startpositions, start_codons = find_all_positions(sequence, startautomaton)
    stoppositions, stop_codons = find_all_positions(sequence, stopautomaton)

    for frame, startpositions in startpositions.items():
        for position in startpositions:
            # Find valid stop codons downstream of start position
            valid_stops = [i for i in stoppositions[frame] if i > position]
            
            if valid_stops:
                stopposition = min(valid_stops)
                stopcodon = stop_codons[stopposition]
            else:
                stopposition = len(sequence)
                stopcodon = sequence[-3:]

            # Create ORF data
            orf_data = {
                "tran_id": tran_id,
                "start": position - 2,  # Adjust for 0-based indexing
                "stop": stopposition - 3 if stopcodon in ["TAA", "TAG", "TGA"] else stopposition,
                "length": stopposition - position,
                "startorf": start_codons[position],
                "stoporf": stopcodon
            }

            # Filter by length
            if minlength < orf_data["length"] < maxlength:
                orf_list.append(orf_data)

    return orf_list


def build_codon_automaton(codons):
    """
    Build an Aho-Corasick automaton for codon pattern matching.

    Args:
        codons (list): List of codon sequences to match

    Returns:
        ahocorasick.Automaton: Compiled automaton for pattern matching
    """
    automaton = ahocorasick.Automaton()
    
    for codon in codons:
        automaton.add_word(codon, codon)
    
    automaton.make_automaton()
    return automaton


def preporfs(sequence_input, start_codons=None, stop_codons=None, minlength=0, maxlength=1000000):
    """
    Predict ORFs from transcript sequences.

    Args:
        sequence_input (str or dict): Either a path to a FASTA file or a dictionary mapping transcript IDs to sequences
        start_codons (list, optional): List of start codon sequences. Defaults to ["ATG"]
        stop_codons (list, optional): List of stop codon sequences. Defaults to ["TAA", "TAG", "TGA"]
        minlength (int, optional): Minimum ORF length. Defaults to 0
        maxlength (int, optional): Maximum ORF length. Defaults to 1000000

    Returns:
        DataFrame: DataFrame containing predicted ORFs with their coordinates
    """
    # Set default codons if not provided
    if start_codons is None:
        start_codons = ["ATG"]
    if stop_codons is None:
        stop_codons = ["TAA", "TAG", "TGA"]
    
    # Handle input type
    if isinstance(sequence_input, str):
        # Input is a file path, read FASTA
        fasta = Fasta(sequence_input)
        sequence_dict = {name: str(seq) for name, seq in fasta.items()}
    elif isinstance(sequence_input, dict):
        # Input is already a dictionary
        sequence_dict = sequence_input
    else:
        raise ValueError("sequence_input must be either a file path (str) or a dictionary")

    # Build automata
    start_automaton = build_codon_automaton(start_codons)
    stop_automaton = build_codon_automaton(stop_codons)

    all_orfs = []
    for tran_id, sequence in sequence_dict.items():
        orfs = find_orfs(
            sequence=sequence,
            tran_id=tran_id,
            startautomaton=start_automaton,
            stopautomaton=stop_automaton,
            minlength=minlength,
            maxlength=maxlength
        )
        all_orfs.extend(orfs)
    
    if not all_orfs:
        return pl.DataFrame(schema={
            'tran_id': pl.Utf8,
            'start': pl.Int64,
            'stop': pl.Int64,
            'length': pl.Int64,
            'sequence': pl.Utf8
        })
    
    return pl.DataFrame(all_orfs) 

def extract_transcript_id(attr_str):
    """
    Extracts transcript ID from a GTF/GFF attribute string.

    This function takes a string representing attributes in GTF/GFF format
    and extracts the transcript ID from it. The transcript ID is typically
    found in attributes such as 'Parent=transcript:', 'ID=transcript:',
    'transcript_id=', or ' transcript_id '.

    Parameters:
        attr_str (str): A string containing attributes in GTF/GFF format.

    Returns:
        str: The extracted transcript ID, or an empty string if not found.

    Example:
        transcript_id = extract_transcript_id('gene_id="ENSG00000223972"; transcript_id="ENST00000456328"; ')
    """
    for attr in attr_str.split(";"):
        if attr.startswith("Parent=transcript:") or attr.startswith("ID=transcript:"):
            return attr.split(":")[1]
        elif attr.startswith("transcript_id="):
            return attr.split("=")[1]
        elif attr.startswith(" transcript_id "):
            return attr.split(" ")[2].replace('"', "")
    return ""


def exontranscriptcoords(df: pl.DataFrame, posstrand=True) -> pl.DataFrame:
    """
    Calculates transcript-level coordinates for exons.

    This function takes a DataFrame containing exon coordinates and calculates
    the transcript-level coordinates for each exon. For exons on the positive strand,
    the transcript coordinates start from 0 and increase. For exons on the negative
    strand, the transcript coordinates are calculated in reverse order.

    Parameters:
        df (polars.DataFrame): DataFrame containing exon coordinates.
        posstrand (bool): Flag indicating whether the exons are on the positive strand.
                          If True, transcript coordinates are calculated assuming exons
                          are on the positive strand. If False, transcript coordinates
                          are calculated assuming exons are on the negative strand.
                          Default is True.

    Returns:
        polars.DataFrame: DataFrame containing transcript-level exon coordinates.

    Example:
        exon_transcript_coords = exontranscriptcoords(exon_df, posstrand=True)
    """
    # Initialize new columns
    new_column_1 = []
    new_column_2 = []
    # Iterate over rows
    start_column = "start"
    end_column = "stop"

    for i in range(len(df)):
        start_values = (
            df[start_column][i]
            if posstrand
            else sorted(df[start_column][i], reverse=True)
        )
        end_values = (
            df[end_column][i] if posstrand else sorted(df[end_column][i], reverse=True)
        )
        new_start_values = []  # Starting value is 0
        new_stop_values = []
        for j in range(len(start_values)):
            if j == 0:
                new_start = 0
            else:
                new_start = new_stop_values[j - 1] + 1
            # Calculate stop coordinate
            stop_coordinate = end_values[j] - start_values[j]
            new_start_values.append(new_start)
            new_stop_values.append(new_start + stop_coordinate)
        new_column_1.append(new_start_values)
        new_column_2.append(new_stop_values)

    # Add new columns to the dataframe
    df = df.with_columns((pl.Series(new_column_1)).alias("tran_start"))
    df = df.with_columns((pl.Series(new_column_2)).alias("tran_stop"))
    return df


def gettranscriptcoords(cds_df, exon_df, posstrand=True):
    """
    Calculates transcript-level coordinates for coding sequences (CDS).

    This function takes a DataFrame containing CDS coordinates and a DataFrame
    containing exon coordinates. It then calculates the transcript-level coordinates
    for the CDS based on the exon coordinates.

    Parameters:
        cds_df (polars.DataFrame): DataFrame containing CDS coordinates.
        exon_df (polars.DataFrame): DataFrame containing exon coordinates.

    Returns:
        polars.DataFrame: DataFrame containing transcript-level CDS coordinates.

    Example:
        transcript_cds_coords = gettranscriptcoords(cds_df, exon_df)
    """
    # Explode exon start and stop lists to individual rows
    exploded_exon_df = exon_df.explode(["start", "stop", "tran_start", "tran_stop"])
    exploded_cds_df = cds_df.explode(["start", "stop"])
    # Join the DataFrames on the transcript ID
    combined_df = exploded_cds_df.join(exploded_exon_df, on="tran_id", how="inner")
    # Calculate transcript-level start and stop coordinates
    if posstrand:
        combined_df = combined_df.with_columns(
            [
                (
                    pl.when(
                        (pl.col("start") >= pl.col("start_right"))
                        & (pl.col("start") <= pl.col("stop_right"))
                    )
                    .then(
                        pl.col("tran_start") + (pl.col("start") - pl.col("start_right"))
                    )
                    .otherwise(None)
                ).alias("tran_start_cd"),
                (
                    pl.when(
                        (pl.col("stop") >= pl.col("start_right"))
                        & (pl.col("stop") <= pl.col("stop_right"))
                    )
                    .then(pl.col("tran_stop") - (pl.col("stop_right") - pl.col("stop")))
                    .otherwise(None)
                ).alias("tran_stop_cd"),
            ]
        )
    else:
        combined_df = combined_df.with_columns(
            [
                (
                    pl.when(
                        (pl.col("start") >= pl.col("start_right"))
                        & (pl.col("start") <= pl.col("stop_right"))
                    )
                    .then(
                        pl.col("tran_stop") - (pl.col("start") - pl.col("start_right"))
                    )
                    .otherwise(None)
                ).alias("tran_stop_cd"),
                (
                    pl.when(
                        (pl.col("stop") >= pl.col("start_right"))
                        & (pl.col("stop") <= pl.col("stop_right"))
                    )
                    .then(
                        pl.col("tran_start") + (pl.col("stop_right") - pl.col("stop"))
                    )
                    .otherwise(None)
                ).alias("tran_start_cd"),
            ]
        )
    # Drop rows with None values in calculated columns
    combined_df = combined_df.drop_nulls(["tran_start_cd", "tran_stop_cd"])
    combined_df = combined_df.group_by("tran_id").agg(
        pl.min("tran_start_cd"), pl.max("tran_stop_cd")
    )
    # Select and rename relevant columns
    result_df = combined_df.select(
        [
            pl.col("tran_id"),
            pl.col("tran_start_cd").alias("tran_start"),
            pl.col("tran_stop_cd").alias("tran_stop"),
        ]
    )

    return result_df


def procesexons(df):
    """
    Processes exon data by separating them based on strand orientation.

    This function takes a DataFrame containing exon data and separates the exons
    based on their strand orientation (positive or negative). It groups the exons
    by transcript ID and aggregates the start, stop, strand, and chromosome information
    for each group.

    Parameters:
        df (polars.DataFrame): DataFrame containing exon data.

    Returns:
        tuple: A tuple containing two polars DataFrames:
               - The first DataFrame contains exons on the positive strand.
               - The second DataFrame contains exons on the negative strand.

    Example:
        pos_exons, neg_exons = procesexons(exon_df)
    """
    exonplus = df.filter((pl.col("strand") == "+"))
    exonneg = df.filter((pl.col("strand") == "-"))

    groupedexonspos = (
        exonplus.group_by("tran_id")
        .agg([
            pl.col("start"),
            pl.col("stop"),
            pl.col("strand"),
            # Take first chromosome as they should all be the same for a transcript
            pl.col("chr").first().alias("chr")
        ])
        .select(["chr", "tran_id", "start", "stop", "strand"])
    )
    groupedexonsneg = (
        exonneg.group_by("tran_id")
        .agg([
            pl.col("start"),
            pl.col("stop"),
            pl.col("strand"),
            # Take first chromosome as they should all be the same for a transcript
            pl.col("chr").first().alias("chr")
        ])
        .select(["chr", "tran_id", "start", "stop", "strand"])
    )
    return groupedexonspos, groupedexonsneg


def getexons_and_cds(annotation_file, tran=[]):
    """
    Extracts CDS and exon coordinates from an annotation file.

    This function reads an annotation file in GTF/GFF format and extracts the
    coordinates of coding sequences (CDS) and exons. It then processes these
    coordinates to obtain transcript-level coordinates for exons and returns
    the results.

    Parameters:
        annotation_file (str): The path to the annotation file in GTF/GFF format.
        tran (list): A list of transcript IDs to filter. Only coordinates corresponding
                     to these transcripts will be extracted if provided. Default is [].

    Returns:
        tuple: A tuple containing two polars DataFrames:
               - The first DataFrame contains CDS coordinates.
               - The second DataFrame contains exon coordinates.

    Notes:
        - This function assumes the annotation file has columns separated by tabs ('\t').
        - The annotation file is expected to have no header, with comment lines starting with '#'.
        - The following columns are expected in the annotation file: 'chr', 'type', 'start', 'stop',
          'strand', 'attributes'.
        - The 'attributes' column is expected to contain transcript IDs.
        - The function 'extract_transcript_id' is used to extract transcript IDs from the 'attributes' column.
        - If 'tran' is provided, only coordinates corresponding to the specified transcripts will be extracted.

    Example:
        cds_coords, exon_coords = getexons_and_cds("annotation.gff", tran=['ENST00000223972', 'ENST00000456328'])
    """
    df = (
        pl.read_csv(
            annotation_file,
            separator="\t",
            ignore_errors=True,
            has_header=False,
            truncate_ragged_lines=True,
            comment_prefix="#",
        )
        .select(
            ["column_1", "column_3", "column_4", "column_5", "column_7", "column_9"]
        )
        .rename(
            {
                "column_1": "chr",
                "column_3": "type",
                "column_4": "start",
                "column_5": "stop",
                "column_7": "strand",
                "column_9": "attributes",
            }
        )
    )
    df = df.with_columns(
        pl.col("attributes")
        .apply(lambda attributes: extract_transcript_id(attributes))
        .alias("tran_id")
    ).select(pl.all().exclude("attributes"))

    if tran:
        df = df.filter((pl.col("tran_id").is_in(tran)))

    # Getting CDS
    coding_regions = df.filter((pl.col("type") == "CDS"))

    groupedcds = (
        coding_regions.group_by("tran_id")
        .agg(pl.col("start"), pl.col("stop"), pl.col("chr"))
        .select(["chr", "tran_id", "start", "stop"])
    )
    # Getting exons
    exon_regions = df.filter((pl.col("type") == "exon"))
    pos_exons, neg_exons = procesexons(exon_regions)

    # column names switched to calculate inverse of positions for negative strands
    exon_coords_pos = exontranscriptcoords(pos_exons, posstrand=True)
    exon_coords_neg = exontranscriptcoords(neg_exons, posstrand=False)

    cds_coords_pos = gettranscriptcoords(groupedcds, exon_coords_pos, posstrand=True)
    cds_coords_neg = gettranscriptcoords(groupedcds, exon_coords_neg, posstrand=False)

    cds_coords = pl.concat([cds_coords_pos, cds_coords_neg])
    exondf = pl.concat([exon_coords_pos, exon_coords_neg]).select(
        pl.all().exclude("strand")
    )
    return cds_coords, exondf

def classify_orf(row):
    """
    Classify an ORF (Open Reading Frame) based on its relative position to transcript start and stop sites.

    Parameters:
    - row (Series or dict-like): A Pandas Series or dictionary-like object containing ORF information,
                                including 'start', 'stop', 'tran_start', and 'tran_stop' values.

    Returns:
    - str: A string indicating the classification of the ORF based on its relative position:
           - "uORF": Upstream ORF (stop < tran_start)
           - "CDS": Coding Sequence (start == tran_start and stop == tran_stop)
           - "dORF": Downstream ORF (start > tran_stop)
           - "uoORF": Upstream Overlapping ORF (start < tran_start and stop >= tran_start)
           - "doORF": Downstream Overlapping ORF (start <= tran_stop and stop > tran_stop)
           - "iORF": Internal ORF (start >= tran_start and stop <= tran_stop)
           - "eoORF": Encapsulated Overlapping ORF (start < tran_start and stop > tran_stop)
           - "extORF": Extended ORF (start < tran_start and stop == tran_stop)
           - "Unexpected": Indicates unexpected conditions where none of the above criteria are met.

    This function categorizes an ORF based on its positional relationship with the transcript start (`tran_start`)
    and stop (`tran_stop`) sites. It evaluates the relative positions of 'start' and 'stop' compared to
    'tran_start' and 'tran_stop' to determine the appropriate classification.

    If the relative position does not match any expected categories, an "Unexpected" classification is printed
    with the values of 'start', 'stop', 'tran_start', and 'tran_stop', and "Unexpected" is returned.

    Note: This function assumes the input `row` contains numerical values for 'start', 'stop', 'tran_start',
    and 'tran_stop', typically retrieved from a Pandas DataFrame or similar data structure.
    """
    print(row)
    if row['tran_start'] == None and row['tran_stop'] == None:
        return "Non Coding"
    elif row["stop"] < row["tran_start"]:
        return "uORF"
    elif row["start"] == row["tran_start"] and row["stop"] == row["tran_stop"]:
        return "CDS"
    elif row["start"] > row["tran_stop"]:
        return "dORF"
    elif row["start"] < row["tran_start"] and row["stop"] >= row["tran_start"]:
        return "uoORF"
    elif row["start"] <= row["tran_stop"] and row["stop"] > row["tran_stop"]:
        return "doORF"
    elif row["start"] >= row["tran_start"] and row["stop"] <= row["tran_stop"]:
        return "iORF"
    elif row["start"] < row["tran_start"] and row["stop"] > row["tran_stop"]:
        return "eoORF"
    elif row["start"] < row["tran_start"] and row["stop"] == row["tran_stop"]:
        return "extORF"
    else:
        print(
            "unexpected", row["start"], row["stop"], row["tran_start"], row["tran_stop"]
        )
        return "Unexpected"
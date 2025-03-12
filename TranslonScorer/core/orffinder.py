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
    
    for i, (order, codon) in automaton.iter(sequence):
        # position of last nucleotide of codon returned
        frames[(i - 2) % 3].append(i)
        codons[i] = codon
        
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

    log_info(f"Found {len(orf_list)} ORFs in transcript {tran_id}")
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
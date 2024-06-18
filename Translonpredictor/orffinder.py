def find_all_positions(sequence, automaton):
    """
    Find positions of all occurrences of patterns from an Aho-Corasick automaton in a given sequence.

    Parameters:
    - sequence (str): The input sequence (e.g., DNA or RNA sequence) to search for patterns.
    - automaton (ahocorasick.Automaton): An Aho-Corasick automaton containing patterns to search for.

    Returns:
    - tuple: A tuple containing:
        - dict: A dictionary where keys are frame indices (0, 1, 2) representing reading frames
                and values are lists of positions where patterns were found in the sequence.
        - dict: A dictionary where keys are positions in the sequence where patterns were found
                and values are the corresponding codons identified by the automaton.

    This function uses an Aho-Corasick automaton (`automaton`) to efficiently search for multiple
    patterns (codons) within the given `sequence`. It iterates over matches found by the automaton
    and organizes them into reading frames (0, 1, 2) based on their position in the sequence.

    The function returns two dictionaries:
    - `frames`: Contains lists of positions (0-based) in the sequence for each reading frame
                where patterns were found.
    - `codons`: Maps positions in the sequence to the specific codons identified by the automaton.

    Note: This function assumes the use of the `ahocorasick` library, which provides an efficient
    implementation of the Aho-Corasick algorithm for pattern matching.
    """
    frames = {0: [], 1: [], 2: []}
    codons = {}
    for i, (order, codon) in automaton.iter(sequence):
        # position of last nucleotide of codon returned
        frames[(i - 2) % 3].append(i)
        codons[i] = codon
    return frames, codons


def find_orfs(
    sequence, tran_id, startautomaton, stopautomaton, minlength=0, maxlength=1000000
):
    """
    Predict Open Reading Frames (ORFs) in a given nucleotide sequence using start and stop codon patterns.

    Parameters:
    - sequence (str): Nucleotide sequence (e.g., DNA or RNA) in which ORFs will be predicted.
    - tran_id (str): Identifier for the transcript sequence.
    - startautomaton (ahocorasick.Automaton): A pre-built Aho-Corasick automaton for start codons (e.g., 'ATG').
    - stopautomaton (ahocorasick.Automaton): A pre-built Aho-Corasick automaton for stop codons (e.g., 'TAA').
    - minlength (int, optional): Minimum length threshold for predicted ORFs. Defaults to 0.
    - maxlength (int, optional): Maximum length threshold for predicted ORFs. Defaults to 1000000 (1 million).

    Returns:
    - list: A list of dictionaries, where each dictionary represents an ORF with the following keys:
        - "tran_id": Identifier of the transcript sequence.
        - "start": Start position of the ORF in 0-based indexing.
        - "stop": Stop position of the ORF in 0-based indexing.
        - "length": Length of the ORF.
        - "startorf": Start codon sequence.
        - "stoporf": Stop codon sequence.

    This function identifies potential ORFs in the provided `sequence` by searching for start and stop codon
    patterns using pre-built Aho-Corasick automata (`startautomaton` and `stopautomaton`). For each identified
    start codon, it searches for a valid stop codon downstream and calculates the ORF properties.

    ORFs are represented as dictionaries containing information such as transcript ID (`tran_id`), start and
    stop positions (adjusted for frame), ORF length, and the sequences of start and stop codons (`startorf`
    and `stoporf`).

    ORFs are filtered based on user-defined `minlength` and `maxlength` thresholds before being added to the
    `orf_list`, which is returned as the output.

    Note: This function assumes the use of the `ahocorasick` library for efficient pattern matching
    with Aho-Corasick automata.
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

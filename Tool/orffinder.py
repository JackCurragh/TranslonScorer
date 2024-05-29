
def find_all_positions(sequence, automaton):
    frames = {0:[], 1:[], 2:[]}
    codons = {}
    for i, (order, codon) in automaton.iter(sequence):
        #position of last nucleotide of codon returned  
        frames[(i-2)%3].append(i)
        codons[i] = codon
    return frames, codons

def find_orfs(
    sequence,
    tran_id,
    startautomaton,
    stopautomaton,
    minlength=0,
    maxlength=1000000
):
    """
    Finds Open Reading Frames (ORFs) within a given sequence based on specified start and stop codons.

    Parameters:
    - sequence (str): Input DNA sequence to search for ORFs.
    - start_codons (list): List of start codons to identify the start of an ORF. Default is ["ATG"].
    - stop_codons (list): List of stop codons to identify the end of an ORF. Default is ["TAA", "TAG", "TGA"].
    - minlength (int): Minimum length of ORF to be considered. Default is 0.
    - maxlength (int): Maximum length of ORF to be considered. Default is 1000000.

    Returns:
    - orfs (list): List of tuples containing the start position, end position, length, and the sequence of each identified ORF.

    This function searches for Open Reading Frames (ORFs) within the given DNA sequence. It iterates over the sequence
    and looks for start and stop codons to identify potential ORFs. The identified ORFs are stored as tuples
    containing the start position, end position, length, and the corresponding sequence. ORFs that meet the specified
    length constraints are included in the output list.
    """
    orf_list=[]
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
            orf_data= {
                "tran_id":tran_id,
                "start":position-2,
                "stop": stopposition-3,
                "length":stopposition - position,
                "startorf":start_codons[position],
                "stoporf": stopcodon
                }
            if orf_data['length'] < maxlength and orf_data['length'] > minlength:
                orf_list.append(orf_data)
    return orf_list

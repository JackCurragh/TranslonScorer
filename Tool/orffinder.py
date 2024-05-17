
def find_orfs(
    sequence,
    tran_id,
    start_codons=["ATG"],
    stop_codons=["TAA", "TAG", "TGA"],
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
    append_list=[]
    in_orf = False
    current_orf_length = 0

    for index in range(1, len(sequence)):
        codon = sequence[index : index + 3]
        if codon in start_codons:
            in_orf = True
            start = index
            startcodon = codon
            current_orf_length = 3
        if in_orf:
            for orf in range(index+3, len(sequence), 3):
                codon = sequence[orf : orf + 3]
                if codon in stop_codons:
                    current_orf_length += 3
                    if current_orf_length > minlength and current_orf_length < maxlength:
                        #Frames
                        if start % 3 == 0:
                            frame = 0
                        elif start % 3 == 1:
                            frame = 1
                        elif start % 3 == 2:
                            frame = 2
                        #add values to dict
                        append_data= {
                            "tran_id":tran_id,
                            "start":start,
                            "end":orf+3,
                            "length":current_orf_length,
                            "startorf":str(startcodon),
                            "stoporf":str(codon),
                            "frame":frame
                            }
                        append_list.append(append_data)
                    index = index + 1
                    in_orf = False
                    start = None
                    current_orf_length = 0
                    break
                elif in_orf:
                    current_orf_length += 3
    
    return append_list

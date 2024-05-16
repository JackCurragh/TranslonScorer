def find_orfs(
    sequence,
    start_codons=["ATG"],
    stop_codons=["TAA", "TAG", "TGA"],
    minlength=0,
    maxlength=1000000,
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
    orfs = []
    in_orf = False
    current_orf_start = None
    current_orf_length = 0

    for index in range(0, len(sequence)):
        codon = sequence[index : index + 3]
        if codon in start_codons:
            in_orf = True
            current_orf_start = index
            current_orf_length = 3
        if in_orf:
            for orf in range(index+3, len(sequence), 3):
                codon = sequence[orf : orf + 3]
                if codon in stop_codons:
                    current_orf_length += 3
                    if current_orf_length > minlength and current_orf_length < maxlength:
                        orfs.append(
                            (
                                current_orf_start+1,
                                orf+3,
                                current_orf_length,
                                sequence[current_orf_start : orf + 3],
                            )
                        )
                    index = current_orf_start + 1
                    in_orf = False
                    current_orf_start = None
                    current_orf_length = 0
                    break
                elif in_orf:
                    current_orf_length += 3
    return orfs

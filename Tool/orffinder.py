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

    for index in range(0, len(sequence), 3):
        codon = sequence[index : index + 3]

        if codon in start_codons:
            if not in_orf:
                in_orf = True
                current_orf_start = index
                current_orf_length = 3
            else:
                current_orf_length += 3
        elif codon in stop_codons:
            if in_orf:
                if current_orf_length > minlength and current_orf_length < maxlength:
                    orfs.append(
                        (
                            current_orf_start,
                            index,
                            current_orf_length,
                            sequence[current_orf_start : index + 3],
                        )
                    )
                in_orf = False
                current_orf_start = None
                current_orf_length = 0
        elif in_orf:
            current_orf_length += 3

    return orfs

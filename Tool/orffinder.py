def find_orfs(sequence, start_codons=['ATG',], stop_codons=['TAA', 'TAG', 'TGA'], minlength=0, maxlength=1000000):
    orfs = []
    in_orf = False
    current_orf_start = None
    current_orf_length = 0
    
    for index in range(0, len(sequence), 3):
        codon = sequence[index:index+3]
        
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
                    orfs.append((current_orf_start, index + 3, current_orf_length, sequence[current_orf_start : index +3]))
                in_orf = False
                current_orf_start = None
                current_orf_length = 0
        elif in_orf:
            current_orf_length += 3
    
    return orfs


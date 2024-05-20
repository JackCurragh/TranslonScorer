import polars as pl
from numpy import log as ln

def sru_score(start, bigwig_df, rng, invert):
    '''
    docstring
    '''
    afterstart=[]
    beforestart=[]
    for i in range(start-rng, start+3+rng):
        if i % 3 == start % 3:
            if i < start:
                before = bigwig_df.filter(pl.col("tran_start") == i)
                if not before.is_empty():
                    before = before["counts"][0]
                    beforestart.append(before)
            elif i == start:
                continue
            else:
                after = bigwig_df.filter(pl.col("tran_start") == i)
                if not after.is_empty():
                    after = after["counts"][0]
                    afterstart.append(after)

    numerator = 1 + sum(afterstart)
    denominator = 1 + sum(beforestart)
    #Check for Start rise up or step down score
    if invert == 0:
        sru = ln(numerator / denominator)
    else:
        sru = ln(denominator / numerator)

    return float(sru)

sru_data = {
"tran_start":[24,25,26,27,29,31,32,33,36,37,39,41,42,46,47,50,2300,2301,2303,2304,2305,2307,2310,2311,2312,2316,2317],
"counts":[1.0,1.0,1.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,1.0,1.0]}

sru_df = pl.DataFrame(sru_data)
sru_score(36, sru_df, 12, 0)

def test_sru():
    assert sru_score(36, sru_df, 12, 0) == -0.2231435513142097
    assert sru_score(2307, sru_df, 12, 1) == 0.28768207245178085

def calculate_scores(start, stop, bigwig_df):
    '''
    Calculate HRF, average, and NZC scores.
    
    Parameters:
        start (int): Start position.
        stop (int): Stop position.
        bigwig_df (polars DataFrame): DataFrame containing relevant data.
        
    Returns:
        dict: A dictionary containing HRF, average, and NZC scores.
    '''
    framecounts = {0: 0, 1: 0, 2: 0}
    framescore = []
    codons_f0 = 0
    total_codons = 0
    frame0 = 0
    frame1 = 0
    frame2 = 0
    for i in range(start, stop + 1):
        counts = bigwig_df.filter(pl.col("tran_start") == i)
        if i % 3 == 0:
            frame0 = 0
            if not counts.is_empty():
                framecounts[0] += counts["counts"][0]
                framescore.append(counts["counts"][0])
                frame0 = counts["counts"][0]
        elif i % 3 == 1:
            frame1 = 0
            if not counts.is_empty():
                framecounts[1] += counts["counts"][0]
                frame1 = counts["counts"][0]
        elif i % 3 == 2:
            frame2 = 0
            if not counts.is_empty():
                framecounts[2] += counts["counts"][0]
                frame2 = counts["counts"][0]
            if frame0 > frame1 and frame0 > frame2:
                codons_f0 += 1
            if frame0 > 0 or frame1 > 0 or frame2 > 0:
                total_codons += 1
    
    hrf = framecounts[0] / max(framecounts[1], framecounts[2]) if framecounts[1] and framecounts[2] != 0 else 0
    avg = sum(framescore) / (len(range(start, stop + 1)) / 3) if framescore else 0
    nzc = codons_f0 / total_codons if total_codons > 0 else 0
    
    return float(hrf), float(avg), float(nzc)


data = {"tran_id":["ENST00000222271.7", "ENST00000222271.7", "ENST00000222271.7","ENST00000222271.7","ENST00000222271.7",
"ENST00000222271.7","ENST00000222271.7"],
"tran_start":[222, 965, 1008, 1130, 1463, 1567, 1598],
"tran_stop":[223,966,1009,1131,1464,1568,1599],
"counts":[1.0,1.0,1.0,2.0,1.0,1.0,1.0]}

df = pl.DataFrame(data)

def test_scores():
    assert calculate_scores(36, 2307, df) == (0.4, 0.002640845070422535, 0.2857142857142857)
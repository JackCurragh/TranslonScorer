'''Script to score the ORF's'''
import polars as pl
from numpy import log as ln

def hrf_score(start, stop, bigwig_df):
    '''
    docstring
    '''
    framecounts = {0:0, 1:0, 2:0}
    for i in range(start, stop+1):
        counts = bigwig_df.filter(pl.col("tran_start") == i)
        if not counts.is_empty():
            framecounts[i % 3] += counts["counts"][0]
    if framecounts[1] and framecounts[2] == 0:
        hrf = framecounts[0] / max(framecounts[1], framecounts[2])
    else: hrf = 0
    return hrf

def sru_score(start, bigwig_df, rng, invert):
    '''
    docstring
    '''
    scoretop=[]
    scorebot=[]
    for i in range(start-rng, start+3+rng):
        if i % 3 == start % 3:
            if i < start:
                bot = bigwig_df.filter(pl.col("tran_start") == i)
                if not bot.is_empty():
                    bot = bot["counts"][0]
                    scorebot.append(bot)
            elif i == start:
                continue
            else:
                top = bigwig_df.filter(pl.col("tran_start") == i)
                if not top.is_empty():
                    top = top["counts"][0]
                    scoretop.append(top)

    numerator = 1 + sum(scoretop)
    denominator = 1 + sum(scorebot)
    #Check for Start rise up or step down score
    if invert == 0:
        sru = ln(numerator / denominator)
    else:
        sru = ln(denominator / numerator)
    return sru

def nzc_score(start, stop, bigwig_df):
    '''
    docstring
    '''
    codons_f0 = 0
    total_codons = 0
    for i in range(start, stop+1):
        counts = bigwig_df.filter(pl.col("tran_start") == i)
        if i % 3 == 0:
            frame0 = 0
            if not counts.is_empty():
                frame0 = counts["counts"][0]
        if i % 3 == 1:
            frame1 =0
            if not counts.is_empty():
                frame1 = counts["counts"][0]
        if i % 3 == 2:
            frame2 = 0
            if not counts.is_empty():
                frame2 = counts["counts"][0]
            if frame0 > frame1 and frame0 > frame2:
                codons_f0 += 1
            if frame0 > 0 or frame1 > 0 or frame2 > 0:
                total_codons += 1
        
    if total_codons > 0:    
        nzc = codons_f0/total_codons
    else: nzc = 0
    return nzc

def avg_score(start, stop, bigwig_df):
    '''
    docstring
    '''
    framescore = []
    for i in range(start, stop+1):
        if i % 3 == 0:
            counts = bigwig_df.filter(pl.col("tran_start") == i)
            if not counts.is_empty():
                framescore.append(counts["counts"][0])
    
    codons = len(range(start, stop+1))/3
    avg =  sum(framescore)/codons
    return avg


def final_score(hrf, nzc, sru1, sru2, avg):
    score = sum([hrf, avg, nzc, sru1, sru2])
    return score
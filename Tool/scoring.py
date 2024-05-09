'''Script to score the ORF's'''
import polars as pl


def hrf_score(start, stop, bigwig_df):
    '''
    docstring
    '''
    framecounts = {0:1, 1:0, 2:0}
    for i in range(start, stop):
        counts = bigwig_df.filter(pl.col("tran_start") == i)
        print(counts)
        if not counts.is_empty():
            framecounts[i % 3] += counts["counts"][0]
    if framecounts[1] or framecounts[2] is not 0:
        hrf = framecounts[0] / max(framecounts[1], framecounts[2])
        print(hrf)
        return hrf

def sru_score():
    '''
    docstring
    '''
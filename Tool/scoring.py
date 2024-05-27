'''Script to score the ORF's'''
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
                    beforestart.append(before["counts"][0])
            elif i == start:
                continue
            else:
                after = bigwig_df.filter(pl.col("tran_start") == i)
                if not after.is_empty():
                    afterstart.append(after["counts"][0])

    numerator = 1 + sum(afterstart)
    denominator = 1 + sum(beforestart)
    #Check for Start rise up or step down score
    if invert == 0:
        sru = ln(numerator / denominator)
    else:
        sru = ln(denominator / numerator)

    return float(sru)



import polars as pl

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
    # Filter the dataframe for the range of interest
    relevant_df = bigwig_df.filter((pl.col("tran_start") >= start) & (pl.col("tran_start") <= stop))
    # Create a new column for frame based on tran_start
    relevant_df = relevant_df.with_columns((pl.col("tran_start") % 3).alias("frame"))
    # Calculate framecounts
    framecounts = relevant_df.groupby("frame").agg(pl.sum("counts").alias("frame_counts")).to_dict(as_series=False)
    # Extract frame counts
    framecounts = {frame: framecounts["frame_counts"][i] if i < len(framecounts["frame_counts"]) else 0 for i, frame in enumerate([0, 1, 2])}
    # Filter the data frame for frame 0
    frame0_df = relevant_df.filter(pl.col("frame") == 0)
    # Calculate framescore
    framescore = frame0_df["counts"].sum()
    # Calculate codons_f0 and total_codons
    codons_df = relevant_df.groupby((pl.col("tran_start") // 3)).agg([
        pl.sum(pl.when(pl.col("frame") == 0).then(pl.col("counts")).otherwise(0)).alias("frame0"),
        pl.sum(pl.when(pl.col("frame") == 1).then(pl.col("counts")).otherwise(0)).alias("frame1"),
        pl.sum(pl.when(pl.col("frame") == 2).then(pl.col("counts")).otherwise(0)).alias("frame2")
    ])
    codons_f0 = codons_df.filter((pl.col("frame0") > pl.col("frame1")) & (pl.col("frame0") > pl.col("frame2"))).height
    total_codons = codons_df.filter((pl.col("frame0") > 0) | (pl.col("frame1") > 0) | (pl.col("frame2") > 0)).height
    
    # Calculate hrf, avg, and nzc
    hrf = framecounts[0] / max(framecounts[1], framecounts[2]) if framecounts[1] and framecounts[2] != 0 else 0
    avg = framescore / (len(range(start, stop + 1)) / 3) if framescore else 0
    nzc = codons_f0 / total_codons if total_codons > 0 else 0
    
    return float(hrf), float(avg), float(nzc)


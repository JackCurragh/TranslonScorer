"""Script to plot file"""

orfs = "data/files/annotated_orfs.csv"
bwfile = "data/files/bw_tran.csv"
range_param = 30

import polars as pl
import plotly.express as px

from scoring import sru_score, calculate_scores

def scoreandplot(orfs, bwfile, range_param, sru_range=12):
    """
    Preprocesses data and generates bar plots for each transcript based on input ORFs and BigWig files.
    
    Parameters:
    - orfs (str): Path to CSV file containing ORF data.
    - bwfile (str): Path to BigWig file.
    - range_param (int): Parameter to define the range for locating the transcripts.
    
    Returns:
    - None
    
    This function reads the ORF data and BigWig file, filters the ORF data for a certain type of ORF,
    sorts the BigWig file by transcript starting positions, selects a range to search for the transcripts.
    For each ORF, it locates the corresponding transcript in the BigWig file, calculates relative start position,
    groups the counts by transcript start position, filters data within the defined range, groups the counts
    by relative location to the ORF start position and then plots a bar graph to visualize
    the count per base data for each transcript.
    """
    # orf df
    orf_df = pl.read_csv(orfs)
    #for typeorf in orf_df['type'].unique():
    orf_df_filtered = orf_df.filter((pl.col("type") == 'uoORF'))
    # bw df
    bigwig_df = pl.read_csv(bwfile).sort(pl.col("tran_start"))
    score_col = []
    for row in range(len(orf_df_filtered)):
        transcript = orf_df_filtered["tran_id"][row]
        relstart = orf_df_filtered["pos"][row]
        end = orf_df_filtered["end"][row]
        
        bigwig_tran = bigwig_df.filter(pl.col("tran_id") == transcript)
        
        if not bigwig_tran.is_empty():
            #passing 0,1 to choose invert of sru score for step down score
            sru1 = sru_score(relstart, bigwig_tran, sru_range, 0)
            sru2 = sru_score(end, bigwig_tran, sru_range, 1)
            hrf, avg, nzc = calculate_scores(relstart, end, bigwig_tran)
            score = sum([hrf, nzc, sru1, sru2, avg])
            score_col.append(score)
        else:
            score_col.append(0.0)

    orf_df_filtered = orf_df_filtered.with_columns(pl.Series(score_col).alias("score"))
    top10_df = orf_df_filtered.sort("score", descending=True).head(10)
        
        #getting counts for top 10 transcripts 
    range_list = list(range(-range_param, range_param + 1))
    relative_df = pl.DataFrame()
    normal_df = pl.DataFrame()

    for row in range(len(top10_df)):
        transcript = top10_df["tran_id"][row]
        relstart = top10_df["pos"][row]
        bigwig_tran2 = bigwig_df.filter(pl.col("tran_id") == transcript)
        
        normalcounts = bigwig_tran2.group_by("tran_id", "tran_start").agg(pl.col("counts").sum())
        normal_df = pl.concat([normal_df, normalcounts])

        relativecounts = normalcounts.with_columns((pl.col("tran_start") - relstart).alias("relativeloc"))
        relative_df = pl.concat([relative_df, relativecounts])

    relative_df = relative_df.filter(pl.col("relativeloc").is_in(range_list))
    relative_df = relative_df.group_by("relativeloc").agg(pl.col("counts").sum())


    print(normal_df)
    fig_separate = px.bar(
        normal_df.to_pandas(), x="tran_start", y="counts", title="Top 10 scoring upstream overlapping ORFs", color='tran_id',
    )
    fig_separate.show()

    # plot metagene
    fig_combined = px.bar(
        relative_df.to_pandas(), x="relativeloc", y="counts", title="typeorf",
    )
    fig_combined.show()


scoreandplot(orfs, bwfile, range_param)

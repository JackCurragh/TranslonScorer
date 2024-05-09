"""Script to plot file"""

orfs = "data/files/annotated_orfs.csv"
bwfile = "data/files/bw_tran.csv"
range_param = 30

import polars as pl
import plotly.express as px

from scoring import hrf_score

def preprocessplot(orfs, bwfile, range_param):
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
    orf_df = orf_df.filter((pl.col("type") == "CDS"))
    # bw df
    bigwig_df = pl.read_csv(bwfile).sort(pl.col("tran_start"))
    # range to search in
    range_list = list(range(-range_param, range_param + 1))
    concat_df = pl.DataFrame()
    for row in range(len(orf_df)):
        transcript = orf_df["tran_id"][row]
        relstart = orf_df["pos"][row]
        end = orf_df["end"][row]
        
        bigwig_tran = bigwig_df.filter(pl.col("tran_id") == transcript)
        if not bigwig_tran.is_empty():
            hrf_score(relstart, end, bigwig_tran)
            concat_df = pl.concat([concat_df, bigwig_tran])
    plot_df = concat_df.group_by("tran_start").agg(pl.col("counts").sum())
    plot_df = plot_df.with_columns((pl.col("tran_start") - relstart).alias("relativeloc"))
    plot_df = plot_df.filter(pl.col("relativeloc").is_in(range_list))
    plot_df = plot_df.group_by("relativeloc").agg(pl.col("counts").sum())
    
    # plot figure
    fig = px.bar(
        plot_df.to_pandas(), x="relativeloc", y="counts", title=transcript
    )
    fig.show()


preprocessplot(orfs, bwfile, range_param)

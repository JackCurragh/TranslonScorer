"""Script to plot file"""

orfs = "data/files/annotated_orfs.csv"
bwfile = "data/files/bw_tran.csv"
range_param = 100

import polars as pl
import plotly.express as px


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
    groups the counts by transcript start position, filters data within the defined range, and then plots a bar graph
    to visualize the count per base data for each transcript.
    """
    # orf df
    orf_df = pl.read_csv(orfs)
    orf_df = orf_df.filter((pl.col("type") == "Non Coding"))
    # bw df
    bw_df = pl.read_csv(bwfile).sort(pl.col("tran_start"))
    # range to search in
    range_list = list(range(-range_param, range_param + 1))
    for i in range(len(orf_df)):
        transcript = orf_df["tran_id"][i]
        relstart = orf_df["pos"][i]

        bw_df_tran = bw_df.filter(pl.col("tran_id") == transcript)
        if not bw_df_tran.is_empty():
            bw_df_tran = bw_df_tran.group_by("tran_start").agg(pl.col("counts").sum())
            bw_df_tran = bw_df_tran.with_columns(
                (pl.col("tran_start") - relstart).alias("relativeloc")
            )
            bw_df_tran = bw_df_tran.filter(pl.col("relativeloc").is_in(range_list))
            print(bw_df_tran)
            # plot figure
            fig = px.bar(
                bw_df_tran.to_pandas(), x="tran_start", y="counts", title=transcript
            )
            fig.show()


preprocessplot(orfs, bwfile, range_param)

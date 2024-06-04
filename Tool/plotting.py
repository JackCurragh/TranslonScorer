"""Script to plot file"""

import polars as pl
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pyBigWig as bw

from bigwigtodf import transcriptreads

def metageneplot(df, bwfile, exon_df, range_list):
    '''
    docstring
    '''
    dflist=[]
    for typeorf in df['type'].unique():
        df_type_filtered = df.filter(pl.col('type') == typeorf)
        df_type_filtered = df_type_filtered.sort('score', descending=True).head(10)
        metagene_list=[]
        for row in range(len(df_type_filtered)):
            tran = df_type_filtered['tran_id'][row]
            start = df_type_filtered['start'][row]

            exons = exon_df.filter(pl.col('tran_id') == tran)
            if exons.is_empty():
                continue
            tran_reads = transcriptreads(bwfile, exons)
            if not tran_reads.is_empty():
                coordinates = tran_reads.group_by("tran_start")\
                                       .agg(pl.col("counts").sum())\
                                       .with_columns(
                                           (pl.col("tran_start") - start)
                                           .alias("relativeloc"))\
                                       .to_dict(as_series=False)
                metagene_list.append(coordinates)

        metagene_df = pl.from_dicts(metagene_list).explode(pl.all())
        metagene_df = metagene_df.filter(pl.col("relativeloc").is_in(range_list))\
                                 .group_by("relativeloc")\
                                 .agg(pl.col("counts").sum())
            
        #plot metagene per type
        fig_combined = px.bar(
        metagene_df.to_pandas(), x="relativeloc", y="counts", title=typeorf,
        )
        fig_combined.update_xaxes(title_text='Relative coordinates')
        fig_combined.update_yaxes(title_text='Counts')
        fig_combined.show()

        df_type_filtered = df_type_filtered.to_dict(as_series=False)
        dflist.append(df_type_filtered)
    
    df = pl.from_dicts(dflist).explode(pl.all())
    return df

def plottop10(df, bigwig, exon, range_param):
    '''
    docstring
    '''
    range_list = list(range(-range_param, range_param + 1))
    df = pl.read_csv(df, has_header=True, separator=",")
    bwfile = bw.open(bigwig)
    exon_df = pl.read_csv(exon, has_header=True, separator=",")
    exon_df = exon_df.with_columns(
            pl.col("start", "stop", "tran_start", "tran_stop").apply(lambda x: x.split(",")),
            pl.col('chr').apply(lambda x: x.split('r')[1]))
    
    top10_df = metageneplot(df, bwfile, exon_df, range_list)
    





    ###############################################################################################
    #getting counts for top 10 transcripts 
    normal_df = pl.DataFrame()
    for row in range(len(top10_df)):
        transcript = top10_df["tran_id"][row]
        relstart = top10_df["start"][row]
        orftype = top10_df["type"][row]
        frame = top10_df["start"][row] % 3

        bigwig_tran2 = bigwig_df.filter((pl.col("tran_id") == transcript))
        
        normalcounts = bigwig_tran2.group_by("tran_id", "start").agg(pl.col("counts").sum())
        normal_df = pl.concat([normal_df, normalcounts])

        relativecounts = normalcounts.with_columns((pl.col("tran_start") - relstart).alias("relativeloc"))
        relative_df = pl.concat([relative_df, relativecounts])
    
    
    relative_df = relative_df.filter(pl.col("relativeloc").is_in(range_list))
    relative_df = relative_df.group_by("relativeloc").agg(pl.col("counts").sum())
    
    
    #PLOTTING
    fig_separate = px.bar(
        normal_df.to_pandas(), x="tran_start", y="counts", title="Top 10 scoring upstream overlapping ORFs", color='type',
    )
    fig_separate.update_xaxes(title_text='Transcriptomic coordinates')
    fig_separate.update_yaxes(title_text='Counts')
    fig_separate.show()

    # plot metagene
    fig_combined = px.bar(
        relative_df.to_pandas(), x="relativeloc", y="counts", title="typeorf",
    )
    fig_combined.update_xaxes(title_text='Relative coordinates around the start coordinate')
    fig_combined.update_yaxes(title_text='Counts')
    fig_combined.show()
    return normal_df, relative_df

    
  
    
    return

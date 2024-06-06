"""Script to plot file"""

import polars as pl
import plotly.express as px
import plotly.subplots as sp
import pyBigWig as bw
import plotly.graph_objects as go

from bigwigtodf import transcriptreads
from report import generate_report

def summaryplot(df, df2):
    '''
    docstring
    '''
    df=df.to_pandas()
    df2=df2.to_pandas()
    
    #Table generation
    fig = go.Figure(data=[go.Table(
        header=dict(values=list(df2.columns),
                fill_color='paleturquoise',
                align='left'),
        cells=dict(values=[df2.tran_id, df2.start, df2.stop,
                           df2.length,df2.startorf,df2.stoporf,df2.type,
                           df2.rise_up,df2.step_down,df2.hrf,df2.avg,
                           df2.nzc,df2.score, df2.frame],
               fill_color='lavender',
               align='left'))
        ])

    figure1 = px.bar(df, x="tran_start", y="counts", color='tran_id')
    figure1.update_xaxes(title_text='Transcript coordinates')
    figure1.update_yaxes(title_text='Counts')

    return figure1, fig

def pertranscriptplot(df):
    '''
    docstring
    '''
    pertranlist=[]
    for tran in df['tran_id'].unique():
        transcriptdf = df.filter(pl.col('tran_id') == tran)
        fig = px.bar(transcriptdf.to_pandas(), x="tran_start", y="counts", title=tran)
        fig.update_xaxes(title_text='Transcript coordinates')
        fig.update_yaxes(title_text='Counts')
        fig = fig.to_html(full_html=False)
        pertranlist.append(fig)

    return pertranlist

def metageneplot(df, bwfile, exon_df, range_list):
    '''
    docstring
    '''
    dflist=[]
    plotlist=[]
    tranlist=[]
    for typeorf in df['type'].unique():
        df_type_filtered = df.filter(pl.col('type') == typeorf)
        df_type_filtered = df_type_filtered.sort('score', descending=True).head(10)
        metagene_start_list=[]
        metagene_stop_list=[]
        for row in range(len(df_type_filtered)):
            tran = df_type_filtered['tran_id'][row]
            start = df_type_filtered['start'][row]
            stop = df_type_filtered['stop'][row]

            exons = exon_df.filter(pl.col('tran_id') == tran)
            if exons.is_empty():
                continue
            tran_reads = transcriptreads(bwfile, exons)
            if not tran_reads.is_empty():
                coordinates = tran_reads.group_by("tran_start")\
                                       .agg(pl.col("counts").sum())
                #for summary plot
                transcriptplot = coordinates.with_columns(
                    (pl.lit(tran).alias('tran_id')),
                    (pl.col('tran_start') % 3).alias('frame'))
                    
                transcriptplot = transcriptplot.to_dict(as_series=False)
                tranlist.append(transcriptplot)
                
                #for stop plot per type
                stopplot = tran_reads.group_by("tran_stop")\
                                       .agg(pl.col("counts").sum())\
                                       .with_columns(
                                           (pl.col("tran_stop") - stop)
                                           .alias("relativeloc"))\
                                       .to_dict(as_series=False)
                metagene_stop_list.append(stopplot)
                #for start plot per type
                coordinates=coordinates.with_columns(
                                           (pl.col("tran_start") - start)
                                           .alias("relativeloc"))\
                                       .to_dict(as_series=False)
                metagene_start_list.append(coordinates)

        metagene_start_df = pl.from_dicts(metagene_start_list).explode(pl.all())
        metagene_start_df = metagene_start_df.filter(pl.col("relativeloc").is_in(range_list))\
                                 .group_by("relativeloc")\
                                 .agg(pl.col("counts").sum())
        metagene_stop_df = pl.from_dicts(metagene_stop_list).explode(pl.all())
        metagene_stop_df = metagene_stop_df.filter(pl.col("relativeloc").is_in(range_list))\
                                 .group_by("relativeloc")\
                                 .agg(pl.col("counts").sum())
        
        #plot metagene per type
        fig_combined_stop = px.bar(
        metagene_stop_df.to_pandas(), x="relativeloc", y="counts", title=' '
        )
        fig_combined_stop.update_xaxes(title_text='Relative coordinates')
        fig_combined_stop.update_yaxes(title_text='Counts')
        metagene_stop = fig_combined_stop.to_html(full_html=False)

        fig_combined = px.bar(
        metagene_start_df.to_pandas(), x="relativeloc", y="counts", title=typeorf,
        )
        fig_combined.update_xaxes(title_text='Relative coordinates')
        fig_combined.update_yaxes(title_text='Counts')
        metagene_start = fig_combined.to_html(full_html=False)
        plotlist.append(metagene_start)
        plotlist.append(metagene_stop)

        df_type_filtered = df_type_filtered.to_dict(as_series=False)
        dflist.append(df_type_filtered)
    
    df_type_filtered = pl.from_dicts(dflist).explode(pl.all())\
                        .with_columns((pl.col('start') % 3).alias('frame'),
                                      (pl.col('avg')).round(2),
                                      (pl.col('nzc')).round(2),
                                      (pl.col('hrf')).round(2),
                                      (pl.col('rise_up')).round(2),
                                      (pl.col('step_down').round(2)),
                                      (pl.col('score').round(2)))
    #One plot containing all transcripts present
    transcriptplotdf = pl.from_dicts(tranlist).explode(pl.all())
    pertranscript = pertranscriptplot(transcriptplotdf)

    #barplot for frame
    summary_plot, table = summaryplot(transcriptplotdf, df_type_filtered)
    table = table.to_html(full_html=False)
    summary_plot = summary_plot.to_html(full_html=False)

    return plotlist, summary_plot, table, pertranscript


def plottop10(df, bigwig, exon, range_param, filename, parameters):
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
    
    plotlist, tranplot, table, pertranscript = metageneplot(df, bwfile, exon_df, range_list)
    
    generate_report(plotlist, tranplot, parameters, table, filename, pertranscript)

    return




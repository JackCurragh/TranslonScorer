"""
Visualization functions for TranslonScorer.

This module contains functions for generating plots and visualizations
of TranslonScorer results, including metagene profiles and transcript-specific plots.
"""

import polars as pl
import plotly.express as px
import plotly.subplots as sp
import pyBigWig as bw
import plotly.graph_objects as go
from pathlib import Path
from ..utils.logging import log_info, log_warning, log_error
from .report import generate_report
from ..file_handlers.bigwig import transcriptreads


def pertranscriptplot(df, exon_df, bwfile):
    """
    Generate plots and tables summarizing features and read counts for top ORFs of each type per transcript.

    Args:
        df (DataFrame): DataFrame containing ORF information
        exon_df (DataFrame): DataFrame containing exon information
        bwfile (str): Path to BigWig file for transcript read counts

    Returns:
        tuple: (summary_plot, table, pertranlist) containing HTML strings for plots and tables
    """
    tranlist = []
    dflist = []

    for typeorf in df["type"].unique():
        df_type_filtered = df.filter(pl.col("type") == typeorf)
        df_type_filtered = df_type_filtered.sort("score", descending=True).head(10)

        for row in range(len(df_type_filtered)):
            tran = df_type_filtered["tran_id"][row]
            exons = exon_df.filter(pl.col("tran_id") == tran)
            if exons.is_empty():
                continue

            tran_reads = transcriptreads(bwfile, exons)
            if not tran_reads.is_empty():
                coordinates = tran_reads.group_by("tran_start").agg(pl.col("counts").sum())
                # for summary plot
                transcriptplot = coordinates.with_columns(
                    pl.lit(tran).alias("tran_id"), (pl.col("tran_start") % 3).alias("frame")
                )
                transcriptplot = transcriptplot.to_dict(as_series=False)
                tranlist.append(transcriptplot)

        df_type_filtered = df_type_filtered.to_dict(as_series=False)
        dflist.append(df_type_filtered)

    df_type_filtered = (
        pl.from_dicts(dflist)
        .explode(pl.all())
        .with_columns(
            (pl.col("start") % 3).alias("frame"),
            pl.col("avg").round(2),
            pl.col("nzc").round(2),
            pl.col("hrf").round(2),
            pl.col("rise_up").round(2),
            pl.col("step_down").round(2),
            pl.col("score").round(2),
        )
        .to_pandas()
    )

    # Table generation
    table = go.Figure(
        data=[
            go.Table(
                header=dict(
                    values=list(df_type_filtered.columns), fill_color="paleturquoise", align="left"
                ),
                cells=dict(
                    values=[df_type_filtered[col] for col in df_type_filtered.columns],
                    fill_color="lavender",
                    align="left",
                ),
            )
        ]
    )
    table = table.to_html(full_html=False)

    # Summary plot
    transcriptplotdf = pl.from_dicts(tranlist).explode(pl.all())
    summary = transcriptplotdf.to_pandas()
    summary_plot = px.bar(summary, x="tran_start", y="counts", color="tran_id")
    summary_plot.update_xaxes(title_text="Transcript coordinates")
    summary_plot.update_yaxes(title_text="Counts")
    summary_plot = summary_plot.to_html(full_html=False)

    # Per-transcript plots
    pertranlist = []
    for tran in transcriptplotdf["tran_id"].unique():
        transcriptdf = transcriptplotdf.filter(pl.col("tran_id") == tran)
        fig = px.bar(transcriptdf.to_pandas(), x="tran_start", y="counts", title=tran)
        fig.update_xaxes(title_text="Transcript coordinates")
        fig.update_yaxes(title_text="Counts")
        fig = fig.to_html(full_html=False)
        pertranlist.append(fig)

    return summary_plot, table, pertranlist


def metageneplot(df, bwfile, exon_df, range_list):
    """
    Generate metagene plots for each type of ORF.

    Args:
        df (DataFrame): DataFrame containing ORF information
        bwfile (str): Path to BigWig file
        exon_df (DataFrame): DataFrame containing exon information
        range_list (list): List of integers for relative coordinates range

    Returns:
        list: HTML strings of metagene plots for each ORF type
    """
    plotlist = []

    for typeorf in df["type"].unique():
        df_type_filtered = df.filter(pl.col("type") == typeorf)
        metagene_start_dict = {i: 0 for i in range_list}
        metagene_stop_dict = {i: 0 for i in range_list}

        for tran in df_type_filtered["tran_id"].unique():
            df_tran = df_type_filtered.filter(pl.col("tran_id") == tran)
            exons = exon_df.filter(pl.col("tran_id") == tran)
            if exons.is_empty():
                continue

            tran_reads = transcriptreads(bwfile, exons)
            if not tran_reads.is_empty():
                starts = df_tran.get_column("start").to_list()
                stops = df_tran.get_column("stop").to_list()

                # Start plot per type - transcriptreads returns tran_start
                startplot = tran_reads.group_by("tran_start").agg(pl.col("counts").sum())
                for start in starts:
                    startplot_final = (
                        startplot.with_columns((pl.col("tran_start") - start).alias("relativeloc"))
                        .filter(pl.col("relativeloc").is_in(range_list))
                        .group_by("relativeloc")
                        .agg(pl.col("counts").sum())
                        .to_dict(as_series=False)
                    )
                    if startplot_final:
                        startplot_final = dict(
                            zip(startplot_final["relativeloc"], startplot_final["counts"])
                        )
                        for i in startplot_final:
                            metagene_start_dict[i] += startplot_final[i]

                # Stop plot per type - use tran_start also for stop analysis
                # Since transcriptreads doesn't return tran_stop column
                stopplot = tran_reads.group_by("tran_start").agg(pl.col("counts").sum())
                for stop in stops:
                    stopplot_final = (
                        stopplot.with_columns((pl.col("tran_start") - stop).alias("relativeloc"))
                        .filter(pl.col("relativeloc").is_in(range_list))
                        .group_by("relativeloc")
                        .agg(pl.col("counts").sum())
                        .to_dict(as_series=False)
                    )
                    if stopplot_final:
                        stopplot_final = dict(
                            zip(stopplot_final["relativeloc"], stopplot_final["counts"])
                        )
                        for i in stopplot_final:
                            metagene_stop_dict[i] += stopplot_final[i]

        # Create plots
        start_dict = {
            "relativeloc": list(metagene_start_dict.keys()),
            "counts": list(metagene_start_dict.values()),
        }
        stop_dict = {
            "relativeloc": list(metagene_stop_dict.keys()),
            "counts": list(metagene_stop_dict.values()),
        }

        # Plot metagene per type
        fig_combined_stop = px.bar(stop_dict, x="relativeloc", y="counts")
        fig_combined_stop.update_xaxes(title_text="Relative coordinates")
        fig_combined_stop.update_yaxes(title_text="Counts")
        metagene_stop = fig_combined_stop.to_html(full_html=False)

        fig_combined = px.bar(start_dict, x="relativeloc", y="counts", title=typeorf)
        fig_combined.update_xaxes(title_text="Relative coordinates")
        fig_combined.update_yaxes(title_text="Counts")
        metagene_start = fig_combined.to_html(full_html=False)

        plotlist.append(metagene_start)
        plotlist.append(metagene_stop)

    return plotlist


def plottop10(df, bigwig, exon, range_param, filename, parameters=None):
    """
    Generate plots and tables summarizing top 10 ORFs per type and metagene profiles.

    Args:
        df (str or DataFrame): Path to CSV file or DataFrame containing ORF information
        bigwig (str): Path to BigWig file
        exon (str or DataFrame): Path to CSV file or DataFrame containing exon information
        range_param (int): Range for relative coordinates
        filename (str): Output filename for report
        parameters (dict, optional): Additional parameters for report generation
    """
    log_info("Generating plots for top 10 ORFs...")

    range_list = list(range(-range_param, range_param + 1))

    # Read input files if paths are provided, otherwise use the DataFrame directly
    if isinstance(df, str):
        df = pl.read_csv(df, has_header=True, separator=",")

    bwfile = bw.open(bigwig)

    if isinstance(exon, str):
        exon_df = pl.read_csv(exon, has_header=True, separator=",")
    else:
        exon_df = exon  # Use the DataFrame directly

    # Ensure exon_df has the correct column formats
    exon_df = exon_df.with_columns(
        pl.col("start", "stop", "tran_start", "tran_stop").apply(
            lambda x: x.split(",") if isinstance(x, str) else x
        )
    )

    # Generate plots
    plotlist = metageneplot(df, bwfile, exon_df, range_list)
    tranplot, table, pertranscript = pertranscriptplot(df, exon_df, bwfile)

    # Generate report
    generate_report(plotlist, tranplot, parameters, table, filename, pertranscript)

    log_info("Plots generated successfully")

"""Script to plot file"""

import polars as pl
import plotly.express as px
import plotly.subplots as sp
import pyBigWig as bw
import plotly.graph_objects as go

from .bigwigtodf import transcriptreads
from .report import generate_report


def pertranscriptplot(df, exon_df, bwfile):
    """
    Generate plots and tables summarizing features and read counts for top ORFs of each type per transcript.

    Parameters:
    - df (DataFrame): Pandas DataFrame containing ORF information, including columns like 'tran_id', 'start', 'stop',
                     'length', 'startorf', 'stoporf', 'type', 'rise_up', 'step_down', 'hrf', 'avg', 'nzc', 'score'.
    - exon_df (DataFrame): Pandas DataFrame containing exon information, with columns including 'tran_id', 'start',
                           'stop', 'tran_start', 'tran_stop'.
    - bwfile (str): Path to the BigWig file used for obtaining transcript read counts.

    Returns:
    - tuple: A tuple containing:
        - list: List of HTML strings representing individual transcript plots for top ORFs of each type.
        - str: HTML string representing a summary plot showing transcript read counts for top ORFs.
        - str: HTML string representing a table summarizing features of top ORFs.

    This function generates plots and tables summarizing features and read counts for top ORFs of each type per transcript.
    It iterates over each ORF type in the input DataFrame `df`, filters the top 10 ORFs based on the 'score' column, and
    retrieves corresponding exon information from `exon_df`.

    For each top ORF, it calculates transcript read counts using the `transcriptreads` function with data from `bwfile`.
    It then creates:
    - A summary plot (`summary_plot`) using Plotly Express (`px.bar`) showing read counts across transcript coordinates.
    - A table (`table`) using Plotly Graph Objects (`go.Table`) summarizing features such as start, stop, length, type,
      scores, and other metrics for each ORF.
    - Individual transcript-specific plots (`pertranlist`) for each transcript showing read counts across transcript
      coordinates for the top ORFs.

    These plots and tables are converted to HTML strings (`summary_plot`, `table`, and `pertranlist`) for integration
    into web applications or reports.

    Note: This function assumes the use of Pandas (`pl`), Plotly (`px`, `go`), and functions like `transcriptreads` for
    data manipulation and visualization.
    """
    for typeorf in df["type"].unique():
        df_type_filtered = df.filter(pl.col("type") == typeorf)
        df_type_filtered = df_type_filtered.sort("score", descending=True).head(10)
        tranlist = []
        dflist = []
        for row in range(len(df_type_filtered)):
            tran = df_type_filtered["tran_id"][row]
            exons = exon_df.filter(pl.col("tran_id") == tran)
            if exons.is_empty():
                continue
            tran_reads = transcriptreads(bwfile, exons)
            if not tran_reads.is_empty():
                coordinates = tran_reads.group_by("tran_start").agg(
                    pl.col("counts").sum()
                )
                # for summary plot
                transcriptplot = coordinates.with_columns(
                    (pl.lit(tran).alias("tran_id")),
                    (pl.col("tran_start") % 3).alias("frame"),
                )
                transcriptplot = transcriptplot.to_dict(as_series=False)
                tranlist.append(transcriptplot)
        df_type_filtered = df_type_filtered.to_dict(as_series=False)
        dflist.append(df_type_filtered)

    df_type_filtered = pl.from_dicts(dflist).explode(pl.all()).with_columns(
        (pl.col("start") % 3).alias("frame"),
        (pl.col("avg")).round(2),
        (pl.col("nzc")).round(2),
        (pl.col("hrf")).round(2),
        (pl.col("rise_up")).round(2),
        (pl.col("step_down").round(2)),
        (pl.col("score").round(2)),
        ).to_pandas()
    # Table generation
    table = go.Figure(
        data=[
            go.Table(
                header=dict(
                    values=list(df_type_filtered.columns),
                    fill_color="paleturquoise",
                    align="left",
                ),
                cells=dict(
                    values=[
                        df_type_filtered.tran_id,
                        df_type_filtered.start,
                        df_type_filtered.stop,
                        df_type_filtered.length,
                        df_type_filtered.startorf,
                        df_type_filtered.stoporf,
                        df_type_filtered.type,
                        df_type_filtered.rise_up,
                        df_type_filtered.step_down,
                        df_type_filtered.hrf,
                        df_type_filtered.avg,
                        df_type_filtered.nzc,
                        df_type_filtered.score,
                        df_type_filtered.frame,
                    ],
                    fill_color="lavender",
                    align="left",
                ),
            )
        ]
    )
    table = table.to_html(full_html=False)

    transcriptplotdf = pl.from_dicts(tranlist).explode(pl.all())
    #summary plot
    summary = transcriptplotdf.to_pandas()
    summary_plot = px.bar(summary, x="tran_start", y="counts", color="tran_id")
    summary_plot.update_xaxes(title_text="Transcript coordinates")
    summary_plot.update_yaxes(title_text="Counts")
    summary_plot = summary_plot.to_html(full_html=False)
    #transcript plots
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
    Generate metagene plots for each type of ORF based on transcript read counts relative to exon coordinates.

    Parameters:
    - df (DataFrame): Pandas DataFrame containing ORF information, including columns like 'tran_id', 'start', 'stop',
                     'type'.
    - bwfile (str): Path to the BigWig file used for obtaining transcript read counts.
    - exon_df (DataFrame): Pandas DataFrame containing exon information, with columns including 'tran_id', 'start',
                           'stop', 'tran_start', 'tran_stop'.
    - range_list (list): List of integers representing the range of relative coordinates around exon boundaries
                         for plotting metagene profiles.

    Returns:
    - list: A list of HTML strings representing metagene plots for each type of ORF.

    This function generates metagene plots for each type of ORF (e.g., uORF, CDS, dORF) based on transcript read counts
    relative to exon coordinates. It iterates over each ORF type in the input DataFrame `df`, filters ORFs by type, and
    calculates metagene profiles for start and stop positions relative to exon boundaries.

    For each ORF type:
    - It initializes dictionaries (`metagene_start_dict` and `metagene_stop_dict`) to accumulate counts of start and stop
      positions relative to exon boundaries within the specified `range_list`.
    - For each transcript associated with the current ORF type, it retrieves exon information from `exon_df` and calculates
      transcript read counts using the `transcriptreads` function with data from `bwfile`.
    - It generates metagene profiles (`start_dict` and `stop_dict`) by aggregating read counts of start and stop positions
      relative to exon boundaries across all transcripts of the current ORF type.
    - It creates Plotly Express bar charts (`fig_combined_stop` and `fig_combined`) for start and stop positions, respectively,
      and converts them to HTML strings (`metagene_stop` and `metagene_start`).

    The function returns a list (`plotlist`) containing HTML strings of metagene plots for each type of ORF, ready for
    integration into web applications or reports.

    Note: This function assumes the use of Pandas (`pl`), Plotly Express (`px`), and functions like `transcriptreads` for
    data manipulation and visualization.
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
                # for start plot per type
                startplot = tran_reads.group_by("tran_start").agg(
                    pl.col("counts").sum()
                )
                for start in starts:
                    startplot_final = (
                        startplot.with_columns(
                            (pl.col("tran_start") - start).alias("relativeloc")
                        )
                        .filter(pl.col("relativeloc").is_in(range_list))
                        .group_by("relativeloc")
                        .agg(pl.col("counts").sum())
                        .to_dict(as_series=False)
                    )
                    if startplot_final:
                        startplot_final = dict(
                            zip(
                                startplot_final["relativeloc"],
                                startplot_final["counts"],
                            )
                        )
                        for i in startplot_final:
                            metagene_start_dict[i] += startplot_final[i]

                # for stop plot per type
                stopplot = tran_reads.group_by("tran_stop").agg(pl.col("counts").sum())
                for stop in stops:
                    stopplot_final = (
                        stopplot.with_columns(
                            (pl.col("tran_stop") - stop).alias("relativeloc")
                        )
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

        start_dict = {
            "relativeloc": list(metagene_start_dict.keys()),
            "counts": list(metagene_start_dict.values()),
        }
        stop_dict = {
            "relativeloc": list(metagene_stop_dict.keys()),
            "counts": list(metagene_stop_dict.values()),
        }
        # plot metagene per type
        fig_combined_stop = px.bar(start_dict, x="relativeloc", y="counts")
        fig_combined_stop.update_xaxes(title_text="Relative coordinates")
        fig_combined_stop.update_yaxes(title_text="Counts")
        metagene_stop = fig_combined_stop.to_html(full_html=False)

        fig_combined = px.bar(
            stop_dict,
            x="relativeloc",
            y="counts",
            title=typeorf,
        )
        fig_combined.update_xaxes(title_text="Relative coordinates")
        fig_combined.update_yaxes(title_text="Counts")
        metagene_start = fig_combined.to_html(full_html=False)
        plotlist.append(metagene_start)
        plotlist.append(metagene_stop)

    return plotlist


def plottop10(df, bigwig, exon, range_param, filename, parameters):
    """
    Generate plots and tables summarizing top 10 ORFs per type and metagene profiles based on transcript read counts.

    Parameters:
    - df (str): Path to a CSV file containing ORF information, with headers and columns like 'tran_id', 'start', 'stop',
               'type', 'score', etc.
    - bigwig (str): Path to the BigWig file used for obtaining transcript read counts.
    - exon (str): Path to a CSV file containing exon information, with headers and columns like 'start', 'stop',
                  'tran_start', 'tran_stop', 'tran_id', etc.
    - range_param (int): Integer representing the range of relative coordinates around exon boundaries for plotting
                         metagene profiles.
    - filename (str): Name of the output file for the generated report.
    - parameters (dict): Dictionary containing additional parameters or settings for generating the report.

    Returns:
    - None

    This function coordinates the generation of plots and tables summarizing top 10 ORFs per type and metagene profiles
    based on transcript read counts relative to exon coordinates. It performs the following steps:

    1. Constructs a `range_list` of integers representing the range of relative coordinates around exon boundaries for
       plotting metagene profiles.
    2. Reads ORF data from the CSV file (`df`) into a Pandas DataFrame (`df`).
    3. Opens the BigWig file (`bigwig`) to obtain transcript read counts (`bwfile`).
    4. Reads exon information from the CSV file (`exon`) into a Pandas DataFrame (`exon_df`). Splits certain columns
       ('start', 'stop', 'tran_start', 'tran_stop') into lists of integers using `apply` and `split`.
    5. Calls the `metageneplot` function to generate metagene profiles (`plotlist`) for each type of ORF based on
       transcript read counts relative to exon coordinates.
    6. Calls the `pertranscriptplot` function to generate individual transcript plots (`tranplot`), a summary table
       (`table`), and transcript-specific plots (`pertranscript`) for the top 10 ORFs per type.
    7. Calls `generate_report` (assumed to be defined elsewhere) to compile all generated plots and tables into a report,
       using provided `parameters` and saving it under `filename`.

    Note: This function assumes the use of Pandas (`pl`), functions like `metageneplot` and `pertranscriptplot` for
    data manipulation and visualization, and a `generate_report` function for compiling the final report.
    """
    range_list = list(range(-range_param, range_param + 1))
    df = pl.read_csv(df, has_header=True, separator=",")
    bwfile = bw.open(bigwig)
    exon_df = pl.read_csv(exon, has_header=True, separator=",")
    exon_df = exon_df.with_columns(
        pl.col("start", "stop", "tran_start", "tran_stop").apply(lambda x: x.split(","))
    )

    plotlist = metageneplot(df, bwfile, exon_df, range_list)

    tranplot, table, pertranscript = pertranscriptplot(df, exon_df, bwfile)
    generate_report(plotlist, tranplot, parameters, table, filename, pertranscript)

    return

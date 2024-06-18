# /usr/bin/python3
"""Script to specify arguments on CLI"""
import os
import click
import polars as pl
import warnings

from .readfiles import readbam
from .fileprocessor import dftobed, bedtobigwig
from .getcandidates import gettranscripts, preporfs, orfrelativeposition
from .filewriter import saveorfsandexons
from .bigwigtodf import scoring
from .plotting import plottop10
from .report import getparameters

warnings.filterwarnings("ignore")


@click.group()
def function():
    pass


@function.command()
@click.option("--bam", "-b", help="Provide a bam file")
@click.option(
    "--chromsize", "-c", help="Provide a file containing the chromosome sizes"
)
@click.option(
    "--seq", "-s", help="Provide a file containing the genomic sequence (.fa)"
)
@click.option(
    "--tran", "-t", help="Provide a file containing the transcript sequences (.fa)"
)
@click.option("--ann", "-a", help="Provide a file containing the annotation (.gtf)")
@click.option("--starts", "-sta", default="ATG", help="Provide a list of start codons")
@click.option(
    "--stops", "-stp", default="TAA,TAG,TGA", help="Provide a list of stop codons"
)
@click.option("--minlen", "-min", default=0, help="Provide the minimum length")
@click.option("--maxlen", "-max", default=1000000, help="Provide the maximum length")
@click.option("--bigwig", "-bw", help="Provide a Bigwig file to convert")
@click.option("--exon", "-ex", help="Provide a file containing exon positions")
@click.option("--bedfile", "-bw", help="Provide a Bigwig file to convert")
@click.option("--orfs", "-of", help="Provide a file containing annotated orfs")
@click.option(
    "--range_param",
    "-rp",
    default=30,
    help="Provide an integer that indicates the range in which a plot will be constructed \
                 around the relative start position",
)
@click.option(
    "--sru_range",
    "-sru",
    default=15,
    help="Provide an integer that indicates the range for \
             the Start Rise Up score. This sets the amount of nucleotides before and after \
             the stop codon will regarded when calculating",
)
@click.option("--offsets", "-ofs", help="Provide a file containing offset parameters")
@click.option(
    "--scoretype",
    "-s",
    default=False,
    help="Select the scoring algorithm, Default = False (old scoring algorithm)",
)
@click.option(
    "--plotfile",
    "-pf",
    help="Provide a '.csv' file containing scored ORFs to use for plotting",
)
@click.option(
    "--outfilename",
    "-ofn",
    help="Provide a name for the files that are written using this tool",
)
def translonpredictor(
    bam,
    bedfile,
    chromsize,
    bigwig,
    seq,
    tran,
    ann,
    starts,
    stops,
    minlen,
    maxlen,
    exon,
    orfs,
    range_param,
    sru_range,
    offsets,
    scoretype,
    plotfile,
    outfilename,
):
    """
    Perform transcript-level analysis pipeline for predicting translons and generating reports.

    Parameters:
    - bam (str): Path to BAM file containing transcriptome sequencing data.
    - bedfile (str): Path to BedGraph file containing genomic coverage information.
    - chromsize (str): Path to file containing chromosome sizes.
    - bigwig (str): Path to BigWig file containing genome-wide read coverage information.
    - seq (str): Path to file containing transcript sequences in FASTA format.
    - tran (str): Specific transcript identifier to analyze.
    - ann (str): Path to annotation file (e.g., GTF format) containing genomic features.
    - starts (str): Comma-separated string of start codons to consider during ORF prediction.
    - stops (str): Comma-separated string of stop codons to consider during ORF prediction.
    - minlen (int): Minimum length threshold for predicted ORFs.
    - maxlen (int): Maximum length threshold for predicted ORFs.
    - exon (str): Path to file containing exon information.
    - orfs (str): Path to file containing pre-annotated ORFs.
    - range_param (int): Parameter for specifying the range around ORFs for metagene analysis.
    - sru_range (int): Range parameter for Start Rise Up (SRU) scoring.
    - offsets (str): Comma-separated string of offsets to apply during transcriptome analysis.
    - scoretype (str): Type of scoring method to apply (e.g., HRF, average, NZC).
    - plotfile (str): Path to file containing data for generating plots.
    - outfilename (str): Output filename prefix for generated files and reports.

    Raises:
    - Exception: If required input parameters are not provided or invalid combinations are used.

    This function orchestrates a pipeline for predicting translons and generating reports based on the provided parameters:
    - Processes BAM, BedGraph, and BigWig files for read coverage and transcript annotations.
    - Extracts transcript sequences or uses specific transcript identifiers.
    - Predicts candidate ORFs based on provided start and stop codons and length thresholds.
    - Scores ORFs using specified scoring methods and parameters.
    - Generates and saves reports, including top-10 ORF plots and metagene plots.

    Example:
    >>> translonpredictor(bam='transcriptome.bam', ann='annotation.gtf', starts='ATG', stops='TAA,TAG,TGA',
                          minlen=50, maxlen=500, exon='exon_info.csv', range_param=50, sru_range=30,
                          scoretype='HRF', plotfile='orf_scores.csv', outfilename='translon_results')

    This example runs the pipeline using a BAM file, GTF annotation, specified start and stop codons, and other parameters,
    and generates ORF scoring reports and plots as 'translon_results_orfs_scored.csv' and associated HTML reports.

    """
    parameters = getparameters(vars())

    if bam or chromsize or bedfile:
        if bam and chromsize and ann and outfilename:
            print("Processing BAM file")
            # READ IN BAM START
            location = os.getcwd() + "/" + bam
            # if file is provided
            if os.path.isfile(location):
                # read in bam file
                df = readbam(location)
                # calculate asite + converting to BedGraph
                print("Calculating and applying offsets")
                beddf, exondf, cdsdf = dftobed(df, ann, offsets)
                print("Writing bed file")
                if not os.path.exists(f"{outfilename}.bedGraph"):
                    beddf.write_csv(
                        f"{outfilename}.bedGraph", separator="\t", include_header=False
                    )

                    # Converting Bedgrapgh to Bigwig format
                print("Writing bigwig file")
                bigwig = bedtobigwig(f"{outfilename}.bedGraph", chromsize, outfilename)

        elif bedfile and chromsize and outfilename:
            print("Writing bigwig file")
            bigwig = bedtobigwig(bedfile, chromsize, outfilename)

        elif bigwig:
            pass
        else:
            raise Exception(
                "Must provide valuable input to . The options are the following:\n \
                1. Bam file (.bam) + file containing chromosome information\n \
                2. BedGraph (.bedGraph) file + file containing chromosome information\n \
                3. Bigwig (.bw) file \n"
                "Please do not forget to always provide a filename for any files that may be written during the process"
            )

    if seq or tran:
        if seq and ann and outfilename:
            print("Extracting transcripts")
            transcript = gettranscripts(seq, ann, outfilename)
        elif tran and ann and outfilename:
            transcript = tran
        print("Getting candidate ORFs")
        orfdf = preporfs(
            transcript, starts.split(","), stops.split(","), minlen, maxlen
        )
        if not "cdsdf" in globals():
            cdsdf = 0
        orf_ann_df, exon_df = orfrelativeposition(ann, orfdf, cdsdf)
        orfs, exon = saveorfsandexons(orf_ann_df, exon_df, outfilename)

        if not bigwig:
            bigwig = f"{outfilename}.bw"
        print("Scoring ORFs")
        scoredorfs = scoring(bigwig, exon, orfs, scoretype, sru_range)
        scoredorfs.write_csv(f"{outfilename}_orfs_scored.csv")

        plotfile = f"{outfilename}_orfs_scored.csv"
        print("Generating report")
        plottop10(plotfile, bigwig, exon, range_param, outfilename, parameters)

    elif orfs and exon and bigwig and outfilename:
        print("Scoring orfs")
        scoredorfs = scoring(bigwig, exon, orfs, scoretype, sru_range)
        scoredorfs.write_csv(f"{outfilename}_orfs_scored.csv")

        plotfile = f"{outfilename}_orfs_scored.csv"

        print("Generating report")
        plottop10(plotfile, bigwig, exon, range_param, outfilename, parameters)

    elif plotfile and bigwig and exon and outfilename:
        print("Generating report")
        plottop10(plotfile, bigwig, exon, range_param, outfilename, parameters)

    else:
        raise Exception(
            "Must provide valuable input for extracting ORFS. The options are the following:\n \
                1. A file containing FASTA sequences (.fa) + annotation file (.gtf)\n \
                2. A file containing transcript sequences (.fa) + annotation file (.gtf)\n \
                3. A file containing annotated ORFS (.csv) + file containg exon information (.csv)\n"
            "Please do not forget to always provide a filename for any files that may be written during the process"
        )


if __name__ == "__main__":
    function()

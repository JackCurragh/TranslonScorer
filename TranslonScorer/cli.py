"""Command-line interface for TranslonScorer."""

import os
import click
import polars as pl
import warnings
from typing import List, Optional

from .readfiles import readbam
from .fileprocessor import dftobed, bedtobigwig
from .getcandidates import gettranscripts, preporfs, orfrelativeposition
from .filewriter import saveorfsandexons
from .bigwigtodf import scoring
from .plotting import plottop10
from .report import getparameters

warnings.filterwarnings("ignore")

@click.group()
def cli():
    """TranslonScorer: A tool for identifying and scoring translational events from Ribo-seq data.
    
    This tool provides several workflows:
    
    1. all: Run the complete pipeline end-to-end
    2. process-bam: Process Ribo-seq BAM files to generate coverage tracks
    3. find-orfs: Identify and score potential ORFs from sequence data
    4. score-orfs: Score existing ORFs using coverage data
    5. plot: Generate visualization reports from scored ORFs
    
    For detailed instructions, use --help with any command:
    translonpredictor all --help
    """
    pass

@cli.command()
@click.option('--bam', '-b', required=True,
              help='Input BAM file from Ribo-seq data. Supports both genomic and transcriptomic alignments (required)')
@click.option('--chromsizes', '-c', required=True,
              help='Chromosome sizes file (required for bigWig conversion)')
@click.option('--sequence', '-s', required=True,
              help='Input FASTA file (genomic or transcriptomic)')
@click.option('--annotation', '-a', required=True,
              help='GTF annotation file (required)')
@click.option('--offsets', '-off',
              help='File containing read length-specific offsets for A-site calculation')
@click.option('--start-codons', default="ATG",
              help='Comma-separated list of start codons (default: ATG)')
@click.option('--stop-codons', default="TAA,TAG,TGA",
              help='Comma-separated list of stop codons (default: TAA,TAG,TGA)')
@click.option('--min-len', type=int, default=0,
              help='Minimum ORF length in nucleotides (default: 0)')
@click.option('--max-len', type=int, default=1000000,
              help='Maximum ORF length in nucleotides (default: 1000000)')
@click.option('--sru-range', type=int, default=15,
              help='Nucleotide range for Start Rise Up score calculation (default: 15)')
@click.option('--scoring-method', type=click.Choice(['classic', 'modern']), 
              default='modern', help='Scoring algorithm to use (default: modern)')
@click.option('--plot-range', type=int, default=30,
              help='Plot range around start position (default: 30)')
@click.option('--outfile', '-o', required=True,
              help='Base name for output files')
def all(bam: str, chromsizes: str, sequence: str, annotation: str,
        outfile: str, offsets: Optional[str] = None,
        start_codons: str = "ATG", stop_codons: str = "TAA,TAG,TGA",
        min_len: int = 0, max_len: int = 1000000,
        sru_range: int = 15, scoring_method: str = 'modern',
        plot_range: int = 30):
    """Run the complete TranslonScorer pipeline end-to-end.
    
    This command runs all steps of the pipeline in sequence:
    1. Process Ribo-seq BAM file to generate coverage tracks
    2. Extract transcripts
    3. Find and score potential ORFs
    4. Generate visualization reports
    
    Required files:
    - BAM file: Ribo-seq reads aligned to either genome or transcriptome
    - Chromosome sizes: Tab-separated file with chr\tsize
    - Sequence: FASTA file (genomic or transcriptomic)
    - Annotation: GTF file with transcript annotations
    
    Example:
    translonpredictor all -b ribo.bam -c chrom.sizes -s genome.fa -a anno.gtf -o output
    
    Note: This command will generate all intermediate files with the specified output prefix.
    The tool automatically detects whether the BAM file contains genomic or transcriptomic alignments.
    """
    click.echo("Starting complete TranslonScorer pipeline...")
    
    # Step 1: Process BAM file
    click.echo("\nStep 1/4: Processing BAM file")
    location = os.path.abspath(bam)
    if not os.path.isfile(location):
        raise click.BadParameter(f"BAM file not found: {bam}")
    
    df = readbam(location)
    click.echo("Calculating and applying offsets")
    beddf, exondf, cdsdf = dftobed(df, annotation, offsets)
    
    bedgraph_path = f"{outfile}.bedGraph"
    click.echo(f"Writing bedGraph file: {bedgraph_path}")
    if not os.path.exists(bedgraph_path):
        beddf.write_csv(bedgraph_path, separator="\t", include_header=False)
    
    bigwig_path = f"{outfile}.bw"
    click.echo(f"Writing bigWig file: {bigwig_path}")
    bedtobigwig(bedgraph_path, chromsizes, outfile)
    
    # Step 2: Extract transcripts
    click.echo("\nStep 2/4: Extracting transcripts")
    transcript = gettranscripts(sequence, annotation, outfile)
    
    # Step 3: Find and score ORFs
    click.echo("\nStep 3/4: Finding and scoring ORFs")
    click.echo("Finding candidate ORFs")
    orfdf = preporfs(
        transcript, 
        start_codons.split(","), 
        stop_codons.split(","), 
        min_len, 
        max_len
    )
    
    click.echo("Mapping ORFs to transcript coordinates")
    orf_ann_df, exon_df = orfrelativeposition(annotation, orfdf, cdsdf)
    orfs, exon = saveorfsandexons(orf_ann_df, exon_df, outfile)
    
    click.echo("Scoring ORFs")
    scoredorfs = scoring(bigwig_path, exon, orfs, scoring_method == 'modern', sru_range)
    scoredorfs.write_csv(f"{outfile}_orfs_scored.csv")
    
    # Step 4: Generate report
    click.echo("\nStep 4/4: Generating visualization report")
    plottop10(f"{outfile}_orfs_scored.csv", bigwig_path, exon, plot_range, outfile, getparameters(locals()))
    
    click.echo("\nPipeline completed successfully!")
    click.echo(f"Output files generated with prefix: {outfile}")
    click.echo("Files generated:")
    click.echo(f"  - {outfile}.bedGraph: Coverage in bedGraph format")
    click.echo(f"  - {outfile}.bw: Coverage in bigWig format")
    click.echo(f"  - {outfile}_orfs_scored.csv: Scored ORFs")
    click.echo(f"  - {outfile}_report.html: Visualization report")

@cli.command()
@click.option('--bam', '-b', required=True,
              help='Input BAM file from Ribo-seq data. Supports both genomic and transcriptomic alignments (required)')
@click.option('--chromsizes', '-c', required=True,
              help='Chromosome sizes file (required for bigWig conversion)')
@click.option('--annotation', '-a', required=True,
              help='GTF annotation file (required)')
@click.option('--offsets', '-off',
              help='File containing read length-specific offsets for A-site calculation')
@click.option('--outfile', '-o', required=True,
              help='Base name for output files')
def process_bam(bam: str, chromsizes: str, annotation: str, 
                outfile: str, offsets: Optional[str] = None):
    """Process Ribo-seq BAM files to generate coverage tracks.
    
    This command processes a Ribo-seq BAM file to:
    1. Calculate A-site positions
    2. Generate coverage tracks (bedGraph and bigWig)
    3. Extract exon and CDS information
    
    Required files:
    - BAM file: Ribo-seq reads aligned to either genome or transcriptome
    - Chromosome sizes: Tab-separated file with chr\tsize
    - Annotation: GTF file with transcript annotations
    
    Example:
    translonpredictor process-bam -b sample.bam -c chrom.sizes -a anno.gtf -o output
    
    Note: The tool automatically detects whether the BAM file contains genomic or 
    transcriptomic alignments and processes it accordingly.
    """
    print("Processing BAM file")
    location = os.path.abspath(bam)
    
    if not os.path.isfile(location):
        raise click.BadParameter(f"BAM file not found: {bam}")
        
    # Read in BAM file
    df = readbam(location)
    
    # Calculate A-site and convert to BedGraph
    print("Calculating and applying offsets")
    beddf, exondf, cdsdf = dftobed(df, annotation, offsets)
    
    # Write bedGraph file
    bedgraph_path = f"{outfile}.bedGraph"
    print(f"Writing bedGraph file: {bedgraph_path}")
    if not os.path.exists(bedgraph_path):
        beddf.write_csv(bedgraph_path, separator="\t", include_header=False)
    
    # Convert to bigWig
    bigwig_path = f"{outfile}.bw"
    print(f"Writing bigWig file: {bigwig_path}")
    bedtobigwig(bedgraph_path, chromsizes, outfile)
    
    print("BAM processing complete!")

@cli.command()
@click.option('--sequence', '-s', required=True,
              help='Input FASTA file (genomic or transcriptomic)')
@click.option('--annotation', '-a', required=True,
              help='GTF annotation file (required)')
@click.option('--bigwig', '-bw', required=True,
              help='BigWig file containing Ribo-seq coverage')
@click.option('--start-codons', default="ATG",
              help='Comma-separated list of start codons (default: ATG)')
@click.option('--stop-codons', default="TAA,TAG,TGA",
              help='Comma-separated list of stop codons (default: TAA,TAG,TGA)')
@click.option('--min-len', type=int, default=0,
              help='Minimum ORF length in nucleotides (default: 0)')
@click.option('--max-len', type=int, default=1000000,
              help='Maximum ORF length in nucleotides (default: 1000000)')
@click.option('--outfile', '-o', required=True,
              help='Base name for output files')
@click.option('--sru-range', type=int, default=15,
              help='Nucleotide range for Start Rise Up score calculation (default: 15)')
@click.option('--scoring-method', type=click.Choice(['classic', 'modern']), 
              default='modern', help='Scoring algorithm to use (default: modern)')
def find_orfs(sequence: str, annotation: str, bigwig: str, outfile: str,
              start_codons: str, stop_codons: str, min_len: int, max_len: int,
              sru_range: int, scoring_method: str):
    """Identify and score potential ORFs from sequence data.
    
    This command takes sequence data and:
    1. Identifies potential ORFs based on start/stop codons
    2. Maps them to transcript coordinates
    3. Scores them using Ribo-seq coverage
    
    Required files:
    - Sequence: FASTA file (genomic or transcriptomic)
    - Annotation: GTF file with transcript annotations
    - BigWig: Ribo-seq coverage from process-bam
    
    Example:
    translonpredictor find-orfs -s genome.fa -a anno.gtf -bw coverage.bw -o output
    """
    print("Extracting transcripts")
    transcript = gettranscripts(sequence, annotation, outfile)
    
    print("Finding candidate ORFs")
    orfdf = preporfs(
        transcript, 
        start_codons.split(","), 
        stop_codons.split(","), 
        min_len, 
        max_len
    )
    
    print("Mapping ORFs to transcript coordinates")
    orf_ann_df, exon_df = orfrelativeposition(annotation, orfdf, None)
    orfs, exon = saveorfsandexons(orf_ann_df, exon_df, outfile)
    
    print("Scoring ORFs")
    scoredorfs = scoring(bigwig, exon, orfs, scoring_method == 'modern', sru_range)
    scoredorfs.write_csv(f"{outfile}_orfs_scored.csv")
    
    print("Generating report")
    plottop10(f"{outfile}_orfs_scored.csv", bigwig, exon, 30, outfile, getparameters(locals()))
    
    print("ORF finding and scoring complete!")

@cli.command()
@click.option('--orfs', '-f', required=True,
              help='CSV file containing pre-annotated ORFs')
@click.option('--bigwig', '-bw', required=True,
              help='BigWig file containing Ribo-seq coverage')
@click.option('--exons', '-e', required=True,
              help='CSV file containing exon positions')
@click.option('--outfile', '-o', required=True,
              help='Base name for output files')
@click.option('--scoring-method', type=click.Choice(['classic', 'modern']), 
              default='modern', help='Scoring algorithm to use (default: modern)')
@click.option('--sru-range', type=int, default=15,
              help='Nucleotide range for Start Rise Up score calculation (default: 15)')
def score_orfs(orfs: str, bigwig: str, exons: str, outfile: str,
               scoring_method: str, sru_range: int):
    """Score existing ORFs using Ribo-seq coverage data.
    
    This command takes pre-annotated ORFs and scores them using:
    1. Ribo-seq coverage from bigWig
    2. Exon position information
    3. Specified scoring parameters
    
    Required files:
    - ORFs: CSV file with ORF annotations
    - BigWig: Ribo-seq coverage
    - Exons: CSV with exon positions
    
    Example:
    translonpredictor score-orfs -f orfs.csv -bw coverage.bw -e exons.csv -o output
    """
    print("Scoring ORFs")
    scoredorfs = scoring(bigwig, exons, orfs, scoring_method == 'modern', sru_range)
    scoredorfs.write_csv(f"{outfile}_orfs_scored.csv")
    
    print("Generating report")
    plottop10(f"{outfile}_orfs_scored.csv", bigwig, exons, 30, outfile, getparameters(locals()))
    
    print("ORF scoring complete!")

@cli.command()
@click.option('--scored-orfs', '-s', required=True,
              help='CSV file containing scored ORFs')
@click.option('--bigwig', '-bw', required=True,
              help='BigWig file containing Ribo-seq coverage')
@click.option('--exons', '-e', required=True,
              help='CSV file containing exon positions')
@click.option('--plot-range', type=int, default=30,
              help='Plot range around start position (default: 30)')
@click.option('--outfile', '-o', required=True,
              help='Base name for output files')
def plot(scored_orfs: str, bigwig: str, exons: str, 
         plot_range: int, outfile: str):
    """Generate visualization reports from scored ORFs.
    
    This command creates visualization reports including:
    1. Top 10 ORF plots
    2. Metagene analysis plots
    3. HTML report with interactive visualizations
    
    Required files:
    - Scored ORFs: CSV from find-orfs or score-orfs
    - BigWig: Ribo-seq coverage
    - Exons: CSV with exon positions
    
    Example:
    translonpredictor plot -s scored_orfs.csv -bw coverage.bw -e exons.csv -o report
    """
    print("Generating visualization report")
    plottop10(scored_orfs, bigwig, exons, plot_range, outfile, getparameters(locals()))
    print("Report generation complete!")

if __name__ == '__main__':
    cli() 
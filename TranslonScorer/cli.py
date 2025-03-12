"""Command-line interface for TranslonScorer."""

import os
import click
import polars as pl
import warnings
import logging
from typing import List, Optional
import pysam
import oxbow as ox
from pathlib import Path
from memory_profiler import profile

from .core import scoring, coordinates
from .core import orffinder
from .file_handlers import bam, bed, bigwig
from .utils.logging import setup_logging, log_info, log_error
from .visualization import plots

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
    translonscorer all --help
    """
    pass

@profile(precision=4)
def process_bam_file(bam_path, annotation_file):
    """Process BAM file and get exons/CDS."""
    log_info('Processing BAM file...')
    bam_df = bam.readbam(bam_path)
    
    # Get exons and CDS
    cds_df, exon_df = bam.getexons_and_cds(annotation_file)
    
    # Process BAM data
    bam_type, _ = bam.detect_bam_type(bam_df, exon_df)
    if bam_type == 'genomic':
        # First convert genomic coordinates to transcript coordinates
        bam_df = bam.bamtranscript(bam_df, exon_df)
        # Then calculate positions relative to CDS
        bam_df = bam.process_transcriptomic_bam(bam_df, cds_df)
    else:
        # For transcriptomic BAM, just calculate CDS positions
        bam_df = bam.process_transcriptomic_bam(bam_df, cds_df)
    
    return bam_df, exon_df

@profile(precision=4)
@cli.command()
@click.option('--bam_path', '-b',
              help='Input BAM file from Ribo-seq data. Required if not providing bigwig')
@click.option('--chromsizes', '-c',
              help='Chromosome sizes file (required if processing BAM)')
@click.option('--sequence', '-s', required=True,
              help='Input FASTA file (genomic or transcriptomic)')
@click.option('--annotation', '-a', required=True,
              help='GTF annotation file (required)')
@click.option('--bigwig', '-bw',
              help='BigWig file containing Ribo-seq coverage. If provided, skips BAM processing')
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
@click.option('--log-file', type=str,
              help='Path to log file. If not provided, logs will only be written to console.')
@click.option('--log-level', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO', help='Set the logging level (default: INFO)')
def all(bam_path: str, chromsizes: str, sequence: str, annotation: str,
        bigwig: str, outfile: str, offsets: Optional[str] = None,
        start_codons: str = "ATG", stop_codons: str = "TAA,TAG,TGA",
        min_len: int = 0, max_len: int = 1000000,
        sru_range: int = 15, scoring_method: str = 'modern',
        plot_range: int = 30, log_file: Optional[str] = None,
        log_level: str = 'INFO'):
    """Run the TranslonScorer pipeline, automatically determining stages based on input.
    
    This command intelligently determines which pipeline stages to run based on provided inputs:
    
    1. If BAM file is provided:
       - Process BAM file to generate coverage tracks
       - Convert to BigWig
       - Find and score ORFs
       - Generate visualizations
       
    2. If BigWig file is provided:
       - Skip BAM processing
       - Find and score ORFs
       - Generate visualizations
       
    3. If both BAM and BigWig are provided:
       - Use the BigWig file directly
       - Find and score ORFs
       - Generate visualizations
    
    Required files:
    - Either BAM file or BigWig file
    - Sequence: FASTA file (genomic or transcriptomic)
    - Annotation: GTF file with transcript annotations
    - Chromosome sizes (only if processing BAM)
    
    Example with BAM:
    translonscorer all -b ribo.bam -c chrom.sizes -s genome.fa -a anno.gtf -o output
    
    Example with BigWig:
    translonscorer all -bw coverage.bw -s genome.fa -a anno.gtf -o output
    """
    # Configure logging
    setup_logging(level=getattr(logging, log_level))
    
    # Validate inputs
    if not bam_path and not bigwig:
        raise click.BadParameter("Either BAM file (-b) or BigWig file (-bw) must be provided")
    
    if bam_path and not chromsizes:
        raise click.BadParameter("Chromosome sizes file (-c) is required when processing BAM files")
    
    # Determine pipeline stages
    bigwig_path = bigwig
    if bam_path:
        if not bigwig:
            log_info("BAM file provided without BigWig. Will process BAM to generate coverage.")
            location = os.path.abspath(bam_path)
            if not os.path.isfile(location):
                raise click.BadParameter(f"BAM file not found: {bam_path}")
            
            # Process BAM and get annotations
            bam_df, exon_df = process_bam_file(location, annotation)
            
            # Calculate A-site positions
            offsets = coordinates.change_point_analysis(bam_df)
            bed_df = bed.asitecalc(bam_df, offsets)
            
            log_info(f"A-site positions calculated for {len(bed_df)} transcripts")
            
            # Convert to BigWig
            bedgraph_path = f"{outfile}.bedGraph"
            bed_df.write_csv(bedgraph_path, separator="\t", include_header=False)
            bigwig_path = f"{outfile}.bw"
            bed.bedtobigwig(bedgraph_path, chromsizes, bigwig_path)
        else:
            log_info("Both BAM and BigWig provided. Using BigWig directly.")
            bigwig_path = bigwig
    
    # Get exons and CDS for ORF finding
    log_info("Finding ORFs...")
    cds_df, exon_df = bam.getexons_and_cds(annotation)
    
    # Find ORFs
    orf_df = orffinder.preporfs(sequence, start_codons.split(","), stop_codons.split(","), min_len, max_len)
    
    # Score ORFs
    log_info("Scoring ORFs...")
    scored_orfs = bigwig.scoring(bigwig_path, exon_df, orf_df, scoring_method == 'classic', sru_range)
    scored_orfs.write_csv(f"{outfile}_orfs_scored.csv")
    
    # Generate plots
    log_info("Generating plots...")
    plots.plottop10(f"{outfile}_orfs_scored.csv", bigwig_path, exon_df, plot_range, outfile)
    
    log_info("Pipeline completed successfully!")

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
    """
    print("Processing BAM file...")
    location = os.path.abspath(bam)
    
    if not os.path.isfile(location):
        raise click.BadParameter(f"BAM file not found: {bam}")
    
    # Read BAM file
    bam_df = pl.read_csv(location)
    
    # Get exons and CDS
    cds_df, exon_df = bam.getexons_and_cds(annotation)
    
    # Process BAM data
    bam_type, _ = bam.detect_bam_type(bam_df, exon_df)
    if bam_type == 'genomic':
        # First convert genomic coordinates to transcript coordinates
        bam_df = bam.bamtranscript(bam_df, exon_df)
        # Then calculate positions relative to CDS
        bam_df = bam.process_transcriptomic_bam(bam_df, cds_df)
    else:
        # For transcriptomic BAM, just calculate CDS positions
        bam_df = bam.process_transcriptomic_bam(bam_df, cds_df)
    
    # Calculate A-site positions
    offsets = coordinates.change_point_analysis(bam_df)
    bed_df = bed.asitecalc(bam_df, offsets)
    
    # Convert to BigWig
    bedgraph_path = f"{outfile}.bedGraph"
    bed_df.write_csv(bedgraph_path, separator="\t", include_header=False)
    bed.bedtobigwig(bedgraph_path, chromsizes, outfile)
    
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
    """
    print("Finding ORFs...")
    orf_df = orffinder.preporfs(sequence, start_codons.split(","), stop_codons.split(","), min_len, max_len)
    
    # Get exons and CDS
    cds_df, exon_df = bam.getexons_and_cds(annotation)
    
    # Score ORFs
    print("Scoring ORFs...")
    scored_orfs = bigwig.scoring(bigwig, exon_df, orf_df, scoring_method == 'classic', sru_range)
    scored_orfs.write_csv(f"{outfile}_orfs_scored.csv")
    
    # Generate plots
    print("Generating plots...")
    plots.plottop10(f"{outfile}_orfs_scored.csv", bigwig, exon_df, 30, outfile)
    
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
    """
    print("Scoring ORFs...")
    scored_orfs = bigwig.scoring(bigwig, exons, orfs, scoring_method == 'classic', sru_range)
    scored_orfs.write_csv(f"{outfile}_orfs_scored.csv")
    
    # Generate plots
    print("Generating plots...")
    plots.plottop10(f"{outfile}_orfs_scored.csv", bigwig, exons, 30, outfile)
    
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
    """
    print("Generating visualization report...")
    plots.plottop10(scored_orfs, bigwig, exons, plot_range, outfile)
    print("Report generation complete!")

@click.command()
@click.option('--bam-file', required=True, help='Path to BAM file')
@click.option('--gtf-file', required=True, help='Path to GTF annotation file')
@click.option('--chrom-sizes', required=True, help='Path to chromosome sizes file')
@click.option('--output-prefix', required=True, help='Prefix for output files')
@click.option('--min-quality', default=50, help='Minimum mapping quality', type=int)
@click.option('--genomic/--transcriptomic', default=True, help='Whether BAM is genomic or transcriptomic')
def main(bam_file, gtf_file, chrom_sizes, output_prefix, min_quality, genomic):
    """
    Main entry point for TranslonScorer analysis.
    """
    try:
        # Create output directory if needed
        output_dir = Path(output_prefix).parent
        output_dir.mkdir(parents=True, exist_ok=True)

        # Process BAM file
        log_info('Processing BAM file...')
        bam_df = bam.readbam(bam_file)
        
        # Get exons and CDS
        cds_df, exon_df = bam.getexons_and_cds(gtf_file)

        # Process BAM data
        if genomic:
            bam_df = bam.bamtranscript(bam_df, exon_df)
        else:
            bam_df = bam.process_transcriptomic_bam(bam_df, cds_df)

        # Calculate A-site positions
        offsets = coordinates.change_point_analysis(bam_df)
        bed_df = bed.asitecalc(bam_df, offsets)

        # Convert to BigWig
        bedgraph = f"{output_prefix}.bedGraph"
        bigwig_out = f"{output_prefix}.bw"
        bed.bedtobigwig(bedgraph, chrom_sizes, bigwig_out)

        # Score ORFs
        orfs_df = bigwig.scoring(bigwig_out, exon_df, cds_df, False, min_quality)

        # Save results
        bed.saveorfsandexons(orfs_df, exon_df, output_prefix)
        
        log_info('Analysis complete!')
        
    except Exception as e:
        log_error(f"Error during analysis: {str(e)}")
        raise

if __name__ == '__main__':
    cli() 
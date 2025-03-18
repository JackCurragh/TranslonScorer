"""Command-line interface for TranslonScorer."""

import os
import click
import polars as pl
import warnings
import logging
from typing import List, Optional
from pathlib import Path
from memory_profiler import profile

from .core import scoring, coordinates
from .core import orffinder
from .file_handlers import bam, bed, bigwig
from .core.scoring import orfrelativeposition
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
# Updates to cli.py

@cli.command()
@click.option('--bam_path', '-b',
              help='Input BAM file from Ribo-seq data. Required if not providing bigwig')
@click.option('--chromsizes', '-c',
              help='Chromosome sizes file (required if processing BAM)')
@click.option('--sequence', '-s', required=True,
              help='Input FASTA file (genomic or transcriptomic)')
@click.option('--annotation', '-a', required=True,
              help='GTF annotation file (required)')
@click.option('--bigwig_path', '-bw',
              help='BigWig file containing Ribo-seq coverage. If provided, skips BAM processing')
@click.option('--forward_bigwig', '-fwd',
              help='Forward strand BigWig file. Used with --reverse_bigwig for strand-specific analysis')
@click.option('--reverse_bigwig', '-rev',
              help='Reverse strand BigWig file. Used with --forward_bigwig for strand-specific analysis')
@click.option('--stranded/--unstranded', default=False,
              help='Process strands separately (default: False)')
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
def all(bam_path, chromsizes, sequence, annotation, bigwig_path, forward_bigwig, reverse_bigwig,
        stranded, outfile, offsets, start_codons="ATG", stop_codons="TAA,TAG,TGA",
        min_len=0, max_len=1000000, sru_range=15, scoring_method='modern',
        plot_range=30, log_file=None, log_level='INFO'):
    """Run the TranslonScorer pipeline, automatically determining stages based on input.
    
    This command intelligently determines which pipeline stages to run based on provided inputs.
    Now supports strand-specific analysis with --stranded flag or by providing separate
    --forward_bigwig and --reverse_bigwig files.
    """
    # Configure logging
    setup_logging(level=getattr(logging, log_level))
    
    # Validate inputs
    if not bam_path and not bigwig_path and not (forward_bigwig and reverse_bigwig):
        raise click.BadParameter(
            "Either BAM file (-b) or BigWig file (-bw) or forward/reverse BigWig files must be provided"
        )
    
    if bam_path and not chromsizes:
        raise click.BadParameter("Chromosome sizes file (-c) is required when processing BAM files")
    
    # Handle strand-specific bigwig inputs
    if forward_bigwig and reverse_bigwig:
        bigwig_paths = {'forward': forward_bigwig, 'reverse': reverse_bigwig}
        stranded = True
    elif bigwig_path:
        bigwig_paths = bigwig_path
    else:
        bigwig_paths = None
    
    # Initialize exon_df
    exon_df = None

    # Determine pipeline stages
    if bam_path:
        if not bigwig_paths:
            log_info("BAM file provided without BigWig. Will process BAM to generate coverage.")
            location = os.path.abspath(bam_path)
            if not os.path.isfile(location):
                raise click.BadParameter(f"BAM file not found: {bam_path}")
            
            # Process BAM and get annotations
            bam_df, exon_df = process_bam_file(location, annotation)
            
            # Calculate A-site positions
            offsets = coordinates.change_point_analysis(bam_df)
            
            # Process by strand if requested
            if stranded:
                # Split BAM by strand
                fwd_bam = bam_df.filter(pl.col("strand") == "+")
                rev_bam = bam_df.filter(pl.col("strand") == "-")
                
                # Process forward strand
                fwd_bed = bed.asitecalc(fwd_bam, offsets)
                fwd_bedgraph = f"{outfile}_fwd.bedGraph"
                fwd_bed.write_csv(fwd_bedgraph, separator="\t", include_header=False)
                fwd_bigwig = f"{outfile}_fwd.bw"
                bed.bedtobigwig(fwd_bedgraph, chromsizes, fwd_bigwig)
                
                # Process reverse strand
                rev_bed = bed.asitecalc(rev_bam, offsets)
                rev_bedgraph = f"{outfile}_rev.bedGraph"
                rev_bed.write_csv(rev_bedgraph, separator="\t", include_header=False)
                rev_bigwig = f"{outfile}_rev.bw"
                bed.bedtobigwig(rev_bedgraph, chromsizes, rev_bigwig)
                
                # Set up paths for scoring
                bigwig_paths = {'forward': fwd_bigwig, 'reverse': rev_bigwig}
            else:
                # Process combined strands (original behavior)
                bed_df = bed.asitecalc(bam_df, offsets)
                bedgraph_path = f"{outfile}.bedGraph"
                bed_df.write_csv(bedgraph_path, separator="\t", include_header=False)
                bigwig_path = f"{outfile}.bw"
                bed.bedtobigwig(bedgraph_path, chromsizes, bigwig_path)
                bigwig_paths = bigwig_path
        else:
            log_info("Both BAM and BigWig provided. Using BigWig directly.")
    
    # Ensure exon_df is set when using BigWig directly
    if exon_df is None:
        log_info("Loading exon data from annotation file...")
        cds_df, exon_df = bam.getexons_and_cds(annotation)
    
    # Process BigWig files to ensure they're in transcript coordinates
    log_info("Processing BigWig files...")
    bigwig_paths = bigwig.get_bigwig_paths_for_strands(
        bigwig_paths, exon_df, annotation, outfile, stranded
    )
    
    # Ensure transcriptomic input for ORF finding
    log_info("Ensuring transcriptomic input for ORF finding...")
    transcript_fasta = sequence
    if not sequence.endswith('_transcripts.fa'):
        log_info("Generating transcript sequences from genomic FASTA and GTF annotation...")
        transcript_fasta = coordinates.gettranscripts(sequence, annotation, outfile)
    
    # Find ORFs
    log_info("Finding ORFs...")
    orf_df = orffinder.preporfs(transcript_fasta, start_codons.split(","), stop_codons.split(","), min_len, max_len)
    
    log_info("Determining relative position of ORFs to CDS...")
    # Determine the relative position of ORFs to CDS
    orf_df, exon_coords = orfrelativeposition(annotation, orf_df, cds_df)

    # Score ORFs with strand-aware scoring
    log_info("Scoring ORFs...")
    scored_orfs = score_orfs_by_strand(
        bigwig_paths, exon_df, orf_df, 
        scoring_method == 'classic', 
        sru_range
    )
    
    scored_orfs.write_csv(f"{outfile}_orfs_scored.csv")
    
    # Generate plots
    log_info("Generating plots...")
    plots.plottop10(f"{outfile}_orfs_scored.csv", bigwig_paths['forward'], exon_df, plot_range, outfile)
    
    log_info("Pipeline completed successfully!")


def score_orfs_by_strand(bigwig_paths, exon_df, orf_df, old_scoring, sru_range):
    """
    Score ORFs based on strand-specific BigWig data.
    
    Parameters:
        bigwig_paths (dict): Dictionary with 'forward' and 'reverse' paths for BigWig files
        exon_df (DataFrame): DataFrame containing exon information
        orf_df (DataFrame): DataFrame containing ORF information
        old_scoring (bool): Whether to use old scoring method
        sru_range (int): Nucleotide range for Start Rise Up score calculation
        
    Returns:
        DataFrame: Scored ORFs
    """
    from ..utils.logging import log_info
    import polars as pl
    
    # Split ORFs by strand
    pos_strand_orfs = orf_df.filter(pl.col("strand") == "+")
    neg_strand_orfs = orf_df.filter(pl.col("strand") == "-")
    
    # Handle ORFs without strand information
    unstrand_orfs = orf_df.filter(~pl.col("strand").is_in(["+", "-"]))
    
    # Score each strand separately
    results = []
    
    if not pos_strand_orfs.is_empty():
        log_info("Scoring forward strand ORFs...")
        forward_results = bigwig.scoring(
            bigwig_paths['forward'], exon_df, pos_strand_orfs, old_scoring, sru_range
        )
        results.append(forward_results)
    
    if not neg_strand_orfs.is_empty():
        log_info("Scoring reverse strand ORFs...")
        reverse_results = bigwig.scoring(
            bigwig_paths['reverse'], exon_df, neg_strand_orfs, old_scoring, sru_range
        )
        results.append(reverse_results)
    
    if not unstrand_orfs.is_empty():
        log_info("Scoring unstranded ORFs using forward data...")
        unstrand_results = bigwig.scoring(
            bigwig_paths['forward'], exon_df, unstrand_orfs, old_scoring, sru_range
        )
        results.append(unstrand_results)
    
    # Combine results
    if results:
        return pl.concat(results)
    else:
        return pl.DataFrame()

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
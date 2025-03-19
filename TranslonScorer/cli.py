"""Command-line interface for TranslonScorer."""

import click
from typing import Optional
from .pipeline.workflow import process_bam_workflow, find_orfs_workflow, score_orfs_workflow, plot_workflow, all_workflow
from .pipeline.config import Config
from .pipeline.validator import validate_config
from .utils.logging import setup_logging, log_info, log_error


def common_options(func):
    """Apply common options to command functions."""
    # Input files
    func = click.option('--bam', '-b',
              help='Input BAM file from Ribo-seq data. Required if not providing bigwig files')(func)
    func = click.option('--chromsizes', '-c',
              help='Chromosome sizes file (required if processing BAM)')(func)
    func = click.option('--sequence', '-s', required=True,
              help='Input FASTA file (genomic or transcriptomic)')(func)
    func = click.option('--annotation', '-a', required=True,
              help='GTF annotation file for identifying exons and transcripts')(func)
    
    # BigWig options
    func = click.option('--bigwig', '-w',
              help='BigWig file containing Ribo-seq coverage. If provided, skips BAM processing')(func)
    func = click.option('--forward_bigwig', '-fw',
              help='Forward strand BigWig file for strand-specific analysis')(func)
    func = click.option('--reverse_bigwig', '-rv',
              help='Reverse strand BigWig file for strand-specific analysis')(func)
    
    # Analysis options
    func = click.option('--stranded/--unstranded', default=False,
              help='Process strands separately (default: False)')(func)
    func = click.option('--offsets', 
              help='File containing read length-specific offsets for A-site calculation')(func)
    func = click.option('--start_codons', '--starts', default="ATG",
              help='Comma-separated list of start codons (default: ATG)')(func)
    func = click.option('--stop_codons', '--stops', default="TAA,TAG,TGA",
              help='Comma-separated list of stop codons (default: TAA,TAG,TGA)')(func)
    func = click.option('--min_length', '--min', type=int, default=0,
              help='Minimum ORF length in nucleotides (default: 0)')(func)
    func = click.option('--max_length', '--max', type=int, default=1000000,
              help='Maximum ORF length in nucleotides (default: 1000000)')(func)
    func = click.option('--sru_range', '--sru', type=int, default=15,
              help='Nucleotide range for Start Rise Up score calculation (default: 15)')(func)
    func = click.option('--plot_range', '--plot', type=int, default=30,
              help='Plot range around start position (default: 30)')(func)
    
    # Output options
    func = click.option('--output', '-o', required=True,
              help='Base name for output files')(func)
    func = click.option('--log_file', '--log',
              help='Path to log file. If not provided, logs will only be written to console.')(func)
    func = click.option('--log_level', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO', help='Set the logging level (default: INFO)')(func)
    
    return func

@click.group(invoke_without_command=True)
@click.pass_context
@common_options
def cli(ctx, **kwargs):
    """TranslonScorer: A tool for identifying and scoring translational events from Ribo-seq data.
    
    This tool provides several workflows:
    
    1. all: Run the complete pipeline end-to-end
    2. process-bam: Process Ribo-seq BAM files to generate coverage tracks
    3. find-orfs: Identify and score potential ORFs from sequence data
    4. score-orfs: Score existing ORFs using coverage data
    5. plot: Generate visualization reports from scored ORFs
    
    If no command is specified, the complete pipeline will be run.
    For detailed instructions, use --help with any command:
    translonscorer --help
    translonscorer process-bam --help
    """
    # If no subcommand is provided, run the 'all' logic
    if ctx.invoked_subcommand is None:
        # Call the all function with the provided parameters
        all_workflow(**kwargs)


@cli.command()
@common_options
def all(**kwargs):
    """Run the TranslonScorer pipeline, automatically determining stages based on input.
    
    This command intelligently determines which pipeline stages to run based on provided inputs.
    """
    # Configure logging
    setup_logging(level=getattr(logging, log_level))

    # Validate configuration
    config = Config(**kwargs)
    validate_config(config)

    # Run the complete pipeline
    all_workflow(config)    

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
@click.option('--output', '-o', required=True,
              help='Base name for output files')
def process_bam(bam: str, chromsizes: str, annotation: str, 
                output: str, offsets: Optional[str] = None):
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
    bedgraph_path = f"{output}.bedGraph"
    bed_df.write_csv(bedgraph_path, separator="\t", include_header=False)
    bed.bedtobigwig(bedgraph_path, chromsizes, output)
    
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
@click.option('--output', '-o', required=True,
              help='Base name for output files')
@click.option('--sru-range', type=int, default=15,
              help='Nucleotide range for Start Rise Up score calculation (default: 15)')
@click.option('--scoring-method', type=click.Choice(['classic', 'modern']), 
              default='modern', help='Scoring algorithm to use (default: modern)')
def find_orfs(sequence: str, annotation: str, bigwig: str, output: str,
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
    scored_orfs.write_csv(f"{output}_orfs_scored.csv")
    
    # Generate plots
    print("Generating plots...")
    plots.plottop10(f"{output}_orfs_scored.csv", bigwig, exon_df, 30, output)
    
    print("ORF finding and scoring complete!")

@cli.command()
@click.option('--orfs', '-f', required=True,
              help='CSV file containing pre-annotated ORFs')
@click.option('--bigwig', '-bw', required=True,
              help='BigWig file containing Ribo-seq coverage')
@click.option('--exons', '-e', required=True,
              help='CSV file containing exon positions')
@click.option('--output', '-o', required=True,
              help='Base name for output files')
@click.option('--scoring-method', type=click.Choice(['classic', 'modern']), 
              default='modern', help='Scoring algorithm to use (default: modern)')
@click.option('--sru-range', type=int, default=15,
              help='Nucleotide range for Start Rise Up score calculation (default: 15)')
def score_orfs(orfs: str, bigwig: str, exons: str, output: str,
               scoring_method: str, sru_range: int):
    """Score existing ORFs using Ribo-seq coverage data.
    
    This command takes pre-annotated ORFs and scores them using:
    1. Ribo-seq coverage from bigWig
    2. Exon position information
    3. Specified scoring parameters
    """
    print("Scoring ORFs...")
    scored_orfs = bigwig.scoring(bigwig, exons, orfs, scoring_method == 'classic', sru_range)
    scored_orfs.write_csv(f"{output}_orfs_scored.csv")
    
    # Generate plots
    print("Generating plots...")
    plots.plottop10(f"{output}_orfs_scored.csv", bigwig, exons, 30, output)
    
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
@click.option('--output', '-o', required=True,
              help='Base name for output files')
def plot(scored_orfs: str, bigwig: str, exons: str, 
         plot_range: int, output: str):
    """Generate visualization reports from scored ORFs.
    
    This command creates visualization reports including:
    1. Top 10 ORF plots
    2. Metagene analysis plots
    3. HTML report with interactive visualizations
    """
    print("Generating visualization report...")
    plots.plottop10(scored_orfs, bigwig, exons, plot_range, output)
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
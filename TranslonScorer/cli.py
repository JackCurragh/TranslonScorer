"""Command-line interface for TranslonScorer."""

import click
import polars as pl
from .pipeline.workflow import (
    process_bam_workflow,
    find_orfs_workflow,
    score_orfs_workflow,
    plot_workflow,
    all_workflow,
)
from .pipeline.config import Config
from .pipeline.validator import validate_config
from .utils.logging import setup_logging, log_info


def common_options(func):
    """Apply common options to command functions."""
    # Input files
    func = click.option('--bam', '-b', help='Input BAM file from Ribo-seq data.')(func)
    func = click.option('--bam-collapsed', is_flag=True, default=False, help='Treat BAM as collapsed (counts in name or tag).')(func)
    func = click.option('--bam-count-from', type=click.Choice(['name','tag']), help='Where to parse counts for collapsed BAM.')(func)
    func = click.option('--bam-count-pattern', help='Regex with named group "count" for name-based collapsed BAM.')(func)
    func = click.option('--bam-count-tag', help='SAM tag for tag-based collapsed BAM (e.g. RC).')(func)
    func = click.option('--chromsizes', '-c', help='Chromosome sizes file (required if processing BAM)')(func)
    func = click.option('--sequence', '-s', required=True, help='Input FASTA file (genomic or transcriptomic)')(func)
    func = click.option('--annotation', '-a', required=True, help='GTF annotation file for identifying exons and transcripts')(func)
    
    # BigWig options
    func = click.option('--bigwig', '-w', help='BigWig file containing Ribo-seq coverage. If provided, skips BAM/Zarr processing')(func)
    func = click.option('--forward-bigwig', 'forward_bigwig', help='Forward strand BigWig file for strand-specific analysis')(func)
    func = click.option('--reverse-bigwig', 'reverse_bigwig', help='Reverse strand BigWig file for strand-specific analysis')(func)
    
    # Analysis options
    func = click.option('--stranded/--unstranded', default=False, help='Process strands separately (default: False)')(func)
    func = click.option('--offsets', help='File containing read length-specific offsets for A-site calculation')(func)
    func = click.option('--start-codons', 'start_codons', default='ATG', help='Comma-separated list of start codons (default: ATG)')(func)
    func = click.option('--stop-codons', 'stop_codons', default='TAA,TAG,TGA', help='Comma-separated list of stop codons (default: TAA,TAG,TGA)')(func)
    func = click.option('--min-length', 'min_length', type=int, default=0, help='Minimum ORF length in nucleotides (default: 0)')(func)
    func = click.option('--max-length', 'max_length', type=int, default=1000000, help='Maximum ORF length in nucleotides (default: 1000000)')(func)
    func = click.option('--sru-range', 'sru_range', type=int, default=15, help='Nucleotide range for Start Rise Up score calculation (default: 15)')(func)
    func = click.option('--plot-range', 'plot_range', type=int, default=30, help='Plot range around start position (default: 30)')(func)
    
    # Output options
    func = click.option('--output', '-o', required=True, help='Base name for output files')(func)
    func = click.option('--log-file', 'log_file', help='Path to log file. If not provided, logs will only be written to console.')(func)
    func = click.option('--log-level', 'log_level', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']), default='INFO', help='Set the logging level (default: INFO)')(func)
    
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
        ctx.invoke(all, **kwargs)


@cli.command()
@common_options
def all(**kwargs):
    """Run the complete pipeline end-to-end based on provided inputs."""
    setup_logging()
    config = Config.from_click_args(**kwargs)
    validate_config(config)
    all_workflow(config)
    log_info("Pipeline completed successfully!")

    
@cli.command("process-bam")
@click.option('--bam', '-b', required=True, help='Input BAM file (genomic or transcriptomic).')
@click.option('--chromsizes', '-c', required=True, help='Chromosome sizes file (required for bigWig conversion).')
@click.option('--annotation', '-a', required=True, help='GTF annotation file.')
@click.option('--output', '-o', required=True, help='Base name for output files.')
@click.option('--stranded/--unstranded', default=False, help='Write strand-specific bigWigs if enabled.')
def process_bam(bam: str, chromsizes: str, annotation: str, output: str, stranded: bool):
    """Process BAM to coverage (bedGraph/bigWig)."""
    setup_logging()
    config = Config(bam=bam, chromsizes=chromsizes, sequence='placeholder.fa', annotation=annotation, output=output, stranded=stranded)
    # Minimal validation for this subcommand
    config._validate_file_exists(config.bam, 'BAM')
    config._validate_file_exists(config.chromsizes, 'Chromosome sizes')
    config._validate_file_exists(config.annotation, 'Annotation')
    process_bam_workflow(config)
    log_info("BAM processing complete!")

@cli.command("find-orfs")
@common_options
def find_orfs(**kwargs):
    """Find ORFs and classify relative to CDS (no scoring)."""
    setup_logging()
    config = Config.from_click_args(**kwargs)
    # minimal validation of inputs for this step
    config._validate_file_exists(config.sequence, 'Sequence')
    config._validate_file_exists(config.annotation, 'Annotation')
    orf_df, _ = find_orfs_workflow(config)
    orf_df.write_csv(f"{config.output}_orfs.csv")
    log_info("ORF finding complete!")


@cli.command("profiles")
@common_options
@click.option('--profiles-out', required=True, help='Output path (Parquet for transcript profiles or Zarr path for locus mode).')
@click.option('--loci-bed', help='BED of loci to build genomic A-site profile matrices (Zarr).')
@click.option('--offsets-file', help='CSV with columns length,offset to override defaults.')
def profiles(**kwargs):
    """Generate transcript-space A-site profiles from BAM (classic/collapsed), Zarr, or BigWig."""
    setup_logging()
    config = Config.from_click_args(**kwargs)
    # Minimal file checks
    config._validate_file_exists(config.annotation, 'Annotation')
    if not (config.bam or config.zarr_root or config.bigwig or (config.forward_bigwig and config.reverse_bigwig)):
        raise click.BadParameter('Provide one of: --bam (classic/collapsed), --zarr-root with --read-index-parquet and --sample, or --bigwig/--forward-bigwig+--reverse-bigwig')

    # Load annotation
    cds_df, exon_df = __import__('TranslonScorer.file_handlers.bam', fromlist=['']).getexons_and_cds(config.annotation)

    from .pipeline.profiles import profiles_from_bam, profiles_from_zarr, profiles_from_bigwig, write_profiles_parquet
    from .pipeline.locus_profiles import build_locus_profiles_zarr

    if config.bam:
        prof = profiles_from_bam(
            config.bam,
            exon_df,
            cds_df,
            collapsed=config.bam_collapsed,
            count_from=config.bam_count_from,
            count_pattern=config.bam_count_pattern,
            count_tag=config.bam_count_tag,
        )
        write_profiles_parquet(prof, kwargs['profiles_out'])
        log_info("Profiles written")
    elif config.zarr_root:
        if not (config.read_index_parquet and config.samples):
            raise click.BadParameter('Zarr mode requires --read-index-parquet and at least one --sample')
        if kwargs.get('loci_bed'):
            # Locus matrices path
            build_locus_profiles_zarr(
                config.zarr_root,
                config.read_index_parquet,
                list(config.samples),
                kwargs['loci_bed'],
                offsets_file=kwargs.get('offsets_file'),
                out_zarr=kwargs['profiles_out'],
            )
            log_info("Locus profile matrices written to Zarr")
        else:
            # Transcript-space tidy profiles
            for sample, prof in profiles_from_zarr(
                config.zarr_root,
                config.read_index_parquet,
                list(config.samples),
                exon_df,
                cds_df,
            ):
                out = kwargs['profiles_out']
                stem, ext = (out.rsplit('.', 1) + ['parquet'])[:2]
                write_profiles_parquet(prof, f"{stem}_{sample}.{ext}", sample=sample)
            log_info("Profiles written for all samples")
    else:
        # BigWig lane: single or stranded inputs
        bw = { 'forward': config.forward_bigwig, 'reverse': config.reverse_bigwig } if (config.forward_bigwig and config.reverse_bigwig) else config.bigwig
        prof = profiles_from_bigwig(bw, exon_df, stranded=config.stranded)
        write_profiles_parquet(prof, kwargs['profiles_out'])
        log_info("Profiles written from BigWig")

@cli.command("score-orfs")
@click.option('--orfs', '-f', required=True, help='CSV file containing ORFs to score')
@click.option('--exons', '-e', required=True, help='CSV file containing exon positions')
@click.option('--bigwig', '-w', required=True, help='BigWig file containing Ribo-seq coverage')
@click.option('--output', '-o', required=True, help='Base name for output files')
@click.option('--scoring-method', type=click.Choice(['classic', 'modern']), default='modern')
@click.option('--sru-range', type=int, default=15)
def score_orfs(orfs: str, exons: str, bigwig: str, output: str, scoring_method: str, sru_range: int):
    """Score ORFs and write results + report."""
    setup_logging()
    config = Config(sequence='placeholder.fa', annotation='placeholder.gtf', bigwig=bigwig, output=output, scoring_method=scoring_method, sru_range=sru_range)
    # load inputs
    orf_df = pl.read_csv(orfs)
    exon_df = pl.read_csv(exons)
    scored = score_orfs_workflow(config, bigwig, exon_df, orf_df)
    scored_path = f"{output}_orfs_scored.csv"
    scored.write_csv(scored_path)
    plots.plottop10(scored_path, bigwig, exons, 30, output)
    log_info("ORF scoring complete!")

@cli.command("plot")
@click.option('--scored-orfs', '-s', required=True, help='CSV file of scored ORFs')
@click.option('--bigwig', '-w', required=True, help='BigWig file containing Ribo-seq coverage')
@click.option('--exons', '-e', required=True, help='CSV file containing exon positions')
@click.option('--plot-range', default=30, type=int, help='Plot range around start position (default: 30)')
@click.option('--output', '-o', required=True, help='Base name for output files')
def plot(scored_orfs: str, bigwig: str, exons: str, plot_range: int, output: str):
    """Generate visualization reports from scored ORFs."""
    setup_logging()
    plots.plottop10(scored_orfs, bigwig, exons, plot_range, output)
    log_info("Report generation complete!")

@cli.command("features")
@click.option('--annotation', '-a', required=True, help='GTF annotation file.')
@click.option('--output-prefix', '-o', required=True, help='Output prefix for feature tables (Parquet).')
def features(annotation: str, output_prefix: str):
    """Emit locus features and transcript→feature mappings (no signals)."""
    setup_logging()
    from .pipeline.locus_features import build_locus_features
    feats, fmap = build_locus_features(annotation)
    feats.write_parquet(f"{output_prefix}.features.parquet")
    fmap.write_parquet(f"{output_prefix}.feature_map.parquet")
    log_info("Feature tables written")


@cli.command("feature-metrics")
@click.option('--profiles', required=True, help='Transcript profiles Parquet (canonical).')
@click.option('--features', required=True, help='Locus feature table Parquet.')
@click.option('--feature-map', 'feature_map', required=True, help='Transcript→feature map Parquet.')
@click.option('--splits-csv', help='Optional CSV of split junction counts with columns chr,donor_pos,acceptor_pos,strand,count')
@click.option('--out', 'out_parquet', required=True, help='Output feature metrics Parquet path.')
def feature_metrics_cmd(profiles: str, features: str, feature_map: str, splits_csv: str | None, out_parquet: str):
    """Compute per-feature metrics including junction LLR scores and SRU for TIS/TTS."""
    setup_logging()
    from .pipeline.feature_metrics import feature_metrics
    splits_df = None
    if splits_csv:
        import polars as pl
        splits_df = pl.read_csv(splits_csv)
    feature_metrics(
        profiles_parquet=profiles,
        feature_parquet=features,
        feature_map_parquet=feature_map,
        genome_bam_splits=splits_df,
        out_parquet=out_parquet,
    )
    log_info("Feature metrics computed")


@cli.command("assemble")
@click.option('--orfs-parquet', required=True, help='ORF candidates with composite scores (Parquet).')
@click.option('--out', 'out_parquet', required=True, help='Output assembled translome Parquet.')
@click.option('--solver', type=click.Choice(['PULP', 'GREEDY']), default='PULP')
@click.option('--timeout', 'timeout_sec', type=int, default=60)
def assemble_cmd(orfs_parquet: str, out_parquet: str, solver: str, timeout_sec: int):
    """Assemble a locus translome by selecting a consistent set of ORFs under soft penalties."""
    setup_logging()
    from .pipeline.assemble import assemble_translome
    assemble_translome(orfs_parquet, out_parquet, solver=solver, timeout_sec=timeout_sec)
    log_info("Translome assembly complete")


@cli.command("orf-composite")
@click.option('--orfs', required=True, help='ORF candidates (Parquet/CSV) with tran_id,start,stop,frame,type,length[,locus_id].')
@click.option('--feature-metrics', required=True, help='Feature metrics Parquet.')
@click.option('--feature-map', 'feature_map', required=True, help='Transcript→feature map Parquet.')
@click.option('--out', 'out_parquet', required=True, help='Output ORFs with composite scores (Parquet).')
def orf_composite_cmd(orfs: str, feature_metrics: str, feature_map: str, out_parquet: str):
    """Aggregate per-feature metrics into composite ORF scores."""
    setup_logging()
    from .pipeline.orf_composite import orf_composite
    orf_composite(orfs, feature_metrics, feature_map, out_parquet)
    log_info("Composite ORF scores written")

if __name__ == '__main__':
    cli()
    # Zarr unique read matrix
    func = click.option('--zarr-root', help='Root path to Zarr unique-read matrix (samples x read_id).')(func)
    func = click.option('--read-index-parquet', help='Parquet index mapping read_id to genomic coords.')(func)
    func = click.option('--sample', 'samples', multiple=True, help='Sample(s) to process from Zarr. Repeatable.')(func)

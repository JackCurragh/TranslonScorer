"""Command-line interface for TranslonScorer."""
from __future__ import annotations

import os
import click
from typing import Optional
import polars as pl
from .pipeline.config import Config
from .pipeline.validator import validate_config
from .utils.logging import setup_logging, log_info
from .file_handlers import bam as bam_handlers


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
    func = click.option('--annotation', '-a', required=False, help='GTF annotation file for identifying exons and transcripts')(func)
    func = click.option('--annotation-dir', required=False, help='Annotation bundle directory (preferred)')(func)
    
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
    # Not all commands require an explicit -o (e.g., profiles when --profiles-out is given).
    # Keep optional here; commands that require -o should validate explicitly.
    func = click.option('--output', '-o', required=False, help='Base name for output files')(func)
    func = click.option('--log-file', 'log_file', help='Path to log file. If not provided, logs will only be written to console.')(func)
    func = click.option('--log-level', 'log_level', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']), default='INFO', help='Set the logging level (default: INFO)')(func)
    
    return func

@click.group(invoke_without_command=True)
@click.pass_context
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
    from .pipeline.workflow import all_workflow
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
    from .pipeline.workflow import process_bam_workflow
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
    from .pipeline.workflow import find_orfs_workflow
    orf_df, _ = find_orfs_workflow(config)
    orf_df.write_csv(f"{config.output}_orfs.csv")
    log_info("ORF finding complete!")


@cli.command("profiles")
@common_options
@click.option(
    '--profiles-out',
    required=False,
    help='Optional output path. Defaults under -o: transcript_profiles.parquet (transcript) or locus_profiles.zarr (locus).'
)
@click.option('--loci-bed', help='BED of loci to build genomic A-site profile matrices (Zarr).')
@click.option('--offsets-file', help='CSV with columns length,offset to override defaults.')
# Zarr + indexing helpers
@click.option('--zarr-root', help='Root of Zarr counts matrix (samples x reads or reads x samples).')
@click.option('--read-index-parquet', help='Parquet mapping read_id (Zarr row) to genomic alignment (chr,start,stop,strand,length).')
@click.option('--splits-index-parquet', help='Optional Parquet of per-read junctions (read_id, chr, donor_pos, acceptor_pos, strand).')
@click.option('--sample', 'samples', multiple=True, help='Sample name(s) to process from the Zarr matrix.')
@click.option('--all-samples', is_flag=True, default=False, help='Process all samples found in the Zarr counts array (uses /samples names if present).')
@click.option('--zarr-metadata-parquet', help='Parquet with read row metadata (e.g., row_id,qname or sequence).')
@click.option('--zarr-reads-fasta', help='FASTA with reads in Zarr row order; headers or sequences used for mapping.')
@click.option('--zarr-read-key', type=click.Choice(['auto','row_id','qname','sequence']), default='auto', help='Key to align Zarr rows to BAM.')
@click.option('--bam-key', type=click.Choice(['auto','qname','sequence']), default='auto', help='Key to align BAM reads to Zarr.')
@click.option('--hash-alg', type=click.Choice(['sha1','md5','xxh64']), default='sha1', help='Hash algorithm for sequence-based mapping.')
# Outputs & policy
@click.option('--junctions-out', help='Path to write junction counts when available.')
@click.option('--offsets-out', help='Path to write discovered offsets per sample/length.')
@click.option('--partitioned/--no-partitioned', default=False, help='Write profiles as a partitioned dataset by sample.')
@click.option('--offsets-mode', type=click.Choice(['auto','global','required']), default='auto', help='Offsets selection policy.')
def profiles(**kwargs):
    """Generate transcript-space A-site profiles from BAM (classic/collapsed), Zarr, or BigWig."""
    setup_logging()
    config = Config.from_click_args(**kwargs)
    # Decide default output path when --profiles-out is omitted (do not require -o)
    if not kwargs.get('profiles_out'):
        out_base = kwargs.get('output')
        if out_base:
            # If -o provided, derive from it
            if out_base.endswith('/') or os.path.isdir(out_base):
                out_dir = out_base
                kwargs['profiles_out'] = os.path.join(
                    out_dir,
                    'locus_profiles.zarr' if kwargs.get('loci_bed') else 'transcript_profiles.parquet'
                )
            else:
                base = out_base.rstrip('/')
                kwargs['profiles_out'] = (
                    f"{base}_locus_profiles.zarr" if kwargs.get('loci_bed')
                    else f"{base}_transcript_profiles.parquet"
                )
        else:
            # Neither --profiles-out nor -o given: write into cwd with sensible default name
            kwargs['profiles_out'] = (
                'locus_profiles.zarr' if kwargs.get('loci_bed') else 'transcript_profiles.parquet'
            )
    # Minimal file checks: accept --annotation-dir (bundle) or --annotation (GTF)
    if not (config.annotation_dir or config.annotation):
        raise click.BadParameter('Provide --annotation-dir (bundle) or --annotation (GTF)')
    if config.annotation_dir:
        # Bundle provided; no need to validate GTF path here
        pass
    else:
        config._validate_file_exists(config.annotation, 'Annotation')
    if not (config.bam or config.zarr_root or config.bigwig or (config.forward_bigwig and config.reverse_bigwig)):
        raise click.BadParameter('Provide one of: --bam (classic/collapsed), --zarr-root with --read-index-parquet and --sample, or --bigwig/--forward-bigwig+--reverse-bigwig')

    log_info("Starting profiles workflow…")
    # Load annotation: prefer bundle if provided
    exon_df = cds_df = None
    if config.annotation_dir:
        from .pipeline.annotation_bundle import load_annotation_bundle
        exon_df, cds_df, feats_df, fmap_df, tx_df, loci_bed, manifest = load_annotation_bundle(config.annotation_dir)
    else:
        cds_df, exon_df = bam_handlers.getexons_and_cds(config.annotation)

    from .pipeline.profiles import profiles_from_bam, profiles_from_zarr, profiles_from_bigwig, write_profiles_parquet
    from .pipeline.locus_profiles import build_locus_profiles_zarr

    log_info("Inputs detected: " + (
        "BAM" if config.bam else ("Zarr" if config.zarr_root else ("BigWig" if (config.bigwig or (config.forward_bigwig and config.reverse_bigwig)) else "Unknown"))
    ))
    # BAM lane
    if config.bam and not config.zarr_root:
        prof, offsets = profiles_from_bam(
            config.bam,
            exon_df,
            cds_df,
            collapsed=config.bam_collapsed,
            count_from=config.bam_count_from,
            count_pattern=config.bam_count_pattern,
            count_tag=config.bam_count_tag,
        )
        log_info(f"Writing profiles to: {kwargs['profiles_out']}")
        write_profiles_parquet(prof, kwargs['profiles_out'])
        if config.offsets_out:
            # Persist discovered offsets
            odf = pl.from_dicts([{'sample': kwargs.get('sample') or 'sample', 'length': int(L), 'offset': int(ofs)} for L, ofs in offsets.items()])
            odf.write_csv(config.offsets_out)
        if config.junctions_out:
            # Aggregate junctions from BAM
            from .pipeline.junctions import aggregate_bam_junctions
            log_info("Aggregating junctions from BAM…")
            j = aggregate_bam_junctions(config.bam)
            if not j.is_empty():
                from .utils.io import write_parquet_safe
                log_info(f"Writing junctions to: {config.junctions_out}")
                write_parquet_safe(j, config.junctions_out)
                log_info("Junction counts written from BAM")
        log_info("Profiles written")
    # Zarr lane (with optional index build from BAM)
    elif config.zarr_root:
        # Ensure samples provided
        if not config.samples or len(config.samples) == 0:
            if kwargs.get('all_samples'):
                # Discover sample names from Zarr store
                import importlib
                zmod = importlib.import_module('zarr')
                store = zmod.open(config.zarr_root, mode='r')
                counts = store['counts']
                # Prefer an existing '/samples' 1D array of names if present
                if 'samples' in store:
                    names_arr = store['samples'][...]
                    try:
                        names = [n.decode() if isinstance(n, (bytes, bytearray)) else str(n) for n in list(names_arr)]
                    except Exception:
                        names = [str(n) for n in list(names_arr)]
                else:
                    n = counts.shape[0] if counts.shape[0] <= counts.shape[1] else counts.shape[1]
                    names = [str(i) for i in range(n)]
                config.samples = names
            else:
                raise click.BadParameter('Zarr mode requires at least one --sample or use --all-samples')

            raise click.BadParameter('Zarr mode requires at least one --sample')
        # If index missing but BAM is present, build index now
        if not config.read_index_parquet and config.bam:
            from .pipeline.index_from_bam import build_read_index_from_bam
            log_info('Building Zarr read index from BAM…')
            out_idx = os.path.join(os.path.dirname(kwargs['profiles_out']), 'read_index.parquet') if kwargs.get('profiles_out') else 'read_index.parquet'
            config.read_index_parquet = build_read_index_from_bam(
                bam_path=config.bam,
                zarr_root=config.zarr_root,
                out_index=out_idx,
                zarr_metadata=config.zarr_metadata_parquet,
                zarr_reads_fasta=config.zarr_reads_fasta,
                zarr_read_key=config.zarr_read_key,
                bam_key=config.bam_key,
                hash_alg=config.hash_alg,
            )
        if not config.read_index_parquet:
            raise click.BadParameter('Zarr mode requires --read-index-parquet (or provide --bam to build it).')
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
                offsets_mode=config.offsets_mode,
                offsets_out=config.offsets_out,
            ):
                out = kwargs['profiles_out']
                stem, ext = (out.rsplit('.', 1) + ['parquet'])[:2]
                out_path = f"{stem}_{sample}.{ext}"
                log_info(f"Writing sample {sample} profiles to: {out_path}")
                write_profiles_parquet(prof, out_path, sample=sample)
            log_info("Profiles written for all samples")
            # Optional junction aggregation from splits index
            if config.junctions_out:
                if config.splits_index_parquet and os.path.isfile(config.splits_index_parquet):
                    # Aggregate per junction over read_ids using counts as weights
                    s = pl.read_parquet(config.splits_index_parquet)
                    # Join with counts per read_id per sample
                    # For Zarr, counts per read_id per sample are in the matrix; here we approximate unweighted counts (1 per split event)
                    j = (
                        s.group_by(['chr','donor_pos','acceptor_pos','strand'])
                        .agg(pl.len().alias('count'))
                    )
                    j.write_parquet(config.junctions_out)
                    log_info("Junction counts written from splits index")
                else:
                    log_info("No splits index provided; skipping junction aggregation for Zarr lane")
    else:
        # BigWig lane: single or stranded inputs
        bw = { 'forward': config.forward_bigwig, 'reverse': config.reverse_bigwig } if (config.forward_bigwig and config.reverse_bigwig) else config.bigwig
        log_info("Computing profiles from BigWig…")
        prof = profiles_from_bigwig(bw, exon_df, stranded=config.stranded)
        log_info(f"Writing profiles to: {kwargs['profiles_out']}")
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
    from .pipeline.workflow import score_orfs_workflow
    scored = score_orfs_workflow(config, bigwig, exon_df, orf_df)
    scored_path = f"{output}_orfs_scored.csv"
    scored.write_csv(scored_path)
    from .visualization import plots
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
    from .visualization import plots
    plots.plottop10(scored_orfs, bigwig, exons, plot_range, output)
    log_info("Report generation complete!")

@cli.command("index-from-bam")
@click.option('--bam', '-b', required=True, help='Input BAM used to align reads.')
@click.option('--zarr-root', required=True, help='Root of Zarr counts matrix (samples x reads or reads x samples).')
@click.option('--out-index', required=True, help='Output Parquet path for read_id→alignment index.')
@click.option('--zarr-metadata-parquet', help='Parquet with read row metadata (e.g., row_id,qname or sequence).')
@click.option('--zarr-reads-fasta', help='FASTA with reads in Zarr row order; headers or sequences used for mapping.')
@click.option('--zarr-read-key', type=click.Choice(['auto','row_id','qname','sequence']), default='auto')
@click.option('--bam-key', type=click.Choice(['auto','qname','sequence']), default='auto')
@click.option('--hash-alg', type=click.Choice(['sha1','md5','xxh64']), default='sha1')
def index_from_bam_cmd(bam: str, zarr_root: str, out_index: str, zarr_metadata_parquet: str | None, zarr_reads_fasta: str | None, zarr_read_key: str, bam_key: str, hash_alg: str):
    """Build read_index.parquet for Zarr counts by aligning rows to BAM reads."""
    setup_logging()
    from .pipeline.index_from_bam import build_read_index_from_bam
    path = build_read_index_from_bam(
        bam_path=bam,
        zarr_root=zarr_root,
        out_index=out_index,
        zarr_metadata=zarr_metadata_parquet,
        zarr_reads_fasta=zarr_reads_fasta,
        zarr_read_key=zarr_read_key,
        bam_key=bam_key,
        hash_alg=hash_alg,
    )
    log_info(f"Read index written: {path}")

@cli.command("orfs-import")
@click.option('--bed12', required=True, help='Input ORFs in BED12 format.')
@click.option('--annotation', '-a', required=True, help='GTF annotation file for exon/transcript models.')
@click.option('--out', 'out_parquet', required=True, help='Output canonical ORFs (Parquet).')
@click.option('--assign-policy', type=click.Choice(['best', 'all']), default='best', help="Assignment policy when multiple transcripts match (default: best).")
@click.option('--require-junction-match/--allow-junction-mismatch', default=True, help='Require all BED12 junctions to match transcript junctions (default: require).')
@click.option('--progress/--no-progress', default=True, help='Show progress bars (default: on).')
def orfs_import_cmd(bed12: str, annotation: str, out_parquet: str, assign_policy: str, require_junction_match: bool, progress: bool):
    """Import ORFs from BED12, map to transcripts, and write canonical ORFs (Parquet)."""
    setup_logging()
    from .pipeline.orfs_import import import_bed12
    import_bed12(
        bed12_path=bed12,
        gtf_path=annotation,
        out_parquet=out_parquet,
        assign_policy=assign_policy,
        require_junction_match=require_junction_match,
        progress=progress,
    )
    log_info("ORFs imported from BED12")

@cli.command("features")
@click.option('-a', '--annotation', required=True, help='GTF annotation file.')
@click.option('--out-dir', required=False, help='Output annotation bundle directory (recommended).')
@click.option('-o', '--output', '--output-prefix', 'output_prefix', required=False, help='Legacy output prefix for feature tables (Parquet).')
@click.option('--progress/--no-progress', default=True, help='Show progress bars (default: on).')
def features(annotation: str, out_dir: str = None, output_prefix: str = None, progress: bool = True):
    """Build an annotation bundle (exons, CDS, features, feature_map, transcripts, loci, manifest)."""
    setup_logging()
    # Backward-compat mode: if only -o is provided, write old files and also emit a bundle next to it
    if not out_dir and output_prefix:
        import os
        base = os.path.basename(output_prefix)
        root = os.path.dirname(output_prefix) or '.'
        out_dir = os.path.join(root, f"{base}.annotation")
    if not out_dir and not output_prefix:
        raise click.BadParameter('Provide either --out-dir for the bundle or -o for legacy outputs')

    # Always build the bundle
    from .pipeline.annotation_bundle import build_annotation_bundle
    bdir, paths = build_annotation_bundle(annotation, out_dir, progress=progress)
    log_info(f"Annotation bundle written at: {bdir}")

    # Write legacy outputs if requested
    if output_prefix:
        import polars as pl
        feats = pl.read_parquet(paths['features'])
        fmap = pl.read_parquet(paths['feature_map'])
        feats.write_parquet(f"{output_prefix}.features.parquet")
        fmap.write_parquet(f"{output_prefix}.feature_map.parquet")
        log_info("Legacy feature tables written")


@cli.command("feature-metrics")
@click.option('--profiles', required=True, help='Transcript profiles Parquet (canonical).')
@click.option('--features', required=True, help='Locus feature table Parquet.')
@click.option('--feature-map', 'feature_map', required=True, help='Transcript→feature map Parquet.')
@click.option('--splits-csv', help='Optional CSV of split junction counts with columns chr,donor_pos,acceptor_pos,strand,count')
@click.option('--out', 'out_parquet', required=True, help='Output feature metrics Parquet path.')
@click.option('--progress/--no-progress', default=True, help='Show progress bars (default: on).')
def feature_metrics_cmd(profiles: str, features: str, feature_map: str, splits_csv: Optional[str], out_parquet: str, progress: bool):
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

@cli.command("map-orfs")
@click.option('--orfs', 'orfs_parquet', required=True, help='Canonical ORFs with tran_id,start_pos_tran,stop_pos_tran (Parquet).')
@click.option('--feature-map', 'feature_map_parquet', required=True, help='Transcript→feature mapping Parquet.')
@click.option('--features', 'features_parquet', required=True, help='Feature table Parquet.')
@click.option('--out', 'out_parquet', required=True, help='Output per-ORF feature chains Parquet.')
@click.option('--progress/--no-progress', default=True, help='Show progress bars (default: on).')
@click.option('--tis-range', default=15, type=int, help='Half-window around TIS/TTS in transcript coords (default: 15).')
@click.option('--flank-nt', default=60, type=int, help='Upstream/downstream chunk window size in nt (default: 60).')
def map_orfs_cmd(orfs_parquet: str, feature_map_parquet: str, features_parquet: str, out_parquet: str, progress: bool, tis_range: int, flank_nt: int):
    """Slice transcript feature chains into per-ORF chains and ranges."""
    setup_logging()
    from .pipeline.map_orfs import map_orfs
    map_orfs(orfs_parquet, feature_map_parquet, features_parquet, out_parquet, progress=progress, tis_range=tis_range, flank_nt=flank_nt)
    log_info("Mapped ORFs to feature chains")

@cli.command("inspect")
@click.option('--parquet', 'parquet_path', required=True, help='Parquet file to inspect.')
@click.option('--limit', default=5, type=int, help='Number of ORFs to preview (default: 5).')
@click.option('--expand/--no-expand', default=True, help='Expand composite parts in the preview (default: on).')
@click.option('--out-csv', type=click.Path(), help='Optional path to write an expanded CSV preview.')
def inspect_cmd(parquet_path: str, limit: int, expand: bool, out_csv: Optional[str]):
    """Inspect a Parquet file (schema, summary, and an optional expanded composite preview)."""
    setup_logging()
    from .pipeline.inspect import inspect_parquet
    inspect_parquet(parquet_path, limit=limit, expand=expand, out_csv=out_csv)

if __name__ == '__main__':
    cli()

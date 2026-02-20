from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Tuple, Union

import polars as pl

from ..utils import log_info
from .config import Config
from ..file_handlers import bam as bam_handlers
from ..file_handlers import bed as bed_handlers
from ..file_handlers import bigwig as bw_handlers
from ..file_handlers import zarr as zarr_handlers
from ..core import coordinates, orffinder
from ..visualization import plots


def process_bam_workflow(config: Config) -> Tuple[Union[str, Dict[str, str]], pl.DataFrame, pl.DataFrame]:
    """Process BAM -> bedGraph/bigWig and return (bigwig_paths, exon_df, cds_df).

    If `config.stranded` is True, returns a dict with forward/reverse bigWig paths.
    Otherwise returns a single bigWig path string.
    """
    log_info("Processing BAM file → coverage tracks…")

    # Read BAM and annotation
    bam_df = bam_handlers.readbam(
        config.bam,
        collapsed=config.bam_collapsed,
        count_from=config.bam_count_from,
        count_pattern=config.bam_count_pattern,
        count_tag=config.bam_count_tag,
        include_qname=bool(config.bam_count_from in ['name', 'tag'])
    )
    cds_df, exon_df = bam_handlers.getexons_and_cds(config.annotation)

    # Detect BAM type and normalize
    bam_type, _ = bam_handlers.detect_bam_type(bam_df, exon_df)
    if bam_type == "genomic":
        bam_df = bam_handlers.bamtranscript(bam_df, exon_df)
    # Now annotate relative to CDS
    bam_df = bam_handlers.process_transcriptomic_bam(bam_df, cds_df)

    # Offsets and A-site positions
    offsets = coordinates.change_point_analysis(bam_df)
    bed_df = bed_handlers.asitecalc(bam_df, offsets)

    # Write bedGraph and bigWig(s)
    bedgraph_path = f"{config.output}.bedGraph"
    bed_df.write_csv(bedgraph_path, separator="\t", include_header=False)

    # Current A-site output is unstranded; write a single bigWig
    bigwig_paths = bed_handlers.bedtobigwig(bedgraph_path, config.chromsizes, config.output)

    return bigwig_paths, exon_df, cds_df


def _default_offsets_from_lengths(lengths: list[int]) -> dict:
    # Conservative default if no offsets provided for Zarr lane
    return {int(L): 15 for L in set(int(x) for x in lengths)}


def process_zarr_workflow(config: Config) -> Dict[str, str]:
    """Process Zarr unique-read matrix into per-sample bigWigs.

    Returns:
        Dict[sample -> bigwig_path]
    """
    log_info("Processing Zarr unique-read matrix → per-sample coverage…")

    if not (config.zarr_root and config.read_index_parquet and config.samples):
        raise ValueError("Zarr lane requires --zarr-root, --read-index-parquet and at least one --sample")

    # Load annotation once
    cds_df, exon_df = bam_handlers.getexons_and_cds(config.annotation)

    # Determine offsets:
    # - If user provided offsets file later (not yet wired), we would load it.
    # - Else derive a conservative default=15 for observed lengths in first chunk.
    offsets = None

    # Prepare temp bedGraph writers per sample
    tmp_paths: Dict[str, str] = {s: f"{config.output}_{s}.bedGraph" for s in config.samples}
    # Clear existing temp files
    for p in tmp_paths.values():
        try:
            os.remove(p)
        except FileNotFoundError:
            pass

    # Iterate Zarr by chunks; function yields (sample, df)
    for sample, chunk_df in zarr_handlers.iter_reads_from_zarr(
        config.zarr_root,
        config.read_index_parquet,
        config.samples,
        chunk_size=1_000_000,
    ):
        if chunk_df.is_empty():
            continue

        # Derive default offsets once from observed lengths if needed
        if offsets is None:
            lengths = chunk_df.get_column("length").unique().to_list()
            offsets = _default_offsets_from_lengths(lengths)

        # Use bed.asitecalc directly on genomic starts with offsets (consistent with BAM lane)
        bed_df = bed_handlers.asitecalc(chunk_df.select(["chr", "start", "length", "count"]), offsets)
        if not bed_df.is_empty():
            # Append to per-sample bedGraph
            bed_df.write_csv(tmp_paths[sample], separator="\t", include_header=False, mode="a")

    # Convert each bedGraph to bigWig
    sample_bw: Dict[str, str] = {}
    for s, bedgraph in tmp_paths.items():
        if not os.path.exists(bedgraph):
            log_info(f"No data for sample {s}; skipping")
            continue
        bw_path = bed_handlers.bedtobigwig(bedgraph, config.chromsizes, f"{config.output}_{s}")
        sample_bw[s] = bw_path

    return sample_bw


def find_orfs_workflow(config: Config) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Find ORFs and classify relative to CDS. Returns (orf_df, exon_df)."""
    log_info("Finding ORFs…")

    # Ensure transcript FASTA if input is genomic
    transcript_fasta = config.sequence
    if not transcript_fasta.endswith("_transcripts.fa"):
        log_info("Generating transcript FASTA from genomic sequence and annotation")
        transcript_fasta = coordinates.gettranscripts(config.sequence, config.annotation, config.output)

    # Predict ORFs
    orf_df = orffinder.preporfs(
        transcript_fasta,
        start_codons=config.start_codons,
        stop_codons=config.stop_codons,
        minlength=config.min_length,
        maxlength=config.max_length,
    )

    # Classify ORFs relative to CDS and get exon coordinates
    orf_df, exon_df = coordinates.orfrelativeposition(config.annotation, orf_df)
    return orf_df, exon_df


def score_orfs_workflow(config: Config, bigwig_paths: Union[str, Dict[str, str]], exon_df: pl.DataFrame, orf_df: pl.DataFrame) -> pl.DataFrame:
    """Score ORFs using bigWig(s) and return scored DataFrame."""
    log_info("Scoring ORFs…")
    old_scoring = (config.scoring_method == "classic")
    scored = bw_handlers.scoring(bigwig_paths, exon_df, orf_df, old_scoring, config.sru_range, stranded=config.stranded)
    return scored


def plot_workflow(config: Config, scored_orfs: pl.DataFrame, exon_df: pl.DataFrame) -> None:
    """Generate plots and HTML report."""
    log_info("Generating plots and report…")
    # Persist scored ORFs to CSV for plotting convenience
    scored_path = f"{config.output}_orfs_scored.csv"
    scored_orfs.write_csv(scored_path)

    # Choose a plotting bigWig: forward if dict provided, else single path
    if isinstance(config.bigwig_paths, dict):
        plot_bw = config.bigwig_paths.get('forward') or next(iter(config.bigwig_paths.values()))
    else:
        plot_bw = config.bigwig_paths or ""

    # Pass exon_df directly; plotting can handle DataFrame input
    plots.plottop10(scored_path, plot_bw, exon_df, config.plot_range, config.output)


def all_workflow(config: Config) -> None:
    """Run the entire pipeline based on provided config."""
    log_info("Starting pipeline…")

    # Determine input lane: BigWig vs BAM (Zarr lane stubbed for future)
    bigwig_paths: Union[str, Dict[str, str], None] = config.bigwig_paths
    exon_df = None

    if config.bigwig or (config.forward_bigwig and config.reverse_bigwig):
        # BigWig lane: ensure exon_df and continue
        _, exon_df = bam_handlers.getexons_and_cds(config.annotation)
        if config.forward_bigwig and config.reverse_bigwig:
            config.bigwig_paths = {"forward": config.forward_bigwig, "reverse": config.reverse_bigwig}
            config.stranded = True
        else:
            config.bigwig_paths = config.bigwig
    elif config.bam:
        # BAM lane (classic or collapsed)
        if not config.bam:
            raise ValueError("Provide either BigWig(s) or a BAM file (classic or collapsed)")
        bigwig_paths, exon_df, _ = process_bam_workflow(config)
        config.bigwig_paths = bigwig_paths
    elif config.zarr_root:
        # Zarr lane (multi-sample)
        sample_to_bw = process_zarr_workflow(config)
        # Find ORFs once
        orf_df, exon_df = find_orfs_workflow(config)
        # Score per sample
        for sample, bw_path in sample_to_bw.items():
            scored = score_orfs_workflow(config, bw_path, exon_df, orf_df)
            out_base = f"{config.output}_{sample}"
            scored.write_csv(f"{out_base}_orfs_scored.csv")
            # Report
            plots.plottop10(f"{out_base}_orfs_scored.csv", bw_path, exon_df, config.plot_range, out_base)
        log_info("Pipeline completed successfully for Zarr samples.")
        return

    # Find ORFs and classify
    orf_df, exon_df2 = find_orfs_workflow(config)
    if exon_df is None:
        exon_df = exon_df2

    # Score ORFs
    scored = score_orfs_workflow(config, config.bigwig_paths, exon_df, orf_df)

    # Save results
    scored_path = f"{config.output}_orfs_scored.csv"
    scored.write_csv(scored_path)

    # Report
    plot_workflow(config, scored, exon_df)

    log_info("Pipeline completed successfully.")

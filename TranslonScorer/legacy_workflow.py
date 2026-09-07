from __future__ import annotations

import os
from typing import Dict, Tuple, Union

import polars as pl

from .config import Config
from .core import coordinates
from .coverage.transcript_coords import cds_to_transcript_space
from .file_handlers import bam as bam_handlers
from .file_handlers import bed as bed_handlers
from .file_handlers import zarr as zarr_handlers
from .utils import log_info


def process_bam_workflow(
    config: Config,
) -> Tuple[Union[str, Dict[str, str]], pl.DataFrame, pl.DataFrame]:
    """Process BAM -> bedGraph/bigWig and return (bigwig_paths, exon_df, cds_df).

    If `config.stranded` is True, returns a dict with forward/reverse bigWig paths.
    Otherwise returns a single bigWig path string.
    """
    log_info("Processing BAM file → coverage tracks…")

    # Read BAM and annotation
    if config.bam is None:
        raise ValueError("process-bam requires --bam")
    bam_df = bam_handlers.readbam(
        config.bam,
        collapsed=config.bam_collapsed,
        count_from=config.bam_count_from,
        count_pattern=config.bam_count_pattern,
        count_tag=config.bam_count_tag,
        include_qname=bool(config.bam_count_from in ["name", "tag"]),
    )
    cds_df, exon_df = bam_handlers.getexons_and_cds(config.annotation)

    # Detect BAM type and normalize
    bam_type, _ = bam_handlers.detect_bam_type(bam_df, exon_df)
    if bam_type == "genomic":
        bam_df = bam_handlers.bamtranscript(bam_df, exon_df)
    # Now annotate relative to CDS
    cds_tran_df = cds_to_transcript_space(cds_df, exon_df)
    bam_df = bam_handlers.process_transcriptomic_bam(bam_df, cds_tran_df)

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
        raise ValueError(
            "Zarr lane requires --zarr-root, --read-index-parquet and at least one --sample"
        )

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
        bed_df = bed_handlers.asitecalc(
            chunk_df.select(["chr", "start", "length", "count"]), offsets
        )
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

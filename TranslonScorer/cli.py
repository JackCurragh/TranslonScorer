"""Command-line interface for TranslonScorer."""

from __future__ import annotations

import os
from typing import Optional

import click
import polars as pl

from .config import Config
from .file_handlers import bam as bam_handlers
from .utils.logging import log_info, setup_logging


def _build_and_write_frame_support(profiles, cds_df, cfg):
    """Shell adapter: Config → FrameSupportParams, build frame support (pure),
    write to cfg.frame_support_out if set. Returns the support DataFrame."""
    from .frame_support import build_frame_support
    from .model import FrameSupportParams

    params = FrameSupportParams(
        frame_method=getattr(cfg, "frame_method", "none"),
        frame_by_length=bool(getattr(cfg, "frame_by_length", False)),
        frame_background=getattr(cfg, "frame_background", "uniform"),
        frame_hmm_lambda=float(getattr(cfg, "frame_hmm_lambda", 1.0)),
    )
    out = build_frame_support(profiles, cds_df, params)
    out_path = getattr(cfg, "frame_support_out", None)
    if out_path and not out.is_empty():
        out.write_parquet(out_path)
    return out


def common_options(func):
    """Apply common options to command functions."""
    # Input files
    func = click.option("--bam", "-b", help="Input BAM file from Ribo-seq data.")(func)
    func = click.option(
        "--bam-collapsed",
        is_flag=True,
        default=False,
        help="Treat BAM as collapsed (counts in name or tag).",
    )(func)
    func = click.option(
        "--bam-count-from",
        type=click.Choice(["name", "tag"]),
        help="Where to parse counts for collapsed BAM.",
    )(func)
    func = click.option(
        "--bam-count-pattern", help='Regex with named group "count" for name-based collapsed BAM.'
    )(func)
    func = click.option("--bam-count-tag", help="SAM tag for tag-based collapsed BAM (e.g. RC).")(
        func
    )
    func = click.option(
        "--chromsizes", "-c", help="Chromosome sizes file (required if processing BAM)"
    )(func)
    func = click.option(
        "--sequence", "-s", required=True, help="Input FASTA file (genomic or transcriptomic)"
    )(func)
    func = click.option(
        "--annotation",
        "-a",
        required=False,
        help="GTF annotation file for identifying exons and transcripts",
    )(func)
    func = click.option(
        "--annotation-dir", required=False, help="Annotation bundle directory (preferred)"
    )(func)

    # BigWig options
    func = click.option(
        "--bigwig",
        "-w",
        help="BigWig file containing Ribo-seq coverage. If provided, skips BAM/Zarr processing",
    )(func)
    func = click.option(
        "--forward-bigwig",
        "forward_bigwig",
        help="Forward strand BigWig file for strand-specific analysis",
    )(func)
    func = click.option(
        "--reverse-bigwig",
        "reverse_bigwig",
        help="Reverse strand BigWig file for strand-specific analysis",
    )(func)

    # Analysis options
    func = click.option(
        "--stranded/--unstranded", default=False, help="Process strands separately (default: False)"
    )(func)
    func = click.option(
        "--offsets", help="File containing read length-specific offsets for A-site calculation"
    )(func)
    func = click.option(
        "--start-codons",
        "start_codons",
        default="ATG",
        help="Comma-separated list of start codons (default: ATG)",
    )(func)
    func = click.option(
        "--stop-codons",
        "stop_codons",
        default="TAA,TAG,TGA",
        help="Comma-separated list of stop codons (default: TAA,TAG,TGA)",
    )(func)
    func = click.option(
        "--min-length",
        "min_length",
        type=int,
        default=0,
        help="Minimum ORF length in nucleotides (default: 0)",
    )(func)
    func = click.option(
        "--max-length",
        "max_length",
        type=int,
        default=1000000,
        help="Maximum ORF length in nucleotides (default: 1000000)",
    )(func)
    func = click.option(
        "--sru-range",
        "sru_range",
        type=int,
        default=15,
        help="Nucleotide range for Start Rise Up score calculation (default: 15)",
    )(func)
    func = click.option(
        "--plot-range",
        "plot_range",
        type=int,
        default=30,
        help="Plot range around start position (default: 30)",
    )(func)

    # Output options
    # Not all commands require an explicit -o (e.g., profiles when --profiles-out is given).
    # Keep optional here; commands that require -o should validate explicitly.
    func = click.option("--output", "-o", required=False, help="Base name for output files")(func)
    func = click.option(
        "--log-file",
        "log_file",
        help="Path to log file. If not provided, logs will only be written to console.",
    )(func)
    func = click.option(
        "--log-level",
        "log_level",
        type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
        default="INFO",
        help="Set the logging level (default: INFO)",
    )(func)

    return func


@click.group()
def cli():
    """TranslonScorer: identify and score translational events from Ribo-seq data.

    Event-scoring workflow (recommended):
      pipeline        one-shot: extract-events → score → report
      extract-events  annotation sqlite → genomic event store
      score-matrix    score events against the sparse annotation-scale matrix
      score-bams      score events against 1-20 genome/transcriptome BAMs
      report          compose per-translon report + consequentiality
      consequential   re-gate an existing report under a policy

    Run `translonscorer <command> --help` for details.
    """


@cli.command("process-bam")
@click.option("--bam", "-b", required=True, help="Input BAM file (genomic or transcriptomic).")
@click.option(
    "--chromsizes",
    "-c",
    required=True,
    help="Chromosome sizes file (required for bigWig conversion).",
)
@click.option("--annotation", "-a", required=True, help="GTF annotation file.")
@click.option("--output", "-o", required=True, help="Base name for output files.")
@click.option(
    "--stranded/--unstranded", default=False, help="Write strand-specific bigWigs if enabled."
)
def process_bam(bam: str, chromsizes: str, annotation: str, output: str, stranded: bool):
    """Process BAM to coverage (bedGraph/bigWig)."""
    setup_logging()
    config = Config(
        bam=bam,
        chromsizes=chromsizes,
        sequence="placeholder.fa",
        annotation=annotation,
        output=output,
        stranded=stranded,
    )
    # Minimal validation for this subcommand
    config._validate_file_exists(config.bam, "BAM")
    config._validate_file_exists(config.chromsizes, "Chromosome sizes")
    config._validate_file_exists(config.annotation, "Annotation")
    from .legacy_workflow import process_bam_workflow

    process_bam_workflow(config)
    log_info("BAM processing complete!")


@cli.command("profiles")
@common_options
@click.option(
    "--profiles-out",
    required=False,
    help="Optional output path. Defaults under -o: transcript_profiles.parquet (transcript) or locus_profiles.zarr (locus).",
)
@click.option("--loci-bed", help="BED of loci to build genomic A-site profile matrices (Zarr).")
@click.option("--offsets-file", help="CSV with columns length,offset to override defaults.")
# Matrix + indexing helpers
@click.option(
    "--zarr-root", help="Root of Zarr counts matrix (samples x reads or reads x samples)."
)
@click.option(
    "--sparse-matrix-manifest", help="Sparse Parquet global matrix manifest or output directory."
)
@click.option(
    "--matrix-scoring-mode",
    type=click.Choice(["aggregate", "per_sample"]),
    default="aggregate",
    show_default=True,
    help="Scoring mode for sparse Parquet matrix: aggregate (all samples) or per_sample.",
)
@click.option(
    "--sample-offsets-path",
    help="CSV/Parquet with per-sample A-site offsets (sample_id, length, offset) for per_sample mode.",
)
@click.option(
    "--read-index-parquet",
    help="Parquet mapping read_id (Zarr row) to genomic alignment (chr,start,stop,strand,length).",
)
@click.option(
    "--splits-index-parquet",
    help="Optional Parquet of per-read junctions (read_id, chr, donor_pos, acceptor_pos, strand).",
)
@click.option(
    "--sample", "samples", multiple=True, help="Sample name(s) to process from the Zarr matrix."
)
@click.option(
    "--all-samples",
    is_flag=True,
    default=False,
    help="Process all samples found in the Zarr counts array (uses /samples names if present).",
)
@click.option(
    "--zarr-metadata-parquet",
    help="Parquet with read row metadata (e.g., row_id,qname or sequence).",
)
@click.option(
    "--zarr-reads-fasta",
    help="FASTA with reads in Zarr row order; headers or sequences used for mapping.",
)
@click.option(
    "--zarr-read-key",
    type=click.Choice(["auto", "row_id", "qname", "sequence"]),
    default="auto",
    help="Key to align Zarr rows to BAM.",
)
@click.option(
    "--bam-key",
    type=click.Choice(["auto", "qname", "sequence"]),
    default="auto",
    help="Key to align BAM reads to Zarr.",
)
@click.option(
    "--hash-alg",
    type=click.Choice(["sha1", "md5", "xxh64"]),
    default="sha1",
    help="Hash algorithm for sequence-based mapping.",
)
# Outputs & policy
@click.option("--junctions-out", help="Path to write junction counts when available.")
@click.option("--offsets-out", help="Path to write discovered offsets per sample/length.")
@click.option(
    "--gene-expression-out", help="Optional wide gene x sample expression matrix Parquet output."
)
@click.option(
    "--gene-expression-long-out",
    help="Optional long gene expression Parquet output with gene_id,sample_id,count.",
)
@click.option(
    "--partitioned/--no-partitioned",
    default=False,
    help="Write profiles as a partitioned dataset by sample.",
)
@click.option(
    "--mapped-index",
    "mapped_index_parquet",
    help="Optional path to a precomputed mapped index (read_id→tran_id/tran_start_bam). If provided and missing, it will be built.",
)
@click.option(
    "--offsets-mode",
    type=click.Choice(["auto", "global", "required"]),
    default="auto",
    help="Offsets selection policy.",
)
# Frame assignment options (opt-in)
@click.option(
    "--frame-method",
    type=click.Choice(["none", "linear", "linear+hmm", "deblur+linear+hmm", "latent"]),
    default="none",
    help="Frame assignment method to run (default: none).",
)
@click.option(
    "--frame-by-length/--no-frame-by-length",
    default=True,
    help="Model frame per read length when available (default: on).",
)
@click.option(
    "--frame-support-out",
    help="Output Parquet for frame posteriors (default: derived from --profiles-out/-o).",
)
@click.option(
    "--profiles-with-length",
    is_flag=True,
    default=False,
    help="Retain read length in profiles (enables length-aware modeling).",
)
def profiles(**kwargs):
    """Generate transcript-space A-site profiles from BAM (classic/collapsed), Zarr, or BigWig."""
    setup_logging()
    config = Config.from_click_args(**kwargs)
    # Decide default output path when --profiles-out is omitted (do not require -o)
    if not kwargs.get("profiles_out"):
        out_base = kwargs.get("output")
        if out_base:
            # If -o provided, derive from it
            if out_base.endswith("/") or os.path.isdir(out_base):
                out_dir = out_base
                kwargs["profiles_out"] = os.path.join(
                    out_dir,
                    (
                        "locus_profiles.zarr"
                        if kwargs.get("loci_bed")
                        else "transcript_profiles.parquet"
                    ),
                )
            else:
                base = out_base.rstrip("/")
                kwargs["profiles_out"] = (
                    f"{base}_locus_profiles.zarr"
                    if kwargs.get("loci_bed")
                    else f"{base}_transcript_profiles.parquet"
                )
        else:
            # Neither --profiles-out nor -o given: write into cwd with sensible default name
            kwargs["profiles_out"] = (
                "locus_profiles.zarr" if kwargs.get("loci_bed") else "transcript_profiles.parquet"
            )
    # Minimal file checks: accept --annotation-dir (bundle) or --annotation (GTF)
    if not (config.annotation_dir or config.annotation):
        raise click.BadParameter("Provide --annotation-dir (bundle) or --annotation (GTF)")
    if config.annotation_dir:
        # Bundle provided; no need to validate GTF path here
        pass
    else:
        config._validate_file_exists(config.annotation, "Annotation")
    if not (
        config.bam
        or config.zarr_root
        or config.sparse_matrix_manifest
        or config.bigwig
        or (config.forward_bigwig and config.reverse_bigwig)
    ):
        raise click.BadParameter(
            "Provide one of: --bam (classic/collapsed), --sparse-matrix-manifest with --bam, --zarr-root with --read-index-parquet and --sample, or --bigwig/--forward-bigwig+--reverse-bigwig"
        )
    if config.sparse_matrix_manifest and not config.bam:
        raise click.BadParameter(
            "Sparse Parquet matrix mode requires --bam pointing at the global unique-read BAM."
        )

    log_info("Starting profiles workflow…")
    # Frame-related options (opt-in)
    frame_method = kwargs.get("frame_method") or "none"
    frame_by_length = bool(kwargs.get("frame_by_length")) if "frame_by_length" in kwargs else True
    frame_support_out = kwargs.get("frame_support_out")
    profiles_with_length = (
        bool(kwargs.get("profiles_with_length")) if "profiles_with_length" in kwargs else False
    )
    # Load annotation: prefer bundle if provided
    exon_df = cds_df = None
    if config.annotation_dir:
        from .io.annotation_bundle import load_annotation_bundle

        exon_df, cds_df, feats_df, fmap_df, tx_df, loci_bed, manifest = load_annotation_bundle(
            config.annotation_dir
        )
    else:
        cds_df, exon_df = bam_handlers.getexons_and_cds(config.annotation)

    from .coverage.locus_profiles import build_locus_profiles_zarr
    from .coverage.profiles import (
        gene_expression_matrix_from_profiles,
        profiles_from_bam,
        profiles_from_bigwig,
        profiles_from_sparse_parquet_matrix,
        profiles_from_zarr,
        write_profiles_parquet,
    )
    from .coverage.transcript_coords import cds_to_transcript_space

    log_info(
        "Inputs detected: "
        + (
            "Sparse Parquet"
            if config.sparse_matrix_manifest
            else (
                "BAM"
                if config.bam
                else (
                    "Zarr"
                    if config.zarr_root
                    else (
                        "BigWig"
                        if (config.bigwig or (config.forward_bigwig and config.reverse_bigwig))
                        else "Unknown"
                    )
                )
            )
        )
    )
    # Sparse Parquet lane: global unique-read BAM + sparse per-run count matrix
    if config.sparse_matrix_manifest:
        prof, offsets, _genomic_counts = profiles_from_sparse_parquet_matrix(
            bam_path=config.bam,
            manifest_path=config.sparse_matrix_manifest,
            exon_df=exon_df,
            cds_df=cds_df,
            samples=list(config.samples) if config.samples else None,
            keep_length=(profiles_with_length or (frame_by_length and (frame_method != "none"))),
        )
        log_info(f"Writing sample-resolved sparse matrix profiles to: {kwargs['profiles_out']}")
        write_profiles_parquet(prof, kwargs["profiles_out"])
        if config.offsets_out and offsets:
            rows = [
                {"sample": "all", "length": int(length), "offset": int(offset)}
                for length, offset in sorted(offsets.items())
            ]
            pl.from_dicts(rows).write_csv(config.offsets_out)
        if config.gene_expression_out or config.gene_expression_long_out:
            long, wide = gene_expression_matrix_from_profiles(
                prof,
                exon_df=exon_df,
                transcripts_df=tx_df if "tx_df" in locals() else None,
                feature_map_df=fmap_df if "fmap_df" in locals() else None,
            )
            if config.gene_expression_out:
                log_info(f"Writing wide gene expression matrix to: {config.gene_expression_out}")
                wide.write_parquet(config.gene_expression_out)
            if config.gene_expression_long_out:
                log_info(
                    f"Writing long gene expression matrix to: {config.gene_expression_long_out}"
                )
                long.write_parquet(config.gene_expression_long_out)
        log_info("Sparse Parquet profiles written")
        if frame_method and frame_method.lower() != "none":
            if not frame_support_out:
                base = kwargs["profiles_out"].rsplit(".", 1)[0]
                frame_support_out = (
                    f"{base.replace('_transcript_profiles','')}_frame_support.parquet"
                )
            cfg = Config.from_click_args(**kwargs)
            cfg.frame_method = frame_method
            cfg.frame_by_length = frame_by_length
            cfg.frame_support_out = frame_support_out
            cds_tran_df = cds_to_transcript_space(cds_df, exon_df)
            log_info(f"Building frame support with method={frame_method}...")
            _build_and_write_frame_support(prof, cds_tran_df, cfg)
            log_info(f"Frame support written to: {cfg.frame_support_out}")
    # BAM lane
    elif config.bam and not config.zarr_root:
        prof, offsets = profiles_from_bam(
            config.bam,
            exon_df,
            cds_df,
            collapsed=config.bam_collapsed,
            count_from=config.bam_count_from,
            count_pattern=config.bam_count_pattern,
            count_tag=config.bam_count_tag,
            keep_length=(profiles_with_length or (frame_by_length and (frame_method != "none"))),
        )
        log_info(f"Writing profiles to: {kwargs['profiles_out']}")
        write_profiles_parquet(prof, kwargs["profiles_out"])
        if config.offsets_out:
            # Persist discovered offsets
            odf = pl.from_dicts(
                [
                    {
                        "sample": kwargs.get("sample") or "sample",
                        "length": int(L),
                        "offset": int(ofs),
                    }
                    for L, ofs in offsets.items()
                ]
            )
            odf.write_csv(config.offsets_out)
        if config.junctions_out:
            # Aggregate junctions from BAM
            from .coverage.junctions import aggregate_bam_junctions

            log_info("Aggregating junctions from BAM…")
            j = aggregate_bam_junctions(config.bam)
            if not j.is_empty():
                from .utils.io import write_parquet_safe

                log_info(f"Writing junctions to: {config.junctions_out}")
                write_parquet_safe(j, config.junctions_out)
                log_info("Junction counts written from BAM")
        log_info("Profiles written")
        # Optional frame support
        if frame_method and frame_method.lower() != "none":
            # Derive default path if not provided
            if not frame_support_out:
                base = kwargs["profiles_out"].rsplit(".", 1)[0]
                frame_support_out = (
                    f"{base.replace('_transcript_profiles','')}_frame_support.parquet"
                )
            cfg = Config.from_click_args(**kwargs)
            cfg.frame_method = frame_method
            cfg.frame_by_length = frame_by_length
            cfg.frame_support_out = frame_support_out
            cds_tran_df = cds_to_transcript_space(cds_df, exon_df)
            log_info(f"Building frame support with method={frame_method}…")
            _build_and_write_frame_support(prof, cds_tran_df, cfg)
            log_info(f"Frame support written to: {cfg.frame_support_out}")
    # Zarr lane (with optional index build from BAM)
    elif config.zarr_root:
        # Ensure samples provided
        if not config.samples or len(config.samples) == 0:
            if kwargs.get("all_samples"):
                # Discover sample names from Zarr store (v2/v3; group/array roots)
                try:
                    # Reuse robust opener to locate counts array
                    from .file_handlers.zarr import (
                        _open_counts as _open_counts_helper,  # type: ignore
                    )
                except Exception:
                    _open_counts_helper = None
                names: list[str] | None = None
                if _open_counts_helper is not None:
                    try:
                        arr = _open_counts_helper(config.zarr_root)
                        # Heuristic: smaller dim = samples
                        n = arr.shape[0] if arr.shape[0] <= arr.shape[1] else arr.shape[1]
                        names = [str(i) for i in range(n)]
                    except Exception:
                        names = None
                if names is None:
                    # Fallback: open root, prefer group attrs['samples'] if present
                    import importlib

                    zmod = importlib.import_module("zarr")
                    store = zmod.open(config.zarr_root, mode="r")
                    # Prefer root attribute 'samples' (Zarr v3/v2 attrs)
                    try:
                        attrs = getattr(store, "attrs", None)
                        if attrs is not None and "samples" in attrs:
                            names_arr = attrs["samples"]
                            names = [
                                n.decode() if isinstance(n, (bytes, bytearray)) else str(n)
                                for n in list(names_arr)
                            ]
                        else:
                            raise KeyError("no attrs[samples]")
                    except Exception:
                        # Otherwise: if group has a dataset named 'samples', use it; else fall back to counts shape
                        grp = store if hasattr(store, "__contains__") else None
                        if grp is not None and "samples" in grp:
                            names_arr = grp["samples"][...]
                            try:
                                names = [
                                    n.decode() if isinstance(n, (bytes, bytearray)) else str(n)
                                    for n in list(names_arr)
                                ]
                            except Exception:
                                names = [str(n) for n in list(names_arr)]
                        else:
                            counts = zmod.open(config.zarr_root.rstrip("/") + "/counts", mode="r")
                            n = (
                                counts.shape[0]
                                if counts.shape[0] <= counts.shape[1]
                                else counts.shape[1]
                            )
                            names = [str(i) for i in range(n)]
                config.samples = names
            else:
                raise click.BadParameter(
                    "Zarr mode requires at least one --sample or use --all-samples"
                )
        # If index missing but BAM is present, build index now
        if not config.read_index_parquet and config.bam:
            from .coverage.index_from_bam import build_read_index_from_bam

            log_info("Building Zarr read index from BAM…")
            out_idx = (
                os.path.join(os.path.dirname(kwargs["profiles_out"]), "read_index.parquet")
                if kwargs.get("profiles_out")
                else "read_index.parquet"
            )
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
            raise click.BadParameter(
                "Zarr mode requires --read-index-parquet (or provide --bam to build it)."
            )
        if kwargs.get("loci_bed"):
            # Locus matrices path
            build_locus_profiles_zarr(
                config.zarr_root,
                config.read_index_parquet,
                list(config.samples),
                kwargs["loci_bed"],
                offsets_file=kwargs.get("offsets_file"),
                out_zarr=kwargs["profiles_out"],
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
                mapped_index_parquet=kwargs.get("mapped_index_parquet"),
            ):
                out = kwargs["profiles_out"]
                stem, ext = (out.rsplit(".", 1) + ["parquet"])[:2]
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
                    j = s.group_by(["chr", "donor_pos", "acceptor_pos", "strand"]).agg(
                        pl.len().alias("count")
                    )
                    j.write_parquet(config.junctions_out)
                    log_info("Junction counts written from splits index")
                else:
                    log_info(
                        "No splits index provided; skipping junction aggregation for Zarr lane"
                    )
    else:
        # BigWig lane: single or stranded inputs
        bw = (
            {"forward": config.forward_bigwig, "reverse": config.reverse_bigwig}
            if (config.forward_bigwig and config.reverse_bigwig)
            else config.bigwig
        )
        log_info("Computing profiles from BigWig…")
        prof = profiles_from_bigwig(bw, exon_df, stranded=config.stranded)
        log_info(f"Writing profiles to: {kwargs['profiles_out']}")
        write_profiles_parquet(prof, kwargs["profiles_out"])
        log_info("Profiles written from BigWig")
        # Optional frame support from BigWig-derived profiles
        if frame_method and frame_method.lower() != "none":
            if not frame_support_out:
                base = kwargs["profiles_out"].rsplit(".", 1)[0]
                frame_support_out = (
                    f"{base.replace('_transcript_profiles','')}_frame_support.parquet"
                )
            cfg = Config.from_click_args(**kwargs)
            cfg.frame_method = frame_method
            cfg.frame_by_length = frame_by_length
            cfg.frame_support_out = frame_support_out
            cds_tran_df = cds_to_transcript_space(cds_df, exon_df)
            log_info(f"Building frame support with method={frame_method} (BigWig lane)…")
            _build_and_write_frame_support(prof, cds_tran_df, cfg)
            log_info(f"Frame support written to: {cfg.frame_support_out}")


@cli.command("score-compare-frame")
@click.option("--orfs", required=True, help="ORFs to score (CSV/TSV/Parquet).")
@click.option("--exons", required=True, help="Transcript exon table (CSV/TSV/Parquet).")
@click.option("--bigwig", required=True, help="BigWig coverage track.")
@click.option(
    "--frame-support",
    "frame_support",
    required=True,
    help="Frame support Parquet from the profiles command.",
)
@click.option("--out-prefix", required=True, help="Output prefix for raw/frame/comparison tables.")
@click.option(
    "--profiles",
    "profiles_path",
    help="Optional transcript profiles Parquet for frame-weighted-count summaries.",
)
@click.option(
    "--panel-manifest", help="Optional frozen panel manifest to merge onto ORFs before scoring."
)
@click.option("--scoring-method", type=click.Choice(["classic", "modern"]), default="modern")
@click.option("--sru-range", type=int, default=15)
@click.option("--label-column", help="Optional truth/category column for score-gap summaries.")
@click.option("--max-workers", type=int, default=None)
def score_compare_frame_cmd(
    orfs: str,
    exons: str,
    bigwig: str,
    frame_support: str,
    out_prefix: str,
    profiles_path: Optional[str],
    panel_manifest: Optional[str],
    scoring_method: str,
    sru_range: int,
    label_column: Optional[str],
    max_workers: Optional[int],
):
    """Gate 1: compare raw ORF scores with frame-weighted scores."""
    setup_logging()
    from .orf.score_gates import compare_raw_frame_scoring

    paths = compare_raw_frame_scoring(
        orfs_path=orfs,
        exons_path=exons,
        bigwig_path=bigwig,
        frame_support_path=frame_support,
        out_prefix=out_prefix,
        scoring_method=scoring_method,
        sru_range=sru_range,
        profiles_path=profiles_path,
        panel_manifest_path=panel_manifest,
        label_column=label_column,
        max_workers=max_workers,
    )
    for label, path in paths.items():
        log_info(f"{label}: {path}")


@cli.command("validate-panel")
@click.option(
    "--panel-manifest", required=True, help="Panel manifest to validate/freeze (CSV/TSV/Parquet)."
)
@click.option("--out", "out_csv", help="Optional CSV path for the validation/freeze report.")
def validate_panel_cmd(panel_manifest: str, out_csv: Optional[str]):
    """Validate a frozen score panel manifest and report its SHA256."""
    setup_logging()
    from .orf.panel_manifest import panel_freeze_report

    report = panel_freeze_report(panel_manifest)
    if out_csv:
        from .utils.io import write_csv_safe

        write_csv_safe(report, out_csv)
        log_info(f"panel report: {out_csv}")
    else:
        click.echo(report)


@cli.command("score-compare-existing")
@click.option("--raw-scores", required=True, help="Existing raw score table (CSV/TSV/Parquet).")
@click.option(
    "--frame-scores", required=True, help="Existing frame-weighted score table (CSV/TSV/Parquet)."
)
@click.option("--out-prefix", required=True, help="Output prefix for comparison tables.")
@click.option("--label-column", help="Optional truth/category column for score-gap summaries.")
def score_compare_existing_cmd(
    raw_scores: str, frame_scores: str, out_prefix: str, label_column: Optional[str]
):
    """Compare already generated raw and frame-weighted score tables."""
    setup_logging()
    from .orf.score_gates import summarize_existing_score_pair

    paths = summarize_existing_score_pair(
        raw_scores_path=raw_scores,
        frame_scores_path=frame_scores,
        out_prefix=out_prefix,
        label_column=label_column,
    )
    for label, path in paths.items():
        log_info(f"{label}: {path}")


@cli.command("compare-profiles")
@click.option("--profile-a", required=True, help="First transcript profile table.")
@click.option("--profile-b", required=True, help="Second transcript profile table.")
@click.option("--out-prefix", required=True, help="Output prefix for profile comparison.")
@click.option("--label-a", default="bigwig", help="Label for first profile table.")
@click.option("--label-b", default="bam", help="Label for second profile table.")
@click.option(
    "--write-deltas/--no-write-deltas", default=False, help="Write per-position delta table."
)
def compare_profiles_cmd(
    profile_a: str, profile_b: str, out_prefix: str, label_a: str, label_b: str, write_deltas: bool
):
    """Gate 2: compare two transcript-space profile tables."""
    setup_logging()
    from .orf.profile_compare import compare_profile_files

    paths = compare_profile_files(
        profile_a_path=profile_a,
        profile_b_path=profile_b,
        out_prefix=out_prefix,
        label_a=label_a,
        label_b=label_b,
        write_deltas=write_deltas,
    )
    for label, path in paths.items():
        log_info(f"{label}: {path}")


@cli.command("compare-frame-methods")
@click.option(
    "--profiles",
    required=True,
    help="Transcript-space profile table with tran_id, pos, count[, length].",
)
@click.option(
    "--cds",
    "cds_path",
    required=True,
    help="Transcript-space CDS table with tran_id,start,stop or tran_start,tran_stop.",
)
@click.option(
    "--out-prefix", required=True, help="Output prefix for frame-support and comparison tables."
)
@click.option(
    "--methods",
    default="linear,latent",
    help="Comma-separated methods to compare; first is baseline.",
)
@click.option(
    "--frame-by-length/--no-frame-by-length",
    default=True,
    help="Model frame leakage per read length when available.",
)
@click.option(
    "--hmm-lambda", type=float, default=2.0, help="HMM smoothing strength for +hmm methods."
)
@click.option(
    "--background",
    type=click.Choice(["flat", "zero"]),
    default="flat",
    help="Background model for latent EM.",
)
@click.option(
    "--trim-nt",
    type=int,
    default=30,
    help="Trim this many nucleotides from CDS ends for validation summaries.",
)
@click.option(
    "--write-validation-rows/--no-write-validation-rows",
    default=False,
    help="Write CDS-interior per-row validation tables.",
)
def compare_frame_methods_cmd(
    profiles: str,
    cds_path: str,
    out_prefix: str,
    methods: str,
    frame_by_length: bool,
    hmm_lambda: float,
    background: str,
    trim_nt: int,
    write_validation_rows: bool,
):
    """Compare frame-only correction methods, e.g. linear versus latent EM."""
    setup_logging()
    from .frame.frame_method_compare import compare_frame_methods

    method_list = [m.strip().lower() for m in methods.split(",") if m.strip()]
    paths = compare_frame_methods(
        profiles_path=profiles,
        cds_path=cds_path,
        out_prefix=out_prefix,
        methods=method_list,
        frame_by_length=frame_by_length,
        hmm_lambda=hmm_lambda,
        background=background,
        trim_nt=trim_nt,
        write_validation_rows=write_validation_rows,
    )
    for label, path in paths.items():
        log_info(f"{label}: {path}")


@cli.command("export-rdg-flux")
@click.option(
    "--profiles",
    required=True,
    help="Transcript profile parquet/CSV with tran_id,pos,count[, sample/length].",
)
@click.option(
    "--out",
    "out_path",
    required=True,
    help="Output Parquet path, or directory when --partition-by-sample is set.",
)
@click.option(
    "--frame-support",
    "frame_support_path",
    help="Existing frame_support parquet/CSV. If omitted, build from --cds/--annotation/--annotation-dir.",
)
@click.option(
    "--sample-id",
    help="Sample/replicate id to add when profiles do not already contain sample_id/sample.",
)
@click.option(
    "--cds",
    "cds_path",
    help="Transcript-space CDS table. Used only when --frame-support is omitted.",
)
@click.option(
    "--annotation",
    "annotation_path",
    help="GTF annotation. Used only when --frame-support is omitted.",
)
@click.option(
    "--annotation-dir",
    help="Annotation bundle directory. Used only when --frame-support is omitted.",
)
@click.option("--transcriptome-fasta", help="FASTA path recorded in RDG-Flux metadata.")
@click.option(
    "--psite-offset-model", default="unknown", help="P-site offset model/path recorded in metadata."
)
@click.option("--annotation-source", help="Annotation source string/path recorded in metadata.")
@click.option(
    "--normalization", default="raw_psite_counts", help="Normalization label recorded in metadata."
)
@click.option(
    "--model-stage", default="stage1_frame_posterior", help="Model stage recorded in metadata."
)
@click.option(
    "--frame-method",
    type=click.Choice(["linear", "linear+hmm", "deblur+linear+hmm", "latent"]),
    default="linear+hmm",
    help="Frame method to build when --frame-support is omitted.",
)
@click.option(
    "--frame-by-length/--no-frame-by-length",
    default=True,
    help="Use read-length-specific frame evidence when length is available.",
)
@click.option(
    "--background-probability",
    type=float,
    default=0.0,
    help="Constant background posterior mass for v1 exports without a learned background model.",
)
@click.option(
    "--background-model", default="none", help="Background model label recorded in metadata."
)
@click.option(
    "--preserve-read-length",
    is_flag=True,
    default=False,
    help="Emit read_length_bin and keep length-specific rows when profiles/frame support contain length.",
)
@click.option(
    "--partition-by-sample/--single-file",
    default=False,
    help="Write Hive-style sample_id=<id>/part.parquet directories.",
)
@click.option(
    "--metadata-out", "metadata_path", help="Sidecar metadata JSON path. Defaults next to output."
)
def export_rdg_flux_cmd(
    profiles: str,
    out_path: str,
    frame_support_path: Optional[str],
    sample_id: Optional[str],
    cds_path: Optional[str],
    annotation_path: Optional[str],
    annotation_dir: Optional[str],
    transcriptome_fasta: Optional[str],
    psite_offset_model: str,
    annotation_source: Optional[str],
    normalization: str,
    model_stage: str,
    frame_method: str,
    frame_by_length: bool,
    background_probability: float,
    background_model: str,
    preserve_read_length: bool,
    partition_by_sample: bool,
    metadata_path: Optional[str],
):
    """Export RDG-Flux v1 per-position frame posterior substrate."""
    setup_logging()
    from .orf.rdg_flux_export import export_rdg_flux_v1

    paths = export_rdg_flux_v1(
        profiles_path=profiles,
        out_path=out_path,
        frame_support_path=frame_support_path,
        sample_id=sample_id,
        annotation_path=annotation_path,
        annotation_dir=annotation_dir,
        cds_path=cds_path,
        transcriptome_fasta=transcriptome_fasta,
        psite_offset_model=psite_offset_model,
        annotation_source=annotation_source,
        normalization=normalization,
        model_stage=model_stage,
        frame_method=frame_method,
        frame_by_length=frame_by_length,
        background_probability=background_probability,
        background_model=background_model,
        preserve_read_length=preserve_read_length,
        partition_by_sample=partition_by_sample,
        metadata_path=metadata_path,
    )
    for label, path in paths.items():
        log_info(f"{label}: {path}")


@cli.command("frame-disambiguation")
@click.option("--bam", help="BAM to project onto transcript candidates.")
@click.option("--annotation", "-a", help="GTF annotation used for transcript models.")
@click.option(
    "--candidates",
    help="Precomputed candidate table with read_key, tran_id, and transcript position.",
)
@click.option("--cds", "cds_path", help="Transcript-space CDS table for --candidates mode.")
@click.option("--out-prefix", required=True, help="Output prefix for disambiguation tables.")
@click.option("--bam-collapsed", is_flag=True, default=False)
@click.option("--bam-count-from", type=click.Choice(["name", "tag"]))
@click.option("--bam-count-pattern")
@click.option("--bam-count-tag")
def frame_disambiguation_cmd(
    bam: Optional[str],
    annotation: Optional[str],
    candidates: Optional[str],
    cds_path: Optional[str],
    out_prefix: str,
    bam_collapsed: bool,
    bam_count_from: Optional[str],
    bam_count_pattern: Optional[str],
    bam_count_tag: Optional[str],
):
    """Gate 3: quantify frame-discordant ambiguous read assignments."""
    setup_logging()
    from .frame.frame_disambiguation import (
        frame_disambiguation_from_bam,
        frame_disambiguation_from_candidates,
    )

    if candidates:
        if not cds_path:
            raise click.BadParameter("--cds is required with --candidates")
        paths = frame_disambiguation_from_candidates(
            candidates_path=candidates,
            cds_path=cds_path,
            out_prefix=out_prefix,
        )
    else:
        if not bam or not annotation:
            raise click.BadParameter("Provide either --candidates + --cds or --bam + --annotation")
        paths = frame_disambiguation_from_bam(
            bam_path=bam,
            annotation_path=annotation,
            out_prefix=out_prefix,
            collapsed=bam_collapsed,
            count_from=bam_count_from,
            count_pattern=bam_count_pattern,
            count_tag=bam_count_tag,
        )
    for label, path in paths.items():
        log_info(f"{label}: {path}")


@cli.command("compare-read-assignment")
@click.option(
    "--candidates",
    required=True,
    help="Candidate table with read_key/read_id, tran_id, transcript position, and optional truth columns.",
)
@click.option("--out-prefix", required=True, help="Output prefix for assignment comparison tables.")
@click.option(
    "--frame-support",
    "frame_support_path",
    help="Optional frame_support table with tran_id,codon,p0,p1,p2.",
)
@click.option(
    "--cds", "cds_path", help="Optional transcript-space CDS table. Used for diagnostics only."
)
@click.option(
    "--methods",
    default="unique,fractional,em,frame_em",
    help="Comma-separated methods: unique,fractional,frame_fractional,em,frame_em,rdg_local_frame_em,rdg_gated_frame_em.",
)
@click.option(
    "--abundance-key",
    default="tran_id",
    help="Column whose abundance is estimated, e.g. tran_id for isoforms or locus_id for genomic origin.",
)
@click.option(
    "--truth-column",
    default="is_true",
    help="Optional boolean truth column for synthetic/benchmark summaries.",
)
@click.option(
    "--write-assignments/--no-write-assignments",
    default=False,
    help="Write per-candidate posterior assignment Parquets.",
)
@click.option("--max-iter", type=int, default=100, help="Maximum EM iterations for em/frame_em.")
@click.option(
    "--tol", type=float, default=1e-7, help="EM convergence tolerance on target abundance L1 delta."
)
@click.option(
    "--abundance-prior",
    type=float,
    default=1e-3,
    help="Small symmetric abundance prior for EM targets.",
)
@click.option(
    "--min-frame-support-count",
    type=float,
    default=0.0,
    help="Minimum frame-support total_count required before frame likelihood can affect assignment.",
)
@click.option(
    "--min-frame-support-evidence",
    type=float,
    default=0.0,
    help="Minimum support_evidence required before frame likelihood can affect assignment.",
)
@click.option(
    "--require-complete-frame-support/--allow-partial-frame-support",
    default=False,
    help="Only use frame likelihood for a read when every assignment target has gated frame support.",
)
@click.option(
    "--frame-gate-min-likelihood-range",
    type=float,
    default=0.25,
    help="Minimum read-local frame likelihood range for rdg_*_frame_em gate.",
)
@click.option(
    "--frame-gate-min-read-fraction",
    type=float,
    default=0.25,
    help="Minimum read-weighted local gate fraction for rdg_gated_frame_em consensus.",
)
def compare_read_assignment_cmd(
    candidates: str,
    out_prefix: str,
    frame_support_path: Optional[str],
    cds_path: Optional[str],
    methods: str,
    abundance_key: str,
    truth_column: str,
    write_assignments: bool,
    max_iter: int,
    tol: float,
    abundance_prior: float,
    min_frame_support_count: float,
    min_frame_support_evidence: float,
    require_complete_frame_support: bool,
    frame_gate_min_likelihood_range: float,
    frame_gate_min_read_fraction: float,
):
    """Compare unique, fractional, EM, and frame-aware read assignment."""
    setup_logging()
    from .frame.read_assignment import compare_read_assignment_methods

    method_list = [m.strip().lower() for m in methods.split(",") if m.strip()]
    paths = compare_read_assignment_methods(
        candidates_path=candidates,
        out_prefix=out_prefix,
        frame_support_path=frame_support_path,
        cds_path=cds_path,
        methods=method_list,
        abundance_key=abundance_key,
        truth_column=truth_column,
        write_assignments=write_assignments,
        max_iter=max_iter,
        tol=tol,
        abundance_prior=abundance_prior,
        min_frame_support_count=min_frame_support_count,
        min_frame_support_evidence=min_frame_support_evidence,
        require_complete_frame_support=require_complete_frame_support,
        frame_gate_min_likelihood_range=frame_gate_min_likelihood_range,
        frame_gate_min_read_fraction=frame_gate_min_read_fraction,
    )
    for label, path in paths.items():
        log_info(f"{label}: {path}")


@cli.command("plot")
@click.option("--scored-orfs", "-s", required=True, help="CSV file of scored ORFs")
@click.option("--bigwig", "-w", required=True, help="BigWig file containing Ribo-seq coverage")
@click.option("--exons", "-e", required=True, help="CSV file containing exon positions")
@click.option(
    "--plot-range", default=30, type=int, help="Plot range around start position (default: 30)"
)
@click.option("--output", "-o", required=True, help="Base name for output files")
def plot(scored_orfs: str, bigwig: str, exons: str, plot_range: int, output: str):
    """Generate visualization reports from scored ORFs."""
    setup_logging()
    from .visualization import plots

    plots.plottop10(scored_orfs, bigwig, exons, plot_range, output)
    log_info("Report generation complete!")


@cli.command("index-from-bam")
@click.option("--bam", "-b", required=True, help="Input BAM used to align reads.")
@click.option(
    "--zarr-root",
    required=True,
    help="Root of Zarr counts matrix (samples x reads or reads x samples).",
)
@click.option("--out-index", required=True, help="Output Parquet path for read_id→alignment index.")
@click.option(
    "--zarr-metadata-parquet",
    help="Parquet with read row metadata (e.g., row_id,qname or sequence).",
)
@click.option(
    "--zarr-reads-fasta",
    help="FASTA with reads in Zarr row order; headers or sequences used for mapping.",
)
@click.option(
    "--zarr-read-key", type=click.Choice(["auto", "row_id", "qname", "sequence"]), default="auto"
)
@click.option("--bam-key", type=click.Choice(["auto", "qname", "sequence"]), default="auto")
@click.option("--hash-alg", type=click.Choice(["sha1", "md5", "xxh64"]), default="sha1")
def index_from_bam_cmd(
    bam: str,
    zarr_root: str,
    out_index: str,
    zarr_metadata_parquet: str | None,
    zarr_reads_fasta: str | None,
    zarr_read_key: str,
    bam_key: str,
    hash_alg: str,
):
    """Build read_index.parquet for Zarr counts by aligning rows to BAM reads."""
    setup_logging()
    from .coverage.index_from_bam import build_read_index_from_bam

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
@click.option("--bed12", required=True, help="Input ORFs in BED12 format.")
@click.option(
    "--annotation", "-a", required=True, help="GTF annotation file for exon/transcript models."
)
@click.option("--out", "out_parquet", required=True, help="Output canonical ORFs (Parquet).")
@click.option(
    "--assign-policy",
    type=click.Choice(["best", "all"]),
    default="best",
    help="Assignment policy when multiple transcripts match (default: best).",
)
@click.option(
    "--require-junction-match/--allow-junction-mismatch",
    default=True,
    help="Require all BED12 junctions to match transcript junctions (default: require).",
)
@click.option("--progress/--no-progress", default=True, help="Show progress bars (default: on).")
def orfs_import_cmd(
    bed12: str,
    annotation: str,
    out_parquet: str,
    assign_policy: str,
    require_junction_match: bool,
    progress: bool,
):
    """Import ORFs from BED12, map to transcripts, and write canonical ORFs (Parquet)."""
    setup_logging()
    from .orf.orfs_import import import_bed12

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
@click.option("-a", "--annotation", required=True, help="GTF annotation file.")
@click.option("--out-dir", required=False, help="Output annotation bundle directory (recommended).")
@click.option(
    "-o",
    "--output",
    "--output-prefix",
    "output_prefix",
    required=False,
    help="Legacy output prefix for feature tables (Parquet).",
)
@click.option("--progress/--no-progress", default=True, help="Show progress bars (default: on).")
def features(
    annotation: str, out_dir: str = None, output_prefix: str = None, progress: bool = True
):
    """Build an annotation bundle (exons, CDS, features, feature_map, transcripts, loci, manifest)."""
    setup_logging()
    # Backward-compat mode: if only -o is provided, write old files and also emit a bundle next to it
    if not out_dir and output_prefix:
        import os

        base = os.path.basename(output_prefix)
        root = os.path.dirname(output_prefix) or "."
        out_dir = os.path.join(root, f"{base}.annotation")
    if not out_dir and not output_prefix:
        raise click.BadParameter("Provide either --out-dir for the bundle or -o for legacy outputs")

    # Always build the bundle
    from .io.annotation_bundle import build_annotation_bundle

    bdir, paths = build_annotation_bundle(annotation, out_dir, progress=progress)
    log_info(f"Annotation bundle written at: {bdir}")

    # Write legacy outputs if requested
    if output_prefix:
        import polars as pl

        feats = pl.read_parquet(paths["features"])
        fmap = pl.read_parquet(paths["feature_map"])
        feats.write_parquet(f"{output_prefix}.features.parquet")
        fmap.write_parquet(f"{output_prefix}.feature_map.parquet")
        log_info("Legacy feature tables written")


@cli.command("assemble")
@click.option(
    "--orfs-parquet", required=True, help="ORF candidates with composite scores (Parquet)."
)
@click.option("--out", "out_parquet", required=True, help="Output assembled translome Parquet.")
@click.option("--solver", type=click.Choice(["PULP", "GREEDY"]), default="PULP")
@click.option("--timeout", "timeout_sec", type=int, default=60)
def assemble_cmd(orfs_parquet: str, out_parquet: str, solver: str, timeout_sec: int):
    """Assemble a locus translome by selecting a consistent set of ORFs under soft penalties."""
    setup_logging()
    from .orf.assemble import assemble_translome

    assemble_translome(orfs_parquet, out_parquet, solver=solver, timeout_sec=timeout_sec)
    log_info("Translome assembly complete")


@cli.command("map-orfs")
@click.option(
    "--orfs",
    "orfs_parquet",
    required=True,
    help="Canonical ORFs with tran_id,start_pos_tran,stop_pos_tran (Parquet).",
)
@click.option(
    "--feature-map",
    "feature_map_parquet",
    required=True,
    help="Transcript→feature mapping Parquet.",
)
@click.option("--features", "features_parquet", required=True, help="Feature table Parquet.")
@click.option("--out", "out_parquet", required=True, help="Output per-ORF feature chains Parquet.")
@click.option("--progress/--no-progress", default=True, help="Show progress bars (default: on).")
@click.option(
    "--tis-range",
    default=15,
    type=int,
    help="Half-window around TIS/TTS in transcript coords (default: 15).",
)
@click.option(
    "--flank-nt",
    default=60,
    type=int,
    help="Upstream/downstream chunk window size in nt (default: 60).",
)
def map_orfs_cmd(
    orfs_parquet: str,
    feature_map_parquet: str,
    features_parquet: str,
    out_parquet: str,
    progress: bool,
    tis_range: int,
    flank_nt: int,
):
    """Slice transcript feature chains into per-ORF chains and ranges."""
    setup_logging()
    from .orf.map_orfs import map_orfs

    map_orfs(
        orfs_parquet,
        feature_map_parquet,
        features_parquet,
        out_parquet,
        progress=progress,
        tis_range=tis_range,
        flank_nt=flank_nt,
    )
    log_info("Mapped ORFs to feature chains")


@cli.command("inspect")
@click.option("--parquet", "parquet_path", required=True, help="Parquet file to inspect.")
@click.option("--limit", default=5, type=int, help="Number of ORFs to preview (default: 5).")
@click.option(
    "--expand/--no-expand",
    default=True,
    help="Expand composite parts in the preview (default: on).",
)
@click.option(
    "--out-csv", type=click.Path(), help="Optional path to write an expanded CSV preview."
)
def inspect_cmd(parquet_path: str, limit: int, expand: bool, out_csv: Optional[str]):
    """Inspect a Parquet file (schema, summary, and an optional expanded composite preview)."""
    setup_logging()
    from .io.inspect import inspect_parquet

    inspect_parquet(parquet_path, limit=limit, expand=expand, out_csv=out_csv)


# ===========================================================================
# Event-scoring workflows (functional-core/imperative-shell architecture).
# These wrap TranslonScorer.workflows and supersede the ORF-composite commands.
# ===========================================================================


@cli.command("extract-events")
@click.option(
    "--out-dir",
    required=True,
    help="Output directory for events/, feature_event/, event_overlap/ Parquet trees.",
)
# --- feature source (exactly one) ---
@click.option("--gtf", "gtf_path", help="GTF/GFF annotation; scores --feature-type features.")
@click.option(
    "--feature-type",
    default="CDS",
    show_default=True,
    help="GTF feature type to score (e.g. CDS, exon).",
)
@click.option(
    "--bed12", "bed12_path", help="BED12 file (blockSizes/blockStarts give exon structure)."
)
@click.option("--bigbed", "bigbed_path", help="bigBed (BED12) file.")
@click.option(
    "--fasta", "fasta_path", help="FASTA for de-novo ORF finding (genome/contig records)."
)
@click.option("--sqlite", "sqlite_path", help="Annotation sqlite (translons + translon_blocks).")
# --- de-novo ORF options (with --fasta) ---
@click.option(
    "--start-codons",
    default="ATG",
    show_default=True,
    help="Comma-separated start codons (--fasta).",
)
@click.option(
    "--stop-codons",
    default="TAA,TAG,TGA",
    show_default=True,
    help="Comma-separated stop codons (--fasta).",
)
@click.option(
    "--min-len", type=int, default=0, show_default=True, help="Minimum ORF length nt (--fasta)."
)
@click.option(
    "--max-len",
    type=int,
    default=1_000_000,
    show_default=True,
    help="Maximum ORF length nt (--fasta).",
)
@click.option(
    "--annotation-version",
    default="",
    help="Annotation version string stamped into event records (for reproducible event_ids).",
)
@click.option(
    "--chrom",
    "chroms",
    multiple=True,
    help="Restrict extraction to these chromosome(s) (repeatable; default: all).",
)
def extract_events_cmd(
    out_dir,
    gtf_path,
    feature_type,
    bed12_path,
    bigbed_path,
    fasta_path,
    sqlite_path,
    start_codons,
    stop_codons,
    min_len,
    max_len,
    annotation_version,
    chroms,
):
    """Extract deduplicated genomic events from a feature source.

    Exactly one of --gtf / --bed12 / --bigbed / --fasta / --sqlite. No annotation
    database is required: score annotated CDSs from a GTF, your own ORFs from a
    BED12/bigBed, or find ORFs de-novo in a FASTA.
    """
    setup_logging()
    from .workflows import extract_events_workflow

    n_src = sum(bool(x) for x in (gtf_path, bed12_path, bigbed_path, fasta_path, sqlite_path))
    if n_src != 1:
        raise click.BadParameter(
            "provide exactly one source: --gtf | --bed12 | --bigbed | --fasta | --sqlite"
        )
    summary = extract_events_workflow(
        out_dir,
        sqlite_path=sqlite_path or None,
        gtf_path=gtf_path or None,
        feature_type=feature_type,
        bed12_path=bed12_path or None,
        bigbed_path=bigbed_path or None,
        fasta_path=fasta_path or None,
        start_codons=[c.strip() for c in start_codons.split(",") if c.strip()],
        stop_codons=[c.strip() for c in stop_codons.split(",") if c.strip()],
        min_len=min_len,
        max_len=max_len,
        annotation_version=annotation_version,
        chroms=list(chroms) or None,
    )
    log_info(
        f"Extracted events: {summary['events']} events, "
        f"{summary['feature_event']} feature links, "
        f"{summary['event_overlap']} overlaps across {summary['chroms']} chromosomes"
    )
    for t, n in sorted(summary.get("by_type", {}).items()):
        log_info(f"  {t}: {n}")


@cli.command("build-matrix-cache")
@click.option("--matrix-dir", required=True, help="Matrix ROOT directory (all partitions).")
@click.option("--n-workers", type=int, default=None, help="Worker processes (default: all cores).")
def build_matrix_cache_cmd(matrix_dir, n_workers):
    """Pre-build the per-partition (read_id, total_count) cache (one-time).

    After this, aggregate score-matrix reads compact per-read totals instead of
    the full per-sample count matrix — independent of cohort size. Set
    TS_MATRIX_CACHE_DIR for a read-only matrix.
    """
    setup_logging()
    from .matrix_rollup import build_read_totals_cache

    summary = build_read_totals_cache(matrix_dir, n_workers=n_workers)
    log_info(
        f"cache: {summary['built']} built, {summary['already_cached']} existing, "
        f"{summary['empty']} empty ({summary['partitions']} partitions)"
    )


@cli.command("score-matrix")
@click.option("--events-dir", required=True, help="Events directory produced by extract-events.")
@click.option(
    "--matrix-dir",
    help="Matrix ROOT directory; ALL partitions under it are scanned together "
    "(the matrix is sharded by read sequence and must always be used in full).",
)
@click.option(
    "--partitions",
    multiple=True,
    help="Explicit partition directories (advanced; normally use --matrix-dir). "
    "If used you must pass every partition — a subset gives sequence-biased coverage.",
)
@click.option("--store-dir", required=True, help="Output fact_event_score store directory.")
@click.option(
    "--data-version", required=True, help="Data version label for the append-only score partition."
)
@click.option("--annotation-version", default="", help="Annotation version stamped into the store.")
@click.option(
    "--ref-offset",
    type=int,
    default=15,
    show_default=True,
    help="P-site offset used when the matrix was built.",
)
@click.option(
    "--sample",
    "sample_names",
    multiple=True,
    help="Restrict to these sample name(s) (default: all).",
)
@click.option(
    "--site",
    type=click.Choice(["A", "P"]),
    default="A",
    show_default=True,
    help="Coverage site to query.",
)
@click.option(
    "--n-workers", type=int, default=None, help="Worker processes for partition scanning."
)
def score_matrix_cmd(
    events_dir,
    matrix_dir,
    partitions,
    store_dir,
    data_version,
    annotation_version,
    ref_offset,
    sample_names,
    site,
    n_workers,
):
    """Score extracted events against the sparse annotation-scale matrix.

    The matrix is sharded by read sequence; point --matrix-dir at the root and
    all partitions are scanned together (there is no valid single-partition use).
    """
    setup_logging()
    from .io.matrix import discover_partitions
    from .workflows import score_matrix_workflow

    if bool(matrix_dir) == bool(partitions):
        raise click.BadParameter(
            "provide exactly one of --matrix-dir (recommended) or --partitions"
        )
    part_dirs = (
        [str(p) for p in discover_partitions(matrix_dir)] if matrix_dir else list(partitions)
    )
    log_info(f"Scoring against {len(part_dirs)} matrix partitions")

    written = score_matrix_workflow(
        events_dir,
        part_dirs,
        store_dir,
        data_version=data_version,
        annotation_version=annotation_version,
        ref_offset=ref_offset,
        sample_names=list(sample_names) or None,
        n_workers=n_workers,
        site=site,
    )
    log_info(f"Scores written: {written or '(no events scored)'}")


@cli.command("score-bams")
@click.option("--events-dir", required=True, help="Events directory produced by extract-events.")
@click.option(
    "--bam",
    "bams",
    required=True,
    multiple=True,
    help="Genome-aligned BAM (repeatable, 1-20 samples).",
)
@click.option("--store-dir", required=True, help="Output fact_event_score store directory.")
@click.option(
    "--data-version", required=True, help="Data version label for the append-only score partition."
)
@click.option("--annotation-version", default="", help="Annotation version stamped into the store.")
@click.option(
    "--sample",
    "sample_names",
    multiple=True,
    help="Sample name(s) aligned with --bam order (default: BAM stems).",
)
@click.option(
    "--offset-method",
    type=click.Choice(["global", "file", "metagene"]),
    default="global",
    show_default=True,
    help="P-site offset calibration method.",
)
@click.option(
    "--global-offset",
    type=int,
    default=12,
    show_default=True,
    help="Fixed P-site offset for --offset-method global.",
)
@click.option("--offsets-file", help="CSV with read_length,offset for --offset-method file.")
@click.option(
    "--multimap",
    type=click.Choice(["unique"]),
    default="unique",
    show_default=True,
    help="Multimapper policy (unique only for now).",
)
@click.option(
    "--site",
    type=click.Choice(["A", "P"]),
    default="A",
    show_default=True,
    help="Coverage site to query.",
)
@click.option(
    "--transcriptome",
    is_flag=True,
    default=False,
    help="BAMs are transcriptome-aligned; project reads to genome via --annotation.",
)
@click.option(
    "--annotation",
    "-a",
    help="GTF annotation for transcriptome→genome projection (required with --transcriptome).",
)
def score_bams_cmd(
    events_dir,
    bams,
    store_dir,
    data_version,
    annotation_version,
    sample_names,
    offset_method,
    global_offset,
    offsets_file,
    multimap,
    site,
    transcriptome,
    annotation,
):
    """Score extracted events against 1-20 genome- or transcriptome-aligned BAMs."""
    setup_logging()
    from .model import OffsetParams
    from .workflows import score_bams_workflow

    exon_df = None
    if transcriptome:
        if not annotation:
            raise click.BadParameter("--transcriptome requires --annotation (GTF) for projection.")
        # Full exon structure (mRNA-origin coords) so UTR reads — 5'UTR uORFs in
        # particular — project to the correct genomic position, not just CDS.
        from .io.annotation import build_exon_blocks

        exon_df = build_exon_blocks(annotation)
    offsets = OffsetParams(
        method=offset_method, global_offset=global_offset, offsets_file=offsets_file
    )
    written = score_bams_workflow(
        events_dir,
        list(bams),
        store_dir,
        data_version=data_version,
        annotation_version=annotation_version,
        offsets=offsets,
        sample_names=list(sample_names) or None,
        multimap=multimap,
        site=site,
        transcriptome=transcriptome,
        exon_df=exon_df,
    )
    log_info(f"Scores written: {written or '(no events scored)'}")


@cli.command("report")
@click.option(
    "--store-dir",
    required=True,
    help="fact_event_score store directory produced by score-matrix/score-bams.",
)
@click.option(
    "--events-dir",
    required=True,
    help="Events directory produced by extract-events (must contain feature_event/).",
)
@click.option(
    "--out",
    "out_path",
    required=True,
    help="Output per-translon report Parquet (with consequentiality columns).",
)
@click.option(
    "--data-version", default=None, help="Restrict to this data_version partition (default: all)."
)
@click.option("--tier", default=None, help="Restrict to this tier (default: all).")
@click.option(
    "--min-tier-confidence",
    type=float,
    default=0.0,
    show_default=True,
    help="Confidence floor for the consequential flag.",
)
@click.option(
    "--min-expression-percentile",
    type=float,
    default=0.0,
    show_default=True,
    help="Expression-percentile floor for the consequential flag.",
)
@click.option(
    "--context-weight",
    type=float,
    default=1.0,
    show_default=True,
    help="Scale applied to the composite consequentiality score.",
)
def report_cmd(
    store_dir,
    events_dir,
    out_path,
    data_version,
    tier,
    min_tier_confidence,
    min_expression_percentile,
    context_weight,
):
    """Compose a per-translon report from scored events and apply the policy."""
    setup_logging()
    from .model import ConsequentialityPolicy
    from .workflows import report_workflow

    policy = ConsequentialityPolicy(
        min_tier_confidence=min_tier_confidence,
        min_expression_percentile=min_expression_percentile,
        context_weight=context_weight,
    )
    rep = report_workflow(
        store_dir,
        events_dir,
        out_path,
        data_version=data_version,
        tier=tier,
        policy=policy,
    )
    n_conseq = (
        int(rep["consequential"].sum()) if "consequential" in rep.columns and rep.height else 0
    )
    log_info(f"Report written: {out_path} ({rep.height} translons, {n_conseq} consequential)")


@cli.command("consequential")
@click.option(
    "--report", "report_path", required=True, help="Per-translon report Parquet/CSV to label."
)
@click.option(
    "--out", "out_path", required=True, help="Output Parquet with consequentiality columns."
)
@click.option(
    "--min-tier-confidence",
    type=float,
    default=0.0,
    show_default=True,
    help="Confidence floor for the consequential flag.",
)
@click.option(
    "--min-expression-percentile",
    type=float,
    default=0.0,
    show_default=True,
    help="Expression-percentile floor for the consequential flag.",
)
@click.option(
    "--context-weight",
    type=float,
    default=1.0,
    show_default=True,
    help="Scale applied to the composite consequentiality score.",
)
def consequential_cmd(
    report_path: str, out_path: str, min_tier_confidence, min_expression_percentile, context_weight
):
    """Apply the consequentiality policy to an existing per-translon report."""
    setup_logging()
    from .model import ConsequentialityPolicy
    from .workflows import consequential_workflow

    policy = ConsequentialityPolicy(
        min_tier_confidence=min_tier_confidence,
        min_expression_percentile=min_expression_percentile,
        context_weight=context_weight,
    )
    report = (
        pl.read_parquet(report_path)
        if report_path.endswith(".parquet")
        else pl.read_csv(report_path)
    )
    out = consequential_workflow(report, policy)
    out.write_parquet(out_path)
    log_info(f"Consequentiality labels written: {out_path}")


@cli.command("pipeline")
@click.option(
    "--out-dir", required=True, help="Output directory; writes events/, scores/, report.parquet."
)
# --- feature source for extract-events (exactly one) ---
@click.option("--gtf", "gtf_path", help="GTF/GFF; scores --feature-type features.")
@click.option("--feature-type", default="CDS", show_default=True, help="GTF feature type to score.")
@click.option("--bed12", "bed12_path", help="BED12 feature file.")
@click.option("--bigbed", "bigbed_path", help="bigBed (BED12) feature file.")
@click.option("--fasta", "fasta_path", help="FASTA for de-novo ORF finding.")
@click.option("--sqlite", "sqlite_path", help="Annotation sqlite (translons + translon_blocks).")
@click.option("--start-codons", default="ATG", show_default=True, help="Start codons for --fasta.")
@click.option(
    "--stop-codons", default="TAA,TAG,TGA", show_default=True, help="Stop codons for --fasta."
)
@click.option(
    "--matrix-dir",
    help="Sparse-matrix ROOT directory → matrix mode (ALL partitions scanned together).",
)
@click.option(
    "--bam", "bams", multiple=True, help="Genome/transcriptome BAM (repeatable) → BAM mode."
)
@click.option(
    "--chrom", "chroms", multiple=True, help="Restrict to chromosome(s) (repeatable; default: all)."
)
@click.option(
    "--data-version", default="run", show_default=True, help="Score-store partition label."
)
@click.option(
    "--annotation-version", default="", help="Annotation version stamped into events/scores."
)
@click.option(
    "--sample",
    "sample_names",
    multiple=True,
    help="Sample name(s) aligned with --bam/--matrix order.",
)
@click.option(
    "--offset-method",
    type=click.Choice(["global", "file", "metagene"]),
    default="global",
    show_default=True,
    help="P-site offset method (BAM mode).",
)
@click.option(
    "--global-offset",
    type=int,
    default=12,
    show_default=True,
    help="Fixed P-site offset for global method (BAM mode).",
)
@click.option("--offsets-file", help="CSV read_length,offset for --offset-method file (BAM mode).")
@click.option(
    "--site",
    type=click.Choice(["A", "P"]),
    default="A",
    show_default=True,
    help="Coverage site to query.",
)
@click.option(
    "--transcriptome",
    is_flag=True,
    default=False,
    help="BAMs are transcriptome-aligned; project via --annotation (BAM mode).",
)
@click.option(
    "--annotation", "-a", help="GTF for transcriptome→genome projection (with --transcriptome)."
)
@click.option(
    "--min-tier-confidence",
    type=float,
    default=0.0,
    show_default=True,
    help="Consequential confidence floor.",
)
@click.option(
    "--min-expression-percentile",
    type=float,
    default=0.0,
    show_default=True,
    help="Consequential expression-percentile floor.",
)
def pipeline_cmd(
    out_dir,
    gtf_path,
    feature_type,
    bed12_path,
    bigbed_path,
    fasta_path,
    sqlite_path,
    start_codons,
    stop_codons,
    matrix_dir,
    bams,
    chroms,
    data_version,
    annotation_version,
    sample_names,
    offset_method,
    global_offset,
    offsets_file,
    site,
    transcriptome,
    annotation,
    min_tier_confidence,
    min_expression_percentile,
):
    """One-shot event-scoring run: extract-events → score (matrix or BAMs) → report.

    Feature source (exactly one): --gtf / --bed12 / --bigbed / --fasta / --sqlite.
    Coverage source (exactly one): --matrix-dir (whole sharded matrix) or --bam
    (1-20 BAMs). Equivalent to running extract-events, score-matrix/score-bams and
    report in sequence, into one output directory.
    """
    setup_logging()
    from .io.matrix import discover_partitions
    from .model import ConsequentialityPolicy, OffsetParams
    from .workflows import pipeline_workflow

    if bool(matrix_dir) == bool(bams):
        raise click.BadParameter(
            "provide exactly one of --matrix-dir (matrix mode) or --bam (BAM mode)"
        )
    n_src = sum(bool(x) for x in (gtf_path, bed12_path, bigbed_path, fasta_path, sqlite_path))
    if n_src != 1:
        raise click.BadParameter(
            "provide exactly one feature source: --gtf | --bed12 | --bigbed | --fasta | --sqlite"
        )
    partition_dirs = [str(p) for p in discover_partitions(matrix_dir)] if matrix_dir else None
    if partition_dirs:
        log_info(f"Scoring against {len(partition_dirs)} matrix partitions")
    exon_df = None
    if transcriptome:
        if not annotation:
            raise click.BadParameter("--transcriptome requires --annotation (GTF).")
        from .io.annotation import build_exon_blocks

        exon_df = build_exon_blocks(annotation)
    policy = ConsequentialityPolicy(
        min_tier_confidence=min_tier_confidence,
        min_expression_percentile=min_expression_percentile,
    )
    paths = pipeline_workflow(
        out_dir,
        partition_dirs=partition_dirs,
        bams=list(bams) or None,
        data_version=data_version,
        chroms=list(chroms) or None,
        annotation_version=annotation_version,
        offsets=OffsetParams(
            method=offset_method, global_offset=global_offset, offsets_file=offsets_file
        ),
        sample_names=list(sample_names) or None,
        site=site,
        transcriptome=transcriptome,
        exon_df=exon_df,
        policy=policy,
        # feature source forwarded to extract-events
        sqlite_path=sqlite_path or None,
        gtf_path=gtf_path or None,
        feature_type=feature_type,
        bed12_path=bed12_path or None,
        bigbed_path=bigbed_path or None,
        fasta_path=fasta_path or None,
        start_codons=[c.strip() for c in start_codons.split(",") if c.strip()],
        stop_codons=[c.strip() for c in stop_codons.split(",") if c.strip()],
    )
    log_info(f"Pipeline complete: {paths['report']}")


if __name__ == "__main__":
    cli()

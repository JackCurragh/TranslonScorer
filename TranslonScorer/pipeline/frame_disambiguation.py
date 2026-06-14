from __future__ import annotations

from typing import Any, Optional

import polars as pl

from ..file_handlers import bam as bam_handlers
from ..utils.io import write_csv_safe, write_parquet_safe
from ..utils.logging import log_info
from .transcript_coords import cds_to_transcript_space


def _load_table(path: str) -> pl.DataFrame:
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    if path.endswith(".tsv"):
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def _ensure_read_key(df: pl.DataFrame) -> pl.DataFrame:
    if "read_key" in df.columns:
        return df
    if "qname" in df.columns:
        return df.with_columns(pl.col("qname").cast(pl.Utf8).alias("read_key"))
    if "read_row_id" in df.columns:
        return df.with_columns(pl.concat_str([pl.lit("row"), pl.col("read_row_id").cast(pl.Utf8)], separator=":").alias("read_key"))
    coord_cols = [c for c in ["chr", "start", "stop", "length", "strand"] if c in df.columns]
    if coord_cols:
        return df.with_columns(pl.concat_str([pl.col(c).cast(pl.Utf8) for c in coord_cols], separator=":").alias("read_key"))
    return (
        df.with_row_index("_read_row")
        .with_columns(pl.concat_str([pl.lit("row"), pl.col("_read_row").cast(pl.Utf8)], separator=":").alias("read_key"))
        .drop("_read_row")
    )


def _position_column(df: pl.DataFrame) -> str:
    for col in ("tran_start_bam", "tran_start", "pos", "start_pos_tran"):
        if col in df.columns:
            return col
    raise ValueError("Candidate table needs one transcript-position column")


def candidate_frame_table(candidates: pl.DataFrame, cds_tran_df: pl.DataFrame) -> pl.DataFrame:
    """Annotate candidate read-to-transcript assignments with CDS-relative frame."""
    if candidates.is_empty():
        return candidates
    pos_col = _position_column(candidates)
    df = _ensure_read_key(candidates)
    if "count" not in df.columns:
        df = df.with_columns(pl.lit(1.0).alias("count"))

    cds = cds_tran_df.select(["tran_id", pl.col("start").alias("cds_start")]).unique(subset=["tran_id"])
    df = df.join(cds, on="tran_id", how="left")
    return (
        df.with_columns(
            pl.when(pl.col("cds_start").is_not_null())
            .then((pl.col(pos_col).cast(pl.Int64) - pl.col("cds_start").cast(pl.Int64)) % 3)
            .otherwise(pl.col(pos_col).cast(pl.Int64) % 3)
            .cast(pl.Int64)
            .alias("candidate_frame")
        )
    )


def summarize_frame_disambiguation(candidate_frames: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Summarize how often ambiguous candidate assignments disagree in frame."""
    if candidate_frames.is_empty():
        empty_classes = pl.DataFrame()
        summary = pl.DataFrame([{
            "n_read_keys": 0,
            "n_ambiguous_read_keys": 0,
            "n_frame_discordant_read_keys": 0,
            "ambiguous_fraction": 0.0,
            "frame_discordant_fraction_of_ambiguous": None,
            "weighted_ambiguous_count": 0.0,
            "weighted_frame_discordant_count": 0.0,
            "weighted_frame_discordant_fraction_of_ambiguous": None,
        }])
        return empty_classes, summary

    grouped = (
        candidate_frames.group_by("read_key")
        .agg(
            pl.len().alias("candidate_rows"),
            pl.col("tran_id").n_unique().alias("n_candidate_transcripts"),
            pl.col("candidate_frame").n_unique().alias("n_candidate_frames"),
            pl.col("count").first().cast(pl.Float64).alias("read_weight"),
            pl.col("tran_id").unique().alias("candidate_transcripts"),
            pl.col("candidate_frame").unique().alias("candidate_frames"),
        )
        .with_columns(
            (pl.col("n_candidate_transcripts") > 1).alias("is_ambiguous"),
            ((pl.col("n_candidate_transcripts") > 1) & (pl.col("n_candidate_frames") > 1)).alias("is_frame_discordant"),
        )
    )

    n_read_keys = grouped.height
    ambiguous = grouped.filter(pl.col("is_ambiguous"))
    discordant = grouped.filter(pl.col("is_frame_discordant"))
    weighted_amb = float(ambiguous["read_weight"].sum() or 0.0) if ambiguous.height else 0.0
    weighted_disc = float(discordant["read_weight"].sum() or 0.0) if discordant.height else 0.0

    summary = pl.DataFrame([{
        "n_read_keys": n_read_keys,
        "n_ambiguous_read_keys": ambiguous.height,
        "n_frame_discordant_read_keys": discordant.height,
        "ambiguous_fraction": (ambiguous.height / n_read_keys) if n_read_keys else 0.0,
        "frame_discordant_fraction_of_ambiguous": (discordant.height / ambiguous.height) if ambiguous.height else None,
        "weighted_ambiguous_count": weighted_amb,
        "weighted_frame_discordant_count": weighted_disc,
        "weighted_frame_discordant_fraction_of_ambiguous": (weighted_disc / weighted_amb) if weighted_amb > 0 else None,
    }])
    return grouped, summary


def frame_disambiguation_from_candidates(
    *,
    candidates_path: str,
    cds_path: str,
    out_prefix: str,
) -> dict[str, str]:
    candidates = _load_table(candidates_path)
    cds_tran = _load_table(cds_path)
    frames = candidate_frame_table(candidates, cds_tran)
    classes, summary = summarize_frame_disambiguation(frames)
    paths = {
        "candidate_frames": f"{out_prefix}.candidate_frames.parquet",
        "classes": f"{out_prefix}.read_classes.parquet",
        "summary": f"{out_prefix}.summary.csv",
    }
    write_parquet_safe(frames, paths["candidate_frames"])
    write_parquet_safe(classes, paths["classes"])
    write_csv_safe(summary, paths["summary"])
    return paths


def frame_disambiguation_from_bam(
    *,
    bam_path: str,
    annotation_path: str,
    out_prefix: str,
    collapsed: bool = False,
    count_from: Optional[str] = None,
    count_pattern: Optional[str] = None,
    count_tag: Optional[str] = None,
) -> dict[str, str]:
    """Run the cheap Decision 2 frame-disambiguation analysis from a BAM."""
    log_info("Loading annotation for frame-disambiguation analysis")
    cds_df, exon_df = bam_handlers.getexons_and_cds(annotation_path)
    cds_tran = cds_to_transcript_space(cds_df, exon_df)

    log_info("Reading BAM with read identity preserved")
    reads = bam_handlers.readbam(
        bam_path,
        collapsed=collapsed,
        count_from=count_from,
        count_pattern=count_pattern,
        count_tag=count_tag,
        unique=True,
        include_qname=True,
    ).with_row_index("read_row_id")

    log_info("Projecting reads to all compatible transcript coordinates")
    mapped = bam_handlers.bamtranscript(reads, exon_df)
    frames = candidate_frame_table(mapped, cds_tran)
    classes, summary = summarize_frame_disambiguation(frames)

    paths = {
        "cds_transcript_space": f"{out_prefix}.cds_transcript_space.parquet",
        "candidate_frames": f"{out_prefix}.candidate_frames.parquet",
        "classes": f"{out_prefix}.read_classes.parquet",
        "summary": f"{out_prefix}.summary.csv",
    }
    write_parquet_safe(cds_tran, paths["cds_transcript_space"])
    write_parquet_safe(frames, paths["candidate_frames"])
    write_parquet_safe(classes, paths["classes"])
    write_csv_safe(summary, paths["summary"])
    return paths


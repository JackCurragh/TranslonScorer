from __future__ import annotations

import json
import os
from datetime import datetime
from typing import Any, Dict, Optional, Tuple

import polars as pl

from ..coverage.locus_features import build_locus_features  # reuse existing parsers
from ..file_handlers import bam as bam_handlers
from ..utils import log_info


def _write_json(path: str, obj: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


def _detect_chr_style(exon_df: pl.DataFrame) -> str:
    try:
        chrs = set(exon_df.get_column("chr").unique().cast(pl.Utf8))
        return "chr" if any(str(c).startswith("chr") for c in chrs) else "nochr"
    except Exception:
        return "unknown"


def _make_transcripts_table(exon_df: pl.DataFrame, fmap_df: Optional[pl.DataFrame]) -> pl.DataFrame:
    # length_tran from last tran_stop per transcript
    # gene_id from feature_map if available
    lens = (
        exon_df.with_columns(
            pl.struct(["tran_start", "tran_stop"])
            .map_elements(
                lambda s: int(s["tran_stop"][-1]) if len(s["tran_stop"]) > 0 else 0,
                return_dtype=pl.Int64,
            )
            .alias("length_tran")
        )
        .select(["tran_id", "chr", "strand", "length_tran"])
        .unique(subset=["tran_id"])
    )
    # Join locus/gene id if present in feature_map; tolerate either 'transcript_id' or 'tran_id'
    if fmap_df is not None and "locus_id" in fmap_df.columns:
        tx_col = (
            "tran_id"
            if "tran_id" in fmap_df.columns
            else ("transcript_id" if "transcript_id" in fmap_df.columns else None)
        )
        if tx_col is not None:
            lut = fmap_df.select([tx_col, "locus_id"]).unique()
            # Normalize key name to 'tran_id' to avoid extra columns
            if tx_col != "tran_id":
                lut = lut.rename({tx_col: "tran_id"})
            lens = lens.join(lut, on="tran_id", how="left").rename({"locus_id": "gene_id"})
    return lens


def build_annotation_bundle(
    gtf_path: str, out_dir: str, *, progress: bool = True
) -> Tuple[str, Dict[str, str]]:
    """Build an annotation bundle directory from a GTF.

    Writes exons.parquet, cds.parquet, features.parquet, feature_map.parquet,
    transcripts.parquet, loci.bed, and manifest.json under out_dir.

    Returns (out_dir, paths_dict).
    """
    os.makedirs(out_dir, exist_ok=True)

    # Features and feature map (heavy gene-structured tables)
    log_info(f"Building features bundle in {out_dir} …")
    feats, fmap = build_locus_features(gtf_path, progress=progress)
    feats_path = os.path.join(out_dir, "features.parquet")
    fmap_path = os.path.join(out_dir, "feature_map.parquet")
    feats.write_parquet(feats_path)
    fmap.write_parquet(fmap_path)

    # Exon/CDS transcript models (projection + CDS offsets)
    cds_df, exon_df = bam_handlers.getexons_and_cds(gtf_path)
    exons_path = os.path.join(out_dir, "exons.parquet")
    cds_path = os.path.join(out_dir, "cds.parquet")
    exon_df.write_parquet(exons_path)
    cds_df.write_parquet(cds_path)

    # Transcripts summary table (tran_id → chr,strand,length_tran[,gene_id])
    tx_path = os.path.join(out_dir, "transcripts.parquet")
    tx = _make_transcripts_table(exon_df, fmap)
    tx.write_parquet(tx_path)

    # Loci BED (optional): min..max span per locus_id from features (exon chunks)
    loci_bed = os.path.join(out_dir, "loci.bed")
    try:
        loci = (
            feats.filter(pl.col("feature_type") == "exon_chunk")
            .group_by(["locus_id", "chr", "strand"])
            .agg([pl.col("start").min().alias("start"), pl.col("end").max().alias("end")])
            .select(["chr", "start", "end", "locus_id"])
            .sort(["chr", "start"])
        )
        loci.write_csv(loci_bed, separator="\t", include_header=False)
    except Exception:
        loci_bed = ""

    # Manifest
    chr_style = _detect_chr_style(exon_df)
    manifest = {
        "bundle_version": 1,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "gtf_path": os.path.abspath(gtf_path),
        "chr_style": chr_style,
        "schemas": {
            "exons": list(exon_df.columns),
            "cds": list(cds_df.columns),
            "features": list(feats.columns),
            "feature_map": list(fmap.columns),
            "transcripts": list(tx.columns),
        },
    }
    manifest_path = os.path.join(out_dir, "manifest.json")
    _write_json(manifest_path, manifest)

    paths = {
        "exons": exons_path,
        "cds": cds_path,
        "features": feats_path,
        "feature_map": fmap_path,
        "transcripts": tx_path,
        "loci_bed": loci_bed,
        "manifest": manifest_path,
    }
    log_info("Annotation bundle written")
    return out_dir, paths


def load_annotation_bundle(
    dir_path: str,
) -> Tuple[
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
    Optional[str],
    Dict[str, Any],
]:
    """Load exons, cds, features, feature_map, transcripts, loci_bed and manifest from bundle dir."""
    exons = pl.read_parquet(os.path.join(dir_path, "exons.parquet"))
    cds = pl.read_parquet(os.path.join(dir_path, "cds.parquet"))
    feats = pl.read_parquet(os.path.join(dir_path, "features.parquet"))
    fmap = pl.read_parquet(os.path.join(dir_path, "feature_map.parquet"))
    tx = pl.read_parquet(os.path.join(dir_path, "transcripts.parquet"))
    loci_bed = os.path.join(dir_path, "loci.bed")
    if not os.path.exists(loci_bed):
        loci_bed = None
    manifest_path = os.path.join(dir_path, "manifest.json")
    manifest = {}
    if os.path.exists(manifest_path):
        with open(manifest_path, "r") as fh:
            manifest = json.load(fh)
    return exons, cds, feats, fmap, tx, loci_bed, manifest

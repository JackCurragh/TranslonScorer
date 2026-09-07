"""RDG-Flux v1 export: position-level frame-posterior tables.

Lives with the frame algorithms because that is what it exports — it is built
on ``frame_support.build_frame_support``, not on any ORF-composite code.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import polars as pl

from TranslonScorer import __version__

from ..coverage.transcript_coords import cds_to_transcript_space
from ..file_handlers import bam as bam_handlers
from ..frame_support import build_frame_support
from ..model import FrameSupportParams
from ..utils.io import ensure_dir_for_file, write_parquet_safe

MAX_FRAME_ENTROPY = float(1.584962500721156)

RDG_FLUX_COLUMNS = [
    "sample_id",
    "transcript_id",
    "pos",
    "count",
    "p_frame0",
    "p_frame1",
    "p_frame2",
    "p_background",
    "frame_entropy",
    "effective_depth",
    "p_translated",
    "local_periodicity_score",
]


def _load_table(path: str) -> pl.DataFrame:
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    if path.endswith(".tsv"):
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def _scan_parquet_with_schema(path: str) -> tuple[pl.LazyFrame, set[str]]:
    lf = pl.scan_parquet(path)
    return lf, set(lf.collect_schema().names())


def _normalise_profile_schema(profiles: pl.DataFrame, sample_id: str | None) -> pl.DataFrame:
    df = profiles
    if "transcript_id" in df.columns and "tran_id" not in df.columns:
        df = df.rename({"transcript_id": "tran_id"})
    if "sample_id" not in df.columns:
        if "sample" in df.columns:
            df = df.rename({"sample": "sample_id"})
        elif sample_id:
            df = df.with_columns(pl.lit(sample_id).alias("sample_id"))
        else:
            raise ValueError(
                "profiles must contain sample_id/sample, or --sample-id must be provided"
            )
    required = {"sample_id", "tran_id", "pos", "count"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"profiles missing required columns: {sorted(missing)}")
    return df


def _normalise_profile_lazy(
    profiles_path: str, sample_id: str | None
) -> tuple[pl.LazyFrame, list[str], set[str]]:
    lf, cols = _scan_parquet_with_schema(profiles_path)
    if "transcript_id" in cols and "tran_id" not in cols:
        lf = lf.rename({"transcript_id": "tran_id"})
        cols = (cols - {"transcript_id"}) | {"tran_id"}
    if "sample_id" not in cols:
        if "sample" in cols:
            lf = lf.rename({"sample": "sample_id"})
            cols = (cols - {"sample"}) | {"sample_id"}
        elif sample_id:
            lf = lf.with_columns(pl.lit(sample_id).alias("sample_id"))
            cols = cols | {"sample_id"}
        else:
            raise ValueError(
                "profiles must contain sample_id/sample, or --sample-id must be provided"
            )
    required = {"sample_id", "tran_id", "pos", "count"}
    missing = required - cols
    if missing:
        raise ValueError(f"profiles missing required columns: {sorted(missing)}")
    sample_ids = sorted(
        str(x)
        for x in lf.select(pl.col("sample_id").unique()).collect().get_column("sample_id").to_list()
    )
    return lf, sample_ids, cols


def _normalise_frame_schema(frame_support: pl.DataFrame, sample_ids: list[str]) -> pl.DataFrame:
    df = frame_support
    if "transcript_id" in df.columns and "tran_id" not in df.columns:
        df = df.rename({"transcript_id": "tran_id"})
    if "sample_id" not in df.columns:
        if "sample" in df.columns:
            df = df.rename({"sample": "sample_id"})
        elif len(sample_ids) == 1:
            df = df.with_columns(pl.lit(sample_ids[0]).alias("sample_id"))
        else:
            raise ValueError(
                "frame support lacks sample_id/sample but profiles contain multiple samples"
            )
    required = {"sample_id", "tran_id", "codon", "p0", "p1", "p2"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"frame support missing required columns: {sorted(missing)}")
    return df


def _normalise_frame_lazy(
    frame_support_path: str, sample_ids: list[str]
) -> tuple[pl.LazyFrame, set[str]]:
    lf, cols = _scan_parquet_with_schema(frame_support_path)
    if "transcript_id" in cols and "tran_id" not in cols:
        lf = lf.rename({"transcript_id": "tran_id"})
        cols = (cols - {"transcript_id"}) | {"tran_id"}
    if "sample_id" not in cols:
        if "sample" in cols:
            lf = lf.rename({"sample": "sample_id"})
            cols = (cols - {"sample"}) | {"sample_id"}
        elif len(sample_ids) == 1:
            lf = lf.with_columns(pl.lit(sample_ids[0]).alias("sample_id"))
            cols = cols | {"sample_id"}
        else:
            raise ValueError(
                "frame support lacks sample_id/sample but profiles contain multiple samples"
            )
    required = {"sample_id", "tran_id", "codon", "p0", "p1", "p2"}
    missing = required - cols
    if missing:
        raise ValueError(f"frame support missing required columns: {sorted(missing)}")
    return lf, cols


def _entropy_expr(p0: pl.Expr, p1: pl.Expr, p2: pl.Expr) -> pl.Expr:
    terms = []
    for p in (p0, p1, p2):
        terms.append(pl.when(p > 0).then(p * p.log(2)).otherwise(0.0))
    return (-(terms[0] + terms[1] + terms[2])).alias("frame_entropy")


def _collapse_profiles(profiles: pl.DataFrame, preserve_read_length: bool) -> pl.DataFrame:
    keep_length = preserve_read_length and "length" in profiles.columns
    groups = ["sample_id", "tran_id", "pos"] + (["length"] if keep_length else [])
    return (
        profiles.with_columns(
            [
                pl.col("pos").cast(pl.Int64),
                pl.col("count").cast(pl.Float64),
                (pl.col("pos").cast(pl.Int64) // 3).alias("codon"),
            ]
        )
        .group_by(groups + ["codon"])
        .agg(pl.col("count").sum().alias("count"))
        .sort(groups)
    )


def _collapse_profiles_lazy(
    profiles: pl.LazyFrame, cols: set[str], preserve_read_length: bool
) -> pl.LazyFrame:
    keep_length = preserve_read_length and "length" in cols
    groups = ["sample_id", "tran_id", "pos"] + (["length"] if keep_length else [])
    return (
        profiles.with_columns(
            [
                pl.col("pos").cast(pl.Int64),
                pl.col("count").cast(pl.Float64),
                (pl.col("pos").cast(pl.Int64) // 3).alias("codon"),
            ]
        )
        .group_by(groups + ["codon"])
        .agg(pl.col("count").sum().alias("count"))
    )


def _collapse_frame_support(
    frame_support: pl.DataFrame, preserve_read_length: bool
) -> pl.DataFrame:
    keep_length = preserve_read_length and "length" in frame_support.columns
    groups = ["sample_id", "tran_id", "codon"] + (["length"] if keep_length else [])
    cols = set(frame_support.columns)

    if {"adjusted_f0", "adjusted_f1", "adjusted_f2"}.issubset(cols):
        agg_exprs = [
            pl.col("adjusted_f0").sum().alias("adjusted_f0"),
            pl.col("adjusted_f1").sum().alias("adjusted_f1"),
            pl.col("adjusted_f2").sum().alias("adjusted_f2"),
            (
                pl.col("total_count").sum().alias("frame_total_count")
                if "total_count" in cols
                else pl.len().cast(pl.Float64).alias("frame_total_count")
            ),
        ]
        if "support_evidence" in cols:
            agg_exprs.append(pl.col("support_evidence").mean().alias("support_evidence"))
        if "frame_periodicity_score" in cols:
            agg_exprs.append(
                pl.col("frame_periodicity_score").mean().alias("frame_periodicity_score")
            )
        collapsed = (
            frame_support.group_by(groups)
            .agg(agg_exprs)
            .with_columns(
                (pl.col("adjusted_f0") + pl.col("adjusted_f1") + pl.col("adjusted_f2")).alias(
                    "_adjusted_total"
                )
            )
            .with_columns(
                [
                    pl.when(pl.col("_adjusted_total") > 0)
                    .then(pl.col("adjusted_f0") / pl.col("_adjusted_total"))
                    .otherwise(1.0 / 3.0)
                    .alias("p0"),
                    pl.when(pl.col("_adjusted_total") > 0)
                    .then(pl.col("adjusted_f1") / pl.col("_adjusted_total"))
                    .otherwise(1.0 / 3.0)
                    .alias("p1"),
                    pl.when(pl.col("_adjusted_total") > 0)
                    .then(pl.col("adjusted_f2") / pl.col("_adjusted_total"))
                    .otherwise(1.0 / 3.0)
                    .alias("p2"),
                ]
            )
            .drop("_adjusted_total")
        )
    else:
        weight_col = "_frame_weight"
        weighted = frame_support.with_columns(
            (
                pl.col("total_count").cast(pl.Float64) if "total_count" in cols else pl.lit(1.0)
            ).alias(weight_col)
        )
        agg_exprs = [
            (pl.col("p0") * pl.col(weight_col)).sum().alias("_p0w"),
            (pl.col("p1") * pl.col(weight_col)).sum().alias("_p1w"),
            (pl.col("p2") * pl.col(weight_col)).sum().alias("_p2w"),
            pl.col(weight_col).sum().alias("frame_total_count"),
        ]
        if "support_evidence" in cols:
            agg_exprs.append(pl.col("support_evidence").mean().alias("support_evidence"))
        if "frame_periodicity_score" in cols:
            agg_exprs.append(
                pl.col("frame_periodicity_score").mean().alias("frame_periodicity_score")
            )
        collapsed = (
            weighted.group_by(groups)
            .agg(agg_exprs)
            .with_columns(
                [
                    pl.when(pl.col("frame_total_count") > 0)
                    .then(pl.col("_p0w") / pl.col("frame_total_count"))
                    .otherwise(1.0 / 3.0)
                    .alias("p0"),
                    pl.when(pl.col("frame_total_count") > 0)
                    .then(pl.col("_p1w") / pl.col("frame_total_count"))
                    .otherwise(1.0 / 3.0)
                    .alias("p1"),
                    pl.when(pl.col("frame_total_count") > 0)
                    .then(pl.col("_p2w") / pl.col("frame_total_count"))
                    .otherwise(1.0 / 3.0)
                    .alias("p2"),
                ]
            )
            .drop(["_p0w", "_p1w", "_p2w"])
        )

    return collapsed.with_columns(_entropy_expr(pl.col("p0"), pl.col("p1"), pl.col("p2"))).sort(
        groups
    )


def _collapse_frame_support_lazy(
    frame_support: pl.LazyFrame, cols: set[str], preserve_read_length: bool
) -> pl.LazyFrame:
    keep_length = preserve_read_length and "length" in cols
    groups = ["sample_id", "tran_id", "codon"] + (["length"] if keep_length else [])

    if {"adjusted_f0", "adjusted_f1", "adjusted_f2"}.issubset(cols):
        frame_total_expr = (
            pl.col("total_count").sum().alias("frame_total_count")
            if "total_count" in cols
            else pl.len().cast(pl.Float64).alias("frame_total_count")
        )
        agg_exprs = [
            pl.col("adjusted_f0").sum().alias("adjusted_f0"),
            pl.col("adjusted_f1").sum().alias("adjusted_f1"),
            pl.col("adjusted_f2").sum().alias("adjusted_f2"),
            frame_total_expr,
        ]
        if "support_evidence" in cols:
            agg_exprs.append(pl.col("support_evidence").mean().alias("support_evidence"))
        if "frame_periodicity_score" in cols:
            agg_exprs.append(
                pl.col("frame_periodicity_score").mean().alias("frame_periodicity_score")
            )
        collapsed = (
            frame_support.group_by(groups)
            .agg(agg_exprs)
            .with_columns(
                (pl.col("adjusted_f0") + pl.col("adjusted_f1") + pl.col("adjusted_f2")).alias(
                    "_adjusted_total"
                )
            )
            .with_columns(
                [
                    pl.when(pl.col("_adjusted_total") > 0)
                    .then(pl.col("adjusted_f0") / pl.col("_adjusted_total"))
                    .otherwise(1.0 / 3.0)
                    .alias("p0"),
                    pl.when(pl.col("_adjusted_total") > 0)
                    .then(pl.col("adjusted_f1") / pl.col("_adjusted_total"))
                    .otherwise(1.0 / 3.0)
                    .alias("p1"),
                    pl.when(pl.col("_adjusted_total") > 0)
                    .then(pl.col("adjusted_f2") / pl.col("_adjusted_total"))
                    .otherwise(1.0 / 3.0)
                    .alias("p2"),
                ]
            )
            .drop("_adjusted_total")
        )
    else:
        weight_col = "_frame_weight"
        weighted = frame_support.with_columns(
            (
                pl.col("total_count").cast(pl.Float64) if "total_count" in cols else pl.lit(1.0)
            ).alias(weight_col)
        )
        agg_exprs = [
            (pl.col("p0") * pl.col(weight_col)).sum().alias("_p0w"),
            (pl.col("p1") * pl.col(weight_col)).sum().alias("_p1w"),
            (pl.col("p2") * pl.col(weight_col)).sum().alias("_p2w"),
            pl.col(weight_col).sum().alias("frame_total_count"),
        ]
        if "support_evidence" in cols:
            agg_exprs.append(pl.col("support_evidence").mean().alias("support_evidence"))
        if "frame_periodicity_score" in cols:
            agg_exprs.append(
                pl.col("frame_periodicity_score").mean().alias("frame_periodicity_score")
            )
        collapsed = (
            weighted.group_by(groups)
            .agg(agg_exprs)
            .with_columns(
                [
                    pl.when(pl.col("frame_total_count") > 0)
                    .then(pl.col("_p0w") / pl.col("frame_total_count"))
                    .otherwise(1.0 / 3.0)
                    .alias("p0"),
                    pl.when(pl.col("frame_total_count") > 0)
                    .then(pl.col("_p1w") / pl.col("frame_total_count"))
                    .otherwise(1.0 / 3.0)
                    .alias("p1"),
                    pl.when(pl.col("frame_total_count") > 0)
                    .then(pl.col("_p2w") / pl.col("frame_total_count"))
                    .otherwise(1.0 / 3.0)
                    .alias("p2"),
                ]
            )
            .drop(["_p0w", "_p1w", "_p2w"])
        )

    return collapsed.with_columns(_entropy_expr(pl.col("p0"), pl.col("p1"), pl.col("p2")))


def _rdg_flux_position_lazy(
    profiles: pl.LazyFrame,
    profile_cols: set[str],
    frame_support: pl.LazyFrame,
    frame_cols: set[str],
    *,
    background_probability: float = 0.0,
    preserve_read_length: bool = False,
) -> pl.LazyFrame:
    if not 0.0 <= background_probability < 1.0:
        raise ValueError("background_probability must be >= 0 and < 1")

    prof = _collapse_profiles_lazy(profiles, profile_cols, preserve_read_length)
    fs = _collapse_frame_support_lazy(frame_support, frame_cols, preserve_read_length)
    has_support_evidence = "support_evidence" in frame_cols
    has_periodicity_score = "frame_periodicity_score" in frame_cols

    join_keys = ["sample_id", "tran_id", "codon"]
    if preserve_read_length and "length" in profile_cols and "length" in frame_cols:
        join_keys.append("length")

    background_floor = float(background_probability)
    translated_floor = 1.0 - background_floor
    support_expr = (
        pl.col("support_evidence").fill_null(0.0).clip(0.0, 1.0)
        if has_support_evidence
        else (pl.col("frame_total_count") > 0).cast(pl.Float64)
    )
    periodicity_expr = (
        pl.col("frame_periodicity_score").fill_null(0.0).clip(0.0, 1.0)
        if has_periodicity_score
        else (1.0 - (pl.col("frame_entropy") / MAX_FRAME_ENTROPY)).clip(0.0, 1.0)
    )
    out = (
        prof.join(fs, on=join_keys, how="left")
        .with_columns(
            [
                pl.col("p0").fill_null(1.0 / 3.0),
                pl.col("p1").fill_null(1.0 / 3.0),
                pl.col("p2").fill_null(1.0 / 3.0),
                pl.col("frame_entropy").fill_null(MAX_FRAME_ENTROPY),
                pl.col("frame_total_count").fill_null(0.0),
            ]
        )
        .with_columns(
            [
                support_expr.alias("_support_evidence_for_export"),
                periodicity_expr.alias("local_periodicity_score"),
            ]
        )
        .with_columns(
            (pl.lit(translated_floor) * pl.col("_support_evidence_for_export")).alias(
                "p_translated"
            )
        )
        .with_columns(
            [
                (pl.col("p0") * pl.col("p_translated")).alias("p_frame0"),
                (pl.col("p1") * pl.col("p_translated")).alias("p_frame1"),
                (pl.col("p2") * pl.col("p_translated")).alias("p_frame2"),
                (1.0 - pl.col("p_translated")).alias("p_background"),
                (pl.col("count") * pl.col("p_translated")).alias("effective_depth"),
            ]
        )
        .rename({"tran_id": "transcript_id"})
    )

    columns = list(RDG_FLUX_COLUMNS)
    if preserve_read_length and "length" in profile_cols:
        out = out.with_columns(pl.col("length").alias("read_length_bin"))
        columns.append("read_length_bin")

    return out.select(columns).sort(
        ["sample_id", "transcript_id", "pos"]
        + (["read_length_bin"] if "read_length_bin" in columns else [])
    )


def rdg_flux_position_table(
    profiles: pl.DataFrame,
    frame_support: pl.DataFrame,
    *,
    sample_id: str | None = None,
    background_probability: float = 0.0,
    preserve_read_length: bool = False,
) -> pl.DataFrame:
    """Build the RDG-Flux v1 position-level frame-posterior table."""
    if not 0.0 <= background_probability < 1.0:
        raise ValueError("background_probability must be >= 0 and < 1")

    profiles_n = _normalise_profile_schema(profiles, sample_id)
    sample_ids = sorted(str(x) for x in profiles_n.get_column("sample_id").unique().to_list())
    frame_n = _normalise_frame_schema(frame_support, sample_ids)

    prof = _collapse_profiles(profiles_n, preserve_read_length)
    fs = _collapse_frame_support(frame_n, preserve_read_length)
    fs_cols = set(fs.columns)
    has_support_evidence = "support_evidence" in fs_cols
    has_periodicity_score = "frame_periodicity_score" in fs_cols

    join_keys = ["sample_id", "tran_id", "codon"]
    if preserve_read_length and "length" in prof.columns and "length" in fs.columns:
        join_keys.append("length")

    background_floor = float(background_probability)
    translated_floor = 1.0 - background_floor
    support_expr = (
        pl.col("support_evidence").fill_null(0.0).clip(0.0, 1.0)
        if has_support_evidence
        else (pl.col("frame_total_count") > 0).cast(pl.Float64)
    )
    periodicity_expr = (
        pl.col("frame_periodicity_score").fill_null(0.0).clip(0.0, 1.0)
        if has_periodicity_score
        else (1.0 - (pl.col("frame_entropy") / MAX_FRAME_ENTROPY)).clip(0.0, 1.0)
    )
    joined = prof.join(fs, on=join_keys, how="left")

    out = (
        joined.with_columns(
            [
                pl.col("p0").fill_null(1.0 / 3.0),
                pl.col("p1").fill_null(1.0 / 3.0),
                pl.col("p2").fill_null(1.0 / 3.0),
                pl.col("frame_entropy").fill_null(MAX_FRAME_ENTROPY),
                pl.col("frame_total_count").fill_null(0.0),
            ]
        )
        .with_columns(
            [
                support_expr.alias("_support_evidence_for_export"),
                periodicity_expr.alias("local_periodicity_score"),
            ]
        )
        .with_columns(
            (pl.lit(translated_floor) * pl.col("_support_evidence_for_export")).alias(
                "p_translated"
            )
        )
        .with_columns(
            [
                (pl.col("p0") * pl.col("p_translated")).alias("p_frame0"),
                (pl.col("p1") * pl.col("p_translated")).alias("p_frame1"),
                (pl.col("p2") * pl.col("p_translated")).alias("p_frame2"),
                (1.0 - pl.col("p_translated")).alias("p_background"),
                (pl.col("count") * pl.col("p_translated")).alias("effective_depth"),
            ]
        )
        .rename({"tran_id": "transcript_id"})
    )

    columns = list(RDG_FLUX_COLUMNS)
    if preserve_read_length and "length" in out.columns:
        out = out.with_columns(pl.col("length").alias("read_length_bin"))
        columns.append("read_length_bin")

    return out.select(columns).sort(
        ["sample_id", "transcript_id", "pos"]
        + (["read_length_bin"] if "read_length_bin" in columns else [])
    )


def make_rdg_flux_metadata(
    *,
    profile_path: str,
    frame_support_path: str | None,
    sample_ids: list[str],
    psite_offset_model: str,
    annotation_source: str,
    transcriptome_fasta: str,
    model_stage: str,
    normalization: str,
    frame_method: str,
    background_model: str,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "coordinate_system": "transcript_0_based_half_open_positions",
        "psite_offset_model": psite_offset_model,
        "annotation_source": annotation_source,
        "transcriptome_fasta": transcriptome_fasta,
        "translonscorer_version": __version__,
        "model_stage": model_stage,
        "normalization": normalization,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "profile_path": profile_path,
        "frame_support_path": frame_support_path,
        "sample_ids": sample_ids,
        "frame_method": frame_method,
        "background_model": background_model,
        "probability_semantics": "p0/p1/p2 in frame-support inputs are conditional frame posteriors. RDG-Flux p_frame0..2 are scaled by p_translated; residual mass is p_background so p_frame0+p_frame1+p_frame2+p_background=1 within floating point tolerance.",
        "translation_support_semantics": "When support_evidence is present in frame support, p_translated=(1-background_probability)*support_evidence. Otherwise p_translated falls back to 1-background_probability only for codons with frame support and 0 for unsupported codons.",
    }
    if extra:
        metadata.update(extra)
    return metadata


def write_metadata_json(metadata: dict[str, Any], metadata_path: str) -> str:
    ensure_dir_for_file(metadata_path)
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return metadata_path


def _default_metadata_path(out_path: str, partition_by_sample: bool) -> str:
    if partition_by_sample:
        return os.path.join(out_path, "metadata.json")
    if out_path.endswith(".parquet"):
        return f"{out_path[:-8]}.metadata.json"
    return f"{out_path}.metadata.json"


def write_rdg_flux_position_table(
    table: pl.DataFrame,
    out_path: str,
    *,
    partition_by_sample: bool = False,
) -> str:
    if not partition_by_sample:
        return write_parquet_safe(table, out_path)

    os.makedirs(out_path, exist_ok=True)
    for sample in table.get_column("sample_id").unique().to_list():
        safe = str(sample).replace("/", "_")
        sample_dir = Path(out_path) / f"sample_id={safe}"
        sample_dir.mkdir(parents=True, exist_ok=True)
        table.filter(pl.col("sample_id") == sample).write_parquet(sample_dir / "part.parquet")
    return out_path


def _export_rdg_flux_v1_lazy_precomputed(
    *,
    profiles_path: str,
    frame_support_path: str,
    out_path: str,
    sample_id: str | None,
    background_probability: float,
    preserve_read_length: bool,
) -> tuple[list[str], int]:
    profiles_lf, sample_ids, profile_cols = _normalise_profile_lazy(profiles_path, sample_id)
    frame_lf, frame_cols = _normalise_frame_lazy(frame_support_path, sample_ids)
    out_lf = _rdg_flux_position_lazy(
        profiles_lf,
        profile_cols,
        frame_lf,
        frame_cols,
        background_probability=background_probability,
        preserve_read_length=preserve_read_length,
    )
    ensure_dir_for_file(out_path)
    out_lf.sink_parquet(out_path)
    rows = int(pl.scan_parquet(out_path).select(pl.len().alias("rows")).collect()["rows"][0])
    return sample_ids, rows


def _build_frame_support_from_annotation(
    profiles: pl.DataFrame,
    *,
    annotation_path: str | None,
    annotation_dir: str | None,
    cds_path: str | None,
    frame_method: str,
    frame_by_length: bool,
) -> pl.DataFrame:
    if cds_path:
        cds_tran = _load_table(cds_path)
    elif annotation_dir:
        from ..io.annotation_bundle import load_annotation_bundle

        exon_df, cds_df, _feats_df, _fmap_df, _tx_df, _loci_bed, _manifest = load_annotation_bundle(
            annotation_dir
        )
        cds_tran = cds_to_transcript_space(cds_df, exon_df)
    elif annotation_path:
        cds_df, exon_df = bam_handlers.getexons_and_cds(annotation_path)
        cds_tran = cds_to_transcript_space(cds_df, exon_df)
    else:
        raise ValueError(
            "Provide --frame-support, or provide --cds/--annotation/--annotation-dir to build frame support"
        )

    params = FrameSupportParams(frame_method=frame_method, frame_by_length=frame_by_length)
    return build_frame_support(profiles, cds_tran, params)


def export_rdg_flux_v1(
    *,
    profiles_path: str,
    out_path: str,
    frame_support_path: str | None = None,
    sample_id: str | None = None,
    annotation_path: str | None = None,
    annotation_dir: str | None = None,
    cds_path: str | None = None,
    transcriptome_fasta: str | None = None,
    psite_offset_model: str = "unknown",
    annotation_source: str | None = None,
    normalization: str = "raw_psite_counts",
    model_stage: str = "stage1_frame_posterior",
    frame_method: str = "linear+hmm",
    frame_by_length: bool = True,
    background_probability: float = 0.0,
    background_model: str = "none",
    preserve_read_length: bool = False,
    partition_by_sample: bool = False,
    metadata_path: str | None = None,
) -> dict[str, str]:
    lazy_precomputed = (
        frame_support_path is not None
        and profiles_path.endswith(".parquet")
        and frame_support_path.endswith(".parquet")
        and not partition_by_sample
    )

    if lazy_precomputed:
        if frame_support_path is None:
            raise ValueError("lazy_precomputed export requires frame_support_path")
        sample_ids, rows = _export_rdg_flux_v1_lazy_precomputed(
            profiles_path=profiles_path,
            frame_support_path=frame_support_path,
            out_path=out_path,
            sample_id=sample_id,
            background_probability=background_probability,
            preserve_read_length=preserve_read_length,
        )
        frame_method_used = "precomputed"
    else:
        profiles = _load_table(profiles_path)
        if frame_support_path:
            frame_support = _load_table(frame_support_path)
            frame_method_used = "precomputed"
        else:
            frame_support = _build_frame_support_from_annotation(
                profiles,
                annotation_path=annotation_path,
                annotation_dir=annotation_dir,
                cds_path=cds_path,
                frame_method=frame_method,
                frame_by_length=frame_by_length,
            )
            frame_method_used = frame_method

        table = rdg_flux_position_table(
            profiles,
            frame_support,
            sample_id=sample_id,
            background_probability=background_probability,
            preserve_read_length=preserve_read_length,
        )
        write_rdg_flux_position_table(table, out_path, partition_by_sample=partition_by_sample)
        sample_ids = sorted(str(x) for x in table.get_column("sample_id").unique().to_list())
        rows = table.height

    metadata = make_rdg_flux_metadata(
        profile_path=profiles_path,
        frame_support_path=frame_support_path,
        sample_ids=sample_ids,
        psite_offset_model=psite_offset_model,
        annotation_source=annotation_source
        or annotation_path
        or annotation_dir
        or cds_path
        or "unknown",
        transcriptome_fasta=transcriptome_fasta or "unknown",
        model_stage=model_stage,
        normalization=normalization,
        frame_method=frame_method_used,
        background_model=background_model,
        extra={
            "partition_by_sample": partition_by_sample,
            "preserve_read_length": preserve_read_length,
            "rows": rows,
            "lazy_precomputed_export": lazy_precomputed,
        },
    )
    meta_path = metadata_path or _default_metadata_path(out_path, partition_by_sample)
    write_metadata_json(metadata, meta_path)
    return {"parquet": out_path, "metadata": meta_path}

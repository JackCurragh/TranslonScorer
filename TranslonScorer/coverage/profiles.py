"""
Profiles pipeline: compute transcript-space A-site coverage profiles
from classic BAMs, collapsed BAMs, or Zarr unique-read matrices.

Outputs tidy profiles for clustering/scoring without requiring bigWig.

Schema (profiles):
  - tran_id: str
  - pos: int (transcript coordinate, 0-based)
  - count: float

Optionally include a 'sample' column for multi-sample sources.
"""
from __future__ import annotations

from typing import Dict, Iterable, Iterator, Optional, Tuple
import os

import polars as pl

from ..file_handlers import bam as bam_handlers
from ..file_handlers import sparse_parquet as sparse_parquet_handlers
from ..file_handlers import zarr as zarr_handlers
from .mapped_index import build_mapped_index
from ..file_handlers import bigwig as bigwig_handlers
from ..core import coordinates
from ..utils import log_info, log_warning
from ..utils.io import write_parquet_safe
from .transcript_coords import cds_to_transcript_space


def _ensure_transcript_coords(reads_df: pl.DataFrame, exon_df: pl.DataFrame) -> pl.DataFrame:
    """Ensure transcript coordinates exist for reads.

    If reads_df already has 'tran_start_bam' and 'tran_id', return as-is.
    If reads_df appears transcriptomic (chr equals transcript ids), map columns accordingly.
    Else, project genomic reads to transcript coords via bamtranscript().
    """
    cols = set(reads_df.columns)
    if {'tran_id', 'tran_start_bam'}.issubset(cols):
        return reads_df

    # Heuristic: if 'chr' overlaps many tran_ids in exon_df, treat as transcriptomic
    if 'chr' in cols:
        # Sample small subset to test overlap
        sample_chr = set(reads_df.get_column('chr').head(1000).to_list())
        exon_tran = set(exon_df.get_column('tran_id').head(1000).to_list())
        if len(sample_chr.intersection(exon_tran)) > 10:
            df = reads_df.rename({'chr': 'tran_id'})
            return df.with_columns(
                pl.col('start').alias('tran_start_bam')
            )

    # Default: genomic reads → transcript projection
    projected = bam_handlers.bamtranscript(reads_df, exon_df)
    return projected


def _compute_asite_profiles(df_with_tran: pl.DataFrame, offsets: Dict[int, int], *, keep_length: bool = False) -> pl.DataFrame:
    """Compute A-site transcript profiles (tran_id, pos, count) from transcript-mapped reads.

    Expects columns: tran_id, tran_start_bam, length, count
    """
    if df_with_tran.is_empty():
        return pl.DataFrame({'tran_id': [], 'pos': [], 'count': []})

    # Map per-length offsets; unknown lengths default to 15
    def _ofs(L: int) -> int:
        try:
            return int(offsets.get(int(L), 15))
        except Exception:
            return 15

    sample_cols = [c for c in ("sample_id", "sample_index", "study_id", "study_id_int") if c in df_with_tran.columns]
    cols = sample_cols + ['tran_id', 'pos'] + (['length'] if keep_length and 'length' in df_with_tran.columns else [])
    out = (
        df_with_tran
        .with_columns(
            pl.col('length').map_elements(_ofs).alias('ofs'),
            pl.col('tran_start_bam').cast(pl.Int64)
        )
        .with_columns((pl.col('tran_start_bam') + pl.col('ofs')).alias('pos'))
        .select(sample_cols + ['tran_id', 'pos', 'count'] + (['length'] if keep_length and 'length' in df_with_tran.columns else []))
        .group_by(cols)
        .agg(pl.col('count').sum())
        .rename({'count': 'count'})
        .sort(cols)
    )
    return out


def profiles_from_bam(
    bam_path: str,
    exon_df: pl.DataFrame,
    cds_df: pl.DataFrame,
    *,
    collapsed: bool = False,
    count_from: Optional[str] = None,
    count_pattern: Optional[str] = None,
    count_tag: Optional[str] = None,
    keep_length: bool = False,
) -> tuple[pl.DataFrame, Dict[int,int]]:
    """Compute transcript A-site profiles from a BAM (classic or collapsed).

    Returns (profiles_df, offsets_dict).
    """
    log_info("Reading BAM for profiles…")
    cds_tran_df = cds_to_transcript_space(cds_df, exon_df)
    reads = bam_handlers.readbam(
        bam_path,
        collapsed=collapsed,
        count_from=count_from,
        count_pattern=count_pattern,
        count_tag=count_tag,
        include_qname=False,
    )

    # Determine BAM type and, if genomic, project to transcripts
    try:
        bam_type, _ = bam_handlers.detect_bam_type(reads, exon_df)
    except Exception:
        bam_type = 'genomic'

    if bam_type == 'genomic':
        reads = bam_handlers.bamtranscript(reads, exon_df)
    else:
        reads = reads.rename({'chr': 'tran_id'}).with_columns(
            pl.col('start').alias('tran_start_bam')
        )

    # Compute offsets from a change-point analysis on transcriptomic reads
    # Reuse existing routine which expects positions relative to CDS; for now,
    # derive a coarse per-length offset by calling change_point_analysis on the
    # transcriptomic read table after joining CDS (as done in process_transcriptomic_bam).
    reads_rel = bam_handlers.process_transcriptomic_bam(reads, cds_tran_df)
    offsets = coordinates.change_point_analysis(reads_rel)

    profiles = _compute_asite_profiles(reads, offsets, keep_length=keep_length)
    return profiles, {int(k): int(v) for k, v in offsets.items()}


def _default_offsets_for_profiles(reads: pl.DataFrame, default_offset: int) -> dict[int, int]:
    if reads.is_empty() or "length" not in reads.columns:
        return {0: int(default_offset)}
    return {
        int(length): int(default_offset)
        for length in reads.get_column("length").drop_nulls().unique().to_list()
    }


def _transcript_gene_map(
    *,
    exon_df: pl.DataFrame,
    transcripts_df: pl.DataFrame | None = None,
    feature_map_df: pl.DataFrame | None = None,
) -> pl.DataFrame:
    if transcripts_df is not None and {"tran_id", "gene_id"}.issubset(set(transcripts_df.columns)):
        return transcripts_df.select(["tran_id", "gene_id"]).drop_nulls().unique()
    if feature_map_df is not None and "locus_id" in feature_map_df.columns:
        tx_col = "tran_id" if "tran_id" in feature_map_df.columns else ("transcript_id" if "transcript_id" in feature_map_df.columns else None)
        if tx_col:
            out = feature_map_df.select([tx_col, "locus_id"]).drop_nulls().unique()
            if tx_col != "tran_id":
                out = out.rename({tx_col: "tran_id"})
            return out.rename({"locus_id": "gene_id"})
    return exon_df.select("tran_id").unique().with_columns(pl.col("tran_id").alias("gene_id"))


def gene_expression_matrix_from_profiles(
    profiles: pl.DataFrame,
    *,
    exon_df: pl.DataFrame,
    transcripts_df: pl.DataFrame | None = None,
    feature_map_df: pl.DataFrame | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Build long and wide gene expression matrices from sample-resolved profiles."""
    if profiles.is_empty():
        empty_long = pl.DataFrame(schema={"sample_id": pl.Utf8, "gene_id": pl.Utf8, "count": pl.Float64})
        empty_wide = pl.DataFrame(schema={"gene_id": pl.Utf8})
        return empty_long, empty_wide
    if "sample_id" not in profiles.columns:
        raise ValueError("Gene expression matrix requires profiles with sample_id")
    tx_gene = _transcript_gene_map(
        exon_df=exon_df,
        transcripts_df=transcripts_df,
        feature_map_df=feature_map_df,
    )
    long = (
        profiles.join(tx_gene, on="tran_id", how="left")
        .with_columns(pl.coalesce([pl.col("gene_id"), pl.col("tran_id")]).alias("gene_id"))
        .group_by(["gene_id", "sample_id"])
        .agg(pl.col("count").sum().alias("count"))
        .sort(["gene_id", "sample_id"])
    )
    wide = (
        long.pivot(index="gene_id", on="sample_id", values="count", aggregate_function="sum")
        .fill_null(0)
        .sort("gene_id")
    )
    return long, wide


def profiles_from_sparse_parquet_matrix(
    *,
    bam_path: str,
    manifest_path: str,
    exon_df: pl.DataFrame,
    cds_df: pl.DataFrame,
    samples: list[str] | None = None,
    default_offset: int = 15,
    keep_length: bool = False,
    regions: list[sparse_parquet_handlers.Region] | None = None,
) -> tuple[pl.DataFrame, Dict[int, int], pl.DataFrame]:
    """Compute sample-resolved transcript A-site profiles from a sparse Parquet matrix.

    Returns ``(profiles, offsets, genomic_counts)``. ``profiles`` keeps one
    long-form row per sample/run, transcript, and position.
    """
    log_info("Joining sparse Parquet matrix counts to global BAM alignments...")
    genomic_counts = sparse_parquet_handlers.sparse_matrix_genomic_counts(
        bam_path=bam_path,
        manifest_path=manifest_path,
        exon_df=exon_df if regions is None else None,
        regions=regions,
        sample_names=samples,
    )
    if genomic_counts.is_empty():
        return (
            pl.DataFrame(schema={"sample_id": pl.Utf8, "tran_id": pl.Utf8, "pos": pl.Int64, "count": pl.Float64}),
            {},
            genomic_counts,
        )
    mapped = bam_handlers.bamtranscript(genomic_counts, exon_df)
    if mapped.is_empty():
        return (
            pl.DataFrame(schema={"sample_id": pl.Utf8, "tran_id": pl.Utf8, "pos": pl.Int64, "count": pl.Float64}),
            {},
            genomic_counts,
        )
    cds_tran_df = cds_to_transcript_space(cds_df, exon_df)
    try:
        rel = bam_handlers.process_transcriptomic_bam(mapped, cds_tran_df)
        offsets = coordinates.change_point_analysis(rel)
        offsets = {int(k): int(v) for k, v in offsets.items()}
    except Exception as exc:
        log_warning(f"Could not infer offsets from sparse matrix profiles; using default offset {default_offset}: {exc}")
        offsets = _default_offsets_for_profiles(mapped, default_offset)
    for length in mapped.get_column("length").drop_nulls().unique().to_list():
        offsets.setdefault(int(length), int(default_offset))
    profiles = _compute_asite_profiles(mapped, offsets, keep_length=keep_length)
    return profiles, offsets, genomic_counts


def profiles_from_zarr(
    zarr_root: str,
    read_index_parquet: str,
    samples: list[str],
    exon_df: pl.DataFrame,
    cds_df: pl.DataFrame,
    *,
    default_offset: int = 15,
    offsets_mode: str = 'auto',
    offsets_out: str | None = None,
    sample_cap_per_len: int = 50000,
    mapped_index_parquet: str | None = None,
) -> Iterator[Tuple[str, pl.DataFrame]]:
    """Compute transcript A-site profiles from Zarr (multi-sample). Yields per-sample profiles.

    If mapped_index_parquet is provided (or path does not exist and can be written),
    we use a precomputed mapping read_id->(tran_id,tran_start_bam,length,strand) to
    avoid re-running genomic->transcript projection for each sample/chunk.
    """
    # zarr_handlers safely defers importing the heavy zarr dependency
    log_info("Streaming Zarr + index for profiles…")
    cds_tran_df = cds_to_transcript_space(cds_df, exon_df)
    offsets_by_sample: Dict[str, Dict[int, int]] = {}
    sampled_by_sample_len: Dict[tuple[str,int], int] = {}

    mapped_path: str | None = None
    if mapped_index_parquet:
        mapped_path = mapped_index_parquet
        if not os.path.exists(mapped_path):
            log_info("Building mapped index (one-time)…")
            build_mapped_index(read_index_parquet, exon_df, out_parquet=mapped_path)

    for sample, chunk in zarr_handlers.iter_reads_from_zarr(
        zarr_root, read_index_parquet, samples, include_read_id=bool(mapped_path)
    ):
        if chunk.is_empty():
            continue

        if mapped_path:
            # Join to mapped index via read_id (faster; no projection)
            # Ensure read_id present; the index provides it
            if 'read_id' not in chunk.columns:
                # Add sequential row ids won't match; enforce safeguard
                raise RuntimeError("read_id missing in chunk; expected from index join")
            mapped = (
                chunk.join(
                    pl.scan_parquet(mapped_path).select(["read_id","tran_id","tran_start_bam","length","strand"]).collect(),
                    on=["read_id","length","strand"],  # use length/strand to reduce collisions
                    how='inner'
                )
            )
            if mapped.is_empty():
                continue
        else:
            # Project to transcript coords on the fly
            mapped = _ensure_transcript_coords(chunk, exon_df)

        # Offsets per sample
        if sample not in offsets_by_sample:
            offsets_by_sample[sample] = {}
        offsets = offsets_by_sample[sample]

        if offsets_mode == 'auto':
            # Accumulate a capped sample per read length for change-point
            need: set[int] = set(int(x) for x in mapped.get_column('length').unique().to_list() if int(x) not in offsets)
            if need:
                # Prepare input for relative-to-CDS shift detection:
                # ensure the 'start' column refers to transcript start.
                if 'tran_start_bam' in mapped.columns:
                    rel_input = mapped.with_columns(pl.col('tran_start_bam').alias('start'))
                elif 'tran_start' in mapped.columns:
                    rel_input = mapped.with_columns(pl.col('tran_start').alias('start'))
                else:
                    rel_input = mapped
                # Join CDS to compute relative to CDS start using existing helper
                rel = bam_handlers.process_transcriptomic_bam(
                    rel_input,
                    cds_tran_df
                )
                # For each needed length, take up to sample_cap_per_len rows
                rows = []
                for L in list(need):
                    sub = rel.filter(pl.col('length') == int(L)).select(['bamcds_start','length','count'])
                    if sub.height == 0:
                        continue
                    # Downsample if necessary
                    take = min(sample_cap_per_len, sub.height)
                    rows.append(sub.head(take))
                if rows:
                    pool = pl.concat(rows)
                    off = coordinates.change_point_analysis(pool)
                    offsets.update({int(k): int(v) for k,v in off.items()})
                # Fill any remaining unseen with default
                for L in need:
                    offsets.setdefault(int(L), default_offset)
        elif offsets_mode == 'required':
            # Expect provided offsets via file (handled by caller); nothing to do
            pass
        else:
            # global default
            lengths = set(int(x) for x in mapped.get_column('length').unique().to_list())
            for L in lengths:
                offsets.setdefault(int(L), default_offset)

        prof = _compute_asite_profiles(mapped, offsets, keep_length=True)
        yield sample, prof

    # Persist offsets if requested
    if offsets_out:
        rows = []
        for s, od in offsets_by_sample.items():
            for L, ofs in od.items():
                rows.append({'sample': s, 'length': int(L), 'offset': int(ofs)})
        if rows:
            pl.from_dicts(rows).write_csv(offsets_out)


def write_profiles_parquet(df: pl.DataFrame, out_path: str, sample: Optional[str] = None) -> str:
    """Write profiles to Parquet; include sample column if provided (safe I/O)."""
    if sample:
        df = df.with_columns(pl.lit(sample).alias('sample'))
    return write_parquet_safe(df, out_path)


def profiles_from_bigwig(
    bigwig_path: str | dict,
    exon_df: pl.DataFrame,
    *,
    stranded: bool = False,
) -> pl.DataFrame:
    """Compute transcript-space profiles from a genomic or transcriptomic BigWig.

    If ``bigwig_path`` is a dict with keys {'forward','reverse'}, strands are summed unless
    ``stranded=True`` is requested (stranded currently not emitted separately; we sum).

    Returns tidy profiles with columns: tran_id, pos, count.
    """
    # bigwig_handlers imports pyBigWig at module import; required for this path

    def _bw_to_profiles(bw) -> pl.DataFrame:
        tran = bigwig_handlers.transcriptreads(bw, exon_df)
        if tran.is_empty():
            return pl.DataFrame({'tran_id': [], 'pos': [], 'count': []})
        return (
            tran.select(['tran_id', 'tran_start', 'counts'])
                .rename({'tran_start': 'pos', 'counts': 'count'})
                .group_by(['tran_id', 'pos']).agg(pl.col('count').sum())
                .sort(['tran_id', 'pos'])
        )

    # Single BigWig path
    if isinstance(bigwig_path, str):
        return _bw_to_profiles(bigwig_path)

    # Strand-specific dict
    fwd = bigwig_path.get('forward')
    rev = bigwig_path.get('reverse')
    f_df = _bw_to_profiles(fwd) if fwd else pl.DataFrame({'tran_id': [], 'pos': [], 'count': []})
    r_df = _bw_to_profiles(rev) if rev else pl.DataFrame({'tran_id': [], 'pos': [], 'count': []})

    if f_df.is_empty() and r_df.is_empty():
        return pl.DataFrame({'tran_id': [], 'pos': [], 'count': []})

    # Sum strands into total profiles for now
    both = pl.concat([f_df, r_df]) if not (f_df.is_empty() or r_df.is_empty()) else (f_df if not f_df.is_empty() else r_df)
    return (
        both.group_by(['tran_id', 'pos']).agg(pl.col('count').sum())
            .sort(['tran_id', 'pos'])
    )

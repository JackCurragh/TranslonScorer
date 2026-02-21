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
from ..file_handlers import zarr as zarr_handlers
from ..file_handlers import bigwig as bigwig_handlers
from ..core import coordinates
from ..utils import log_info, log_warning


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


def _compute_asite_profiles(df_with_tran: pl.DataFrame, offsets: Dict[int, int]) -> pl.DataFrame:
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

    out = (
        df_with_tran
        .with_columns(
            pl.col('length').map_elements(_ofs).alias('ofs'),
            pl.col('tran_start_bam').cast(pl.Int64)
        )
        .with_columns((pl.col('tran_start_bam') + pl.col('ofs')).alias('pos'))
        .select(['tran_id', 'pos', 'count'])
        .group_by(['tran_id', 'pos'])
        .agg(pl.col('count').sum())
        .rename({'count': 'count'})
        .sort(['tran_id', 'pos'])
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
) -> tuple[pl.DataFrame, Dict[int,int]]:
    """Compute transcript A-site profiles from a BAM (classic or collapsed).

    Returns (profiles_df, offsets_dict).
    """
    log_info("Reading BAM for profiles…")
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
    reads_rel = bam_handlers.process_transcriptomic_bam(reads, cds_df)
    offsets = coordinates.change_point_analysis(reads_rel)

    profiles = _compute_asite_profiles(reads, offsets)
    return profiles, {int(k): int(v) for k, v in offsets.items()}


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
) -> Iterator[Tuple[str, pl.DataFrame]]:
    """Compute transcript A-site profiles from Zarr (multi-sample). Yields per-sample profiles."""
    # zarr_handlers safely defers importing the heavy zarr dependency
    log_info("Streaming Zarr + index for profiles…")
    offsets_by_sample: Dict[str, Dict[int, int]] = {}
    sampled_by_sample_len: Dict[tuple[str,int], int] = {}

    for sample, chunk in zarr_handlers.iter_reads_from_zarr(
        zarr_root, read_index_parquet, samples
    ):
        if chunk.is_empty():
            continue

        # Project to transcript coords
        mapped = _ensure_transcript_coords(chunk, exon_df)

        # Offsets per sample
        if sample not in offsets_by_sample:
            offsets_by_sample[sample] = {}
        offsets = offsets_by_sample[sample]

        if offsets_mode == 'auto':
            # Accumulate a capped sample per read length for change-point
            need: set[int] = set(int(x) for x in mapped.get_column('length').unique().to_list() if int(x) not in offsets)
            if need:
                # Join CDS to compute relative to CDS start using existing helper
                rel = bam_handlers.process_transcriptomic_bam(
                    mapped.rename({'tran_start_bam': 'start'}) if 'tran_start_bam' in mapped.columns else mapped.rename({'tran_start': 'start'}),
                    cds_df
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

        prof = _compute_asite_profiles(mapped, offsets)
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
    """Write profiles to Parquet; include sample column if provided."""
    if sample:
        df = df.with_columns(pl.lit(sample).alias('sample'))
    # Ensure parent directory exists (Polars does not create it)
    parent = os.path.dirname(out_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    df.write_parquet(out_path)
    return out_path


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

"""P-site index: build-once, query-many cohort cache.

Design
------
Build (once per cohort):
  Scan each partition BAM, join with count parquets, emit
  (chrom, pos5, strand, length, sample_id, count) for every in-CDS read.
  Write to Hive-partitioned Parquet sorted by pos5 within each chrom.

  ``pos5`` is the raw 5' end of the read (reference_start for +,
  reference_end-1 for -).  No P-site offset is applied at build time —
  the offset is a calibration artefact applied analytically at query time.

Query (O(log N + hits) per feature, all features on a chrom in one pass):
  1. Polars range scan on chrom shard with row-group pushdown.
  2. Vectorised join on calibrated offset table → p_site = pos5 + offset.
  3. Sweep-line assignment to (feature_id, tx_pos).
  4. Group-by → FrameRollup or full CoverageIndex.

Storage layout::

    psite_index/
      chrom=chr1/data.parquet
      chrom=chr22/data.parquet
      samples.parquet            # (sample_id UInt16, sample_name Utf8)

Schema per chrom shard (pos5-sorted):
    pos5      Int32
    strand    Boolean     (True = '+')
    length    UInt8
    sample_id UInt16
    count     Float32
"""

from __future__ import annotations

import io
import multiprocessing as mp
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
import pysam

from .io.annotation import build_gene_spans
from .io.bam import normalise_chrom as _normalise_chrom
from .io.matrix import _count_parquets, _discover_bam, _manifest, _parse_read_id, _samples_df
from .matrix_qc import _bam_chroms, _build_frame_intervals
from .utils.logging import log_info

# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

_PSITE_SCHEMA = {
    "pos5": pl.Int32,
    "strand": pl.Boolean,
    "length": pl.UInt8,
    "sample_id": pl.UInt16,
    "count": pl.Float32,
}

# Schema returned by query functions
_ROLLUP_SCHEMA = {
    "sample_id": pl.UInt16,
    "feature_id": pl.Utf8,
    "length": pl.UInt8,
    "frame0": pl.Float32,
    "frame1": pl.Float32,
    "frame2": pl.Float32,
    "n_reads": pl.Float32,
}

# ---------------------------------------------------------------------------
# Build: worker
# ---------------------------------------------------------------------------

_BCTX: Dict = {}


def _build_worker_init(frame_ivs, bam_refs, sample_map) -> None:
    _BCTX.update(frame_ivs=frame_ivs, bam_refs=bam_refs, sample_map=sample_map)


def _build_worker(pdir_str: str) -> bytes:
    """Scan one partition → IPC bytes of (chrom, pos5, strand, length, sample_id, count)."""
    pdir = Path(pdir_str)
    bam = _discover_bam(pdir)
    empty_schema = {"chrom": pl.Utf8, **_PSITE_SCHEMA}
    empty = pl.DataFrame(schema=empty_schema)

    if bam is None:
        buf = io.BytesIO()
        empty.write_ipc(buf)
        return buf.getvalue()

    frame_ivs = _BCTX["frame_ivs"]
    bam_refs = _BCTX["bam_refs"]
    sample_map: pl.DataFrame = _BCTX["sample_map"]

    # --- scan BAM for in-CDS reads, keep raw position ---
    read_id_list: List[int] = []
    pos5_list: List[int] = []
    length_list: List[int] = []
    strand_list: List[bool] = []
    chrom_list: List[str] = []
    rec_count: Dict[int, int] = {}

    with pysam.AlignmentFile(str(bam), "rb") as bf:
        for rec in bf.fetch(until_eof=True):
            if rec.is_unmapped or rec.reference_name is None:
                continue
            rid = _parse_read_id(rec.query_name)
            if rid is None:
                continue
            length = int(rec.query_length or 0)
            if length == 0:
                continue
            rec_count[rid] = rec_count.get(rid, 0) + 1
            is_rev = rec.is_reverse
            p5 = int(rec.reference_end) - 1 if is_rev else int(rec.reference_start)
            chrom = _normalise_chrom(rec.reference_name, bam_refs) or rec.reference_name
            read_id_list.append(rid)
            pos5_list.append(p5)
            length_list.append(length)
            strand_list.append(not is_rev)
            chrom_list.append(chrom)

    if not read_id_list:
        buf = io.BytesIO()
        empty.write_ipc(buf)
        return buf.getvalue()

    # filter to in-CDS reads using frame_ivs membership (any CDS frame ≥ 0)
    # reuse _assign_frames_sweep just for the membership test
    from .matrix_qc import _assign_frames_sweep

    n = len(read_id_list)
    in_cds = np.zeros(n, dtype=bool)

    # group by (chrom, strand) for sweep-line
    from collections import defaultdict
    grp: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for i in range(n):
        s = "+" if strand_list[i] else "-"
        grp[(chrom_list[i], s)].append(i)

    for (chrom, strand), idxs in grp.items():
        fiv = frame_ivs.get((chrom, strand))
        if not fiv:
            continue
        asites = np.array([pos5_list[i] for i in idxs], dtype=np.int64)
        hits = _assign_frames_sweep(idxs, asites, strand, fiv)
        for idx in hits:
            in_cds[idx] = True

    in_cds_idx = np.where(in_cds)[0]
    if len(in_cds_idx) == 0:
        buf = io.BytesIO()
        empty.write_ipc(buf)
        return buf.getvalue()

    aln = pl.DataFrame(
        {
            "read_id": pl.Series([read_id_list[i] for i in in_cds_idx], dtype=pl.UInt64),
            "chrom": pl.Series([chrom_list[i] for i in in_cds_idx], dtype=pl.Utf8),
            "pos5": pl.Series([pos5_list[i] for i in in_cds_idx], dtype=pl.Int32),
            "strand": pl.Series([strand_list[i] for i in in_cds_idx], dtype=pl.Boolean),
            "length": pl.Series([length_list[i] for i in in_cds_idx], dtype=pl.UInt8),
        }
    )

    # join counts
    try:
        mp_, manifest = _manifest(pdir)
        cfiles = _count_parquets(mp_, manifest)
    except FileNotFoundError:
        cfiles = []
    if not cfiles:
        buf = io.BytesIO()
        empty.write_ipc(buf)
        return buf.getvalue()

    counts = pl.read_parquet(cfiles).select(["read_id", "sample_id", "count"])
    out = (
        counts.join(aln, on="read_id", how="inner")
        .join(sample_map.select("sample_id"), on="sample_id", how="inner")
        .with_columns(pl.col("count").cast(pl.Float32))
        .group_by(["chrom", "pos5", "strand", "length", "sample_id"])
        .agg(pl.col("count").sum())
        .with_columns(
            pl.col("pos5").cast(pl.Int32),
            pl.col("length").cast(pl.UInt8),
            pl.col("sample_id").cast(pl.UInt16),
        )
    )

    buf = io.BytesIO()
    out.write_ipc(buf)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Build: orchestrator
# ---------------------------------------------------------------------------


def build_psite_index(
    partition_dirs,
    cds_df: pl.DataFrame,
    out_dir: str | Path,
    *,
    n_workers: Optional[int] = None,
    chrom: Optional[str] = None,
) -> None:
    """Build the P-site index from partition BAMs.

    Scans every partition once, joins read counts, and writes
    chrom-partitioned Parquet to ``out_dir``.  Subsequent calls with a
    different annotation or ``chrom`` filter append to / overwrite only the
    affected shards.

    Args:
        partition_dirs: path to global_partitioned directory, or list of dirs.
        cds_df:         CDS blocks (from build_cds_blocks or build_cds_blocks_from_bigbed)
                        used to restrict the index to in-CDS reads.
        out_dir:        where to write the index.
        n_workers:      worker processes (default: cpu_count).
        chrom:          restrict to one chromosome (useful for testing).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir() and _discover_bam(d))
    else:
        dirs = [Path(d) for d in partition_dirs if _discover_bam(Path(d))]

    if not dirs:
        log_info("build_psite_index: no partitions found")
        return

    bam_refs = _bam_chroms(_discover_bam(dirs[0]))
    frame_ivs = _build_frame_intervals(cds_df, cds_df, bam_refs)
    if chrom:
        frame_ivs = {k: v for k, v in frame_ivs.items() if k[0] == chrom or k[0] == f"chr{chrom}"}

    mp0, manifest0 = _manifest(dirs[0])
    sample_map = _samples_df(mp0, manifest0).select(
        pl.col("sample_id").cast(pl.UInt16), pl.col("sample_name")
    )

    if n_workers is None:
        n_workers = min(len(dirs), os.cpu_count() or 1)

    log_info(
        f"build_psite_index: {len(dirs)} partitions, "
        f"{sum(len(v) for v in frame_ivs.values()):,} CDS intervals, "
        f"{n_workers} worker(s) → {out_dir}"
    )

    parts: List[pl.DataFrame] = []
    dir_strs = [str(d) for d in dirs]

    if n_workers <= 1:
        _build_worker_init(frame_ivs, bam_refs, sample_map)
        for ds in dir_strs:
            df = pl.read_ipc(io.BytesIO(_build_worker(ds)))
            if not df.is_empty():
                parts.append(df)
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            n_workers,
            initializer=_build_worker_init,
            initargs=(frame_ivs, bam_refs, sample_map),
        ) as pool:
            for i, ipc in enumerate(pool.imap_unordered(_build_worker, dir_strs, chunksize=1)):
                df = pl.read_ipc(io.BytesIO(ipc))
                if not df.is_empty():
                    parts.append(df)
                if (i + 1) % 64 == 0:
                    log_info(f"  scanned {i + 1}/{len(dirs)} partitions")

    if not parts:
        log_info("build_psite_index: no in-CDS reads found")
        return

    combined = pl.concat(parts)
    log_info(f"  {combined.height:,} raw rows; grouping and writing per-chrom shards …")

    # aggregate duplicates that arise from merging partition outputs
    combined = (
        combined.group_by(["chrom", "pos5", "strand", "length", "sample_id"])
        .agg(pl.col("count").sum())
        .with_columns(
            pl.col("pos5").cast(pl.Int32),
            pl.col("length").cast(pl.UInt8),
            pl.col("sample_id").cast(pl.UInt16),
            pl.col("count").cast(pl.Float32),
        )
    )

    for chrom_val, shard in combined.group_by("chrom"):
        chrom_str = str(chrom_val[0]) if isinstance(chrom_val, tuple) else str(chrom_val)
        shard_dir = out_dir / f"chrom={chrom_str}"
        shard_dir.mkdir(exist_ok=True)
        (
            shard.drop("chrom")
            .sort("pos5")
            .write_parquet(shard_dir / "data.parquet", statistics=True, compression="zstd")
        )

    sample_map.write_parquet(out_dir / "samples.parquet")
    log_info(f"build_psite_index: done — {out_dir}")


# ---------------------------------------------------------------------------
# Query helpers
# ---------------------------------------------------------------------------


def load_samples(index_dir: str | Path) -> pl.DataFrame:
    """Return (sample_id UInt16, sample_name Utf8) table."""
    return pl.read_parquet(Path(index_dir) / "samples.parquet")


def available_chroms(index_dir: str | Path) -> List[str]:
    return [
        d.name.split("=", 1)[1]
        for d in sorted(Path(index_dir).iterdir())
        if d.is_dir() and d.name.startswith("chrom=")
    ]


def _load_shard(index_dir: Path, chrom: str) -> Optional[pl.DataFrame]:
    shard = index_dir / f"chrom={chrom}" / "data.parquet"
    if not shard.exists():
        return None
    return pl.read_parquet(shard)


# ---------------------------------------------------------------------------
# Query: FrameRollup from P-site index
# ---------------------------------------------------------------------------


def _exon_intervals_for_chrom(
    cds_df: pl.DataFrame, chrom: str
) -> List[Tuple[str, str, int, int, int, List[Tuple[int, int, int]]]]:
    """Return per-feature exon block info for a single chrom.

    Returns list of (feature_id, strand, genomic_min, genomic_max, n_exons,
                     [(exon_start, exon_stop, tran_start), ...])
    """
    sub = cds_df.filter(pl.col("chr") == chrom)
    if sub.is_empty():
        return []
    out = []
    for row in sub.iter_rows(named=True):
        starts = row["start"]
        stops = row["stop"]
        tran_starts = row["tran_start"]
        exons = list(zip(starts, stops, tran_starts))
        g_min = min(s for s, _, _ in exons)
        g_max = max(e for _, e, _ in exons)
        out.append((str(row["tran_id"]), str(row["strand"]), g_min, g_max, len(exons), exons))
    return out


def _apply_offsets_vectorised(
    pos5_arr: np.ndarray,
    strand_arr: np.ndarray,
    length_arr: np.ndarray,
    sample_id_arr: np.ndarray,
    offsets_by_sl: Dict[Tuple[int, int], int],
    default_offset: int,
) -> np.ndarray:
    """Vectorised P-site computation: p_site = pos5 ± offset[sample_id, length]."""
    # build offset lookup as two arrays (unique (sample,length) pairs)
    unique_sl = set(zip(sample_id_arr.tolist(), length_arr.tolist()))
    sl_to_off = {(s, l): offsets_by_sl.get((s, l), default_offset) for s, l in unique_sl}

    offsets_arr = np.array(
        [sl_to_off[(int(s), int(l))] for s, l in zip(sample_id_arr, length_arr)],
        dtype=np.int32,
    )
    # + strand: p_site = pos5 + offset; - strand: p_site = pos5 - offset
    return np.where(strand_arr, pos5_arr.astype(np.int64) + offsets_arr,
                    pos5_arr.astype(np.int64) - offsets_arr)


def _assign_psite_to_features(
    pos5_arr: np.ndarray,
    strand_arr: np.ndarray,
    features: List[Tuple[str, str, int, int, int, List[Tuple[int, int, int]]]],
    offsets_by_sl: Dict[Tuple[int, int], int],
    length_arr: np.ndarray,
    sample_id_arr: np.ndarray,
    count_arr: np.ndarray,
    default_offset: int = 15,
) -> pl.DataFrame:
    """Map each P-site read to (feature_id, tx_pos) using exon blocks.

    Applies offsets vectorised per unique (sample_id, length) pair, then
    does one sweep per feature (sorted-shard filter + exon lookup).

    Returns (feature_id, sample_id, length, tx_pos, count).
    """
    empty = pl.DataFrame(schema={
        "feature_id": pl.Utf8, "sample_id": pl.UInt16,
        "length": pl.UInt8, "tx_pos": pl.Int32, "count": pl.Float32,
    })
    if len(pos5_arr) == 0 or not features:
        return empty

    p_site = _apply_offsets_vectorised(
        pos5_arr, strand_arr, length_arr, sample_id_arr, offsets_by_sl, default_offset
    )

    feat_ids: List[str] = []
    samp_ids: List[int] = []
    lens: List[int] = []
    tx_positions: List[int] = []
    counts: List[float] = []

    # pre-separate + / - indices for fast strand filter
    plus_mask = strand_arr.astype(bool)
    minus_mask = ~plus_mask

    for fid, strand, g_min, g_max, _, exons in features:
        strand_mask = plus_mask if strand == "+" else minus_mask
        in_span = strand_mask & (p_site >= g_min) & (p_site < g_max)
        idxs = np.where(in_span)[0]
        if len(idxs) == 0:
            continue

        # vectorised exon lookup: build exon boundary arrays for this feature
        ex_starts = np.array([e[0] for e in exons], dtype=np.int64)
        ex_stops = np.array([e[1] for e in exons], dtype=np.int64)
        ex_tran = np.array([e[2] for e in exons], dtype=np.int32)

        ps_sub = p_site[idxs]  # (k,)
        # for each read, find first exon containing it
        # broadcast: (k, n_exons)
        in_exon = (ps_sub[:, None] >= ex_starts[None, :]) & (ps_sub[:, None] < ex_stops[None, :])
        exon_idx = np.argmax(in_exon, axis=1)  # (k,) index of first True per row
        hit_mask = in_exon[np.arange(len(idxs)), exon_idx]  # True if any exon matched

        hit_idxs = idxs[hit_mask]
        hit_exon = exon_idx[hit_mask]
        hit_ps = p_site[hit_idxs]

        tx_pos_arr = ex_tran[hit_exon] + (hit_ps - ex_starts[hit_exon]).astype(np.int32)

        feat_ids.extend([fid] * len(hit_idxs))
        samp_ids.extend(sample_id_arr[hit_idxs].tolist())
        lens.extend(length_arr[hit_idxs].tolist())
        tx_positions.extend(tx_pos_arr.tolist())
        counts.extend(count_arr[hit_idxs].tolist())

    if not feat_ids:
        return empty
    return pl.DataFrame({
        "feature_id": pl.Series(feat_ids, dtype=pl.Utf8),
        "sample_id": pl.Series(samp_ids, dtype=pl.UInt16),
        "length": pl.Series(lens, dtype=pl.UInt8),
        "tx_pos": pl.Series(tx_positions, dtype=pl.Int32),
        "count": pl.Series(counts, dtype=pl.Float32),
    })


def query_frame_rollup(
    index_dir: str | Path,
    cds_df: pl.DataFrame,
    offsets: Dict[Tuple[int, int], int],
    *,
    chroms: Optional[List[str]] = None,
    default_offset: int = 15,
) -> pl.DataFrame:
    """Query the P-site index → FrameRollup.

    Args:
        index_dir:      path to built P-site index.
        cds_df:         CDS blocks for the features to score.
        offsets:        {(sample_id, length): offset} calibrated offsets.
        chroms:         restrict to these chromosomes (default: all).
        default_offset: fallback offset for uncalibrated (sample, length) pairs.

    Returns:
        DataFrame with columns (sample_id, feature_id, length, n_reads,
        frame0, frame1, frame2).
    """
    index_dir = Path(index_dir)
    target_chroms = chroms or available_chroms(index_dir)

    parts: List[pl.DataFrame] = []
    for chrom in target_chroms:
        shard = _load_shard(index_dir, chrom)
        if shard is None:
            continue
        features = _exon_intervals_for_chrom(cds_df, chrom)
        if not features:
            continue

        # find max possible offset for pre-filter padding
        max_off = max(offsets.values(), default=default_offset) + 5
        g_min_all = min(g_min for _, _, g_min, _, _, _ in features) - max_off
        g_max_all = max(g_max for _, _, _, g_max, _, _ in features) + max_off

        # fast range filter using sorted pos5
        window = shard.filter(pl.col("pos5").is_between(g_min_all, g_max_all))
        if window.is_empty():
            continue

        assigned = _assign_psite_to_features(
            window["pos5"].to_numpy(),
            window["strand"].to_numpy(),
            features,
            offsets,
            window["length"].to_numpy(),
            window["sample_id"].to_numpy(),
            window["count"].to_numpy(),
            default_offset=default_offset,
        )
        if not assigned.is_empty():
            parts.append(assigned)

    if not parts:
        return pl.DataFrame(schema=_ROLLUP_SCHEMA)

    assigned_all = pl.concat(parts)
    rollup = (
        assigned_all.with_columns((pl.col("tx_pos") % 3).cast(pl.Int8).alias("frame"))
        .group_by(["feature_id", "sample_id", "length", "frame"])
        .agg(pl.col("count").sum())
        .pivot(on="frame", index=["feature_id", "sample_id", "length"], values="count",
               aggregate_function="sum")
    )

    for col in ["0", "1", "2"]:
        if col not in rollup.columns:
            rollup = rollup.with_columns(pl.lit(0.0).cast(pl.Float32).alias(col))

    rollup = (
        rollup.rename({c: f"frame{c}" for c in ["0", "1", "2"] if c in rollup.columns})
        .with_columns(
            (pl.col("frame0").fill_null(0.0) +
             pl.col("frame1").fill_null(0.0) +
             pl.col("frame2").fill_null(0.0)).alias("n_reads")
        )
    )
    return rollup.select(["sample_id", "feature_id", "length", "n_reads",
                          "frame0", "frame1", "frame2"])


def query_coverage_index(
    index_dir: str | Path,
    cds_df: pl.DataFrame,
    offsets: Dict[Tuple[int, int], int],
    feature_ids: Optional[List[str]] = None,
    *,
    chroms: Optional[List[str]] = None,
    default_offset: int = 15,
) -> pl.DataFrame:
    """Query the P-site index → per-position coverage.

    Returns (feature_id, sample_id, length, tx_pos, count).
    """
    index_dir = Path(index_dir)
    target_chroms = chroms or available_chroms(index_dir)
    fid_set = set(feature_ids) if feature_ids else None

    parts: List[pl.DataFrame] = []
    for chrom in target_chroms:
        shard = _load_shard(index_dir, chrom)
        if shard is None:
            continue
        features = _exon_intervals_for_chrom(cds_df, chrom)
        if fid_set:
            features = [f for f in features if f[0] in fid_set]
        if not features:
            continue

        max_off = max(offsets.values(), default=default_offset) + 5
        g_min_all = min(g_min for _, _, g_min, _, _, _ in features) - max_off
        g_max_all = max(g_max for _, _, _, g_max, _, _ in features) + max_off
        window = shard.filter(pl.col("pos5").is_between(g_min_all, g_max_all))
        if window.is_empty():
            continue

        assigned = _assign_psite_to_features(
            window["pos5"].to_numpy(), window["strand"].to_numpy(), features, offsets,
            window["length"].to_numpy(), window["sample_id"].to_numpy(),
            window["count"].to_numpy(), default_offset=default_offset,
        )
        if not assigned.is_empty():
            parts.append(assigned)

    if not parts:
        return pl.DataFrame(schema={
            "feature_id": pl.Utf8, "sample_id": pl.UInt16,
            "length": pl.UInt8, "tx_pos": pl.Int32, "count": pl.Float32,
        })

    return (
        pl.concat(parts)
        .group_by(["feature_id", "sample_id", "length", "tx_pos"])
        .agg(pl.col("count").sum())
    )

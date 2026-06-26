"""P-site index: build-once, query-many cohort cache.

Two-phase design
----------------
Phase 1 — calibrate (minutes, reuses existing map-reduce rollup):
    build_matrix_rollup() → (sample, length, phase0) rollup
    calibrate_offsets(target_frame=0) → per-(sample, length) P-site offsets
    rollup_to_periodicity() → per-sample QC (periodicity score, f0/f1/f2, …)
    Writes: index_dir/offsets.parquet, index_dir/qc_per_sample.parquet

Phase 2 — build index (hours, embarrassingly parallel):
    Reads calibrated offsets from Phase 1.
    Scans each partition BAM, joins count parquets.
    For every in-CDS read: p_site = pos5 ± offset[sample_id, length]  (baked in)
    Writes: index_dir/chrom=X/data.parquet  sorted by p_site

Query (O(log N + hits) per chrom, all features in one pass):
    Range-filter on p_site (uses Parquet row-group pushdown, no offset join)
    Sweep-line assign to (feature_id, tx_pos)
    Group → FrameRollup or full CoverageIndex

Storage layout::

    index_dir/
      offsets.parquet         (sample_id UInt16, length UInt8, offset Int8)
      qc_per_sample.parquet   (sample_id, periodicity_score, f0, f1, f2, …)
      samples.parquet         (sample_id UInt16, sample_name Utf8)
      chrom=chr1/data.parquet
      chrom=chr22/data.parquet
      …

Chrom shard schema (p_site-sorted):
    p_site     Int32    ← calibrated P-site position (offset already applied)
    strand     Boolean  (True = '+')
    length     UInt8    (useful for length-stratified profiles)
    sample_id  UInt16
    count      Float32
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
# Shard schema
# ---------------------------------------------------------------------------

_SHARD_SCHEMA = {
    "p_site": pl.Int32,
    "strand": pl.Boolean,
    "length": pl.UInt8,
    "sample_id": pl.UInt16,
    "count": pl.Float32,
}

_ROLLUP_SCHEMA = {
    "sample_id": pl.UInt16,
    "feature_id": pl.Utf8,
    "length": pl.UInt8,
    "n_reads": pl.Float32,
    "frame0": pl.Float32,
    "frame1": pl.Float32,
    "frame2": pl.Float32,
}

# ---------------------------------------------------------------------------
# Phase 1: calibrate offsets + QC from rollup
# ---------------------------------------------------------------------------


def calibrate_cohort(
    partition_dirs,
    cds_df: pl.DataFrame,
    *,
    n_workers: Optional[int] = None,
    chrom: Optional[str] = None,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Phase 1: compute per-sample calibrated offsets and QC stats.

    Runs build_matrix_rollup (fast, existing map-reduce) then calibrate_offsets
    and rollup_to_periodicity on the compact rollup — equivalent to running
    RiboMetric on every sample without touching individual BAM files.

    Returns:
        offsets_df: (sample_name, length, offset) — one row per (sample, length)
        qc_df:      per-sample QC table (periodicity_score, f0, f1, f2, …)
    """
    from .matrix_rollup import (
        build_matrix_rollup,
        calibrate_offsets,
        rollup_to_periodicity,
    )

    # If chrom-restricted, we still need a representative CDS set for calibration.
    # Use the supplied cds_df as-is (caller may have pre-filtered to target chrom,
    # or may pass the full annotation — both are valid).
    cal_cds = cds_df
    if chrom:
        cal_cds = cds_df.filter(pl.col("chr") == chrom)
        if cal_cds.is_empty():
            # try without "chr" prefix
            alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            cal_cds = cds_df.filter(pl.col("chr") == alt)
        if cal_cds.is_empty():
            log_info(
                f"calibrate_cohort: no CDS entries for chrom={chrom}; "
                "using full annotation for calibration"
            )
            cal_cds = cds_df

    log_info("Phase 1: building matrix rollup for offset calibration …")
    rollup = build_matrix_rollup(
        partition_dirs, cal_cds, n_workers=n_workers
    )

    log_info("Phase 1: calibrating per-(sample, length) P-site offsets …")
    offsets_dict = calibrate_offsets(rollup, target_frame=0)

    log_info("Phase 1: computing per-sample QC stats …")
    qc_df = rollup_to_periodicity(rollup, offsets_dict)

    # Build offsets_df: (sample_name Utf8, length Int64, offset Int64)
    if offsets_dict:
        offsets_df = pl.DataFrame(
            {
                "sample_name": pl.Series(
                    [s for s, _ in offsets_dict.keys()], dtype=pl.Utf8
                ),
                "length": pl.Series(
                    [int(L) for _, L in offsets_dict.keys()], dtype=pl.Int64
                ),
                "offset": pl.Series(list(offsets_dict.values()), dtype=pl.Int64),
            }
        )
    else:
        offsets_df = pl.DataFrame(
            schema={"sample_name": pl.Utf8, "length": pl.Int64, "offset": pl.Int64}
        )

    return offsets_df, qc_df


# ---------------------------------------------------------------------------
# Phase 2: build worker — bakes calibrated offset into p_site at build time
# ---------------------------------------------------------------------------

_BCTX: Dict = {}


def _build_worker_init(frame_ivs, bam_refs, sample_map, offset_lookup) -> None:
    """offset_lookup: {(sample_id: int, length: int): offset: int}"""
    _BCTX.update(
        frame_ivs=frame_ivs,
        bam_refs=bam_refs,
        sample_map=sample_map,
        offset_lookup=offset_lookup,
    )


def _build_worker(pdir_str: str) -> bytes:
    """Scan one partition → IPC bytes of (chrom, p_site, strand, length, sample_id, count)."""
    from collections import defaultdict

    from .matrix_qc import _assign_frames_sweep

    pdir = Path(pdir_str)
    bam = _discover_bam(pdir)
    empty_schema = {"chrom": pl.Utf8, **_SHARD_SCHEMA}
    empty = pl.DataFrame(schema=empty_schema)

    if bam is None:
        buf = io.BytesIO()
        empty.write_ipc(buf)
        return buf.getvalue()

    frame_ivs = _BCTX["frame_ivs"]
    bam_refs = _BCTX["bam_refs"]
    sample_map: pl.DataFrame = _BCTX["sample_map"]
    offset_lookup: Dict[Tuple[int, int], int] = _BCTX["offset_lookup"]
    default_offset: int = 15

    # scan BAM for all mapped reads, keep raw 5' end + metadata
    read_id_list: List[int] = []
    pos5_list: List[int] = []
    length_list: List[int] = []
    strand_list: List[bool] = []
    chrom_list: List[str] = []

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

    # filter to in-CDS reads
    n = len(read_id_list)
    in_cds = np.zeros(n, dtype=bool)
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
            "read_id": pl.Series(
                [read_id_list[i] for i in in_cds_idx], dtype=pl.UInt64
            ),
            "chrom": pl.Series(
                [chrom_list[i] for i in in_cds_idx], dtype=pl.Utf8
            ),
            "pos5": pl.Series(
                [pos5_list[i] for i in in_cds_idx], dtype=pl.Int32
            ),
            "strand": pl.Series(
                [strand_list[i] for i in in_cds_idx], dtype=pl.Boolean
            ),
            "length": pl.Series(
                [length_list[i] for i in in_cds_idx], dtype=pl.UInt8
            ),
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
    joined = (
        counts.join(aln, on="read_id", how="inner")
        .join(sample_map.select("sample_id"), on="sample_id", how="inner")
        .with_columns(pl.col("count").cast(pl.Float32))
    )

    if joined.is_empty():
        buf = io.BytesIO()
        empty.write_ipc(buf)
        return buf.getvalue()

    # bake calibrated offset into p_site at build time
    sid_arr = joined["sample_id"].to_numpy()
    len_arr = joined["length"].to_numpy()
    pos5_arr = joined["pos5"].to_numpy()
    strand_arr = joined["strand"].to_numpy()

    offset_arr = np.array(
        [offset_lookup.get((int(s), int(l)), default_offset)
         for s, l in zip(sid_arr, len_arr)],
        dtype=np.int32,
    )
    # + strand: p_site = pos5 + offset; - strand: p_site = pos5 - offset
    p_site_arr = np.where(strand_arr, pos5_arr.astype(np.int64) + offset_arr,
                          pos5_arr.astype(np.int64) - offset_arr).astype(np.int32)

    out = (
        joined.with_columns(pl.Series("p_site", p_site_arr, dtype=pl.Int32))
        .group_by(["chrom", "p_site", "strand", "length", "sample_id"])
        .agg(pl.col("count").sum())
        .with_columns(
            pl.col("p_site").cast(pl.Int32),
            pl.col("length").cast(pl.UInt8),
            pl.col("sample_id").cast(pl.UInt16),
            pl.col("count").cast(pl.Float32),
        )
    )

    buf = io.BytesIO()
    out.write_ipc(buf)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Phase 2: orchestrator
# ---------------------------------------------------------------------------


def build_psite_index(
    partition_dirs,
    cds_df: pl.DataFrame,
    out_dir: str | Path,
    offsets_df: Optional[pl.DataFrame] = None,
    *,
    n_workers: Optional[int] = None,
    chrom: Optional[str] = None,
    default_offset: int = 15,
) -> None:
    """Phase 2: build the P-site index using calibrated offsets.

    If offsets_df is None, reads from out_dir/offsets.parquet (written by
    calibrate_cohort / Phase 1).  Raises if neither is available.

    Args:
        partition_dirs: path to global_partitioned directory, or list of dirs.
        cds_df:         CDS blocks used to restrict the index to in-CDS reads.
        out_dir:        where to write (and read Phase 1 outputs from).
        offsets_df:     (sample_name, length, offset) table from Phase 1.
                        If None, loaded from out_dir/offsets.parquet.
        n_workers:      worker processes (default: cpu_count).
        chrom:          restrict output to one chromosome (for testing).
        default_offset: fallback for (sample, length) pairs not in offsets_df.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # load offsets if not supplied
    if offsets_df is None:
        off_path = out_dir / "offsets.parquet"
        if not off_path.exists():
            raise FileNotFoundError(
                f"No offsets.parquet in {out_dir}. Run calibrate_cohort first "
                "(or pass offsets_df directly)."
            )
        offsets_df = pl.read_parquet(off_path)

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
        # restrict CDS intervals to target chrom — reads on other chroms are skipped
        frame_ivs = {
            k: v for k, v in frame_ivs.items()
            if k[0] == chrom or k[0] == f"chr{chrom}"
        }

    mp0, manifest0 = _manifest(dirs[0])
    sample_map = _samples_df(mp0, manifest0).select(
        pl.col("sample_id").cast(pl.UInt16), pl.col("sample_name")
    )
    sample_map.write_parquet(out_dir / "samples.parquet")

    # build offset_lookup keyed by (sample_id: int, length: int)
    # join offsets_df (sample_name, length, offset) with sample_map to get sample_id
    offset_lookup: Dict[Tuple[int, int], int] = {}
    if not offsets_df.is_empty():
        merged = offsets_df.join(sample_map, on="sample_name", how="inner")
        for row in merged.iter_rows(named=True):
            offset_lookup[(int(row["sample_id"]), int(row["length"]))] = int(row["offset"])

    # write numeric offsets table for reference / querying
    if offset_lookup:
        pl.DataFrame(
            {
                "sample_id": pl.Series(
                    [k[0] for k in offset_lookup], dtype=pl.UInt16
                ),
                "length": pl.Series(
                    [k[1] for k in offset_lookup], dtype=pl.UInt8
                ),
                "offset": pl.Series(
                    list(offset_lookup.values()), dtype=pl.Int8
                ),
            }
        ).write_parquet(out_dir / "offsets_numeric.parquet")

    if n_workers is None:
        n_workers = min(len(dirs), os.cpu_count() or 1)

    log_info(
        f"Phase 2 build_psite_index: {len(dirs)} partitions, "
        f"{sum(len(v) for v in frame_ivs.values()):,} CDS intervals, "
        f"{len(offset_lookup):,} (sample, length) offsets, "
        f"{n_workers} worker(s) → {out_dir}"
    )

    parts: List[pl.DataFrame] = []
    dir_strs = [str(d) for d in dirs]

    if n_workers <= 1:
        _build_worker_init(frame_ivs, bam_refs, sample_map, offset_lookup)
        for ds in dir_strs:
            df = pl.read_ipc(io.BytesIO(_build_worker(ds)))
            if not df.is_empty():
                parts.append(df)
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            n_workers,
            initializer=_build_worker_init,
            initargs=(frame_ivs, bam_refs, sample_map, offset_lookup),
        ) as pool:
            for i, ipc in enumerate(
                pool.imap_unordered(_build_worker, dir_strs, chunksize=1)
            ):
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

    combined = (
        combined
        .group_by(["chrom", "p_site", "strand", "length", "sample_id"])
        .agg(pl.col("count").sum())
        .with_columns(
            pl.col("p_site").cast(pl.Int32),
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
            .sort("p_site")
            .write_parquet(
                shard_dir / "data.parquet", statistics=True, compression="zstd"
            )
        )

    log_info(f"build_psite_index: done — {out_dir}")


# ---------------------------------------------------------------------------
# Query helpers
# ---------------------------------------------------------------------------


def load_samples(index_dir: str | Path) -> pl.DataFrame:
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


def _exon_intervals_for_chrom(
    cds_df: pl.DataFrame, chrom: str
) -> List[Tuple[str, str, int, int, int, List[Tuple[int, int, int]]]]:
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
        out.append(
            (str(row["tran_id"]), str(row["strand"]), g_min, g_max, len(exons), exons)
        )
    return out


# ---------------------------------------------------------------------------
# Query: assign p_site reads to features (vectorised)
# ---------------------------------------------------------------------------


def _assign_psite_to_features(
    p_site_arr: np.ndarray,
    strand_arr: np.ndarray,
    features: List[Tuple[str, str, int, int, int, List[Tuple[int, int, int]]]],
    length_arr: np.ndarray,
    sample_id_arr: np.ndarray,
    count_arr: np.ndarray,
) -> pl.DataFrame:
    """Map pre-offset-corrected p_site reads to (feature_id, tx_pos).

    Since offsets are baked into p_site at index build time, no offset
    table is needed here — just range filter + exon projection.

    Returns (feature_id, sample_id, length, tx_pos, count).
    """
    empty = pl.DataFrame(
        schema={
            "feature_id": pl.Utf8,
            "sample_id": pl.UInt16,
            "length": pl.UInt8,
            "tx_pos": pl.Int32,
            "count": pl.Float32,
        }
    )
    if len(p_site_arr) == 0 or not features:
        return empty

    feat_ids: List[str] = []
    samp_ids: List[int] = []
    lens: List[int] = []
    tx_positions: List[int] = []
    counts: List[float] = []

    plus_mask = strand_arr.astype(bool)
    minus_mask = ~plus_mask

    for fid, strand, g_min, g_max, _, exons in features:
        strand_mask = plus_mask if strand == "+" else minus_mask
        in_span = strand_mask & (p_site_arr >= g_min) & (p_site_arr < g_max)
        idxs = np.where(in_span)[0]
        if len(idxs) == 0:
            continue

        ex_starts = np.array([e[0] for e in exons], dtype=np.int64)
        ex_stops = np.array([e[1] for e in exons], dtype=np.int64)
        ex_tran = np.array([e[2] for e in exons], dtype=np.int32)

        ps_sub = p_site_arr[idxs]
        in_exon = (ps_sub[:, None] >= ex_starts[None, :]) & (
            ps_sub[:, None] < ex_stops[None, :]
        )
        exon_idx = np.argmax(in_exon, axis=1)
        hit_mask = in_exon[np.arange(len(idxs)), exon_idx]

        hit_idxs = idxs[hit_mask]
        hit_exon = exon_idx[hit_mask]
        hit_ps = p_site_arr[hit_idxs]

        tx_pos_arr = ex_tran[hit_exon] + (hit_ps - ex_starts[hit_exon]).astype(
            np.int32
        )

        feat_ids.extend([fid] * len(hit_idxs))
        samp_ids.extend(sample_id_arr[hit_idxs].tolist())
        lens.extend(length_arr[hit_idxs].tolist())
        tx_positions.extend(tx_pos_arr.tolist())
        counts.extend(count_arr[hit_idxs].tolist())

    if not feat_ids:
        return empty
    return pl.DataFrame(
        {
            "feature_id": pl.Series(feat_ids, dtype=pl.Utf8),
            "sample_id": pl.Series(samp_ids, dtype=pl.UInt16),
            "length": pl.Series(lens, dtype=pl.UInt8),
            "tx_pos": pl.Series(tx_positions, dtype=pl.Int32),
            "count": pl.Series(counts, dtype=pl.Float32),
        }
    )


# ---------------------------------------------------------------------------
# Public query API
# ---------------------------------------------------------------------------


def _query_chrom(
    index_dir: Path,
    chrom: str,
    features: List[Tuple],
) -> Optional[pl.DataFrame]:
    shard = _load_shard(index_dir, chrom)
    if shard is None or not features:
        return None

    # p_site is already offset-corrected — range filter is a simple BETWEEN
    g_min_all = min(g_min for _, _, g_min, _, _, _ in features)
    g_max_all = max(g_max for _, _, _, g_max, _, _ in features)
    window = shard.filter(pl.col("p_site").is_between(g_min_all, g_max_all))
    if window.is_empty():
        return None

    return _assign_psite_to_features(
        window["p_site"].to_numpy(),
        window["strand"].to_numpy(),
        features,
        window["length"].to_numpy(),
        window["sample_id"].to_numpy(),
        window["count"].to_numpy(),
    )


def query_frame_rollup(
    index_dir: str | Path,
    cds_df: pl.DataFrame,
    *,
    chroms: Optional[List[str]] = None,
) -> pl.DataFrame:
    """Query the P-site index → FrameRollup.

    Offsets are already baked into the index — no offset table needed.

    Returns (sample_id, feature_id, length, n_reads, frame0, frame1, frame2).
    """
    index_dir = Path(index_dir)
    target_chroms = chroms or available_chroms(index_dir)

    parts: List[pl.DataFrame] = []
    for chrom in target_chroms:
        features = _exon_intervals_for_chrom(cds_df, chrom)
        result = _query_chrom(index_dir, chrom, features)
        if result is not None and not result.is_empty():
            parts.append(result)

    if not parts:
        return pl.DataFrame(schema=_ROLLUP_SCHEMA)

    assigned_all = pl.concat(parts)
    rollup = (
        assigned_all.with_columns(
            (pl.col("tx_pos") % 3).cast(pl.Int8).alias("frame")
        )
        .group_by(["feature_id", "sample_id", "length", "frame"])
        .agg(pl.col("count").sum())
        .pivot(
            on="frame",
            index=["feature_id", "sample_id", "length"],
            values="count",
            aggregate_function="sum",
        )
    )

    for col in ["0", "1", "2"]:
        if col not in rollup.columns:
            rollup = rollup.with_columns(pl.lit(0.0).cast(pl.Float32).alias(col))

    rollup = rollup.rename(
        {c: f"frame{c}" for c in ["0", "1", "2"] if c in rollup.columns}
    ).with_columns(
        (
            pl.col("frame0").fill_null(0.0)
            + pl.col("frame1").fill_null(0.0)
            + pl.col("frame2").fill_null(0.0)
        ).alias("n_reads")
    )
    return rollup.select(
        ["sample_id", "feature_id", "length", "n_reads", "frame0", "frame1", "frame2"]
    )


def query_coverage_index(
    index_dir: str | Path,
    cds_df: pl.DataFrame,
    feature_ids: Optional[List[str]] = None,
    *,
    chroms: Optional[List[str]] = None,
) -> pl.DataFrame:
    """Query the P-site index → per-position coverage.

    Returns (feature_id, sample_id, length, tx_pos, count).
    """
    index_dir = Path(index_dir)
    target_chroms = chroms or available_chroms(index_dir)
    fid_set = set(feature_ids) if feature_ids else None

    parts: List[pl.DataFrame] = []
    for chrom in target_chroms:
        features = _exon_intervals_for_chrom(cds_df, chrom)
        if fid_set:
            features = [f for f in features if f[0] in fid_set]
        result = _query_chrom(index_dir, chrom, features)
        if result is not None and not result.is_empty():
            parts.append(result)

    if not parts:
        return pl.DataFrame(
            schema={
                "feature_id": pl.Utf8,
                "sample_id": pl.UInt16,
                "length": pl.UInt8,
                "tx_pos": pl.Int32,
                "count": pl.Float32,
            }
        )
    return (
        pl.concat(parts)
        .group_by(["feature_id", "sample_id", "length", "tx_pos"])
        .agg(pl.col("count").sum())
    )

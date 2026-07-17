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

from ..io.annotation import build_gene_spans
from ..io.bam import normalise_chrom as _normalise_chrom
from ..io.matrix import _count_parquets, _discover_bam, _manifest, _parse_read_id, _samples_df
from .qc import _bam_chroms, _build_frame_intervals, fast_reads_qc
from ..utils.logging import log_info

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
    """Phase 1: compute per-sample calibrated offsets and comprehensive QC.

    Runs three parallel streams and merges results:
      1. Matrix rollup → per-(sample,length) frame counts → periodicity QC + offsets
      2. fast_reads_qc across all partitions → length distribution metrics
      3. fast_reads_qc (same scan) → ligation bias metrics

    The combined QC mirrors RiboMetric's per-sample output without BAMs:
      periodicity_score, f0/f1/f2, n_cds_reads, recommended_offsets
      total_reads, peak_length, mean_length, rpf_28_32_prop
      rld_IQR_metric, rld_CV_metric, rld_normality_metric,
      rld_max_prop_metric, rld_bimodality
      ligation_bias_KL_5p/3p, ligation_bias_score_5p/3p,
      ligation_bias_max_abs_5p/3p
      prop_cds, recommended_lengths, n_recommended,
      recommended_read_proportion, library_type

    Returns:
        offsets_df: (sample_name, length, offset) — one row per (sample, length)
        qc_df:      comprehensive per-sample QC table
    """
    import json as _json

    from .rollup import (
        build_matrix_rollup,
        calibrate_offsets,
        rollup_to_periodicity,
    )
    from ..qc import (
        _classify_library,
        _length_distribution_metrics,
        _ligation_bias_metrics,
        _recommend_read_lengths,
    )

    # Resolve partition_dirs to a list of Path objects
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        pdirs: List[Path] = sorted(d for d in parent.iterdir() if d.is_dir())
    else:
        pdirs = [Path(d) for d in partition_dirs]

    # Phase 1a: offset calibration via matrix rollup
    cal_cds = cds_df
    if chrom:
        cal_cds = cds_df.filter(pl.col("chr") == chrom)
        if cal_cds.is_empty():
            alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            cal_cds = cds_df.filter(pl.col("chr") == alt)
        if cal_cds.is_empty():
            log_info(
                f"calibrate_cohort: no CDS entries for chrom={chrom}; "
                "using full annotation for calibration"
            )
            cal_cds = cds_df

    log_info("Phase 1: building matrix rollup for offset calibration …")
    rollup = build_matrix_rollup(pdirs, cal_cds, n_workers=n_workers)

    log_info("Phase 1: calibrating per-(sample, length) P-site offsets …")
    offsets_dict = calibrate_offsets(rollup, target_frame=0)

    log_info("Phase 1: computing periodicity QC …")
    qc_df = rollup_to_periodicity(rollup, offsets_dict)

    # Phase 1b: length distribution + ligation bias across all partitions
    log_info(f"Phase 1: scanning {len(pdirs)} partition(s) for reads QC …")
    len_parts: List[pl.DataFrame] = []
    dinuc_parts: List[pl.DataFrame] = []
    for pdir in pdirs:
        try:
            ld, dd = fast_reads_qc(pdir)
            if not ld.is_empty():
                len_parts.append(ld)
            if not dd.is_empty():
                dinuc_parts.append(dd)
        except Exception as e:
            log_info(f"  skipping {pdir.name}: {e}")

    reads_qc_rows: Dict[str, dict] = {}

    if len_parts:
        length_all = pl.concat(len_parts)
        # Aggregate length histogram per sample across partitions
        length_hist = length_all.group_by(["sample_name", "length"]).agg(pl.col("count").sum())
        for (sname,), grp in length_hist.group_by("sample_name"):
            lc = {int(row["length"]): float(row["count"]) for row in grp.iter_rows(named=True)}
            reads_qc_rows.setdefault(str(sname), {}).update(_length_distribution_metrics(lc))
            reads_qc_rows[str(sname)]["_rld"] = lc

    if dinuc_parts:
        dinuc_all = pl.concat(dinuc_parts)
        # Cohort-wide background: aggregate across all samples per (end, dinuc)
        cohort_bg = dinuc_all.group_by(["end", "dinuc"]).agg(pl.col("count").sum())

        def _to_freq(df: pl.DataFrame) -> Dict[str, float]:
            total = float(df["count"].sum()) or 1.0
            return {row["dinuc"]: row["count"] / total for row in df.iter_rows(named=True)}

        bg5 = _to_freq(cohort_bg.filter(pl.col("end") == "5p"))
        bg3 = _to_freq(cohort_bg.filter(pl.col("end") == "3p"))

        # Per-sample dinucleotide aggregation
        dinuc_agg = dinuc_all.group_by(["sample_name", "end", "dinuc"]).agg(pl.col("count").sum())
        for (sname,), grp in dinuc_agg.group_by(["sample_name"]):
            five = {
                row["dinuc"]: float(row["count"])
                for row in grp.filter(pl.col("end") == "5p").iter_rows(named=True)
            }
            three = {
                row["dinuc"]: float(row["count"])
                for row in grp.filter(pl.col("end") == "3p").iter_rows(named=True)
            }
            reads_qc_rows.setdefault(str(sname), {}).update(
                _ligation_bias_metrics(five, three, bg5, bg3)
            )

    # Phase 1c: derive composite metrics and merge everything
    if reads_qc_rows:
        perio_lookup: Dict[str, dict] = {}
        for row in qc_df.iter_rows(named=True):
            sid = str(row["sample_id"])
            rfd = _json.loads(row.get("read_frame_distribution") or "{}")
            # rfd format: {length_str: [f0_count, f1_count, f2_count]}
            rfd_int: Dict[int, Dict[int, float]] = {
                int(L): {i: float(c) for i, c in enumerate(fd)} for L, fd in rfd.items()
            }
            offsets_for_sample = {
                int(L): int(off) for (sn, L), off in offsets_dict.items() if sn == sid
            }
            perio_lookup[sid] = {
                "rfd": rfd_int,
                "n_cds_reads": row.get("n_cds_reads", 0),
                "offsets": offsets_for_sample,
            }

        extra_rows = []
        for sname, stats in reads_qc_rows.items():
            rld = stats.pop("_rld", {})
            total = stats.get("total_reads", 0.0)
            pinfo = perio_lookup.get(sname, {})
            n_cds = float(pinfo.get("n_cds_reads", 0))
            prop_cds = n_cds / total if total > 0 else 0.0

            rec = _recommend_read_lengths(
                pinfo.get("rfd", {}),
                rld,
                offsets=pinfo.get("offsets"),
            )
            stats["prop_cds"] = prop_cds
            stats["recommended_lengths"] = _json.dumps(rec["recommended_lengths"])
            stats["n_recommended"] = rec["n_recommended"]
            stats["recommended_read_proportion"] = rec["recommended_read_proportion"]
            extra_rows.append({"sample_id": sname, **stats})

        extra_df = pl.from_dicts(extra_rows)
        qc_df = qc_df.join(extra_df, on="sample_id", how="left")

        # library_type requires periodicity_score + prop_cds
        if "periodicity_score" in qc_df.columns and "prop_cds" in qc_df.columns:
            lib_types = [
                _classify_library(float(r["periodicity_score"] or 0.0), float(r["prop_cds"] or 0.0))
                for r in qc_df.select(["periodicity_score", "prop_cds"]).iter_rows(named=True)
            ]
            qc_df = qc_df.with_columns(pl.Series("library_type", lib_types, dtype=pl.Utf8))

    # Build offsets_df
    if offsets_dict:
        offsets_df = pl.DataFrame(
            {
                "sample_name": pl.Series([s for s, _ in offsets_dict.keys()], dtype=pl.Utf8),
                "length": pl.Series([int(L) for _, L in offsets_dict.keys()], dtype=pl.Int64),
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

    from .qc import _assign_frames_sweep

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
        [offset_lookup.get((int(s), int(l)), default_offset) for s, l in zip(sid_arr, len_arr)],
        dtype=np.int32,
    )
    # + strand: p_site = pos5 + offset; - strand: p_site = pos5 - offset
    p_site_arr = np.where(
        strand_arr, pos5_arr.astype(np.int64) + offset_arr, pos5_arr.astype(np.int64) - offset_arr
    ).astype(np.int32)

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
        frame_ivs = {k: v for k, v in frame_ivs.items() if k[0] == chrom or k[0] == f"chr{chrom}"}

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
                "sample_id": pl.Series([k[0] for k in offset_lookup], dtype=pl.UInt16),
                "length": pl.Series([k[1] for k in offset_lookup], dtype=pl.UInt8),
                "offset": pl.Series(list(offset_lookup.values()), dtype=pl.Int8),
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

    combined = (
        combined.group_by(["chrom", "p_site", "strand", "length", "sample_id"])
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
            .write_parquet(shard_dir / "data.parquet", statistics=True, compression="zstd")
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
        out.append((str(row["tran_id"]), str(row["strand"]), g_min, g_max, len(exons), exons))
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
        in_exon = (ps_sub[:, None] >= ex_starts[None, :]) & (ps_sub[:, None] < ex_stops[None, :])
        exon_idx = np.argmax(in_exon, axis=1)
        hit_mask = in_exon[np.arange(len(idxs)), exon_idx]

        hit_idxs = idxs[hit_mask]
        hit_exon = exon_idx[hit_mask]
        hit_ps = p_site_arr[hit_idxs]

        tx_pos_arr = ex_tran[hit_exon] + (hit_ps - ex_starts[hit_exon]).astype(np.int32)

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
        assigned_all.with_columns((pl.col("tx_pos") % 3).cast(pl.Int8).alias("frame"))
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


# ---------------------------------------------------------------------------
# Query: ordinary GENOMIC coverage (for MatrixProvider / CoverageProvider)
# ---------------------------------------------------------------------------


def query_genomic_coverage(
    index_dir: str | Path,
    regions: List[Tuple[str, int, int]],
    *,
    site: str = "P",
    sample_names: Optional[List[str]] = None,
    group_level: str = "aggregate",
) -> pl.DataFrame:
    """Query the P-site index as ordinary genomic coverage (pos, count[, strand]).

    Unlike query_frame_rollup/query_coverage_index above (which re-project
    into transcript/feature space), this returns coverage keyed by GENOMIC
    position — the CoverageProvider.coverage() contract (coverage/base.py) —
    so the index can be read like any other coverage source and scored
    through workflows._score_events_over_provider (the same pipeline used
    for BAM/bigwig sources: splice-aware flanks, junction wiring, mappability
    annotation, all unmodified).

    Offsets are already baked into `p_site` at index-build time (see module
    docstring); site="A" applies the same +3 nt shift MatrixProvider already
    applies for its flat-ref_offset mode. site="P" needs no shift (the index
    IS P-site).

    If `index_dir/usable_sample_lengths.parquet` exists, (sample, length)
    pairs not listed there are excluded — the same QC-gate semantics
    workflows.score_matrix_rollup_workflow already applies for the
    FrameRollup path (its `psite_index_dir` branch), kept here for parity so
    switching a caller from FrameRollup to this path doesn't silently
    reintroduce known-bad (sample, length) combinations.

    Parameters
    ----------
    regions      : (chrom, start, end) tuples — same shape as
                   matrix_rollup.region_coverage takes.
    site         : "P" (initiation, default) or "A" (elongation/termination).
    sample_names : restrict to these sample(s) (None = all).
    group_level  : "aggregate" (sum across samples) or "sample" (per-sample rows).

    Returns
    -------
    DataFrame: pos (Int64), count (Float64), strand (Int64, ±1)
    [, sample_name (Utf8) if group_level="sample"].
    """
    if site not in {"P", "A"}:
        raise ValueError(f"site must be 'P' or 'A', got {site!r}")
    if group_level not in {"aggregate", "sample"}:
        raise ValueError(f"group_level must be 'aggregate' or 'sample', got {group_level!r}")

    index_dir = Path(index_dir)
    a_shift = 3 if site == "A" else 0

    sample_map = load_samples(index_dir)
    allowed_ids: Optional[set] = None
    if sample_names:
        allowed_ids = set(
            sample_map.filter(pl.col("sample_name").is_in(sample_names))["sample_id"].to_list()
        )

    usable_df: Optional[pl.DataFrame] = None
    usable_path = index_dir / "usable_sample_lengths.parquet"
    if usable_path.exists():
        usable_df = (
            pl.read_parquet(usable_path)
            .rename({"sample_id": "sample_name"})
            .select(["sample_name", pl.col("length").cast(pl.Int64)])
            .join(sample_map, on="sample_name", how="inner")
            .select(["sample_id", "length"])
            .with_columns(pl.lit(True).alias("_usable"))
        )

    by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, start, end in regions:
        by_chrom.setdefault(chrom, []).append((start - a_shift, end - a_shift))

    schema: Dict[str, type] = {"pos": pl.Int64, "count": pl.Float64, "strand": pl.Int64}
    if group_level == "sample":
        schema["sample_name"] = pl.Utf8

    parts: List[pl.DataFrame] = []
    for chrom, wins in by_chrom.items():
        shard = _load_shard(index_dir, chrom)
        if shard is None:
            continue
        if allowed_ids is not None:
            shard = shard.filter(pl.col("sample_id").is_in(list(allowed_ids)))
        if usable_df is not None:
            shard = (
                shard.join(usable_df, on=["sample_id", "length"], how="left")
                .filter(pl.col("_usable").fill_null(False))
                .drop("_usable")
            )
        if shard.is_empty():
            continue

        mask = pl.lit(False)
        for s, e in wins:
            mask = mask | ((pl.col("p_site") >= s) & (pl.col("p_site") < e))
        window = shard.filter(mask)
        if window.is_empty():
            continue

        window = window.with_columns(
            (pl.col("p_site") + a_shift).cast(pl.Int64).alias("pos"),
            pl.when(pl.col("strand")).then(1).otherwise(-1).cast(pl.Int64).alias("strand"),
        )
        group_cols = ["pos", "strand"]
        if group_level == "sample":
            window = window.join(sample_map, on="sample_id", how="left")
            group_cols.append("sample_name")
        agg = window.group_by(group_cols).agg(pl.col("count").sum().cast(pl.Float64).alias("count"))
        parts.append(agg.select(list(schema.keys())))

    if not parts:
        return pl.DataFrame(schema=schema)
    return pl.concat(parts).sort("pos")

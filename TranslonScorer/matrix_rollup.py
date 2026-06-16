"""Scalable matrix periodicity: offset-independent alignment index + map-reduce
frame rollup, with multimapping strategies (unique / guilty-by-association).

Design (see translonscorer/docs/matrix_query_engine.md):

  Phase A — build_alignment_index (parallel, offset-free)
      Scan every partition BAM once.  For each *in-CDS alignment* store
      (read_id, length, strand, phase0, gene_code, n_align), where
        phase0 = frame at offset 0   (sample- and offset-independent)
        n_align = total alignment records for the read (NH) → unique iff == 1
      All candidate loci of a multimapper are retained (no last-wins).

  Phase B — tabulate_rollup (vectorised, associative)
      Join counts × index × samples and reduce to a small rollup keyed
      (sample, length, strand, phase0).  multimap_mode:
        "unique" → keep n_align == 1 reads only.
        "gba"    → unique reads at weight 1 build per-(sample, gene) density;
                   each multimapper's count is split across its CDS candidates
                   ∝ that sample's unique density in the candidate's gene.

  Phase C — calibrate_offsets + rollup_to_periodicity (central, tiny)
      frame = (phase0 ± offset) % 3 is analytic, so the per-(sample, length)
      offset is calibrated on the reduced rollup (which already sums all
      partitions) and applied without rescanning.

Per-query cost is O(rollup size) = O(groups × lengths × phases), independent of
nnz — the property needed for 6k-sample cohorts.
"""

from __future__ import annotations

import collections
import heapq
import io
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
import pysam

from .matrix_qc import (
    _assign_frames_sweep,
    _bam_chroms,
    _build_frame_intervals,
    _count_parquets,
    _discover_bam,
    _manifest,
    _normalise_chrom,
    _parse_read_id,
    _ribometric_frame_scores,
    _samples_df,
)
from .utils.logging import log_info

_INDEX_SCHEMA = {
    "read_id": pl.UInt64,
    "length": pl.Int64,
    "strand": pl.Utf8,
    "phase0": pl.Int8,
    "gene_code": pl.Int32,
    "n_align": pl.Int32,
}


# ---------------------------------------------------------------------------
# Annotation: CDS blocks (with gene_id) + gene spans  [shim — moved to io/annotation.py]
# ---------------------------------------------------------------------------
from TranslonScorer.io.annotation import build_cds_blocks, build_gene_spans  # noqa: F401, E402


def _assign_region_sweep(
    ids: List[int], asites: np.ndarray, intervals: List[Tuple[int, int, int]]
) -> Dict[int, int]:
    """First-containing-interval payload per id (sweep-line). intervals: (start, stop, payload)."""
    if not intervals or len(ids) == 0:
        return {}
    order = np.argsort(asites)
    sa = asites[order].tolist()
    si = [ids[i] for i in order]
    res: Dict[int, int] = {}
    active: list = []
    ptr = 0
    n = len(intervals)
    for asite, rid in zip(sa, si):
        while ptr < n and intervals[ptr][0] <= asite:
            es, ee, pay = intervals[ptr]
            heapq.heappush(active, (ee, pay))
            ptr += 1
        while active and active[0][0] <= asite:
            heapq.heappop(active)
        if active:
            res[rid] = active[0][1]
    return res


# ---------------------------------------------------------------------------
# Phase A: per-partition alignment scan (all candidate loci, offset 0 phase)
# ---------------------------------------------------------------------------


def _scan_partition_alignments(
    bam_path, frame_ivs, gene_ivs, bam_refs, ref_offset: int
) -> pl.DataFrame:
    aln_rid: List[int] = []
    aln_len: List[int] = []
    aln_strand: List[str] = []
    aln_chrom: List[str] = []
    aln_asite: List[int] = []
    rec_count: collections.Counter = collections.Counter()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.reference_name is None:
                continue
            rid = _parse_read_id(rec.query_name)
            if rid is None:
                continue
            length = int(rec.query_length or 0)
            if length == 0:
                continue
            rec_count[rid] += 1
            strand = "-" if rec.is_reverse else "+"
            asite = (
                int(rec.reference_end) - 1 - ref_offset
                if rec.is_reverse
                else int(rec.reference_start) + ref_offset
            )
            aln_rid.append(rid)
            aln_len.append(length)
            aln_strand.append(strand)
            aln_chrom.append(_normalise_chrom(rec.reference_name, bam_refs) or rec.reference_name)
            aln_asite.append(asite)

    n = len(aln_rid)
    if n == 0:
        return pl.DataFrame(schema=_INDEX_SCHEMA)

    frame: List[Optional[int]] = [None] * n
    gene: List[int] = [-1] * n
    groups: Dict[Tuple[str, str], Tuple[List[int], List[int]]] = collections.defaultdict(
        lambda: ([], [])
    )
    for i in range(n):
        groups[(aln_chrom[i], aln_strand[i])][0].append(i)
        groups[(aln_chrom[i], aln_strand[i])][1].append(aln_asite[i])

    for (chrom, strand), (idxs, asites) in groups.items():
        arr = np.array(asites, dtype=np.int64)
        fiv = frame_ivs.get((chrom, strand))
        if fiv:
            for k, v in _assign_frames_sweep(idxs, arr, strand, fiv).items():
                frame[k] = v
        giv = gene_ivs.get((chrom, strand))
        if giv:
            for k, v in _assign_region_sweep(idxs, arr, giv).items():
                gene[k] = v

    out_rid, out_len, out_strand, out_phase0, out_gene, out_nalign = [], [], [], [], [], []
    for i in range(n):
        fr = frame[i]
        if fr is None:  # not in CDS — not a frame candidate
            continue
        strand = aln_strand[i]
        phase0 = (fr - ref_offset) % 3 if strand == "+" else (fr + ref_offset) % 3
        out_rid.append(aln_rid[i])
        out_len.append(aln_len[i])
        out_strand.append(strand)
        out_phase0.append(phase0)
        out_gene.append(gene[i])
        out_nalign.append(rec_count[aln_rid[i]])

    return pl.DataFrame(
        {
            "read_id": pl.Series(out_rid, dtype=pl.UInt64),
            "length": pl.Series(out_len, dtype=pl.Int64),
            "strand": pl.Series(out_strand, dtype=pl.Utf8),
            "phase0": pl.Series(out_phase0, dtype=pl.Int8),
            "gene_code": pl.Series(out_gene, dtype=pl.Int32),
            "n_align": pl.Series(out_nalign, dtype=pl.Int32),
        }
    )


def _scan_partition_alignments_oxbow(
    bam_path, frame_ivs, gene_ivs, bam_refs, ref_offset: int
) -> pl.DataFrame:
    """Vectorised oxbow reader — replaces the per-record pysam Python loop.

    Coordinate mapping (verified vs pysam): reference_start = pos - 1,
    reference_end = end (oxbow's ``end`` already accounts for N/D in the CIGAR,
    so no CIGAR parsing is needed). length = len(seq). Produces the same output
    schema as ``_scan_partition_alignments``.
    """
    import oxbow as ox

    df = pl.read_ipc(io.BytesIO(ox.read_bam(str(bam_path))))
    if df.is_empty():
        return pl.DataFrame(schema=_INDEX_SCHEMA)

    df = (
        df.filter((pl.col("flag") & 4) == 0)  # mapped
        .with_columns(
            [
                pl.col("qname").str.extract(r"read_(\d+)").cast(pl.UInt64).alias("read_id"),
                pl.col("seq").str.len_chars().cast(pl.Int64).alias("length"),
                pl.when((pl.col("flag") & 16) != 0)
                .then(pl.lit("-"))
                .otherwise(pl.lit("+"))
                .alias("strand"),
                pl.col("rname").cast(pl.Utf8).alias("chrom"),
            ]
        )
        .filter(pl.col("read_id").is_not_null() & (pl.col("length") > 0))
        .with_columns(
            pl.when(pl.col("strand") == "-")
            .then(pl.col("end") - 1 - ref_offset)
            .otherwise(pl.col("pos") - 1 + ref_offset)
            .cast(pl.Int64)
            .alias("asite")
        )
        .with_columns(pl.len().over("read_id").cast(pl.Int32).alias("n_align"))
        .with_row_index("aln_idx")
        .select(["aln_idx", "read_id", "length", "strand", "chrom", "asite", "n_align"])
    )

    n = df.height
    frame = np.full(n, -128, dtype=np.int16)  # sentinel = not in CDS
    gene = np.full(n, -1, dtype=np.int32)
    for (chrom, strand), grp in df.group_by(["chrom", "strand"]):
        idxs = grp.get_column("aln_idx").to_list()
        asites = grp.get_column("asite").to_numpy()
        fiv = frame_ivs.get((str(chrom), str(strand)))
        if fiv:
            for k, v in _assign_frames_sweep(idxs, asites, str(strand), fiv).items():
                frame[k] = v
        giv = gene_ivs.get((str(chrom), str(strand)))
        if giv:
            for k, v in _assign_region_sweep(idxs, asites, giv).items():
                gene[k] = v

    out = (
        df.with_columns(
            [
                pl.Series("frame", frame),
                pl.Series("gene_code", gene),
            ]
        )
        .filter(pl.col("frame") != -128)
        .with_columns(
            pl.when(pl.col("strand") == "+")
            .then((pl.col("frame") - ref_offset) % 3)
            .otherwise((pl.col("frame") + ref_offset) % 3)
            .cast(pl.Int8)
            .alias("phase0")
        )
        .select(["read_id", "length", "strand", "phase0", "gene_code", "n_align"])
    )
    return out


def _scan_dispatch(bam_path, frame_ivs, gene_ivs, bam_refs, ref_offset, reader):
    if reader == "oxbow":
        return _scan_partition_alignments_oxbow(bam_path, frame_ivs, gene_ivs, bam_refs, ref_offset)
    return _scan_partition_alignments(bam_path, frame_ivs, gene_ivs, bam_refs, ref_offset)


_CTX: Dict = {}


def _worker_init(frame_ivs, gene_ivs, bam_refs, ref_offset, reader) -> None:
    _CTX.update(
        frame_ivs=frame_ivs,
        gene_ivs=gene_ivs,
        bam_refs=bam_refs,
        ref_offset=ref_offset,
        reader=reader,
    )


def _worker_scan(bam_path_str: str) -> bytes:
    df = _scan_dispatch(
        bam_path_str,
        _CTX["frame_ivs"],
        _CTX["gene_ivs"],
        _CTX["bam_refs"],
        _CTX["ref_offset"],
        _CTX["reader"],
    )
    buf = io.BytesIO()
    df.write_ipc(buf)
    return buf.getvalue()


def build_alignment_index(
    partition_dirs,
    cds_df: pl.DataFrame,
    *,
    ref_offset: int = 15,
    n_workers: Optional[int] = None,
    out_path=None,
    reader: str = "pysam",
) -> pl.DataFrame:
    """Phase A, once: per in-CDS alignment (read_id, length, strand, phase0,
    gene_code, n_align).  ref_offset sets CDS membership (kept near the global
    default so start-proximal reads survive); phase0 is offset-independent."""
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir() and _discover_bam(d))
    else:
        dirs = [Path(d) for d in partition_dirs if _discover_bam(Path(d))]
    if not dirs:
        return pl.DataFrame(schema=_INDEX_SCHEMA)

    bam_refs = _bam_chroms(_discover_bam(dirs[0]))
    frame_ivs = _build_frame_intervals(cds_df, cds_df, bam_refs)
    gene_ivs, _gene_id_of = build_gene_spans(cds_df)
    if n_workers is None:
        n_workers = min(len(dirs), os.cpu_count() or 1)
    log_info(
        f"alignment index: {len(dirs)} partitions, "
        f"{sum(len(v) for v in frame_ivs.values()):,} CDS intervals, "
        f"{sum(len(v) for v in gene_ivs.values()):,} gene spans, {n_workers} worker(s)"
    )

    bam_paths = [str(_discover_bam(d)) for d in dirs]
    parts: List[pl.DataFrame] = []
    if n_workers <= 1:
        for i, bp in enumerate(bam_paths):
            df = _scan_dispatch(bp, frame_ivs, gene_ivs, bam_refs, ref_offset, reader)
            if not df.is_empty():
                parts.append(df)
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            n_workers,
            initializer=_worker_init,
            initargs=(frame_ivs, gene_ivs, bam_refs, ref_offset, reader),
        ) as pool:
            for i, ipc in enumerate(pool.imap_unordered(_worker_scan, bam_paths, chunksize=1)):
                df = pl.read_ipc(io.BytesIO(ipc))
                if not df.is_empty():
                    parts.append(df)
                if (i + 1) % 32 == 0:
                    log_info(f"  scanned {i + 1}/{len(dirs)} partitions")

    idx = pl.concat(parts) if parts else pl.DataFrame(schema=_INDEX_SCHEMA)
    if out_path is not None:
        idx.write_parquet(str(out_path))
        log_info(f"alignment index written: {out_path} ({idx.height:,} rows)")
    return idx


# ---------------------------------------------------------------------------
# Fold-in: scan + count-join inside the workers (map-reduce), for 6k cohorts.
# Each worker reduces its own partition's nnz to small partials; the full nnz
# scan is sharded across workers and never materialised in one place.
# ---------------------------------------------------------------------------

_ROLL_KEYS = ["sample_name", "length", "strand", "phase0"]
_RCTX: Dict = {}
_G2CTX: Dict = {}


def _gba_weight_shard(multi: "pl.LazyFrame", density: "pl.LazyFrame") -> "pl.LazyFrame":
    """Distribute each multimapper's count across its CDS candidates ∝ that
    sample's unique density in the candidate's gene (per-read normalised).

    A read's candidate loci are all co-located in one partition/shard, so the
    (read_id, sample_name) window is correct whether applied per-shard or globally.
    """
    return (
        multi.join(density, on=["sample_name", "gene_code"], how="left")
        .with_columns(pl.col("dens").fill_null(0.0))
        .with_columns(
            [
                pl.col("dens").sum().over(["read_id", "sample_name"]).alias("tot_dens"),
                pl.len().over(["read_id", "sample_name"]).alias("n_cand"),
            ]
        )
        .with_columns(
            pl.when(pl.col("tot_dens") > 0)
            .then(pl.col("dens") / pl.col("tot_dens"))
            .otherwise(1.0 / pl.col("n_cand"))
            .alias("w")
        )
        .with_columns((pl.col("count") * pl.col("w")).alias("count"))
        .group_by(_ROLL_KEYS)
        .agg(pl.col("count").sum().alias("count"))
    )


def _gba2_worker_init(density_path: str) -> None:
    _G2CTX["density"] = pl.read_parquet(density_path)


def _gba2_worker(shard_path: str) -> bytes:
    out = _gba_weight_shard(pl.scan_parquet(shard_path), _G2CTX["density"].lazy()).collect(
        engine="streaming"
    )
    buf = io.BytesIO()
    out.write_ipc(buf)
    return buf.getvalue()


def _roll_worker_init(
    frame_ivs, gene_ivs, bam_refs, ref_offset, reader, sample_map, mode, tmpdir
) -> None:
    _RCTX.update(
        frame_ivs=frame_ivs,
        gene_ivs=gene_ivs,
        bam_refs=bam_refs,
        ref_offset=ref_offset,
        reader=reader,
        sample_map=sample_map,
        mode=mode,
        tmpdir=tmpdir,
    )


def _roll_worker(pdir_str: str):
    """Scan one partition, join its counts, return (unique_roll, density, multi_shard_path)."""
    pdir = Path(pdir_str)
    bam = _discover_bam(pdir)
    empty_u = pl.DataFrame(
        schema={**{k: _INDEX_SCHEMA.get(k, pl.Utf8) for k in _ROLL_KEYS}, "count": pl.Float64}
    )
    if bam is None:
        b = io.BytesIO()
        empty_u.write_ipc(b)
        b2 = io.BytesIO()
        pl.DataFrame(
            schema={"sample_name": pl.Utf8, "gene_code": pl.Int32, "dens": pl.Float64}
        ).write_ipc(b2)
        return (b.getvalue(), b2.getvalue(), None)

    aln = _scan_dispatch(
        bam,
        _RCTX["frame_ivs"],
        _RCTX["gene_ivs"],
        _RCTX["bam_refs"],
        _RCTX["ref_offset"],
        _RCTX["reader"],
    )
    try:
        mp_, manifest = _manifest(pdir)
        cfiles = _count_parquets(mp_, manifest)
    except FileNotFoundError:
        cfiles = []

    if aln.is_empty() or not cfiles:
        b = io.BytesIO()
        empty_u.write_ipc(b)
        b2 = io.BytesIO()
        pl.DataFrame(
            schema={"sample_name": pl.Utf8, "gene_code": pl.Int32, "dens": pl.Float64}
        ).write_ipc(b2)
        return (b.getvalue(), b2.getvalue(), None)

    counts = pl.read_parquet(cfiles).select(["read_id", "sample_id", "count"])
    base = (
        counts.join(aln, on="read_id", how="inner")
        .join(_RCTX["sample_map"], on="sample_id", how="inner")
        .with_columns(pl.col("count").cast(pl.Float64))
    )
    uniq = base.filter(pl.col("n_align") == 1)
    uroll = uniq.group_by(_ROLL_KEYS).agg(pl.col("count").sum().alias("count"))
    dens = (
        uniq.filter(pl.col("gene_code") >= 0)
        .group_by(["sample_name", "gene_code"])
        .agg(pl.col("count").sum().alias("dens"))
    )

    multi_path = None
    if _RCTX["mode"] == "gba":
        multi = base.filter(pl.col("n_align") > 1).select(
            ["read_id", "sample_name", "length", "strand", "phase0", "gene_code", "count"]
        )
        if multi.height:
            multi_path = str(Path(_RCTX["tmpdir"]) / f"multi_{pdir.name}.parquet")
            multi.write_parquet(multi_path)

    bu = io.BytesIO()
    uroll.write_ipc(bu)
    bd = io.BytesIO()
    dens.write_ipc(bd)
    return (bu.getvalue(), bd.getvalue(), multi_path)


def build_matrix_rollup(
    partition_dirs,
    cds_df: pl.DataFrame,
    *,
    multimap_mode: str = "unique",
    ref_offset: int = 15,
    reader: str = "pysam",
    sample_names: Optional[List[str]] = None,
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """Single-pass map-reduce: scan + count-join folded into the workers.

    Returns the rollup (sample_name, length, strand, phase0, count) directly,
    without ever materialising the full nnz join — each worker reduces its own
    partition's counts to small partials. This is the 6k-sample path.
    """
    import shutil
    import tempfile

    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir() and _discover_bam(d))
    else:
        dirs = [Path(d) for d in partition_dirs if _discover_bam(Path(d))]
    out_schema = {
        "sample_name": pl.Utf8,
        "length": pl.Int64,
        "strand": pl.Utf8,
        "phase0": pl.Int8,
        "count": pl.Float64,
    }
    if not dirs:
        return pl.DataFrame(schema=out_schema)

    bam_refs = _bam_chroms(_discover_bam(dirs[0]))
    frame_ivs = _build_frame_intervals(cds_df, cds_df, bam_refs)
    gene_ivs, _ = build_gene_spans(cds_df)
    mp0, manifest0 = _manifest(dirs[0])
    sample_map = _samples_df(mp0, manifest0).select(["sample_id", "sample_name"])
    if sample_names:
        sample_map = sample_map.filter(pl.col("sample_name").is_in(sample_names))
    if n_workers is None:
        n_workers = min(len(dirs), os.cpu_count() or 1)

    tmpdir = tempfile.mkdtemp(prefix="ts_rollup_")
    log_info(f"matrix rollup ({multimap_mode}): {len(dirs)} partitions, {n_workers} worker(s)")

    u_parts: List[pl.DataFrame] = []
    d_parts: List[pl.DataFrame] = []
    multi_paths: List[str] = []
    try:
        dir_strs = [str(d) for d in dirs]
        if n_workers <= 1:
            _roll_worker_init(
                frame_ivs, gene_ivs, bam_refs, ref_offset, reader, sample_map, multimap_mode, tmpdir
            )
            results = [_roll_worker(s) for s in dir_strs]
        else:
            ctx = mp.get_context("spawn")
            with ctx.Pool(
                n_workers,
                initializer=_roll_worker_init,
                initargs=(
                    frame_ivs,
                    gene_ivs,
                    bam_refs,
                    ref_offset,
                    reader,
                    sample_map,
                    multimap_mode,
                    tmpdir,
                ),
            ) as pool:
                results = list(pool.imap_unordered(_roll_worker, dir_strs, chunksize=1))

        for bu, bd, mpath in results:
            u_parts.append(pl.read_ipc(io.BytesIO(bu)))
            d_parts.append(pl.read_ipc(io.BytesIO(bd)))
            if mpath:
                multi_paths.append(mpath)

        uroll = (
            pl.concat(u_parts).group_by(_ROLL_KEYS).agg(pl.col("count").sum().alias("count"))
            if u_parts
            else pl.DataFrame(schema=out_schema)
        )

        if multimap_mode == "unique":
            return uroll

        density = (
            pl.concat(d_parts)
            .group_by(["sample_name", "gene_code"])
            .agg(pl.col("dens").sum().alias("dens"))
            if d_parts
            else pl.DataFrame(
                schema={"sample_name": pl.Utf8, "gene_code": pl.Int32, "dens": pl.Float64}
            )
        )
        if not multi_paths:
            return uroll

        # GBA pass-2: weight multimapper shards against the reduced global
        # density.  Each shard is self-contained per read, so this parallelises
        # across workers (broadcast density via a Parquet file) — otherwise a
        # single streaming pass.
        if n_workers <= 1:
            multi = _gba_weight_shard(pl.scan_parquet(multi_paths), density.lazy()).collect(
                engine="streaming"
            )
        else:
            density_path = str(Path(tmpdir) / "density.parquet")
            density.write_parquet(density_path)
            ctx = mp.get_context("spawn")
            with ctx.Pool(
                min(n_workers, len(multi_paths)),
                initializer=_gba2_worker_init,
                initargs=(density_path,),
            ) as pool:
                m_parts = [
                    pl.read_ipc(io.BytesIO(b))
                    for b in pool.imap_unordered(_gba2_worker, multi_paths, chunksize=1)
                ]
            multi = (
                pl.concat(m_parts).group_by(_ROLL_KEYS).agg(pl.col("count").sum().alias("count"))
                if m_parts
                else pl.DataFrame(schema=out_schema)
            )
        return (
            pl.concat([uroll, multi]).group_by(_ROLL_KEYS).agg(pl.col("count").sum().alias("count"))
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Positional rollup (group by genomic A-site position) for profile clustering.
# Region-scoped: only the requested genes' reads are fetched (indexed BAM
# fetch), so nnz is bounded by the loci of interest, not the whole matrix.
# Same map-reduce shape as build_matrix_rollup, keyed by position not phase0.
# ---------------------------------------------------------------------------

_PCTX: Dict = {}


def gene_regions(cds_df: pl.DataFrame, gene_ids: List[str]) -> List[Tuple[str, str, int, int, int]]:
    """(chrom, strand, start, stop, gene_code) CDS span per requested gene."""
    sub = (
        cds_df.filter(pl.col("gene_id").is_in(gene_ids))
        .with_columns(
            [
                pl.col("start").list.min().alias("g0"),
                pl.col("stop").list.max().alias("g1"),
            ]
        )
        .group_by(["gene_id", "chr", "strand"])
        .agg([pl.col("g0").min(), pl.col("g1").max()])
    )
    code = {g: i for i, g in enumerate(sorted(gene_ids))}
    out = []
    for r in sub.iter_rows(named=True):
        out.append(
            (str(r["chr"]), str(r["strand"]), int(r["g0"]), int(r["g1"]), code[r["gene_id"]])
        )
    return out


def _prof_worker_init(
    regions, bam_refs, ref_offset, sample_map, group_level, cluster_labels
) -> None:
    _PCTX.update(
        regions=regions,
        bam_refs=bam_refs,
        ref_offset=ref_offset,
        sample_map=sample_map,
        group_level=group_level,
        cluster_labels=cluster_labels,
    )


def _prof_worker(pdir_str: str) -> bytes:
    pdir = Path(pdir_str)
    bam = _discover_bam(pdir)
    empty = pl.DataFrame(
        schema={"group": pl.Utf8, "gene_code": pl.Int32, "pos": pl.Int64, "count": pl.Float64}
    )
    if bam is None:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()
    bam_refs = _PCTX["bam_refs"]
    ref_offset = _PCTX["ref_offset"]

    rid, gcode, pos = [], [], []
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        for chrom, strand, start, stop, gc in _PCTX["regions"]:
            fetch_chrom = (
                chrom
                if chrom in bam_refs
                else (
                    f"chr{chrom}"
                    if f"chr{chrom}" in bam_refs
                    else (chrom[3:] if chrom[3:] in bam_refs else None)
                )
            )
            if fetch_chrom is None:
                continue
            want_rev = strand == "-"
            for rec in bf.fetch(fetch_chrom, max(0, start), stop):
                if rec.is_unmapped or rec.is_reverse != want_rev:
                    continue
                r = _parse_read_id(rec.query_name)
                if r is None:
                    continue
                a = (
                    int(rec.reference_end) - 1 - ref_offset
                    if rec.is_reverse
                    else int(rec.reference_start) + ref_offset
                )
                if a < start or a >= stop:
                    continue
                rid.append(r)
                gcode.append(gc)
                pos.append(a)
    if not rid:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()

    aln = pl.DataFrame(
        {
            "read_id": pl.Series(rid, dtype=pl.UInt64),
            "gene_code": pl.Series(gcode, dtype=pl.Int32),
            "pos": pl.Series(pos, dtype=pl.Int64),
        }
    ).unique()
    try:
        mp_, manifest = _manifest(pdir)
        cfiles = _count_parquets(mp_, manifest)
    except FileNotFoundError:
        cfiles = []
    if not cfiles:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()

    counts = pl.read_parquet(cfiles).select(["read_id", "sample_id", "count"])
    base = (
        counts.join(aln, on="read_id", how="inner")
        .join(_PCTX["sample_map"], on="sample_id", how="inner")
        .with_columns(pl.col("count").cast(pl.Float64))
    )

    level = _PCTX["group_level"]
    if level == "aggregate":
        base = base.with_columns(pl.lit("aggregate").alias("group"))
    elif level == "cluster":
        base = base.join(_PCTX["cluster_labels"], on="sample_name", how="inner").with_columns(
            pl.col("cluster_id").cast(pl.Utf8).alias("group")
        )
    else:
        base = base.with_columns(pl.col("sample_name").alias("group"))

    out = base.group_by(["group", "gene_code", "pos"]).agg(pl.col("count").sum().alias("count"))
    b = io.BytesIO()
    out.write_ipc(b)
    return b.getvalue()


def tabulate_profiles(
    partition_dirs,
    cds_df: pl.DataFrame,
    gene_ids: List[str],
    *,
    group_level: str = "sample",
    sample_names: Optional[List[str]] = None,
    cluster_labels: Optional[pl.DataFrame] = None,
    ref_offset: int = 15,
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """Per-(group, gene, A-site position) coverage for the requested genes.

    Region-scoped indexed fetch → the same aggregate/cluster/sample map-reduce
    as build_matrix_rollup, but keyed by genomic position. Pivot with
    ``profile_matrix`` to a [group × position] array for profile clustering /
    score_clustered.
    """
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir() and _discover_bam(d))
    else:
        dirs = [Path(d) for d in partition_dirs if _discover_bam(Path(d))]
    out_schema = {"group": pl.Utf8, "gene_code": pl.Int32, "pos": pl.Int64, "count": pl.Float64}
    if not dirs:
        return pl.DataFrame(schema=out_schema)

    regions = gene_regions(cds_df, gene_ids)
    code_to_gene = {c: g for g, c in {g: i for i, g in enumerate(sorted(gene_ids))}.items()}
    bam_refs = _bam_chroms(_discover_bam(dirs[0]))
    mp0, manifest0 = _manifest(dirs[0])
    sample_map = _samples_df(mp0, manifest0).select(["sample_id", "sample_name"])
    if sample_names:
        sample_map = sample_map.filter(pl.col("sample_name").is_in(sample_names))
    if n_workers is None:
        n_workers = min(len(dirs), os.cpu_count() or 1)

    dir_strs = [str(d) for d in dirs]
    if n_workers <= 1:
        _prof_worker_init(regions, bam_refs, ref_offset, sample_map, group_level, cluster_labels)
        parts = [pl.read_ipc(io.BytesIO(_prof_worker(s))) for s in dir_strs]
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            n_workers,
            initializer=_prof_worker_init,
            initargs=(regions, bam_refs, ref_offset, sample_map, group_level, cluster_labels),
        ) as pool:
            parts = [
                pl.read_ipc(io.BytesIO(b))
                for b in pool.imap_unordered(_prof_worker, dir_strs, chunksize=1)
            ]

    parts = [p for p in parts if not p.is_empty()]
    if not parts:
        return pl.DataFrame(schema=out_schema)
    rolled = (
        pl.concat(parts)
        .group_by(["group", "gene_code", "pos"])
        .agg(pl.col("count").sum().alias("count"))
    )
    return rolled.with_columns(
        pl.col("gene_code").replace_strict(code_to_gene, default=None).alias("gene_id")
    ).sort(["gene_id", "pos"])


def _read_totals_path(pdir: Path) -> Path:
    """Cache path for per-partition (read_id, total_count) summed over samples.

    Aggregate-tier coverage never needs the per-sample breakdown, so a read's
    total collapses 6k sample rows to one — making the aggregate path
    independent of cohort size. Cache lives next to the partition; override the
    directory with TS_MATRIX_CACHE_DIR for read-only matrices.
    """
    cache_root = os.environ.get("TS_MATRIX_CACHE_DIR")
    if cache_root:
        d = Path(cache_root) / pdir.name
        d.mkdir(parents=True, exist_ok=True)
        return d / "read_totals.parquet"
    return pdir / "read_totals.parquet"


def _read_totals(pdir: Path, cfiles) -> pl.DataFrame:
    """(read_id, count) summed over samples — cached on first build, reused after.

    Streaming aggregation keeps memory bounded even at 6k samples. If the cache
    cannot be written (read-only matrix, no TS_MATRIX_CACHE_DIR), it falls back
    to computing in-memory each call (still avoids the per-sample materialisation).
    """
    path = _read_totals_path(pdir)
    if path.exists():
        return pl.read_parquet(path).select(["read_id", "count"])
    agg = pl.scan_parquet(cfiles).group_by("read_id").agg(pl.col("count").sum().alias("count"))
    try:
        agg.sink_parquet(str(path))  # streaming: memory-safe even at 6k samples
        return pl.read_parquet(path).select(["read_id", "count"])
    except Exception:
        return agg.collect().select(["read_id", "count"])


def _cov_worker_init(regions, bam_refs, ref_offset, sample_map, group_level="aggregate") -> None:
    _PCTX.update(
        regions=regions,
        bam_refs=bam_refs,
        ref_offset=ref_offset,
        sample_map=sample_map,
        group_level=group_level,
    )


def _cov_worker(pdir_str: str) -> bytes:
    pdir = Path(pdir_str)
    bam = _discover_bam(pdir)
    per_sample = _PCTX.get("group_level") == "sample"
    schema = (
        {"sample_name": pl.Utf8, "strand": pl.Int64, "pos": pl.Int64, "count": pl.Float64}
        if per_sample
        else {"strand": pl.Int64, "pos": pl.Int64, "count": pl.Float64}
    )
    empty = pl.DataFrame(schema=schema)
    if bam is None:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()
    bam_refs = _PCTX["bam_refs"]
    off = _PCTX["ref_offset"]
    _timing = os.environ.get("TS_COV_TIMING")
    rid, pos, strand = [], [], []
    _t0 = time.perf_counter()
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        for chrom, start, stop in _PCTX["regions"]:
            fc = (
                chrom
                if chrom in bam_refs
                else (
                    f"chr{chrom}"
                    if f"chr{chrom}" in bam_refs
                    else (chrom[3:] if chrom[3:] in bam_refs else None)
                )
            )
            if fc is None:
                continue
            for rec in bf.fetch(fc, max(0, start), stop):
                if rec.is_unmapped:
                    continue
                r = _parse_read_id(rec.query_name)
                if r is None:
                    continue
                if rec.is_reverse:
                    rid.append(r)
                    pos.append(int(rec.reference_end) - 1 - off)
                    strand.append(-1)
                else:
                    rid.append(r)
                    pos.append(int(rec.reference_start) + off)
                    strand.append(1)
    t_bam = time.perf_counter() - _t0
    if not rid:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()
    aln = pl.DataFrame(
        {
            "read_id": pl.Series(rid, dtype=pl.UInt64),
            "pos": pl.Series(pos, dtype=pl.Int64),
            "strand": pl.Series(strand, dtype=pl.Int64),
        }
    )
    try:
        mp_, manifest = _manifest(pdir)
        cfiles = _count_parquets(mp_, manifest)
    except FileNotFoundError:
        cfiles = []
    if not cfiles:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()
    if per_sample:
        # Per-sample / cluster path: needs the sample dimension — unchanged.
        _t0 = time.perf_counter()
        counts = pl.read_parquet(cfiles).select(["read_id", "sample_id", "count"])
        t_read = time.perf_counter() - _t0
        _t0 = time.perf_counter()
        base = counts.join(aln, on="read_id", how="inner")
        t_join_read = time.perf_counter() - _t0
        _t0 = time.perf_counter()
        base = base.join(_PCTX["sample_map"], on="sample_id", how="inner").with_columns(
            pl.col("count").cast(pl.Float64)
        )
        t_join_sample = time.perf_counter() - _t0
        _t0 = time.perf_counter()
        out = base.group_by(["sample_name", "strand", "pos"]).agg(
            pl.col("count").sum().alias("count")
        )
        t_group = time.perf_counter() - _t0
        n_count_rows = counts.height
    else:
        # Aggregate fast path: per-read totals (summed over samples) — independent
        # of cohort size; no sample-level read, no sample_map join.
        _t0 = time.perf_counter()
        totals = _read_totals(pdir, cfiles)
        t_read = time.perf_counter() - _t0
        _t0 = time.perf_counter()
        base = totals.join(aln, on="read_id", how="inner").with_columns(
            pl.col("count").cast(pl.Float64)
        )
        t_join_read = time.perf_counter() - _t0
        t_join_sample = 0.0
        _t0 = time.perf_counter()
        out = base.group_by(["strand", "pos"]).agg(pl.col("count").sum().alias("count"))
        t_group = time.perf_counter() - _t0
        n_count_rows = totals.height
    if _timing:
        print(
            f"[cov-timing] {pdir.name} n_reads={len(rid)} count_rows={n_count_rows} "
            f"bam={t_bam:.3f} read={t_read:.3f} join_read={t_join_read:.3f} "
            f"join_sample={t_join_sample:.3f} group={t_group:.3f}",
            file=sys.stderr,
            flush=True,
        )
    b = io.BytesIO()
    out.write_ipc(b)
    return b.getvalue()


def region_coverage(
    partition_dirs,
    regions,
    *,
    ref_offset: int = 15,
    group_level: str = "aggregate",
    sample_names: Optional[List[str]] = None,
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """A-site coverage over genomic regions, strand-aware. group_level='aggregate'
    → (strand, pos, count); 'sample' → (sample_name, strand, pos, count)."""
    per_sample = group_level == "sample"
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir() and _discover_bam(d))
    else:
        dirs = [Path(d) for d in partition_dirs if _discover_bam(Path(d))]
    out_schema = (
        {"sample_name": pl.Utf8, "strand": pl.Int64, "pos": pl.Int64, "count": pl.Float64}
        if per_sample
        else {"strand": pl.Int64, "pos": pl.Int64, "count": pl.Float64}
    )
    if not dirs or not regions:
        return pl.DataFrame(schema=out_schema)
    bam_refs = _bam_chroms(_discover_bam(dirs[0]))
    mp0, manifest0 = _manifest(dirs[0])
    sample_map = _samples_df(mp0, manifest0).select(["sample_id", "sample_name"])
    if sample_names:
        sample_map = sample_map.filter(pl.col("sample_name").is_in(sample_names))
    if n_workers is None:
        n_workers = min(len(dirs), os.cpu_count() or 1)
    dir_strs = [str(d) for d in dirs]
    if n_workers <= 1:
        _cov_worker_init(regions, bam_refs, ref_offset, sample_map, group_level)
        parts = [pl.read_ipc(io.BytesIO(_cov_worker(s))) for s in dir_strs]
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            n_workers,
            initializer=_cov_worker_init,
            initargs=(regions, bam_refs, ref_offset, sample_map, group_level),
        ) as pool:
            parts = [
                pl.read_ipc(io.BytesIO(b))
                for b in pool.imap_unordered(_cov_worker, dir_strs, chunksize=1)
            ]
    parts = [p for p in parts if not p.is_empty()]
    if not parts:
        return pl.DataFrame(schema=out_schema)
    keys = ["sample_name", "strand", "pos"] if per_sample else ["strand", "pos"]
    return pl.concat(parts).group_by(keys).agg(pl.col("count").sum().alias("count")).sort("pos")


_JUNC_OVERHANG = 6  # min flanking aligned length on BOTH sides → confident span


def _junc_worker_init(junctions, bam_refs, sample_map, group_level, cluster_labels) -> None:
    _PCTX.update(
        junctions=junctions,
        bam_refs=bam_refs,
        sample_map=sample_map,
        group_level=group_level,
        cluster_labels=cluster_labels,
    )


def _junc_worker(pdir_str: str) -> bytes:
    pdir = Path(pdir_str)
    bam = _discover_bam(pdir)
    empty = pl.DataFrame(
        schema={"group": pl.Utf8, "junction_id": pl.UInt64, "kind": pl.Utf8, "count": pl.Float64}
    )
    if bam is None:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()
    bam_refs = _PCTX["bam_refs"]
    # Classify reads at each donor: span_conf (intron matches, overhang≥6 both
    # sides), span_short (matches but short overhang → ambiguous mapping),
    # unspliced (aligned straight through the donor → intron retention).
    rid, jid, kind = [], [], []
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        for chrom, donor, acceptor, strand, jjid in _PCTX["junctions"]:
            fc = (
                chrom
                if chrom in bam_refs
                else (
                    f"chr{chrom}"
                    if f"chr{chrom}" in bam_refs
                    else (chrom[3:] if chrom[3:] in bam_refs else None)
                )
            )
            if fc is None:
                continue
            for rec in bf.fetch(fc, max(0, donor), donor + 1):
                if rec.is_unmapped:
                    continue
                r = _parse_read_id(rec.query_name)
                if r is None:
                    continue
                blocks = rec.get_blocks()
                k = None
                for i in range(len(blocks) - 1):
                    if blocks[i][1] == donor and blocks[i + 1][0] == acceptor:
                        oh = min(blocks[i][1] - blocks[i][0], blocks[i + 1][1] - blocks[i + 1][0])
                        k = "span_conf" if oh >= _JUNC_OVERHANG else "span_short"
                        break
                if k is None:  # aligned straight through the donor (no intron there)?
                    for bs, be in blocks:
                        if bs <= donor - 1 and be >= donor + 1:
                            k = "unspliced"
                            break
                if k is not None:
                    rid.append(r)
                    jid.append(jjid)
                    kind.append(k)
    if not rid:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()
    aln = pl.DataFrame(
        {
            "read_id": pl.Series(rid, dtype=pl.UInt64),
            "junction_id": pl.Series(jid, dtype=pl.UInt64),
            "kind": pl.Series(kind, dtype=pl.Utf8),
        }
    ).unique()
    try:
        mp_, manifest = _manifest(pdir)
        cfiles = _count_parquets(mp_, manifest)
    except FileNotFoundError:
        cfiles = []
    if not cfiles:
        b = io.BytesIO()
        empty.write_ipc(b)
        return b.getvalue()
    counts = pl.read_parquet(cfiles).select(["read_id", "sample_id", "count"])
    base = (
        counts.join(aln, on="read_id", how="inner")
        .join(_PCTX["sample_map"], on="sample_id", how="inner")
        .with_columns(pl.col("count").cast(pl.Float64))
    )
    level = _PCTX["group_level"]
    if level == "aggregate":
        base = base.with_columns(pl.lit("aggregate").alias("group"))
    elif level == "cluster":
        base = base.join(_PCTX["cluster_labels"], on="sample_name", how="inner").with_columns(
            pl.col("cluster_id").cast(pl.Utf8).alias("group")
        )
    else:
        base = base.with_columns(pl.col("sample_name").alias("group"))
    out = base.group_by(["group", "junction_id", "kind"]).agg(pl.col("count").sum().alias("count"))
    b = io.BytesIO()
    out.write_ipc(b)
    return b.getvalue()


def tabulate_junctions(
    partition_dirs,
    junctions: List[Tuple[str, int, int, int, int]],
    *,
    group_level: str = "sample",
    sample_names: Optional[List[str]] = None,
    cluster_labels: Optional[pl.DataFrame] = None,
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """Spanning-read support per junction event. junctions: list of
    (chrom, donor, acceptor, strand, junction_id). A read supports a junction if
    a CIGAR intron (gap between aligned blocks) matches (donor, acceptor).

    NB: short overhangs across a junction are ambiguous (~30 nt reads) — an
    identifiability layer for junction mapping is deferred; this is raw support.
    """
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir() and _discover_bam(d))
    else:
        dirs = [Path(d) for d in partition_dirs if _discover_bam(Path(d))]
    out_schema = {"group": pl.Utf8, "junction_id": pl.UInt64, "kind": pl.Utf8, "count": pl.Float64}
    if not dirs or not junctions:
        return pl.DataFrame(schema=out_schema)
    bam_refs = _bam_chroms(_discover_bam(dirs[0]))
    mp0, manifest0 = _manifest(dirs[0])
    sample_map = _samples_df(mp0, manifest0).select(["sample_id", "sample_name"])
    if sample_names:
        sample_map = sample_map.filter(pl.col("sample_name").is_in(sample_names))
    if n_workers is None:
        n_workers = min(len(dirs), os.cpu_count() or 1)
    dir_strs = [str(d) for d in dirs]
    if n_workers <= 1:
        _junc_worker_init(junctions, bam_refs, sample_map, group_level, cluster_labels)
        parts = [pl.read_ipc(io.BytesIO(_junc_worker(s))) for s in dir_strs]
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            n_workers,
            initializer=_junc_worker_init,
            initargs=(junctions, bam_refs, sample_map, group_level, cluster_labels),
        ) as pool:
            parts = [
                pl.read_ipc(io.BytesIO(b))
                for b in pool.imap_unordered(_junc_worker, dir_strs, chunksize=1)
            ]
    parts = [p for p in parts if not p.is_empty()]
    if not parts:
        return pl.DataFrame(schema=out_schema)
    return (
        pl.concat(parts)
        .group_by(["group", "junction_id", "kind"])
        .agg(pl.col("count").sum().alias("count"))
    )


def junction_support_dict(jdf: pl.DataFrame, group: str = "aggregate") -> Dict[int, dict]:
    """Pivot tabulate_junctions output → {junction_id: {span_conf, span_short, unspliced}}."""
    out: Dict[int, dict] = {}
    if jdf.is_empty():
        return out
    sub = jdf.filter(pl.col("group") == group)
    for r in sub.iter_rows(named=True):
        out.setdefault(int(r["junction_id"]), {})[r["kind"]] = float(r["count"])
    return out


def profile_matrix(
    profiles: pl.DataFrame, gene_id: str
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Pivot tabulate_profiles output for one gene → ([group × position] matrix,
    group labels, position vector).  Ready for profile_clustering / score_clustered."""
    sub = profiles.filter(pl.col("gene_id") == gene_id)
    if sub.is_empty():
        return np.zeros((0, 0)), [], np.array([], dtype=np.int64)
    wide = sub.pivot(values="count", index="group", on="pos").fill_null(0.0).sort("group")
    groups = wide.get_column("group").to_list()
    pos_cols = sorted((int(c) for c in wide.columns if c != "group"))
    mat = wide.select([str(p) for p in pos_cols]).to_numpy()
    return mat, groups, np.array(pos_cols, dtype=np.int64)


# ---------------------------------------------------------------------------
# Phase B: vectorised rollup with multimap strategy
# ---------------------------------------------------------------------------


def tabulate_rollup(
    partition_dirs,
    index: "pl.DataFrame | str | Path",
    *,
    multimap_mode: str = "unique",
    sample_names: Optional[List[str]] = None,
) -> pl.DataFrame:
    """Reduce counts × index → rollup (sample_name, length, strand, phase0, count).

    multimap_mode: "unique" (NH==1 only) | "gba" (guilty-by-association rescue).
    """
    if isinstance(index, (str, Path)):
        index = pl.read_parquet(str(index))
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir())
    else:
        dirs = [Path(d) for d in partition_dirs]

    mp0, manifest0 = _manifest(dirs[0])
    samples = _samples_df(mp0, manifest0).select(["sample_id", "sample_name"])
    if sample_names:
        samples = samples.filter(pl.col("sample_name").is_in(sample_names))

    count_files: List[str] = []
    for pdir in dirs:
        try:
            mp_, manifest = _manifest(pdir)
        except FileNotFoundError:
            continue
        count_files.extend(_count_parquets(mp_, manifest))

    out_schema = {
        "sample_name": pl.Utf8,
        "length": pl.Int64,
        "strand": pl.Utf8,
        "phase0": pl.Int8,
        "count": pl.Float64,
    }
    if not count_files:
        return pl.DataFrame(schema=out_schema)

    counts = pl.scan_parquet(count_files).select(["read_id", "sample_id", "count"])
    base = (
        counts.join(index.lazy(), on="read_id", how="inner")
        .join(samples.lazy(), on="sample_id", how="inner")
        .with_columns(pl.col("count").cast(pl.Float64))
    )

    keys = ["sample_name", "length", "strand", "phase0"]

    if multimap_mode == "unique":
        return (
            base.filter(pl.col("n_align") == 1)
            .group_by(keys)
            .agg(pl.col("count").sum().alias("count"))
            .collect(engine="streaming")
        )

    if multimap_mode != "gba":
        raise ValueError(f"unknown multimap_mode {multimap_mode!r}")

    # --- GBA: unique reads weight 1 + build density; multimappers rescued ---
    base = base.collect(engine="streaming").lazy()  # materialise once, reuse

    uniq = base.filter(pl.col("n_align") == 1)
    uniq_roll = uniq.group_by(keys).agg(pl.col("count").sum().alias("count"))

    # per-(sample, gene) unique-read density
    density = (
        uniq.filter(pl.col("gene_code") >= 0)
        .group_by(["sample_name", "gene_code"])
        .agg(pl.col("count").sum().alias("dens"))
    )

    multi = base.filter(pl.col("n_align") > 1)
    multi = (
        multi.join(density, on=["sample_name", "gene_code"], how="left")
        .with_columns(pl.col("dens").fill_null(0.0))
        .with_columns(
            [
                pl.col("dens").sum().over(["read_id", "sample_name"]).alias("tot_dens"),
                pl.len().over(["read_id", "sample_name"]).alias("n_cand"),
            ]
        )
        .with_columns(
            pl.when(pl.col("tot_dens") > 0)
            .then(pl.col("dens") / pl.col("tot_dens"))
            .otherwise(1.0 / pl.col("n_cand"))
            .alias("w")
        )
        .with_columns((pl.col("count") * pl.col("w")).alias("count"))
    )
    multi_roll = multi.group_by(keys).agg(pl.col("count").sum().alias("count"))

    return (
        pl.concat([uniq_roll, multi_roll])
        .group_by(keys)
        .agg(pl.col("count").sum().alias("count"))
        .collect(engine="streaming")
    )


# ---------------------------------------------------------------------------
# Phase C: central offset calibration + scoring (operates on the tiny rollup)
# ---------------------------------------------------------------------------


def _frame_at(phase0: int, strand: str, offset: int) -> int:
    return (phase0 + offset) % 3 if strand == "+" else (phase0 - offset) % 3


def calibrate_offsets(
    rollup: pl.DataFrame,
    *,
    offset_lo: int = 10,
    offset_hi: int = 18,
    min_reads: int = 50,
) -> Dict[Tuple[str, int], int]:
    """Per-(sample, length) offset maximising dominant-frame fraction."""
    offsets: Dict[Tuple[str, int], int] = {}
    if rollup.is_empty():
        return offsets
    for (sample, length), grp in rollup.group_by(["sample_name", "length"]):
        rows = list(grp.select(["strand", "phase0", "count"]).iter_rows())
        total = sum(c for _, _, c in rows)
        if total < min_reads:
            continue
        best_o, best_dom = offset_lo, -1.0
        for o in range(offset_lo, offset_hi + 1):
            fc = [0.0, 0.0, 0.0]
            for strand, phase0, c in rows:
                fc[_frame_at(int(phase0), strand, o)] += c
            dom = max(fc) / total
            if dom > best_dom:
                best_dom, best_o = dom, o
        offsets[(str(sample), int(length))] = best_o
    return offsets


def rollup_to_periodicity(
    rollup: pl.DataFrame,
    offsets: Dict[Tuple[str, int], int],
    *,
    default_offset: int = 15,
) -> pl.DataFrame:
    """Apply per-(sample, length) offsets to the rollup → periodicity_qc schema."""
    from .matrix_qc import _empty_periodicity_schema

    if rollup.is_empty():
        return _empty_periodicity_schema()

    # sample -> {length: {0,1,2}}
    rfd: Dict[str, Dict[int, Dict[int, float]]] = {}
    used_off: Dict[str, Dict[int, int]] = {}
    for row in rollup.iter_rows(named=True):
        s = str(row["sample_name"])
        L = int(row["length"])
        o = offsets.get((s, L), default_offset)
        fr = _frame_at(int(row["phase0"]), str(row["strand"]), o)
        rfd.setdefault(s, {}).setdefault(L, {0: 0.0, 1: 0.0, 2: 0.0})[fr] += float(row["count"])
        used_off.setdefault(s, {})[L] = o

    rows = []
    for sname, lfd in rfd.items():
        scores = _ribometric_frame_scores(lfd)
        rfd_json, dom_json = {}, {}
        for length, fd in lfd.items():
            counts = [fd.get(0, 0.0), fd.get(1, 0.0), fd.get(2, 0.0)]
            rfd_json[str(length)] = counts
            tot = sum(counts)
            dom_json[str(length)] = (max(counts) / tot) if tot > 0 else 0.0
        rows.append(
            {
                "sample_id": sname,
                "periodicity_score": scores["periodicity_score"],
                "f0": scores["f0"],
                "f1": scores["f1"],
                "f2": scores["f2"],
                "n_cds_reads": scores["n_reads"],
                "recommended_offsets": json.dumps({str(k): v for k, v in used_off[sname].items()}),
                "read_frame_distribution": json.dumps(rfd_json),
                "per_length_periodicity": json.dumps(
                    {str(k): v for k, v in scores["per_length"].items()}
                ),
                "per_length_dominance": json.dumps(dom_json),
            }
        )
    return pl.from_dicts(rows).sort("sample_id") if rows else _empty_periodicity_schema()


def periodicity_qc_scalable(
    partition_dirs,
    cds_df: pl.DataFrame,
    *,
    index: "pl.DataFrame | str | Path | None" = None,
    multimap_mode: str = "unique",
    offset_mode: str = "calibrate",  # "calibrate" | "fixed"
    fixed_offset: int = 15,
    ref_offset: int = 15,
    sample_names: Optional[List[str]] = None,
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """End-to-end scalable periodicity QC: index (or reuse) → rollup → calibrate → score."""
    if index is None:
        index = build_alignment_index(
            partition_dirs, cds_df, ref_offset=ref_offset, n_workers=n_workers
        )
    rollup = tabulate_rollup(
        partition_dirs, index, multimap_mode=multimap_mode, sample_names=sample_names
    )
    if offset_mode == "fixed":
        offsets = {
            (s, L): fixed_offset
            for s, L in {
                (str(r["sample_name"]), int(r["length"])) for r in rollup.iter_rows(named=True)
            }
        }
    else:
        offsets = calibrate_offsets(rollup)
    return rollup_to_periodicity(rollup, offsets, default_offset=fixed_offset)

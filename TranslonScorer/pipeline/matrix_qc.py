"""Per-sample QC metrics computed directly from the sparse Parquet matrix.

Two tiers — designed to mirror RiboMetric's per-sample QC but operating on
the global multi-sample matrix without materialising per-sample BAM files.

Tier 1 — fast_length_qc (no BAM required)
    Joins reads.parquet (read_id → length) with counts.parquet
    (read_id → sample_id → count) to produce weighted per-sample length
    distributions. Runs in seconds regardless of cohort size.

    Output per sample:
        total_reads, unique_reads, mean_length, peak_length,
        rpf_28_32_prop, length_cv, compression_ratio

Tier 2 — periodicity_qc (all partitions + annotation)
    Samples a configurable number of read buckets from the global BAM.
    For each sampled read, looks up which samples contain it (and with
    what count) then maps its genomic position to a CDS-relative frame.
    The per-sample per-length frame distribution yields:
        periodicity_score  — entropy reduction (matches RiboMetric formula)
        f0, f1, f2         — frame fractions
        recommended_offsets — {length: offset} nudged to put dominant
                               frame back to frame 0

    n_sample_buckets controls speed/accuracy: 5 buckets is enough for
    offset calibration; 20 gives robust periodicity estimates.

Entry point: matrix_qc() runs both tiers (periodicity only when bam_path
and exon_df are supplied) and returns a joined DataFrame.
"""
from __future__ import annotations

import collections
import heapq
import json
import math
import multiprocessing as mp
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
import pysam

from ..utils.logging import log_info, log_warning


# ---------------------------------------------------------------------------
# Manifest + lookup helpers  [shim — moved to io/matrix.py]
# ---------------------------------------------------------------------------
from TranslonScorer.io.matrix import (  # noqa: F401, E402
    _count_parquets,
    _discover_bam,
    _manifest,
    _parse_read_id,
    _reads_parquet_path,
    _samples_df,
)


# ---------------------------------------------------------------------------
# Tier 1: length distribution QC (no BAM)
# ---------------------------------------------------------------------------

def fast_length_qc(
    partition_dir: str | Path,
    *,
    sample_names: Optional[List[str]] = None,
) -> pl.DataFrame:
    """Per-sample length distribution metrics from reads + counts Parquets.

    No BAM required. Returns one row per sample with columns:
        sample_id, study_id, total_reads, unique_reads, compression_ratio,
        mean_length, peak_length, rpf_28_32_prop, length_cv
    """
    mp, manifest = _manifest(partition_dir)
    count_files  = _count_parquets(mp, manifest)
    reads_path   = _reads_parquet_path(mp, manifest)
    samples      = _samples_df(mp, manifest)

    if sample_names:
        samples = samples.filter(pl.col("sample_name").is_in(sample_names))

    log_info(f"Length QC: {len(count_files)} count file(s), {samples.height} samples")

    counts_lf = pl.scan_parquet(count_files)
    reads_lf  = pl.scan_parquet(reads_path).select(["read_id", "length"])
    sample_lf = samples.lazy().select(["sample_id", "sample_name", "study_id"])

    per_sl = (
        counts_lf
        .join(reads_lf, on="read_id", how="inner")
        .join(sample_lf, on="sample_id", how="inner")
        .group_by(["sample_name", "study_id", "length"])
        .agg([
            pl.col("count").sum().cast(pl.Float64).alias("total_count"),
            pl.col("read_id").n_unique().alias("unique_at_length"),
        ])
        .collect()
    )

    if per_sl.is_empty():
        return _empty_length_schema()

    rows = []
    for (sname, study_id), grp in per_sl.group_by(["sample_name", "study_id"]):
        total  = float(grp["total_count"].sum())
        unique = int(grp["unique_at_length"].sum())
        if total == 0:
            continue
        lengths = grp["length"].to_numpy().astype(float)
        counts  = grp["total_count"].to_numpy()
        mean_l  = float(np.average(lengths, weights=counts))
        std_l   = float(np.sqrt(np.average((lengths - mean_l) ** 2, weights=counts)))
        peak_l  = int(lengths[np.argmax(counts)])
        rpf_m   = (lengths >= 28) & (lengths <= 32)
        rows.append({
            "sample_id":         str(sname),
            "study_id":          str(study_id),
            "total_reads":       total,
            "unique_reads":      unique,
            "compression_ratio": unique / total,
            "mean_length":       mean_l,
            "peak_length":       peak_l,
            "rpf_28_32_prop":    float(counts[rpf_m].sum() / total),
            "length_cv":         float(std_l / mean_l) if mean_l > 0 else 0.0,
        })

    return pl.from_dicts(rows).sort("sample_id") if rows else _empty_length_schema()


def _empty_length_schema() -> pl.DataFrame:
    return pl.DataFrame(schema={
        "sample_id": pl.Utf8, "study_id": pl.Utf8,
        "total_reads": pl.Float64, "unique_reads": pl.UInt32,
        "compression_ratio": pl.Float64, "mean_length": pl.Float64,
        "peak_length": pl.UInt32, "rpf_28_32_prop": pl.Float64,
        "length_cv": pl.Float64,
    })


# ---------------------------------------------------------------------------
# Tier 2: periodicity and offset calibration (sampled BAM)
# ---------------------------------------------------------------------------

def _bam_chroms(bam_path: str | Path) -> set:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return set(bam.references)


from TranslonScorer.io.bam import normalise_chrom as _normalise_chrom  # noqa: F401


def _compute_periodicity_from_frames(frames: Dict[int, float]) -> Dict:
    """Entropy-reduction periodicity score — exact RiboMetric formula (sqrt of info ratio)."""
    total = sum(frames.values())
    if total == 0:
        return {"periodicity_score": 0.0, "f0": 0.0, "f1": 0.0, "f2": 0.0}
    pseudocount = 1e-100
    probs = [(frames.get(f, 0.0) + pseudocount) / (total + pseudocount) for f in range(3)]
    max_ent = math.log2(3)
    entropy = -sum(p * math.log2(p) for p in probs if p > 0)
    score = float(math.sqrt(max(0.0, (max_ent - entropy) / max_ent)))
    return {"periodicity_score": score, "f0": probs[0], "f1": probs[1], "f2": probs[2]}


def _ribometric_frame_scores(
    rfd: Dict[int, Dict[int, float]],
    min_count_frac: float = 0.05,
) -> Dict:
    """Apply RiboMetric's read_frame_information_content + information_metric_cutoff.

    rfd = {length: {0: count, 1: count, 2: count}} — the read_frame_distribution.
    Returns {periodicity_score, f0, f1, f2, n_reads, per_length: {length: score}}.
    """
    pseudocount = 1e-100
    log2_3 = math.log2(3)
    total_reads = sum(sum(fd.values()) for fd in rfd.values())

    per_length: Dict[int, Tuple[float, int]] = {}
    for length, fd in rfd.items():
        n = sum(fd.values())
        probs = [(fd.get(f, 0.0) + pseudocount) / (n + pseudocount) for f in range(3)]
        entropy = -sum(p * math.log2(p) for p in probs if p > 0)
        per_length[length] = (math.sqrt(max(0.0, (log2_3 - entropy) / log2_3)), int(n))

    valid = {l: (s, n) for l, (s, n) in per_length.items() if n > total_reads * min_count_frac}
    if valid:
        w = sum(n for _, n in valid.values())
        global_score = sum(s * n for s, n in valid.values()) / w
    else:
        global_score = 0.0

    agg: Dict[int, float] = {0: 0.0, 1: 0.0, 2: 0.0}
    for fd in rfd.values():
        for f, c in fd.items():
            agg[f] += c
    t = sum(agg.values()) or 1.0

    return {
        "periodicity_score": global_score,
        "f0": agg[0] / t, "f1": agg[1] / t, "f2": agg[2] / t,
        "n_reads": int(total_reads),
        "per_length": {l: s for l, (s, _) in per_length.items()},
    }


def _nudge_to_frame0(
    length_frame_dist: Dict[int, Dict[int, float]],
    base_offsets: Dict[int, int],
    *,
    min_reads: int = 50,
    min_fraction: float = 0.60,
) -> Dict[int, int]:
    """Nudge per-length offsets ±1 to put dominant frame back to 0.

    Matches RiboMetric._frame_calibrated_offsets() logic.
    """
    out = dict(base_offsets)
    for length, fd in length_frame_dist.items():
        total = sum(fd.values())
        if total < min_reads:
            continue
        dom = max(fd, key=fd.get)
        if dom == 0 or fd[dom] / total < min_fraction:
            continue
        shift = -1 if dom == 1 else 1
        if int(length) in out:
            out[int(length)] += shift
    return out


# Re-export shims — implementations live in TranslonScorer.events
from TranslonScorer.events import (  # noqa: F401, E402
    _build_frame_intervals,
    _deconflict_intervals,
)


def _assign_frames_sweep(
    read_ids: List[int],
    asites: "np.ndarray",
    strand: str,
    intervals: List[Tuple[int, int, int]],
) -> Dict[int, int]:
    """Sweep-line O((m+n) log n) frame assignment.

    Iterates reads sorted by A-site position and maintains an active heap of
    intervals sorted by stop. Returns {read_id: frame} for CDS-overlapping reads.
    """
    if not intervals or len(read_ids) == 0:
        return {}

    order = np.argsort(asites)
    sorted_asites = asites[order].tolist()
    sorted_rids   = [read_ids[i] for i in order]

    result: Dict[int, int] = {}
    active: list = []   # min-heap: (ee, phase)
    iv_ptr = 0
    n_ivs = len(intervals)

    for asite, rid in zip(sorted_asites, sorted_rids):
        # Advance intervals whose start ≤ asite
        while iv_ptr < n_ivs and intervals[iv_ptr][0] <= asite:
            es, ee, phase = intervals[iv_ptr]
            heapq.heappush(active, (ee, phase))
            iv_ptr += 1
        # Evict expired intervals (ee ≤ asite means asite ≥ ee → outside)
        while active and active[0][0] <= asite:
            heapq.heappop(active)
        # First remaining active interval contains asite
        if active:
            _, phase_top = active[0]
            if strand == "+":
                result[rid] = (phase_top + asite) % 3
            else:
                result[rid] = (phase_top - asite) % 3

    return result




def periodicity_qc(
    partition_dirs: "str | Path | List[str | Path]",
    exon_df: pl.DataFrame,
    *,
    cds_df: Optional[pl.DataFrame] = None,
    sample_names: Optional[List[str]] = None,
    default_offset: int = 15,
) -> pl.DataFrame:
    """Per-sample periodicity score and recommended A-site offsets.

    The matrix partitions unique reads by sequence prefix — each of the 256+
    partitions holds ~1/256 of all reads, distributed by the first 4 nt of the
    read sequence.  A single codon position can be sampled by reads from any
    partition, so building a representative read_frame_distribution requires
    aggregating across ALL partitions (equivalent to RiboMetric's per-sample
    annotated_read_df, but assembled from genome BAMs rather than one
    transcriptome BAM).

    Scores use RiboMetric's exact formula:
        per_length_score = sqrt((log2(3) - H) / log2(3))
        global_score     = count-weighted mean across lengths with >5% of reads
    """
    # Resolve partition list
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs: List[Path] = sorted(
            d for d in parent.iterdir() if d.is_dir() and _discover_bam(d)
        )
    else:
        dirs = [Path(d) for d in partition_dirs]

    if not dirs:
        log_warning("No partition directories found — skipping periodicity QC")
        return _empty_periodicity_schema()

    log_info(f"Periodicity QC: processing all {len(dirs)} partitions")

    # Get BAM chroms from the first available BAM (same genome across all partitions)
    first_bam = _discover_bam(dirs[0])
    if first_bam is None:
        log_warning("No BAM found in first partition")
        return _empty_periodicity_schema()
    bam_refs = _bam_chroms(first_bam)

    # Build exon-aware CDS frame intervals once
    frame_ivs: Dict[Tuple[str, str], List[Tuple[int, int, int]]] = {}
    if cds_df is not None:
        log_info("Building exon-aware CDS frame intervals…")
        frame_ivs = _build_frame_intervals(exon_df, cds_df, bam_refs)
        n_ivs = sum(len(v) for v in frame_ivs.values())
        log_info(f"  {n_ivs:,} CDS exon intervals across {len(frame_ivs)} (chrom, strand) groups")

    # read_frame_distribution per sample — equivalent to RiboMetric's annotated_read_df
    # grouped by (sample, read_length, read_frame), aggregated across all partitions.
    # {sample_name: {length: {0: count, 1: count, 2: count}}}
    rfd: Dict[str, Dict[int, Dict[int, float]]] = {}
    n_cds_reads: Dict[str, int] = {}
    total_cds_hit = total_reads_seen = 0

    # Build sample_map from the first partition (same samples across all)
    mp0, manifest0 = _manifest(dirs[0])
    samples_df0 = _samples_df(mp0, manifest0)
    if sample_names:
        samples_df0 = samples_df0.filter(pl.col("sample_name").is_in(sample_names))
    sample_map: Dict[int, str] = {
        int(r[0]): str(r[1])
        for r in samples_df0.select(["sample_id", "sample_name"]).iter_rows()
    }

    for pdir in dirs:
        bam_path = _discover_bam(pdir)
        if bam_path is None:
            continue
        try:
            mp, manifest = _manifest(pdir)
        except FileNotFoundError:
            continue
        count_files = _count_parquets(mp, manifest)

        # Read all BAM alignments → (read_id → asite, chrom, length, strand)
        # One position per read_id (last alignment wins for multimappers —
        # acceptable since we're accumulating frame signal not exact positions)
        read_pos: Dict[int, Tuple[int, str, int, str]] = {}
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
                strand = "-" if rec.is_reverse else "+"
                asite = (
                    int(rec.reference_end) - 1 - default_offset
                    if rec.is_reverse
                    else int(rec.reference_start) + default_offset
                )
                read_pos[rid] = (asite, rec.reference_name, length, strand)

        if not read_pos:
            continue
        total_reads_seen += len(read_pos)

        # Batch-assign CDS frames via sweep-line
        read_frames: Dict[int, int] = {}
        if frame_ivs:
            groups: Dict[Tuple[str, str], Tuple[List[int], List[int]]] = {}
            for rid, (asite, chrom, length, strand) in read_pos.items():
                norm = _normalise_chrom(chrom, bam_refs) or chrom
                key = (norm, strand)
                if key not in groups:
                    groups[key] = ([], [])
                groups[key][0].append(rid)
                groups[key][1].append(asite)
            for key, (rids, asite_list) in groups.items():
                ivs = frame_ivs.get(key)
                if ivs:
                    frames = _assign_frames_sweep(
                        rids, np.array(asite_list, dtype=np.int64), key[1], ivs
                    )
                    read_frames.update(frames)
            total_cds_hit += len(read_frames)

        # Load counts for this partition's reads and accumulate
        sampled_ids = pl.Series(list(read_pos.keys()), dtype=pl.UInt64)
        if not count_files:
            continue
        counts_df = (
            pl.concat([pl.read_parquet(f) for f in count_files])
            .filter(pl.col("read_id").is_in(sampled_ids))
        )

        for row in counts_df.iter_rows(named=True):
            rid   = int(row["read_id"])
            sname = sample_map.get(int(row["sample_id"]))
            if sname is None:
                continue
            info = read_pos.get(rid)
            if info is None:
                continue
            # Only include reads that were assigned a CDS frame — excludes
            # poly-A, intronic, and repeat reads whose genomic position % 3
            # introduces a systematic bias unrelated to translation.
            frame = read_frames.get(rid)
            if frame is None:
                continue
            cnt = float(row["count"])
            asite_g, chrom, length, strand = info

            sample_rfd = rfd.setdefault(sname, {})
            sample_rfd.setdefault(length, {0: 0.0, 1: 0.0, 2: 0.0})[frame] += cnt
            n_cds_reads[sname] = n_cds_reads.get(sname, 0) + 1

    if total_reads_seen:
        log_info(
            f"Frame assignment across {len(dirs)} partitions: "
            f"{total_cds_hit:,}/{total_reads_seen:,} reads ({total_cds_hit/total_reads_seen:.1%}) hit CDS exons"
        )

    if not rfd:
        log_warning("No CDS reads accumulated — check annotation and BAM coverage")
        return _empty_periodicity_schema()

    base_offsets = {length: default_offset for lfd in rfd.values() for length in lfd}

    rows = []
    for sname, lfd in rfd.items():
        scores  = _ribometric_frame_scores(lfd)
        offsets = _nudge_to_frame0(lfd, base_offsets)
        # Per-read-length frame distribution and dominance — the read-length×frame
        # signal RiboMetric surfaces, kept here instead of being collapsed into the
        # aggregate f0/f1/f2 (which cancels across lengths with different offsets).
        rfd_json: Dict[str, List[float]] = {}
        dominance_json: Dict[str, float] = {}
        for length, fd in lfd.items():
            counts = [float(fd.get(0, 0.0)), float(fd.get(1, 0.0)), float(fd.get(2, 0.0))]
            rfd_json[str(length)] = counts
            tot = sum(counts)
            dominance_json[str(length)] = (max(counts) / tot) if tot > 0 else 0.0
        rows.append({
            "sample_id":         sname,
            "periodicity_score": scores["periodicity_score"],
            "f0": scores["f0"], "f1": scores["f1"], "f2": scores["f2"],
            "n_cds_reads":       scores["n_reads"],
            "recommended_offsets": json.dumps({str(k): v for k, v in offsets.items()}),
            # length -> [count_f0, count_f1, count_f2]
            "read_frame_distribution": json.dumps(rfd_json),
            # length -> sqrt entropy-reduction periodicity score
            "per_length_periodicity": json.dumps({str(k): v for k, v in scores["per_length"].items()}),
            # length -> dominant-frame fraction (max frame proportion), offset-invariant
            "per_length_dominance": json.dumps(dominance_json),
        })

    return pl.from_dicts(rows).sort("sample_id") if rows else _empty_periodicity_schema()


def _empty_periodicity_schema() -> pl.DataFrame:
    return pl.DataFrame(schema={
        "sample_id": pl.Utf8,
        "periodicity_score": pl.Float64,
        "f0": pl.Float64, "f1": pl.Float64, "f2": pl.Float64,
        "n_cds_reads": pl.Int64,
        "recommended_offsets": pl.Utf8,
        "read_frame_distribution": pl.Utf8,
        "per_length_periodicity": pl.Utf8,
        "per_length_dominance": pl.Utf8,
    })


# ---------------------------------------------------------------------------
# Read-length × sample frame-dominance matrix
# ---------------------------------------------------------------------------

def frame_dominance_matrix(
    perio_df: pl.DataFrame,
    *,
    metric: str = "dominance",
    min_reads: int = 50,
) -> pl.DataFrame:
    """Pivot per-sample periodicity output into a read_length × sample matrix.

    perio_df : output of ``periodicity_qc`` (must carry the per-length JSON
               columns ``read_frame_distribution`` and either
               ``per_length_dominance`` or ``per_length_periodicity``).
    metric   : ``"dominance"`` → dominant-frame fraction (max frame proportion,
               in [1/3, 1]); ``"periodicity"`` → sqrt entropy-reduction score
               (in [0, 1]).

    Returns a wide DataFrame: one row per ``read_length`` (sorted), one column
    per sample, cells holding the chosen metric.  Cells with fewer than
    ``min_reads`` reads at that (sample, length) are left null.
    """
    if perio_df.is_empty() or "read_frame_distribution" not in perio_df.columns:
        return pl.DataFrame(schema={"read_length": pl.Int64})

    src_col = "per_length_dominance" if metric == "dominance" else "per_length_periodicity"
    long_rows: List[dict] = []
    for row in perio_df.iter_rows(named=True):
        sample = str(row["sample_id"])
        rfd = json.loads(row.get("read_frame_distribution") or "{}")
        vals = json.loads(row.get(src_col) or "{}")
        for length_str, counts in rfd.items():
            n = sum(counts)
            if n < min_reads:
                continue
            long_rows.append({
                "read_length": int(length_str),
                "sample_id": sample,
                "value": float(vals.get(length_str, 0.0)),
            })

    if not long_rows:
        return pl.DataFrame(schema={"read_length": pl.Int64})

    long = pl.from_dicts(long_rows)
    wide = long.pivot(values="value", index="read_length", on="sample_id").sort("read_length")
    return wide


# ===========================================================================
# Two-phase fast path: cached read-locus index  +  vectorised tabulation
# ---------------------------------------------------------------------------
# The BAM scan + CDS frame assignment is sample- and query-independent: it is a
# property of the unique reads and the annotation only.  Splitting it out lets
# us pay it ONCE (build_read_loci_index) and then answer any number of
# aggregate / cluster / per-sample queries as pure-Polars joins+group-bys
# (tabulate_frames), with no pysam and no Python row loop.
# ===========================================================================

def _scan_partition_loci(
    bam_path: str | Path,
    frame_ivs: Dict[Tuple[str, str], List[Tuple[int, int, int]]],
    bam_refs: set,
    default_offset: int,
) -> pl.DataFrame:
    """One partition BAM → DataFrame[read_id, length, frame] for CDS-frame reads.

    Frame assignment is identical to ``periodicity_qc`` (same A-site formula,
    same sweep-line); only CDS-overlapping reads are emitted.
    """
    read_pos: Dict[int, Tuple[int, str, str]] = {}
    read_len: Dict[int, int] = {}
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
            asite = (
                int(rec.reference_end) - 1 - default_offset
                if rec.is_reverse
                else int(rec.reference_start) + default_offset
            )
            strand = "-" if rec.is_reverse else "+"
            read_pos[rid] = (asite, rec.reference_name, strand)
            read_len[rid] = length

    if not read_pos or not frame_ivs:
        return pl.DataFrame(schema={"read_id": pl.UInt64, "length": pl.Int64, "frame": pl.Int8})

    groups: Dict[Tuple[str, str], Tuple[List[int], List[int]]] = {}
    for rid, (asite, chrom, strand) in read_pos.items():
        key = (_normalise_chrom(chrom, bam_refs) or chrom, strand)
        groups.setdefault(key, ([], []))
        groups[key][0].append(rid)
        groups[key][1].append(asite)

    read_frames: Dict[int, int] = {}
    for key, (rids, asite_list) in groups.items():
        ivs = frame_ivs.get(key)
        if ivs:
            read_frames.update(
                _assign_frames_sweep(rids, np.array(asite_list, dtype=np.int64), key[1], ivs)
            )

    if not read_frames:
        return pl.DataFrame(schema={"read_id": pl.UInt64, "length": pl.Int64, "frame": pl.Int8})

    rids = list(read_frames.keys())
    return pl.DataFrame({
        "read_id": pl.Series(rids, dtype=pl.UInt64),
        "length": pl.Series([read_len[r] for r in rids], dtype=pl.Int64),
        "frame": pl.Series([read_frames[r] for r in rids], dtype=pl.Int8),
    })


# Worker globals for the parallel index build (set once per spawned process via
# the Pool initializer, so frame_ivs is pickled once per worker, not per task).
_BUILD_CTX: Dict = {}


def _build_worker_init(frame_ivs, bam_refs, default_offset) -> None:
    _BUILD_CTX["frame_ivs"] = frame_ivs
    _BUILD_CTX["bam_refs"] = bam_refs
    _BUILD_CTX["default_offset"] = default_offset


def _build_worker_scan(bam_path_str: str) -> bytes:
    """Scan one partition; return the loci DF as Arrow IPC bytes (cheap to ship)."""
    import io
    df = _scan_partition_loci(
        bam_path_str,
        _BUILD_CTX["frame_ivs"],
        _BUILD_CTX["bam_refs"],
        _BUILD_CTX["default_offset"],
    )
    buf = io.BytesIO()
    df.write_ipc(buf)
    return buf.getvalue()


def build_read_loci_index(
    partition_dirs: "str | Path | List[str | Path]",
    cds_df: pl.DataFrame,
    *,
    default_offset: int = 15,
    out_path: Optional[str | Path] = None,
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """Phase A, paid once: scan every partition BAM, assign CDS frames, and
    return/persist a ``read_id → (length, frame)`` index.

    This is the sample-independent half of periodicity.  Persist it next to the
    matrix (``out_path``) and every later query is a cheap join via
    ``tabulate_frames`` — no BAM rescans.

    The 256 partitions are independent, so the scan parallelises across
    ``n_workers`` processes (default: all cores).  Per-partition work is the
    exact serial ``_scan_partition_loci``, so the result is identical regardless
    of worker count.
    """
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir() and _discover_bam(d))
    else:
        dirs = [Path(d) for d in partition_dirs if _discover_bam(Path(d))]
    if not dirs:
        return pl.DataFrame(schema={"read_id": pl.UInt64, "length": pl.Int64, "frame": pl.Int8})

    bam_refs = _bam_chroms(_discover_bam(dirs[0]))
    frame_ivs = _build_frame_intervals(cds_df, cds_df, bam_refs)
    if n_workers is None:
        n_workers = min(len(dirs), os.cpu_count() or 1)
    log_info(f"read_loci index: {len(dirs)} partitions, "
             f"{sum(len(v) for v in frame_ivs.values()):,} CDS intervals, "
             f"{n_workers} worker(s)")

    bam_paths = [str(_discover_bam(d)) for d in dirs]
    parts: List[pl.DataFrame] = []
    empty_schema = {"read_id": pl.UInt64, "length": pl.Int64, "frame": pl.Int8}

    if n_workers <= 1:
        for i, bp in enumerate(bam_paths):
            df = _scan_partition_loci(bp, frame_ivs, bam_refs, default_offset)
            if not df.is_empty():
                parts.append(df)
            if (i + 1) % 32 == 0:
                log_info(f"  scanned {i + 1}/{len(dirs)} partitions")
    else:
        import io
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            processes=n_workers,
            initializer=_build_worker_init,
            initargs=(frame_ivs, bam_refs, default_offset),
        ) as pool:
            for i, ipc in enumerate(pool.imap_unordered(_build_worker_scan, bam_paths, chunksize=1)):
                df = pl.read_ipc(io.BytesIO(ipc))
                if not df.is_empty():
                    parts.append(df)
                if (i + 1) % 32 == 0:
                    log_info(f"  scanned {i + 1}/{len(dirs)} partitions")

    loci = pl.concat(parts) if parts else pl.DataFrame(schema=empty_schema)
    if out_path is not None:
        loci.write_parquet(str(out_path))
        log_info(f"read_loci index written: {out_path} ({loci.height:,} reads)")
    return loci


def tabulate_frames(
    partition_dirs: "str | Path | List[str | Path]",
    read_loci: "pl.DataFrame | str | Path",
    *,
    group_level: str = "sample",
    sample_names: Optional[List[str]] = None,
    cluster_labels: Optional[pl.DataFrame] = None,
) -> pl.DataFrame:
    """Phase B, vectorised: join the cached read-locus index against the count
    Parquets and group to a per-(group, length, frame) count table.

    group_level:
        "sample"    → one group per sample_name  (per-sample QC)
        "aggregate" → single collapsed pseudo-sample
        "cluster"   → groups from ``cluster_labels`` (sample_name → cluster_id)

    Returns long-form: [<group col>,] length, frame, count.  The same scan
    serves all three levels — this is the aggregate/cluster/individual
    "same pass" primitive.
    """
    if isinstance(read_loci, (str, Path)):
        read_loci = pl.read_parquet(str(read_loci))
    if isinstance(partition_dirs, (str, Path)):
        parent = Path(partition_dirs)
        dirs = sorted(d for d in parent.iterdir() if d.is_dir())
    else:
        dirs = [Path(d) for d in partition_dirs]

    # sample lookup (built once from the first partition, as periodicity_qc does)
    mp0, manifest0 = _manifest(dirs[0])
    samples = _samples_df(mp0, manifest0).select(["sample_id", "sample_name"])
    if sample_names:
        samples = samples.filter(pl.col("sample_name").is_in(sample_names))

    count_files: List[str] = []
    for pdir in dirs:
        try:
            mp, manifest = _manifest(pdir)
        except FileNotFoundError:
            continue
        count_files.extend(_count_parquets(mp, manifest))
    if not count_files:
        return pl.DataFrame(schema={"sample_name": pl.Utf8, "length": pl.Int64, "frame": pl.Int8, "count": pl.Float64})

    counts = pl.scan_parquet(count_files).select(["read_id", "sample_id", "count"])
    joined = (
        counts
        .join(read_loci.lazy(), on="read_id", how="inner")
        .join(samples.lazy(), on="sample_id", how="inner")
        .with_columns(pl.col("count").cast(pl.Float64))
    )

    if group_level == "aggregate":
        keys = ["length", "frame"]
    elif group_level == "cluster":
        if cluster_labels is None:
            raise ValueError("group_level='cluster' requires cluster_labels (sample_name, cluster_id)")
        joined = joined.join(cluster_labels.lazy(), on="sample_name", how="inner")
        keys = ["cluster_id", "length", "frame"]
    else:  # sample
        keys = ["sample_name", "length", "frame"]

    out = (
        joined.group_by(keys)
        .agg(pl.col("count").sum().alias("count"))
        .collect(engine="streaming")
    )
    return out


def periodicity_qc_fast(
    partition_dirs: "str | Path | List[str | Path]",
    cds_df: pl.DataFrame,
    *,
    read_loci: "pl.DataFrame | str | Path | None" = None,
    sample_names: Optional[List[str]] = None,
    default_offset: int = 15,
) -> pl.DataFrame:
    """Drop-in equivalent of ``periodicity_qc`` using the two-phase fast path.

    If ``read_loci`` is None the index is built on the fly (one BAM pass);
    pass a prebuilt index (DataFrame or Parquet path) to skip the BAM scan
    entirely.  Output schema matches ``periodicity_qc`` so ``frame_dominance_matrix``
    works unchanged.
    """
    if read_loci is None:
        read_loci = build_read_loci_index(partition_dirs, cds_df, default_offset=default_offset)

    tab = tabulate_frames(partition_dirs, read_loci, group_level="sample", sample_names=sample_names)
    if tab.is_empty():
        return _empty_periodicity_schema()

    # Rebuild {sample: {length: {0,1,2}}} from the tiny grouped table, then reuse
    # the existing scoring/emission so output is identical to periodicity_qc.
    rfd: Dict[str, Dict[int, Dict[int, float]]] = {}
    for row in tab.iter_rows(named=True):
        s = str(row["sample_name"])
        L = int(row["length"])
        rfd.setdefault(s, {}).setdefault(L, {0: 0.0, 1: 0.0, 2: 0.0})[int(row["frame"])] += float(row["count"])

    base_offsets = {length: default_offset for lfd in rfd.values() for length in lfd}
    rows = []
    for sname, lfd in rfd.items():
        scores = _ribometric_frame_scores(lfd)
        offsets = _nudge_to_frame0(lfd, base_offsets)
        rfd_json, dominance_json = {}, {}
        for length, fd in lfd.items():
            counts = [float(fd.get(0, 0.0)), float(fd.get(1, 0.0)), float(fd.get(2, 0.0))]
            rfd_json[str(length)] = counts
            tot = sum(counts)
            dominance_json[str(length)] = (max(counts) / tot) if tot > 0 else 0.0
        rows.append({
            "sample_id": sname,
            "periodicity_score": scores["periodicity_score"],
            "f0": scores["f0"], "f1": scores["f1"], "f2": scores["f2"],
            "n_cds_reads": scores["n_reads"],
            "recommended_offsets": json.dumps({str(k): v for k, v in offsets.items()}),
            "read_frame_distribution": json.dumps(rfd_json),
            "per_length_periodicity": json.dumps({str(k): v for k, v in scores["per_length"].items()}),
            "per_length_dominance": json.dumps(dominance_json),
        })
    return pl.from_dicts(rows).sort("sample_id") if rows else _empty_periodicity_schema()


# ---------------------------------------------------------------------------
# Combined entry point
# ---------------------------------------------------------------------------

def matrix_qc(
    partitions: "str | Path | List[str | Path]",
    *,
    exon_df: Optional[pl.DataFrame] = None,
    cds_df: Optional[pl.DataFrame] = None,
    sample_names: Optional[List[str]] = None,
    default_offset: int = 15,
) -> pl.DataFrame:
    """Run length QC and (optionally) periodicity QC across all partitions.

    partitions: parent directory whose subdirectories are partitions, or an
                explicit list of partition directories.

    Always runs fast_length_qc (aggregated across all partitions).
    Runs periodicity_qc when exon_df is provided — uses all partitions to
    build a representative read_frame_distribution per sample.

    cds_df (from getexons_and_cds()) enables accurate exon-aware frame
    assignment; without it frame assignment falls back to genomic position % 3.
    """
    # Resolve list of partition directories
    if isinstance(partitions, (str, Path)):
        parent = Path(partitions)
        dirs: List[Path] = sorted(d for d in parent.iterdir() if d.is_dir())
    else:
        dirs = [Path(d) for d in partitions]

    if not dirs:
        raise FileNotFoundError(f"No partition directories found in {partitions}")

    log_info(f"Length QC across {len(dirs)} partitions…")
    length_parts = [fast_length_qc(d, sample_names=sample_names) for d in dirs]
    length_parts = [p for p in length_parts if not p.is_empty()]
    if not length_parts:
        return pl.DataFrame()

    # Aggregate per-sample stats across partitions
    combined = pl.concat(length_parts)
    base = (
        combined
        .group_by(["sample_id", "study_id"])
        .agg([
            pl.col("total_reads").sum(),
            pl.col("unique_reads").sum(),
            (pl.col("total_reads").sum() / pl.col("unique_reads").sum()).alias("compression_ratio"),
            # Weighted mean/peak/rpf from the per-partition rows
            (
                (pl.col("mean_length") * pl.col("total_reads")).sum()
                / pl.col("total_reads").sum()
            ).alias("mean_length"),
            pl.col("peak_length").mode().first(),
            (
                (pl.col("rpf_28_32_prop") * pl.col("total_reads")).sum()
                / pl.col("total_reads").sum()
            ).alias("rpf_28_32_prop"),
            pl.col("length_cv").mean(),
        ])
        .sort("sample_id")
    )

    if exon_df is not None:
        log_info("Running periodicity QC across all partitions…")
        perio = periodicity_qc(
            dirs,
            exon_df,
            cds_df=cds_df,
            sample_names=sample_names,
            default_offset=default_offset,
        )
        if not perio.is_empty():
            base = base.join(perio, on="sample_id", how="left")

    return base

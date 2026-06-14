"""Extract and deduplicate genomic events from translon annotations.

Translons (features) decompose into shared genomic events; scoring each event
once is ~50× less work and keeps shared evidence consistent.

Public API
----------
extract_events  — blocks + translons → (events, feature_event, event_overlap)
_contention      — elongation events that overlap in a different frame
_deconflict_intervals  — remove CDS positions with conflicting frame phases
_build_frame_intervals — per-(chrom,strand) CDS exon intervals with phase
"""
from __future__ import annotations

import collections
import heapq
from typing import Dict, List, Optional, Set, Tuple

import polars as pl

from TranslonScorer.io.bam import normalise_chrom as _normalise_chrom


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mod3(expr: pl.Expr) -> pl.Expr:
    return ((expr % 3) + 3) % 3


def _eid(key: pl.Expr) -> pl.Expr:
    """Deterministic u64 event id from a canonical key string."""
    return key.hash(seed=0).alias("event_id")


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------

def extract_events(
    blocks: pl.DataFrame,
    translons: pl.DataFrame,
    *,
    annotation_version: str = "",
) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Return (events, feature_event, event_overlap).

    blocks columns:    translon_id, translation_block_rank, bed_chrom,
                       bed_start, bed_end, seq_region_strand, block_length_nt
    translons columns: translon_id, bed_chrom, bed_start, bed_end, seq_region_strand
    """
    b = blocks.rename({"seq_region_strand": "strand", "bed_chrom": "chrom"}).with_columns([
        pl.col("bed_start").cast(pl.Int64),
        pl.col("bed_end").cast(pl.Int64),
        pl.col("strand").cast(pl.Int64),
        pl.col("block_length_nt").cast(pl.Int64),
        pl.col("translation_block_rank").cast(pl.Int64),
    ])

    # Phase per block (CDS-relative 5' offset = cumulative prior length)
    b = b.sort(["translon_id", "translation_block_rank"]).with_columns(
        (pl.col("block_length_nt").cum_sum().over("translon_id") - pl.col("block_length_nt")).alias("ts")
    ).with_columns(
        pl.when(pl.col("strand") > 0)
        .then(_mod3(pl.col("ts") - pl.col("bed_start")))
        .otherwise(_mod3(pl.col("ts") + pl.col("bed_end") - 1))
        .alias("phase")
    )

    # Elongation events
    b = b.with_columns(
        pl.concat_str([
            pl.lit("E"), pl.col("chrom"), pl.col("strand").cast(pl.Utf8),
            pl.col("bed_start").cast(pl.Utf8), pl.col("bed_end").cast(pl.Utf8),
            pl.col("phase").cast(pl.Utf8),
        ], separator="|").alias("_k")
    ).with_columns(_eid(pl.col("_k")))
    elong_events = b.unique("event_id").select(
        "event_id", pl.lit("elongation").alias("type"), "chrom", "strand",
        pl.col("bed_start").alias("start"), pl.col("bed_end").alias("end"), "phase",
    )
    fe_elong = b.select(
        pl.col("translon_id").alias("feature_id"), "event_id",
        pl.lit("elongation").alias("role"),
        pl.col("translation_block_rank").alias("rank"), "phase",
    )

    # Junction events (genomic adjacencies within a translon)
    bj = b.sort(["translon_id", "bed_start"]).with_columns(
        pl.col("bed_start").shift(-1).over("translon_id").alias("next_start")
    ).filter(pl.col("next_start").is_not_null() & (pl.col("next_start") > pl.col("bed_end")))
    if bj.height:
        bj = bj.with_columns(
            pl.concat_str([
                pl.lit("J"), pl.col("chrom"), pl.col("strand").cast(pl.Utf8),
                pl.col("bed_end").cast(pl.Utf8), pl.col("next_start").cast(pl.Utf8),
            ], separator="|").alias("_k")
        ).with_columns(_eid(pl.col("_k")))
        junc_events = bj.unique("event_id").select(
            "event_id", pl.lit("junction").alias("type"), "chrom", "strand",
            pl.col("bed_end").alias("start"), pl.col("next_start").alias("end"),
            pl.lit(None, dtype=pl.Int64).alias("phase"),
        )
        fe_junc = bj.select(
            pl.col("translon_id").alias("feature_id"), "event_id",
            pl.lit("junction").alias("role"), pl.col("translation_block_rank").alias("rank"),
            pl.lit(None, dtype=pl.Int64).alias("phase"),
        )
    else:
        junc_events = elong_events.clear()
        fe_junc = fe_elong.clear()

    # Init / term events (translon 5' / 3' ends)
    t = translons.rename({"seq_region_strand": "strand", "bed_chrom": "chrom"}).with_columns([
        pl.col("bed_start").cast(pl.Int64), pl.col("bed_end").cast(pl.Int64), pl.col("strand").cast(pl.Int64),
    ]).with_columns([
        pl.when(pl.col("strand") > 0).then(pl.col("bed_start")).otherwise(pl.col("bed_end") - 1).alias("init_pos"),
        pl.when(pl.col("strand") > 0).then(pl.col("bed_end") - 1).otherwise(pl.col("bed_start")).alias("term_pos"),
    ])

    def _point_events(pos_col: str, typ: str, prefix: str):
        tt = t.with_columns(
            pl.concat_str([pl.lit(prefix), pl.col("chrom"), pl.col("strand").cast(pl.Utf8),
                           pl.col(pos_col).cast(pl.Utf8)], separator="|").alias("_k")
        ).with_columns(_eid(pl.col("_k")))
        ev = tt.unique("event_id").select(
            "event_id", pl.lit(typ).alias("type"), "chrom", "strand",
            pl.col(pos_col).alias("start"), (pl.col(pos_col) + 1).alias("end"),
            pl.lit(None, dtype=pl.Int64).alias("phase"),
        )
        fe = tt.select(pl.col("translon_id").alias("feature_id"), "event_id",
                       pl.lit(typ).alias("role"), pl.lit(0, dtype=pl.Int64).alias("rank"),
                       pl.lit(None, dtype=pl.Int64).alias("phase"))
        return ev, fe

    init_events, fe_init = _point_events("init_pos", "init", "I")
    term_events, fe_term = _point_events("term_pos", "term", "T")

    events = pl.concat([elong_events, junc_events, init_events, term_events]).unique("event_id")
    if annotation_version:
        events = events.with_columns(pl.lit(annotation_version).alias("annotation_version"))
    feature_event = pl.concat([fe_elong, fe_junc, fe_init, fe_term])

    overlap = _contention(elong_events)
    return events, feature_event, overlap


def run_extract(sqlite_path: str, out_dir: str, *, annotation_version: str = "") -> dict:
    """Genome-wide driver: extract events per chromosome, write Parquet."""
    import sqlite3
    from pathlib import Path

    out = Path(out_dir)
    (out / "events").mkdir(parents=True, exist_ok=True)
    (out / "feature_event").mkdir(parents=True, exist_ok=True)
    (out / "event_overlap").mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(sqlite_path)
    chroms = [r[0] for r in con.execute(
        "SELECT DISTINCT bed_chrom FROM translons WHERE bed_chrom IS NOT NULL").fetchall()]

    n_ev = n_fe = n_ov = 0
    by_type: dict = {}
    for chrom in chroms:
        blocks = pl.read_database(
            "SELECT translon_id, translation_block_rank, bed_chrom, bed_start, bed_end, "
            "seq_region_strand, block_length_nt FROM translon_blocks WHERE bed_chrom = ?",
            con, execute_options={"parameters": [chrom]})
        trans = pl.read_database(
            "SELECT translon_id, bed_chrom, bed_start, bed_end, seq_region_strand "
            "FROM translons WHERE bed_chrom = ?",
            con, execute_options={"parameters": [chrom]})
        if blocks.is_empty():
            continue
        events, fe, overlap = extract_events(blocks, trans, annotation_version=annotation_version)
        safe = str(chrom).replace("/", "_")
        events.write_parquet(out / "events" / f"{safe}.parquet")
        fe.write_parquet(out / "feature_event" / f"{safe}.parquet")
        overlap.write_parquet(out / "event_overlap" / f"{safe}.parquet")
        n_ev += events.height; n_fe += fe.height; n_ov += overlap.height
        for t_type, c in events.group_by("type").len().iter_rows():
            by_type[t_type] = by_type.get(t_type, 0) + c
    con.close()
    return {"events": n_ev, "feature_event": n_fe, "event_overlap": n_ov,
            "by_type": by_type, "chroms": len(chroms)}


# ---------------------------------------------------------------------------
# Contention: elongation events overlapping in a different frame
# ---------------------------------------------------------------------------

def _contention(elong: pl.DataFrame) -> pl.DataFrame:
    """Elongation events that overlap another in a DIFFERENT reading-frame
    register → convoluted signal. Returns (event_id, other_event_id,
    overlap_start, overlap_end). Sweep-line per (chrom, strand)."""
    if elong.is_empty():
        return pl.DataFrame(schema={"event_id": pl.UInt64, "other_event_id": pl.UInt64,
                                    "overlap_start": pl.Int64, "overlap_end": pl.Int64})
    e = elong.with_columns(pl.col("phase").alias("reg"))
    rows = []
    for (_chrom, _strand), grp in e.group_by(["chrom", "strand"]):
        recs = grp.select(["event_id", "start", "end", "reg"]).sort("start").rows()
        active: list = []
        for eid, s, en, reg in recs:
            active = [a for a in active if a[0] > s]
            for aend, aid, astart, areg in active:
                if areg != reg:
                    os, oe = max(s, astart), min(en, aend)
                    if oe > os:
                        rows.append((eid, aid, os, oe))
                        rows.append((aid, eid, os, oe))
            active.append((en, eid, s, reg))
    if not rows:
        return pl.DataFrame(schema={"event_id": pl.UInt64, "other_event_id": pl.UInt64,
                                    "overlap_start": pl.Int64, "overlap_end": pl.Int64})
    return pl.DataFrame(rows, schema=["event_id", "other_event_id", "overlap_start", "overlap_end"],
                        orient="row")


# ---------------------------------------------------------------------------
# CDS frame interval utilities (annotation-derived, pure)
# ---------------------------------------------------------------------------

def _deconflict_intervals(
    intervals: List[Tuple[int, int, int]],
) -> List[Tuple[int, int, int]]:
    """Remove positions where overlapping CDS intervals disagree on phase.

    Returns a sorted list of non-overlapping intervals where every overlapping
    input interval shares the same phase % 3.
    """
    if not intervals:
        return []

    events: List[Tuple[int, int, int]] = []
    for es, ee, ph in intervals:
        events.append((es, 0, ph))
        events.append((ee, 1, ph))
    events.sort()

    result: List[Tuple[int, int, int]] = []
    active: collections.Counter = collections.Counter()
    prev_pos: Optional[int] = None

    for pos, typ, phase in events:
        if prev_pos is not None and prev_pos < pos and active:
            distinct = set(active.keys())
            if len(distinct) == 1:
                p = next(iter(distinct))
                result.append((prev_pos, pos, p))
        if typ == 0:
            active[phase] += 1
        else:
            active[phase] -= 1
            if active[phase] == 0:
                del active[phase]
        prev_pos = pos

    merged: List[List[int]] = []
    for es, ee, ph in result:
        if merged and merged[-1][2] == ph and merged[-1][1] == es:
            merged[-1][1] = ee
        else:
            merged.append([es, ee, ph])

    return [(a, b, c) for a, b, c in merged]


def _build_frame_intervals(
    exon_df: pl.DataFrame,   # unused — cds_df already carries CDS exon structure
    cds_df: pl.DataFrame,
    bam_refs: set,
) -> Dict[Tuple[str, str], List[Tuple[int, int, int]]]:
    """Per-(chrom,strand) CDS exon intervals with pre-computed frame phases.

    Returns {(chrom, strand): [(es, ee, phase), ...]} sorted by es, after
    deconfliction of positions where isoforms disagree on reading frame.
    """
    raw: Dict[Tuple[str, str], List[Tuple[int, int, int]]] = {}

    for row in cds_df.iter_rows(named=True):
        chrom = str(row["chr"])
        norm = _normalise_chrom(chrom, bam_refs)
        if norm is None:
            continue
        strand  = str(row["strand"])
        starts  = row["start"]
        stops   = row["stop"]
        ts_list = row["tran_start"]
        key = (norm, strand)
        bucket = raw.setdefault(key, [])
        for es, ee, ts in zip(starts, stops, ts_list):
            if strand == "+":
                phase = int(ts - es) % 3
            else:
                phase = int(ts + ee - 1) % 3
            bucket.append((int(es), int(ee), phase))

    return {
        k: _deconflict_intervals(sorted(v, key=lambda x: x[0]))
        for k, v in raw.items()
    }

from __future__ import annotations

from typing import Any

import polars as pl

from ..utils.logging import log_info, log_warning


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    return [value]


def _map_genomic_to_transcript(
    exons: list[tuple[int, int, int, int]],
    strand: str,
    genomic_pos: int,
) -> int | None:
    """Map one 0-based genomic coordinate to transcript coordinate."""
    for exon_start, exon_stop, tran_start, _tran_stop in exons:
        if exon_start <= genomic_pos < exon_stop:
            if strand == "-":
                return int(tran_start + (exon_stop - genomic_pos - 1))
            return int(tran_start + (genomic_pos - exon_start))
    return None


def cds_to_transcript_space(cds_df: pl.DataFrame, exon_df: pl.DataFrame) -> pl.DataFrame:
    """Convert CDS intervals from genomic coordinates to transcript coordinates.

    `getexons_and_cds()` stores CDS `start`/`stop` as genomic half-open intervals,
    while frame support and transcript profiles use transcript positions. This
    helper returns `tran_id,start,stop[,chr,strand]` in transcript coordinates.

    If coordinates cannot be mapped, the original `cds_df` is returned unchanged.
    That fallback keeps synthetic/unit inputs usable when they already provide
    transcript-space CDS intervals.
    """
    if cds_df.is_empty() or exon_df.is_empty():
        return cds_df

    required = {"tran_id", "start", "stop"}
    exon_required = {"tran_id", "start", "stop", "tran_start", "tran_stop"}
    if not required.issubset(set(cds_df.columns)) or not exon_required.issubset(
        set(exon_df.columns)
    ):
        log_warning("CDS/exon schema lacks coordinate columns; using CDS coordinates unchanged")
        return cds_df

    exon_lookup: dict[str, tuple[str, list[tuple[int, int, int, int]]]] = {}
    for row in exon_df.iter_rows(named=True):
        tid = str(row["tran_id"])
        strand = str(row.get("strand", "+") or "+")
        starts = [int(x) for x in _as_list(row["start"])]
        stops = [int(x) for x in _as_list(row["stop"])]
        tran_starts = [int(x) for x in _as_list(row["tran_start"])]
        tran_stops = [int(x) for x in _as_list(row["tran_stop"])]
        if not (len(starts) == len(stops) == len(tran_starts) == len(tran_stops)):
            continue
        exon_lookup[tid] = (
            strand,
            list(zip(starts, stops, tran_starts, tran_stops)),
        )

    mapped_rows: list[dict[str, Any]] = []
    for row in cds_df.iter_rows(named=True):
        tid = str(row["tran_id"])
        hit = exon_lookup.get(tid)
        if hit is None:
            continue

        strand, exons = hit
        cds_start = int(row["start"])
        cds_stop = int(row["stop"])
        if cds_stop <= cds_start:
            continue

        if strand == "-":
            five_prime_g = cds_stop - 1
            three_prime_g = cds_start
        else:
            five_prime_g = cds_start
            three_prime_g = cds_stop - 1

        tran_start = _map_genomic_to_transcript(exons, strand, five_prime_g)
        tran_stop_last = _map_genomic_to_transcript(exons, strand, three_prime_g)
        if tran_start is None or tran_stop_last is None:
            continue

        mapped = dict(row)
        mapped["start"] = int(tran_start)
        mapped["stop"] = int(tran_stop_last) + 1
        if mapped["start"] > mapped["stop"]:
            mapped["start"], mapped["stop"] = mapped["stop"], mapped["start"]
        mapped_rows.append(mapped)

    if not mapped_rows:
        log_warning("No CDS intervals mapped to transcript space; using CDS coordinates unchanged")
        return cds_df

    mapped_df = pl.from_dicts(mapped_rows)
    if mapped_df.height < max(1, int(0.5 * cds_df.height)):
        log_warning(
            "Only a minority of CDS intervals mapped to transcript space; "
            "using original CDS coordinates to avoid partial-coordinate mixing"
        )
        return cds_df

    log_info(f"Mapped {mapped_df.height} CDS intervals to transcript coordinates")
    return mapped_df

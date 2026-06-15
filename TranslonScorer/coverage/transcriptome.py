"""Transcript→genome coordinate projection for transcriptome-aligned reads.

A read aligned to the spliced transcript (mRNA) carries a *transcript*
coordinate; to score events in genome space we must walk the transcript's exon
structure and emit the *genomic* position of each P/A-site.  Reads that align to
several isoforms (isoform multimappers) are resolved by projecting every
alignment and checking whether they converge on a single genomic site (shared
exon) before counting.

All functions here are pure (no I/O).  The exon structure is supplied as a
DataFrame with one row per transcript:

    tran_id : str
    chr     : str
    strand  : "+" | "-"
    start   : list[int]   genomic exon starts (0-based), transcription order
    stop    : list[int]   genomic exon ends   (exclusive), transcription order
    tran_start : list[int]  cumulative transcript offset at each exon's 5' end

(start/stop are ordered 5'→3' in transcription direction: left→right for + ,
right→left for - , matching io.annotation.build_cds_blocks.)

Public API
----------
build_exon_index   — exon_df → {tran_id: ExonModel}
project_to_genome  — (ExonModel, transcript_pos) → (strand∈{1,-1}, genomic_pos) | None
"""

from __future__ import annotations

from typing import Dict, List, NamedTuple, Optional, Tuple

import polars as pl


class ExonModel(NamedTuple):
    chrom: str
    strand: int  # +1 / -1
    # per exon, transcription order: (tran_start, tran_end, g_start, g_stop)
    blocks: Tuple[Tuple[int, int, int, int], ...]


def build_exon_index(exon_df: pl.DataFrame) -> Dict[str, ExonModel]:
    """Build {tran_id: ExonModel} from a per-transcript exon DataFrame."""
    required = {"tran_id", "chr", "strand", "start", "stop", "tran_start"}
    missing = required - set(exon_df.columns)
    if missing:
        raise ValueError(f"exon_df missing columns for projection: {sorted(missing)}")

    index: Dict[str, ExonModel] = {}
    for row in exon_df.iter_rows(named=True):
        strand = 1 if str(row["strand"]) in ("+", "1") else -1
        starts = list(row["start"])
        stops = list(row["stop"])
        tstarts = list(row["tran_start"])
        blocks: List[Tuple[int, int, int, int]] = []
        for gs, ge, ts in zip(starts, stops, tstarts):
            blocks.append((int(ts), int(ts) + (int(ge) - int(gs)), int(gs), int(ge)))
        index[str(row["tran_id"])] = ExonModel(str(row["chr"]), strand, tuple(blocks))
    return index


def project_to_genome(model: ExonModel, transcript_pos: int) -> Optional[Tuple[int, int]]:
    """Map a transcript coordinate to (strand, genomic_pos).

    Returns None when the coordinate falls outside the transcript's exons.
    """
    for ts, te, gs, ge in model.blocks:
        if ts <= transcript_pos < te:
            offset = transcript_pos - ts
            if model.strand > 0:
                return 1, gs + offset
            return -1, (ge - 1) - offset
    return None

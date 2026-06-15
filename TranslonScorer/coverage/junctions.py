from __future__ import annotations

"""Utilities to aggregate splice junction counts from BAM or splits index."""

from typing import Dict, Tuple

import polars as pl


def aggregate_bam_junctions(bam_path: str) -> pl.DataFrame:
    """
    Parse a BAM and aggregate junction counts (CIGAR 'N') per (chr, donor, acceptor, strand).

    Returns a Polars DataFrame with columns: chr, donor_pos, acceptor_pos, strand, count
    """
    import pysam

    # CIGAR operation codes
    OP_M = 0  # M
    OP_I = 1  # I
    OP_D = 2  # D
    OP_N = 3  # N (splice)
    OP_S = 4  # S
    OP_H = 5  # H
    OP_P = 6  # P
    OP_EQ = 7  # =
    OP_X = 8  # X

    counts: Dict[Tuple[str, int, int, str], int] = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.cigartuples is None:
                continue
            chr_name = bam.get_reference_name(aln.reference_id)
            pos = int(aln.reference_start)
            strand = "-" if aln.is_reverse else "+"
            for op, length in aln.cigartuples:
                if op in (OP_M, OP_D, OP_EQ, OP_X):
                    pos += int(length)
                elif op == OP_N:
                    donor = pos
                    acceptor = pos + int(length)
                    key = (chr_name, donor, acceptor, strand)
                    counts[key] = counts.get(key, 0) + 1
                    pos += int(length)
                elif op in (OP_I, OP_S, OP_H, OP_P):
                    # Does not consume reference
                    continue

    if not counts:
        return pl.DataFrame(
            {"chr": [], "donor_pos": [], "acceptor_pos": [], "strand": [], "count": []}
        )

    rows = [
        {"chr": c, "donor_pos": d, "acceptor_pos": a, "strand": s, "count": cnt}
        for (c, d, a, s), cnt in counts.items()
    ]
    return pl.from_dicts(rows)

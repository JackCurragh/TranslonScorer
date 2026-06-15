"""BAM I/O adapter: reads via oxbow/pysam, CIGAR block extraction, read_id parsing.

Pure I/O: open → read → close, returning DataFrames / basic Python structures.
No state, no offset computation, no profile generation.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

import polars as pl

try:
    import oxbow as ox  # optional fast BAM reader
except Exception:
    ox = None

from ..utils.logging import log_info, log_warning


# ---------------------------------------------------------------------------
# Read-name parsing (matrix-mode unique-read BAMs: name = "read_<N>[...]")
# ---------------------------------------------------------------------------

_READ_ID_RE = re.compile(r"read_(\d+)")


def parse_read_id(qname: str) -> Optional[int]:
    """Extract the integer read_id embedded in a matrix-BAM query name."""
    m = _READ_ID_RE.search(qname)
    return int(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# Chromosome name normalisation
# ---------------------------------------------------------------------------


def normalise_chrom(chrom: str, bam_refs: set) -> Optional[str]:
    """Return the BAM-reference name matching `chrom`, handling chr-prefix mismatches."""
    if chrom in bam_refs:
        return chrom
    alt = f"chr{chrom}" if not chrom.startswith("chr") else chrom[3:]
    return alt if alt in bam_refs else None


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

# SAM CIGAR operations that consume reference bases
_REF_CONSUMING = {0, 2, 3, 7, 8}  # M, D, N, =, X


def cigar_blocks(ref_start: int, cigar_tuples: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Aligned reference blocks from a pysam cigartuples list.

    Returns list of (block_start, block_end) intervals (0-based, half-open)
    that are contiguous aligned segments (no N-skips i.e. no introns).
    """
    blocks: List[Tuple[int, int]] = []
    pos = ref_start
    block_start = pos
    in_block = False
    for op, length in cigar_tuples:
        if op == 3:  # N — intron skip
            if in_block:
                blocks.append((block_start, pos))
                in_block = False
            pos += length
            block_start = pos
        elif op in _REF_CONSUMING:
            if not in_block:
                block_start = pos
                in_block = True
            pos += length
        # I, S, H, P do not consume reference
    if in_block:
        blocks.append((block_start, pos))
    return blocks


def aggregate_junctions(bam_path: str) -> pl.DataFrame:
    """Per-(chr, donor, acceptor, strand) junction counts via CIGAR N-ops."""
    import pysam

    counts: Dict[Tuple[str, int, int, str], int] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.cigartuples is None:
                continue
            chrom = bam.get_reference_name(aln.reference_id)
            pos = int(aln.reference_start)
            strand = "-" if aln.is_reverse else "+"
            for op, length in aln.cigartuples:
                if op in (0, 2, 7, 8):  # M, D, =, X
                    pos += length
                elif op == 3:  # N — intron
                    key = (chrom, pos, pos + length, strand)
                    counts[key] = counts.get(key, 0) + 1
                    pos += length
    if not counts:
        return pl.DataFrame(
            schema={
                "chr": pl.Utf8,
                "donor_pos": pl.Int64,
                "acceptor_pos": pl.Int64,
                "strand": pl.Utf8,
                "count": pl.Int64,
            }
        )
    rows = [
        {"chr": c, "donor_pos": d, "acceptor_pos": a, "strand": s, "count": cnt}
        for (c, d, a, s), cnt in counts.items()
    ]
    return pl.from_dicts(rows)


# ---------------------------------------------------------------------------
# Primary BAM reader (oxbow fast-path; pysam fallback)
# ---------------------------------------------------------------------------


def readbam(
    bampath: str,
    *,
    collapsed: bool = False,
    count_from: Optional[str] = None,
    count_pattern: Optional[str] = None,
    count_tag: Optional[str] = None,
    unique: bool = False,
    include_qname: bool = False,
) -> pl.DataFrame:
    """Read a BAM file into a normalised positional table.

    Returns DataFrame with columns: chr, start, stop, length, strand, count[, qname]
    """
    import os
    import pysam as _pysam

    if not (os.path.exists(f"{bampath}.bai") or os.path.exists(bampath.replace(".bam", ".bai"))):
        log_info("BAM index not found, creating index…")
        _pysam.index(bampath)

    if ox is not None:
        bamfile = ox.read_bam(bampath)
        log_info("BAM read via oxbow")
        df = pl.read_ipc(bamfile)
    else:
        log_info("oxbow unavailable; using pysam")
        rows = []
        with _pysam.AlignmentFile(bampath, "rb") as bam:
            for aln in bam.fetch(until_eof=True):
                if aln.is_unmapped:
                    continue
                rows.append(
                    {
                        "rname": bam.get_reference_name(aln.reference_id),
                        "pos": int(aln.reference_start),
                        "end": int(aln.reference_end),
                        "seq": aln.query_sequence or "",
                        "flag": int(aln.flag),
                        "qname": aln.query_name,
                    }
                )
        if not rows:
            return pl.DataFrame(
                schema={
                    "chr": pl.Utf8,
                    "start": pl.Int64,
                    "stop": pl.Int64,
                    "length": pl.Int64,
                    "strand": pl.Utf8,
                    "count": pl.Int64,
                }
            )
        df = pl.from_dicts(rows)

    # Normalise column names
    renames = {"rname": "chr", "pos": "start", "end": "stop"}
    for old, new in renames.items():
        if old in df.columns:
            df = df.rename({old: new})

    # Read length
    if "length" not in df.columns:
        if "seq" in df.columns:
            df = df.with_columns(
                pl.when(pl.col("seq").is_not_null())
                .then(pl.col("seq").cast(pl.Utf8).str.len_chars())
                .otherwise(0)
                .alias("length")
            )
        elif {"start", "stop"}.issubset(set(df.columns)):
            df = df.with_columns((pl.col("stop") - pl.col("start")).cast(pl.Int64).alias("length"))

    # Strand from flag
    if "flag" in df.columns:
        df = df.with_columns(
            pl.when(pl.col("flag") & 0x10 > 0)
            .then(pl.lit("-"))
            .otherwise(pl.lit("+"))
            .alias("strand")
        ).drop("flag")

    keep = ["chr", "start", "stop", "length", "strand"]
    if include_qname and "qname" in df.columns:
        keep.append("qname")
    df = df.select(keep)

    # Count column
    if collapsed and count_from == "name":
        import re as _re

        pat = _re.compile(count_pattern or r".*_x(?P<count>\d+)$")
        if "qname" in df.columns:
            df = df.with_columns(
                pl.col("qname")
                .map_elements(
                    lambda s: int(pat.match(s).group("count")) if pat.match(s) else 1,
                    return_dtype=pl.Int64,
                )
                .alias("count")
            )
        else:
            log_warning("count_from=name but qname not available; count=1")
            df = df.with_columns(pl.lit(1).alias("count"))
    elif collapsed and count_from == "tag":
        if count_tag:
            try:
                q2c: Dict[str, int] = {}
                with _pysam.AlignmentFile(bampath, "rb") as bam:
                    for aln in bam.fetch(until_eof=True):
                        try:
                            q2c[aln.query_name] = int(aln.get_tag(count_tag))
                        except Exception:
                            q2c[aln.query_name] = 1
                if "qname" in df.columns:
                    df = df.with_columns(
                        pl.col("qname")
                        .map_elements(lambda s: q2c.get(s, 1), return_dtype=pl.Int64)
                        .alias("count")
                    )
                else:
                    log_warning("qname not present; count=1")
                    df = df.with_columns(pl.lit(1).alias("count"))
            except Exception as e:
                log_warning(f"Tag count parse failed: {e}; count=1")
                df = df.with_columns(pl.lit(1).alias("count"))
        else:
            log_warning("count_from=tag but count_tag not provided; count=1")
            df = df.with_columns(pl.lit(1).alias("count"))
    else:
        df = df.with_columns(pl.lit(1).alias("count"))

    if not unique:
        df = (
            df.group_by(["chr", "start", "stop", "length", "strand"])
            .agg(pl.col("count").sum())
            .sort(["chr", "start"])
        )

    log_info(f"readbam: {len(df)} positions")
    return df

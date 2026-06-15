"""Feature sources → (blocks, translons) for event extraction.

The event scorer is format-agnostic: :func:`TranslonScorer.events.extract_events`
consumes two DataFrames and never touches a file. This module builds those two
DataFrames from whatever the user has — a GTF/GFF, a BED12/bigBed, a FASTA (for
de-novo ORF finding), or the annotation sqlite — so you do **not** need a
pre-built annotation database to score features.

Canonical output schema (matches the sqlite ``translon_blocks``/``translons``):

    blocks    : translon_id, translation_block_rank, bed_chrom, bed_start,
                bed_end, seq_region_strand (+1/-1), block_length_nt
    translons : translon_id, bed_chrom, bed_start, bed_end, seq_region_strand

Block order is 5'→3' in *transcription* direction (left→right on +, right→left
on -), so phase is computed consistently regardless of source.

Public API
----------
from_gtf       — GTF/GFF feature blocks (default feature_type="CDS")
from_bed12     — BED12 file (blockSizes/blockStarts give exon structure)
from_bigbed    — bigBed file (via pyBigWig)
from_fasta     — de-novo ORFs found in FASTA records (single-block, both strands)
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import polars as pl

_BLOCKS_SCHEMA = {
    "translon_id": pl.Utf8,
    "translation_block_rank": pl.Int64,
    "bed_chrom": pl.Utf8,
    "bed_start": pl.Int64,
    "bed_end": pl.Int64,
    "seq_region_strand": pl.Int64,
    "block_length_nt": pl.Int64,
}
_TRANSLONS_SCHEMA = {
    "translon_id": pl.Utf8,
    "bed_chrom": pl.Utf8,
    "bed_start": pl.Int64,
    "bed_end": pl.Int64,
    "seq_region_strand": pl.Int64,
}


def _explode_per_tx(per_tx: pl.DataFrame) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Per-transcript list-column blocks (tran_id, chr, strand, start[], stop[])
    in transcription order → (blocks, translons) in the canonical schema."""
    if per_tx.is_empty():
        return pl.DataFrame(schema=_BLOCKS_SCHEMA), pl.DataFrame(schema=_TRANSLONS_SCHEMA)

    strand_int = pl.when(pl.col("strand") == "-").then(-1).otherwise(1)

    ex = per_tx.select("tran_id", "chr", "strand", "start", "stop").explode(["start", "stop"])
    blocks = ex.with_columns(
        pl.col("start").cum_count().over("tran_id").cast(pl.Int64).alias("translation_block_rank"),
        (pl.col("stop") - pl.col("start")).cast(pl.Int64).alias("block_length_nt"),
        strand_int.alias("seq_region_strand"),
    ).select(
        pl.col("tran_id").alias("translon_id"),
        "translation_block_rank",
        pl.col("chr").alias("bed_chrom"),
        pl.col("start").cast(pl.Int64).alias("bed_start"),
        pl.col("stop").cast(pl.Int64).alias("bed_end"),
        "seq_region_strand",
        "block_length_nt",
    )

    translons = per_tx.select(
        pl.col("tran_id").alias("translon_id"),
        pl.col("chr").alias("bed_chrom"),
        pl.col("start").list.min().cast(pl.Int64).alias("bed_start"),
        pl.col("stop").list.max().cast(pl.Int64).alias("bed_end"),
        strand_int.alias("seq_region_strand"),
    )
    return blocks, translons


# ---------------------------------------------------------------------------
# GTF / GFF
# ---------------------------------------------------------------------------


def from_gtf(gtf_path: str, feature_type: str = "CDS") -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Build events input from a GTF/GFF, scoring the given feature type.

    ``feature_type="CDS"`` scores annotated coding regions; ``"exon"`` scores
    whole transcripts. Exon structure (splicing) is taken from the file.
    """
    from TranslonScorer.io.annotation import _blocks_from_gtf

    per_tx = _blocks_from_gtf(gtf_path, feature_type)
    return _explode_per_tx(per_tx)


# ---------------------------------------------------------------------------
# BED12 / bigBed
# ---------------------------------------------------------------------------


def _bed12_rows_to_per_tx(rows: List[dict]) -> pl.DataFrame:
    """rows: dicts with chrom, chromStart, name, strand, blockSizes, blockStarts.
    Returns per-tx list-column frame in transcription (5'→3') order."""
    recs = []
    for r in rows:
        cs = int(r["chromStart"])
        sizes = [int(x) for x in str(r["blockSizes"]).strip(",").split(",") if x != ""]
        starts = [int(x) for x in str(r["blockStarts"]).strip(",").split(",") if x != ""]
        if not sizes or len(sizes) != len(starts):
            continue
        blocks = [(cs + st, cs + st + sz) for st, sz in zip(starts, sizes)]
        blocks.sort()
        if str(r["strand"]) == "-":  # transcription order = right→left
            blocks.reverse()
        recs.append(
            {
                "tran_id": str(r["name"]),
                "chr": str(r["chrom"]),
                "strand": str(r["strand"]),
                "start": [b[0] for b in blocks],
                "stop": [b[1] for b in blocks],
            }
        )
    if not recs:
        return pl.DataFrame(
            schema={
                "tran_id": pl.Utf8,
                "chr": pl.Utf8,
                "strand": pl.Utf8,
                "start": pl.List(pl.Int64),
                "stop": pl.List(pl.Int64),
            }
        )
    return pl.from_dicts(recs)


def from_bed12(bed_path: str) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Build events input from a BED12 file (12 columns; blockSizes/blockStarts
    define exon structure). Names (col 4) become translon ids."""
    rows: List[dict] = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 12:
                raise ValueError(f"Not a BED12 line (need 12 columns, got {len(f)}): {line[:60]}")
            rows.append(
                {
                    "chrom": f[0],
                    "chromStart": f[1],
                    "name": f[3] or f"{f[0]}:{f[1]}",
                    "strand": f[5],
                    "blockSizes": f[10],
                    "blockStarts": f[11],
                }
            )
    return _explode_per_tx(_bed12_rows_to_per_tx(rows))


def from_bigbed(bigbed_path: str) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Build events input from a bigBed (BED12) file via pyBigWig."""
    try:
        import pyBigWig
    except ImportError as exc:
        raise ImportError(
            "from_bigbed requires pyBigWig (pip install 'TranslonScorer[bigwig]')"
        ) from exc
    bb = pyBigWig.open(bigbed_path)
    if not bb.isBigBed():
        raise ValueError(f"{bigbed_path} is not a bigBed file")
    rows: List[dict] = []
    for chrom, length in bb.chroms().items():
        for start, end, rest in bb.entries(chrom, 0, int(length)) or []:
            f = rest.split("\t")  # BED fields from column 4 onward
            # name, score, strand, thickStart, thickEnd, itemRgb, blockCount, sizes, starts
            if len(f) < 9:
                continue
            rows.append(
                {
                    "chrom": chrom,
                    "chromStart": start,
                    "name": f[0] or f"{chrom}:{start}",
                    "strand": f[2],
                    "blockSizes": f[7],
                    "blockStarts": f[8],
                }
            )
    bb.close()
    return _explode_per_tx(_bed12_rows_to_per_tx(rows))


# ---------------------------------------------------------------------------
# De-novo from sequence (FASTA)
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def from_fasta(
    fasta_path: str,
    *,
    start_codons: Optional[List[str]] = None,
    stop_codons: Optional[List[str]] = None,
    min_len: int = 0,
    max_len: int = 1_000_000,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """De-novo ORF finding in FASTA records, both strands, as single-block
    genomic features (coordinates relative to each record = bed_chrom).

    Intended for genome / contig FASTAs (records are the sequences to score). For
    spliced features from a transcriptome, use a GTF/BED instead so exon
    structure is known.
    """
    import pyfaidx

    from TranslonScorer.core.orffinder import build_codon_automaton, find_orfs

    start_codons = start_codons or ["ATG"]
    stop_codons = stop_codons or ["TAA", "TAG", "TGA"]
    sa = build_codon_automaton(start_codons)
    pa = build_codon_automaton(stop_codons)

    recs = []
    fa = pyfaidx.Fasta(fasta_path)
    for name in fa.keys():
        seq = str(fa[name][:]).upper()
        n = len(seq)
        for strand, s in (("+", seq), ("-", _revcomp(seq))):
            for orf in find_orfs(s, str(name), sa, pa, minlength=min_len, maxlength=max_len):
                a, b = int(orf["start"]), int(orf["stop"])
                if strand == "-":  # map back to forward-strand genomic coords
                    a, b = n - b, n - a
                if b <= a:
                    continue
                recs.append(
                    {
                        "tran_id": f"{name}:{a}-{b}:{strand}",
                        "chr": str(name),
                        "strand": strand,
                        "start": [a],
                        "stop": [b],
                    }
                )
    if not recs:
        return pl.DataFrame(schema=_BLOCKS_SCHEMA), pl.DataFrame(schema=_TRANSLONS_SCHEMA)
    per_tx = pl.from_dicts(recs)
    return _explode_per_tx(per_tx)

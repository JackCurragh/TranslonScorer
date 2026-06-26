"""GTF/BED → CDS blocks and gene spans.

Canonical annotation readers.  Pure I/O adapters: open → parse → close,
return DataFrames / dicts with no side effects.
"""

from __future__ import annotations

import collections
from typing import Dict, List, Tuple

import polars as pl


def _blocks_from_gtf(gtf_path: str, feature_type: str) -> pl.DataFrame:
    """Per-transcript blocks of one GTF feature type (5'->3') with cumulative
    transcript-relative ``tran_start`` measured from the 5' end of the first
    block of that feature type.

    For ``feature_type="exon"`` the origin is the mRNA 5' end (UTRs included);
    for ``"CDS"`` it is the first CDS base. Returns:
    tran_id, gene_id, chr, strand, start[list], stop[list], tran_start[list].
    """
    feats = (
        pl.scan_csv(
            gtf_path,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            schema_overrides={"column_1": pl.Utf8},
        )
        .select(
            pl.col("column_1").alias("chr"),
            pl.col("column_3").alias("type"),
            pl.col("column_4").alias("start"),
            pl.col("column_5").alias("stop"),
            pl.col("column_7").alias("strand"),
            pl.col("column_9").alias("attributes"),
        )
        .filter(pl.col("type") == feature_type)
        .with_columns(
            (pl.col("start").cast(pl.Int64) - 1).alias("start"),
            pl.col("stop").cast(pl.Int64).alias("stop"),
            pl.col("attributes").str.extract(r'transcript_id "([^"]*)"').alias("tran_id"),
            pl.col("attributes").str.extract(r'gene_id "([^"]*)"').alias("gene_id"),
        )
        .collect()
    )
    grouped = (
        feats.group_by("tran_id")
        .agg(
            [
                pl.col("gene_id").first(),
                pl.col("chr").first(),
                pl.col("strand").first(),
                pl.col("start").sort(),
                pl.col("stop").sort(),
            ]
        )
        .with_columns(
            [
                pl.when(pl.col("strand") == "-")
                .then(pl.col("start").list.reverse())
                .otherwise(pl.col("start"))
                .alias("start"),
                pl.when(pl.col("strand") == "-")
                .then(pl.col("stop").list.reverse())
                .otherwise(pl.col("stop"))
                .alias("stop"),
            ]
        )
    )

    def _cumstarts(s) -> list:
        lens = [int(b) - int(a) for a, b in zip(s["start"], s["stop"])]
        acc, out = 0, []
        for L in lens:
            out.append(acc)
            acc += L
        return out

    return grouped.with_columns(
        pl.struct(["start", "stop"])
        .map_elements(_cumstarts, return_dtype=pl.List(pl.Int64))
        .alias("tran_start")
    )


def build_cds_blocks_from_bigbed(bigbed_path: str) -> pl.DataFrame:
    """Per-ORF CDS exon blocks from a BED12/bigBed file (5'->3') with tran_start.

    BED12 block columns (blockSizes, blockStarts relative to chromStart) define
    the exon structure of each ORF.  tran_start is computed as cumulative exon
    lengths in transcription order — identical to the convention used by
    build_cds_blocks() for GTF data, so the result is directly usable by
    build_frame_rollup() / build_coverage_index().

    Returns: tran_id, gene_id, chr, strand, start[list], stop[list], tran_start[list]
    """
    try:
        import pyBigWig
    except ImportError as exc:
        raise ImportError(
            "build_cds_blocks_from_bigbed requires pyBigWig "
            "(pip install 'TranslonScorer[bigwig]')"
        ) from exc

    bb = pyBigWig.open(bigbed_path)
    if not bb.isBigBed():
        raise ValueError(f"{bigbed_path} is not a bigBed file")

    recs = []
    for chrom, length in bb.chroms().items():
        for cs, _ce, rest in bb.entries(chrom, 0, int(length)) or []:
            f = rest.split("\t")
            if len(f) < 9:
                continue
            name = f[0] or f"{chrom}:{cs}"
            strand = f[2]
            sizes = [int(x) for x in f[7].strip(",").split(",") if x]
            starts = [int(x) for x in f[8].strip(",").split(",") if x]
            if not sizes or len(sizes) != len(starts):
                continue
            blocks = sorted((cs + st, cs + st + sz) for st, sz in zip(starts, sizes))
            if strand == "-":
                blocks = list(reversed(blocks))
            # cumulative tran_start from block sizes in transcription order
            acc, tran_starts = 0, []
            for b_start, b_stop in blocks:
                tran_starts.append(acc)
                acc += b_stop - b_start
            recs.append(
                {
                    "tran_id": name,
                    "gene_id": name,
                    "chr": chrom,
                    "strand": strand,
                    "start": [b[0] for b in blocks],
                    "stop": [b[1] for b in blocks],
                    "tran_start": tran_starts,
                }
            )
    bb.close()

    if not recs:
        return pl.DataFrame(
            schema={
                "tran_id": pl.Utf8,
                "gene_id": pl.Utf8,
                "chr": pl.Utf8,
                "strand": pl.Utf8,
                "start": pl.List(pl.Int64),
                "stop": pl.List(pl.Int64),
                "tran_start": pl.List(pl.Int64),
            }
        )
    return pl.from_dicts(recs)


def build_cds_blocks(gtf_path: str) -> pl.DataFrame:
    """Per-transcript CDS exon blocks (5'->3') with CDS-relative tran_start."""
    return _blocks_from_gtf(gtf_path, "CDS")


def build_exon_blocks(gtf_path: str) -> pl.DataFrame:
    """Per-transcript full exon blocks (5'->3') with mRNA-relative tran_start.

    Use this (not build_cds_blocks) for transcriptome→genome read projection:
    transcriptome alignments are positioned from the mRNA 5' end, so UTR reads
    (e.g. 5'UTR uORFs) must project correctly — CDS-relative coordinates would
    place the wrong origin.
    """
    return _blocks_from_gtf(gtf_path, "exon")


def build_gene_spans(
    cds_df: pl.DataFrame,
) -> Tuple[Dict[Tuple[str, str], List[Tuple[int, int, int]]], Dict[int, str]]:
    """Per-(chrom, strand) gene spans (min CDS start → max CDS stop), gene→int code.

    Returns ({(chrom, strand): sorted [(start, stop, gene_code)]}, {gene_code: gene_id}).
    """
    per_tx = (
        cds_df.with_columns(
            [
                pl.col("start").list.min().alias("g_start"),
                pl.col("stop").list.max().alias("g_stop"),
            ]
        )
        .group_by(["chr", "strand", "gene_id"])
        .agg(
            [
                pl.col("g_start").min(),
                pl.col("g_stop").max(),
            ]
        )
    )
    gene_ids = per_tx.get_column("gene_id").to_list()
    code_of = {g: i for i, g in enumerate(sorted(set(gene_ids)))}
    id_of = {i: g for g, i in code_of.items()}

    spans: Dict[Tuple[str, str], List[Tuple[int, int, int]]] = collections.defaultdict(list)
    for row in per_tx.iter_rows(named=True):
        spans[(str(row["chr"]), str(row["strand"]))].append(
            (int(row["g_start"]), int(row["g_stop"]), code_of[row["gene_id"]])
        )
    return {k: sorted(v, key=lambda x: x[0]) for k, v in spans.items()}, id_of

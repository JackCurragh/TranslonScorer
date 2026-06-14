"""GTF/BED → CDS blocks and gene spans.

Canonical annotation readers.  Pure I/O adapters: open → parse → close,
return DataFrames / dicts with no side effects.
"""
from __future__ import annotations

import collections
from typing import Dict, List, Tuple

import polars as pl


def build_cds_blocks(gtf_path: str) -> pl.DataFrame:
    """Per-transcript CDS exon blocks (5'->3') with CDS-relative tran_start.

    Returns: tran_id, gene_id, chr, strand, start[list], stop[list], tran_start[list].
    """
    cds = (
        pl.scan_csv(gtf_path, separator="\t", has_header=False, comment_prefix="#",
                    schema_overrides={"column_1": pl.Utf8})
        .select(
            pl.col("column_1").alias("chr"),
            pl.col("column_3").alias("type"),
            pl.col("column_4").alias("start"),
            pl.col("column_5").alias("stop"),
            pl.col("column_7").alias("strand"),
            pl.col("column_9").alias("attributes"),
        )
        .filter(pl.col("type") == "CDS")
        .with_columns(
            (pl.col("start").cast(pl.Int64) - 1).alias("start"),
            pl.col("stop").cast(pl.Int64).alias("stop"),
            pl.col("attributes").str.extract(r'transcript_id "([^"]*)"').alias("tran_id"),
            pl.col("attributes").str.extract(r'gene_id "([^"]*)"').alias("gene_id"),
        )
        .collect()
    )
    grouped = (
        cds.group_by("tran_id")
        .agg([
            pl.col("gene_id").first(),
            pl.col("chr").first(),
            pl.col("strand").first(),
            pl.col("start").sort(),
            pl.col("stop").sort(),
        ])
        .with_columns([
            pl.when(pl.col("strand") == "-").then(pl.col("start").list.reverse()).otherwise(pl.col("start")).alias("start"),
            pl.when(pl.col("strand") == "-").then(pl.col("stop").list.reverse()).otherwise(pl.col("stop")).alias("stop"),
        ])
    )

    def _cumstarts(s) -> list:
        lens = [int(b) - int(a) for a, b in zip(s["start"], s["stop"])]
        acc, out = 0, []
        for L in lens:
            out.append(acc)
            acc += L
        return out

    return grouped.with_columns(
        pl.struct(["start", "stop"]).map_elements(_cumstarts, return_dtype=pl.List(pl.Int64)).alias("tran_start")
    )


def build_gene_spans(
    cds_df: pl.DataFrame,
) -> Tuple[Dict[Tuple[str, str], List[Tuple[int, int, int]]], Dict[int, str]]:
    """Per-(chrom, strand) gene spans (min CDS start → max CDS stop), gene→int code.

    Returns ({(chrom, strand): sorted [(start, stop, gene_code)]}, {gene_code: gene_id}).
    """
    per_tx = cds_df.with_columns([
        pl.col("start").list.min().alias("g_start"),
        pl.col("stop").list.max().alias("g_stop"),
    ]).group_by(["chr", "strand", "gene_id"]).agg([
        pl.col("g_start").min(),
        pl.col("g_stop").max(),
    ])
    gene_ids = per_tx.get_column("gene_id").to_list()
    code_of = {g: i for i, g in enumerate(sorted(set(gene_ids)))}
    id_of = {i: g for g, i in code_of.items()}

    spans: Dict[Tuple[str, str], List[Tuple[int, int, int]]] = collections.defaultdict(list)
    for row in per_tx.iter_rows(named=True):
        spans[(str(row["chr"]), str(row["strand"]))].append(
            (int(row["g_start"]), int(row["g_stop"]), code_of[row["gene_id"]])
        )
    return {k: sorted(v, key=lambda x: x[0]) for k, v in spans.items()}, id_of

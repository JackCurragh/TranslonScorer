from __future__ import annotations

"""
Build a transcript-mapped index from the genomic read index.

Input: Parquet with columns [read_id, chr, start, stop, strand, length]
Output: Parquet with columns [read_id, tran_id, tran_start_bam, length, strand]

The mapping is sample-independent and can be reused across runs.
"""

import os

import polars as pl

from ..file_handlers import bam as bam_handlers
from ..utils.logging import log_info


def build_mapped_index(
    read_index_parquet: str,
    exon_df: pl.DataFrame,
    *,
    out_parquet: str,
    chunk_size: int = 1_000_000,
) -> str:
    """Stream the read index, map to transcript coords, and write a Parquet index.

    Returns path to out_parquet.
    """
    os.makedirs(os.path.dirname(out_parquet) or ".", exist_ok=True)

    # Determine total rows
    n_rows = pl.scan_parquet(read_index_parquet).select(pl.len()).collect().item()
    log_info(f"Mapped-index build: {n_rows:,} reads; chunk_size={chunk_size}")

    writer = None
    total_mapped = 0
    try:
        for start in range(0, n_rows, chunk_size):
            end = min(start + chunk_size, n_rows)
            idx_df = (
                pl.scan_parquet(read_index_parquet)
                .slice(start, end - start)
                .select(["read_id", "chr", "start", "stop", "strand", "length"])  # noqa: F401
                .collect()
            )
            if idx_df.is_empty():
                continue

            # Add a dummy count to satisfy bamtranscript's expected schema
            idx_df = idx_df.with_columns(pl.lit(1).alias("count"))

            # Map genomic -> transcript coords (sample-independent)
            mapped = bam_handlers.bamtranscript(
                idx_df.select(["chr", "start", "stop", "length", "strand", "count"]).unique(),
                exon_df,
            )
            if mapped.is_empty():
                continue

            # Normalize dtypes for a reliable join (avoid cat vs str mismatches)
            mapped = mapped.with_columns(
                [
                    pl.col("chr").cast(pl.Utf8),
                    pl.col("strand").cast(pl.Utf8),
                    pl.col("start").cast(pl.Int64),
                    pl.col("stop").cast(pl.Int64),
                    pl.col("length").cast(pl.Int64),
                ]
            )
            idx_keys = idx_df.select(
                ["chr", "start", "stop", "length", "strand", "read_id"]
            ).with_columns(
                [
                    pl.col("chr").cast(pl.Utf8),
                    pl.col("strand").cast(pl.Utf8),
                    pl.col("start").cast(pl.Int64),
                    pl.col("stop").cast(pl.Int64),
                    pl.col("length").cast(pl.Int64),
                ]
            )

            # Align chromosome naming: add 'chr' prefix to index if absent
            mapped = mapped.with_columns(pl.col("chr").alias("chr_join"))
            idx_keys = idx_keys.with_columns(
                pl.when(pl.col("chr").str.starts_with("chr"))
                .then(pl.col("chr"))
                .otherwise(pl.lit("chr") + pl.col("chr"))
                .alias("chr_join")
            )

            # Re-attach read_id via genomic keys
            joined = mapped.join(
                idx_keys,
                on=["chr_join", "start", "stop", "length", "strand"],
                how="inner",
            ).select(
                ["read_id", "tran_id", "tran_start_bam", "length", "strand"]
            )  # order cols

            # Append to Parquet
            total_mapped += joined.height
            tbl = joined.to_arrow()
            if writer is None:
                from pyarrow import parquet as pq

                writer = pq.ParquetWriter(out_parquet, schema=tbl.schema)
            writer.write_table(tbl)  # type: ignore

            log_info(f"  chunk {start:,}-{end:,} -> {joined.height:,} mapped rows")
    finally:
        if writer is not None:
            writer.close()  # type: ignore

    log_info(f"Mapped-index written: {out_parquet} ({total_mapped:,} rows)")
    return out_parquet

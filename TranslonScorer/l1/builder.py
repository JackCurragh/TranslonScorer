"""L1 builder: join L0b loci (unique reads) × L0a counts → per-chrom L1 Parquet.

No BAM access — purely Parquet join over:
  - L0b per-chrom shards  (read_id, pos5, strand, length, nh, is_secondary, ...)
  - L0a count parquets    (read_id, sample_id, count)

Two output files per chromosome:
  {chrom}_pos.parquet   — positional stream  (pos5, strand, length, sample_id) → count
  {chrom}_junc.parquet  — junction stream    (donor, acceptor, strand, length, sample_id) → count

Memory model
------------
For each chromosome we:
  1. Scan L0b shard, filter unique reads (nh==1, not secondary), pull minimal columns.
  2. Collect the junction rows (spliced reads only) — small relative to all reads.
  3. Scan all count parquets lazily via glob; join lazily to L0b (positional path)
     and collected-then-joined for the junction path.
  4. Group-by + aggregate with Polars streaming engine (positional) or in-memory
     (junction — explode is not supported in streaming).

Peak memory per worker ≈ (L0b_chrom_unique_pos_rows × 4 int cols)
                         + (all_counts in memory for junction path)
                         + output aggregates.
At full cohort scale the counts glob grows; switch to partition-chunked joins then.
"""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import polars as pl

from ..l0b.contracts import VersionKey, write_meta, read_meta
from .schema import L1_JUNC_SCHEMA, L1_POS_SCHEMA

log = logging.getLogger(__name__)

# Column names used by the unique-read filter on L0b
_L0B_UNIQUE_FILTER = (pl.col("nh") == 1) & (~pl.col("is_secondary"))


# ---------------------------------------------------------------------------
# Per-chrom build helpers
# ---------------------------------------------------------------------------


def _counts_glob(partition_dir: Path, counts_subdir: str) -> str:
    """Return a glob string matching all count Parquet files across partitions."""
    return str(partition_dir / "*" / f"global.*_{counts_subdir}" / "**" / "*.parquet")


def _build_pos_stream(
    l0b_shard: Path,
    counts_lf: pl.LazyFrame,
) -> pl.LazyFrame:
    """Lazy plan for the positional L1 stream of one chromosome.

    Returned as a LazyFrame so the caller can ``sink_parquet`` it — streaming the
    result straight to disk instead of materialising it (and the join hash table)
    in RAM, which OOMs at full cohort scale.
    """
    l0b = (
        pl.scan_parquet(str(l0b_shard))
        .filter(_L0B_UNIQUE_FILTER)
        .select(
            pl.col("read_id").cast(pl.UInt64),
            "pos5",
            "strand",
            "length",
        )
    )

    return (
        l0b.join(counts_lf, on="read_id", how="inner")
        .group_by("pos5", "strand", "length", "sample_id")
        .agg(pl.col("count").sum())
        .sort("pos5")
        .cast({"strand": pl.Int8, "length": pl.Int16, "sample_id": pl.UInt32, "count": pl.UInt32})
    )


def _build_junc_stream(
    l0b_shard: Path,
    counts_lf: pl.LazyFrame,
) -> pl.DataFrame:
    """Build the junction L1 stream for one chromosome.

    Collects the spliced L0b rows first (small — only junction-spanning reads),
    then pulls *only those reads'* counts via a semi-join, instead of
    materialising the entire cohort count matrix (which OOMs at scale).
    """
    spliced = (
        pl.scan_parquet(str(l0b_shard))
        .filter(_L0B_UNIQUE_FILTER & pl.col("junctions_crossed").is_not_null())
        .select(
            pl.col("read_id").cast(pl.UInt64),
            "length",
            "strand",
            "junctions_crossed",
        )
        .collect()
    )

    if spliced.is_empty():
        return pl.DataFrame(schema={f.name: _pa_to_pl(f.type) for f in L1_JUNC_SCHEMA})

    exploded = spliced.explode("junctions_crossed").select(
        "read_id",
        "length",
        "strand",
        pl.col("junctions_crossed").struct.field("donor"),
        pl.col("junctions_crossed").struct.field("acceptor"),
    )

    # Only the spliced reads' counts — semi-join keeps the materialised counts
    # tiny (junction-spanning reads are a small fraction of the cohort).
    read_ids = exploded.select(pl.col("read_id").unique())
    counts_df = counts_lf.join(read_ids.lazy(), on="read_id", how="semi").collect(
        engine="streaming"
    )

    result = (
        exploded.join(counts_df, on="read_id", how="inner")
        .group_by("donor", "acceptor", "strand", "length", "sample_id")
        .agg(pl.col("count").sum())
        .sort("donor", "acceptor")
        .cast({"strand": pl.Int8, "length": pl.Int16, "sample_id": pl.UInt32, "count": pl.UInt32})
    )
    return result


def _pa_to_pl(pa_type) -> pl.PolarsDataType:
    """Minimal PyArrow→Polars type mapping for schema construction."""
    import pyarrow as pa

    mapping = {
        pa.int64(): pl.Int64,
        pa.int8(): pl.Int8,
        pa.int16(): pl.Int16,
        pa.uint32(): pl.UInt32,
    }
    return mapping.get(pa_type, pl.Utf8)


# ---------------------------------------------------------------------------
# Per-chrom worker (runs in subprocess via ProcessPoolExecutor)
# ---------------------------------------------------------------------------


def _chrom_worker(
    chrom: str,
    l0b_shard: str,
    out_dir: str,
    partition_dir: str,
    counts_subdir: str,
) -> tuple[str, int, int]:
    """Build L1 pos + junc shards for one chromosome.

    Returns (chrom, n_pos_rows, n_junc_rows).
    """
    import logging as _logging

    _logging.basicConfig(level=logging.INFO)

    pos_path = Path(out_dir) / f"{chrom}_pos.parquet"
    junc_path = Path(out_dir) / f"{chrom}_junc.parquet"

    if pos_path.exists() and junc_path.exists():
        n_pos = len(pl.read_parquet(str(pos_path), columns=["pos5"]))
        n_junc = len(pl.read_parquet(str(junc_path), columns=["donor"]))
        log.info("checkpoint %s: pos=%d junc=%d", chrom, n_pos, n_junc)
        return chrom, n_pos, n_junc

    log.info("building L1 %s ...", chrom)

    shard = Path(l0b_shard)
    glob_str = _counts_glob(Path(partition_dir), counts_subdir)

    # ------------------------------------------------------------------
    # Positional stream — stream result straight to disk (no materialisation)
    # ------------------------------------------------------------------
    counts_lf = pl.scan_parquet(glob_str).select("read_id", "sample_id", "count")
    _build_pos_stream(shard, counts_lf).sink_parquet(str(pos_path), compression="zstd")
    n_pos = pl.scan_parquet(str(pos_path)).select(pl.len()).collect().item()

    # ------------------------------------------------------------------
    # Junction stream — only spliced reads' counts (semi-join, not full matrix)
    # ------------------------------------------------------------------
    counts_lf = pl.scan_parquet(glob_str).select("read_id", "sample_id", "count")
    junc_df = _build_junc_stream(shard, counts_lf)
    junc_df.write_parquet(str(junc_path), compression="zstd", statistics=True)
    n_junc = len(junc_df)

    log.info("done L1 %s: pos=%d junc=%d", chrom, n_pos, n_junc)
    return chrom, n_pos, n_junc


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def build_l1(
    l0b_dir: Path,
    partition_dir: Path,
    out_dir: Path,
    version_key: VersionKey,
    *,
    chroms: list[str] | None = None,
    workers: int = 4,
    counts_subdir: str = "counts",
) -> Path:
    """Build the L1 raw 5′ cache from L0b loci and L0a counts.

    Parameters
    ----------
    l0b_dir:
        Directory of L0b per-chrom Parquet shards (``{chrom}.parquet``).
    partition_dir:
        Root of the ``global_partitioned`` directory (contains per-prefix subdirs).
    out_dir:
        Where to write L1 shards.  Created if absent.
    version_key:
        The 0.4 version key; written to ``out_dir/_meta.json``.
    chroms:
        Subset of chromosomes to build (default: all ``*_pos.parquet``-able chroms
        found in *l0b_dir*).
    workers:
        Parallel chromosome workers.
    counts_subdir:
        Glob fragment matching the counts subdirectory inside each partition prefix
        (default: ``"counts"``).  The actual glob used is
        ``{partition_dir}/**/global.*_{counts_subdir}/**/*.parquet``.

    Returns
    -------
    Path to the L1 output directory.
    """
    # Validate L0b version key
    l0b_meta = read_meta(l0b_dir)
    if l0b_meta.get("version_key") != version_key.as_dict():
        raise ValueError(
            f"L0b version key mismatch.\n"
            f"  requested: {version_key.as_dict()}\n"
            f"  found:     {l0b_meta.get('version_key')}"
        )

    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Discover chroms from L0b shards
    available = {p.stem: p for p in sorted(l0b_dir.glob("*.parquet")) if not p.stem.startswith("_")}
    target_chroms = chroms if chroms is not None else list(available.keys())
    missing = [c for c in target_chroms if c not in available]
    if missing:
        raise FileNotFoundError(f"L0b shards not found for: {missing}")

    log.info("Building L1 for %d chromosomes", len(target_chroms))
    write_meta(out_dir, version_key, {"status": "building", "n_chroms": len(target_chroms)})

    total_pos = total_junc = 0
    errors: list[str] = []

    # Resolve the counts glob pattern to check it matches something
    glob_check = list((partition_dir).glob(f"*/global.*_{counts_subdir}"))
    if not glob_check:
        raise FileNotFoundError(
            f"No counts directories matching 'global.*_{counts_subdir}' under {partition_dir}"
        )

    # Recycle workers after each chrom so per-process RSS can't accumulate.
    pool_kwargs: dict = {"max_workers": workers}
    try:
        import inspect

        if "max_tasks_per_child" in inspect.signature(ProcessPoolExecutor).parameters:
            pool_kwargs["max_tasks_per_child"] = 1
    except (ValueError, TypeError):
        pass

    with ProcessPoolExecutor(**pool_kwargs) as pool:
        futures = {
            pool.submit(
                _chrom_worker,
                chrom,
                str(available[chrom]),
                str(out_dir),
                str(partition_dir),
                counts_subdir,
            ): chrom
            for chrom in target_chroms
        }
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                _, n_pos, n_junc = fut.result()
                total_pos += n_pos
                total_junc += n_junc
                log.info("✓ %s  pos=%d junc=%d", chrom, n_pos, n_junc)
            except Exception as exc:
                log.error("✗ %s  %s", chrom, exc)
                errors.append(f"{chrom}: {exc}")

    write_meta(
        out_dir,
        version_key,
        {
            "status": "complete" if not errors else "partial",
            "n_chroms_built": len(target_chroms) - len(errors),
            "n_chroms_failed": len(errors),
            "total_pos_rows": total_pos,
            "total_junc_rows": total_junc,
            "errors": errors,
        },
    )

    if errors:
        raise RuntimeError(
            f"L1 build finished with {len(errors)} failed chromosomes:\n" + "\n".join(errors)
        )

    log.info("L1 build complete: %d pos rows, %d junc rows → %s", total_pos, total_junc, out_dir)
    return out_dir

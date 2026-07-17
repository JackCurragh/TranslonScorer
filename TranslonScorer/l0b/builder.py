"""L0b builder: merge partitioned prefix BAMs → per-chrom L0b Parquet shards.

Pipeline
--------
257 prefix BAMs
  → samtools merge/sort (per-chrom, coordinate-sorted)
  → single streaming pysam pass per chrom
  → L0b Parquet shards (one file per chrom, zstd, pos5-sorted)

All chromosomes run in parallel (``--workers``); per-chrom shard presence is
the checkpoint so the job is fully resumable.
"""

from __future__ import annotations

import logging
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Generator, Sequence

import polars as pl
import pysam
import pyarrow as pa
import pyarrow.parquet as pq

from .contracts import L0B_SCHEMA, VersionKey, write_meta

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# CIGAR parsing helpers
# ---------------------------------------------------------------------------

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _is_trivial_cigar(cigar_str: str) -> bool:
    """True for reads whose CIGAR is a single M-only operation (no splicing/clipping)."""
    ops = _CIGAR_RE.findall(cigar_str)
    return len(ops) == 1 and ops[0][1] == "M"


def _junctions_from_cigar(cigar_str: str, ref_start: int) -> list[dict[str, int]]:
    """Return list of {donor, acceptor} dicts for N (intron-skip) operations."""
    junctions: list[dict[str, int]] = []
    pos = ref_start
    for length_str, op in _CIGAR_RE.findall(cigar_str):
        n = int(length_str)
        if op in ("M", "D", "=", "X"):
            pos += n
        elif op == "N":
            donor = pos
            acceptor = pos + n
            junctions.append({"donor": donor, "acceptor": acceptor})
            pos += n
        # S, H, I do not advance the reference
    return junctions


# ---------------------------------------------------------------------------
# Per-record extraction
# ---------------------------------------------------------------------------


def _extract_record(
    rec: pysam.AlignedSegment,
    chrom: str,
) -> dict:
    """Convert one pysam record to a row dict matching L0B_SCHEMA."""
    flag = rec.flag
    is_reverse = bool(flag & 0x10)
    strand: int = -1 if is_reverse else 1

    ref_start: int = rec.reference_start  # 0-based
    ref_end: int = rec.reference_end or (ref_start + rec.query_length)
    read_len: int = rec.query_length or 0

    pos5: int = (ref_end - 1) if is_reverse else ref_start

    cigar_str: str | None = rec.cigarstring
    trivial = cigar_str is None or _is_trivial_cigar(cigar_str)
    stored_cigar: str | None = None if trivial else cigar_str

    junctions: list[dict[str, int]] | None = None
    if cigar_str and not trivial:
        j = _junctions_from_cigar(cigar_str, ref_start)
        if j:
            junctions = j

    # Read NH from tag; fall back to 1 only if tag absent (should not happen).
    try:
        nh: int = rec.get_tag("NH")
    except KeyError:
        nh = 1

    try:
        aln_score: int | None = rec.get_tag("AS")
    except KeyError:
        aln_score = None

    try:
        mismatches: int | None = rec.get_tag("NM")
    except KeyError:
        mismatches = None

    # Read ID: strip leading "read_" prefix, parse as integer.
    qname = rec.query_name or ""
    if qname.startswith("read_"):
        try:
            read_id = int(qname[5:])
        except ValueError:
            read_id = hash(qname) & 0x7FFFFFFFFFFFFFFF
    else:
        read_id = hash(qname) & 0x7FFFFFFFFFFFFFFF

    return {
        "read_id": read_id,
        "chrom": chrom,
        "pos5": pos5,
        "end": ref_end,
        "strand": strand,
        "length": read_len,
        "cigar": stored_cigar,
        "mapq": rec.mapping_quality or 0,
        "nh": nh,
        "is_secondary": bool(flag & 0x100),
        "aln_score": aln_score,
        "mismatches": mismatches,
        "junctions_crossed": junctions,
        "weight": None,
    }


# ---------------------------------------------------------------------------
# Per-chrom streaming scan → Parquet
# ---------------------------------------------------------------------------

_BATCH_SIZE = 100_000  # rows per in-memory batch before flushing to Parquet


def _scan_chrom_bam_to_parquet(
    bam_path: Path,
    chrom: str,
    out_path: Path,
    mask_regions: list[tuple[int, int]] | None = None,
) -> int:
    """Stream one chrom from *bam_path* → write L0b Parquet shard.

    Returns the number of rows written.

    ``mask_regions`` is a sorted list of (start, end) 0-based half-open
    intervals to exclude.  Pass ``None`` (the current default) to disable
    masking.
    """
    rows: list[dict] = []
    total = 0

    writer: pq.ParquetWriter | None = None

    def _flush(final: bool = False) -> None:
        nonlocal writer, total
        if not rows:
            return
        batch = _rows_to_batch(rows)
        if writer is None:
            writer = pq.ParquetWriter(
                str(out_path),
                L0B_SCHEMA,
                compression="zstd",
                write_statistics=True,
            )
        writer.write_batch(batch)
        total += len(rows)
        rows.clear()

    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            # --- contamination mask hook (currently disabled) ---
            # if mask_regions and _in_mask(rec.reference_start, mask_regions):
            #     continue
            rows.append(_extract_record(rec, chrom))
            if len(rows) >= _BATCH_SIZE:
                rows.sort(key=lambda r: r["pos5"])
                _flush()

    if rows:
        rows.sort(key=lambda r: r["pos5"])
    _flush(final=True)

    if writer is not None:
        writer.close()
    elif not out_path.exists():
        # Write an empty file so the checkpoint exists.
        _write_empty_shard(out_path)

    return total


def _rows_to_batch(rows: list[dict]) -> pa.RecordBatch:
    """Convert a list of row dicts to a PyArrow RecordBatch matching L0B_SCHEMA."""
    arrays: list[pa.Array] = []
    for field in L0B_SCHEMA:
        col: list = [r[field.name] for r in rows]
        if field.name == "junctions_crossed":
            arrays.append(_build_junction_array(col))
        else:
            arrays.append(pa.array(col, type=field.type))
    return pa.record_batch(arrays, schema=L0B_SCHEMA)


def _build_junction_array(col: list) -> pa.Array:
    """Build a list<struct<donor,acceptor>> array from a list of lists-of-dicts."""
    struct_type = pa.struct([("donor", pa.int64()), ("acceptor", pa.int64())])
    list_type = pa.list_(struct_type)

    offsets = [0]
    donors: list[int] = []
    acceptors: list[int] = []
    validity: list[bool] = []

    for jlist in col:
        if jlist is None:
            validity.append(False)
            offsets.append(offsets[-1])
        else:
            validity.append(True)
            for j in jlist:
                donors.append(j["donor"])
                acceptors.append(j["acceptor"])
            offsets.append(offsets[-1] + len(jlist))

    struct_arr = pa.StructArray.from_arrays(
        [pa.array(donors, type=pa.int64()), pa.array(acceptors, type=pa.int64())],
        fields=[pa.field("donor", pa.int64()), pa.field("acceptor", pa.int64())],
    )
    # mask=True means null; invert validity
    null_mask = pa.array([not v for v in validity], type=pa.bool_())
    return pa.ListArray.from_arrays(
        pa.array(offsets, type=pa.int32()),
        struct_arr,
        mask=null_mask,
    )


def _write_empty_shard(out_path: Path) -> None:
    writer = pq.ParquetWriter(str(out_path), L0B_SCHEMA, compression="zstd")
    empty = pa.record_batch([pa.array([], type=f.type) for f in L0B_SCHEMA], schema=L0B_SCHEMA)
    writer.write_batch(empty)
    writer.close()


# ---------------------------------------------------------------------------
# samtools merge/sort → per-chrom BAMs
# ---------------------------------------------------------------------------


def _get_chroms(bam_paths: list[Path]) -> list[str]:
    """Return chromosomes present in *bam_paths[0]*'s header (all BAMs share the same header)."""
    with pysam.AlignmentFile(str(bam_paths[0]), "rb") as bam:
        return [sq["SN"] for sq in bam.header.to_dict().get("SQ", [])]


def _ensure_indexed(bam_paths: list[Path], threads: int = 4) -> None:
    """Ensure every input BAM has a coordinate index (.bai/.csi).

    ``samtools merge -R <chrom>`` does region-based random retrieval and
    requires each input to be indexed.  Indexing is idempotent and the
    indices are reusable across runs, so we create any that are missing.
    Inputs must be coordinate-sorted (they are, from the matrix pipeline);
    ``samtools index`` will fail loudly if one is not.
    """
    missing = [
        p
        for p in bam_paths
        if not (
            p.with_suffix(p.suffix + ".bai").exists() or p.with_suffix(p.suffix + ".csi").exists()
        )
    ]
    if not missing:
        return
    log.info("Indexing %d/%d partition BAMs missing an index", len(missing), len(bam_paths))

    def _index_one(p: Path) -> None:
        proc = subprocess.run(
            ["samtools", "index", f"-@{threads}", str(p)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"samtools index failed for {p} (is it coordinate-sorted?): "
                f"{proc.stderr.decode()}"
            )

    # Index in parallel, but bound concurrency to the CPU budget so we don't
    # spawn one samtools per BAM (each index uses `threads` threads).
    cpu = os.cpu_count() or 4
    max_parallel = max(1, min(len(missing), cpu // max(1, threads)))
    with ThreadPoolExecutor(max_workers=max_parallel) as pool:
        list(pool.map(_index_one, missing))


def _merge_sort_chrom(
    bam_paths: list[Path],
    chrom: str,
    work_dir: Path,
    threads: int = 2,
) -> Path:
    """Run ``samtools merge -R <chrom> | samtools sort`` → per-chrom BAM.

    Returns path to the sorted per-chrom BAM.
    """
    out_bam = work_dir / f"{chrom}.bam"
    if out_bam.exists() and out_bam.stat().st_size > 0:
        return out_bam  # checkpoint hit

    # samtools merge -R restricts to one reference sequence
    inputs = [str(p) for p in bam_paths]
    merge_cmd = [
        "samtools",
        "merge",
        "-f",  # overwrite output
        "-R",
        chrom,  # region filter
        f"--threads={threads}",
        "-",  # output to stdout
    ] + inputs

    sort_cmd = [
        "samtools",
        "sort",
        f"--threads={threads}",
        "-m",
        "768M",  # cap per-thread RAM; spill beyond it
        "-T",
        str(work_dir / f".sort_{chrom}"),  # spill to work_dir (disk), not /tmp (often tmpfs/RAM)
        "-o",
        str(out_bam),
        "-",
    ]

    log.debug("merge+sort %s", chrom)
    merge_proc = subprocess.Popen(merge_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sort_proc = subprocess.Popen(
        sort_cmd,
        stdin=merge_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if merge_proc.stdout:
        merge_proc.stdout.close()

    _, sort_err = sort_proc.communicate()
    merge_proc.wait()

    if merge_proc.returncode != 0:
        _, merge_err = merge_proc.communicate()
        raise RuntimeError(f"samtools merge failed for {chrom}: {merge_err.decode()}")
    if sort_proc.returncode != 0:
        raise RuntimeError(f"samtools sort failed for {chrom}: {sort_err.decode()}")

    return out_bam


# ---------------------------------------------------------------------------
# Per-chrom worker (runs in subprocess)
# ---------------------------------------------------------------------------


def _chrom_worker(
    bam_paths: list[str],
    chrom: str,
    work_dir: str,
    out_dir: str,
    samtools_threads: int,
) -> tuple[str, int]:
    """Top-level function for ProcessPoolExecutor: build one chrom's L0b shard.

    Returns (chrom, n_rows).
    """
    import logging as _logging

    _logging.basicConfig(level=logging.INFO)
    _log = _logging.getLogger(__name__)

    shard_path = Path(out_dir) / f"{chrom}.parquet"
    if shard_path.exists():
        _log.info("checkpoint hit: %s", chrom)
        pf = pq.read_table(str(shard_path))
        return chrom, len(pf)

    _log.info("building %s ...", chrom)
    chrom_bam = _merge_sort_chrom(
        [Path(p) for p in bam_paths],
        chrom,
        Path(work_dir),
        threads=samtools_threads,
    )
    n = _scan_chrom_bam_to_parquet(chrom_bam, chrom, shard_path)
    _log.info("done %s: %d rows", chrom, n)
    return chrom, n


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def build_l0b(
    partition_dir: Path,
    out_dir: Path,
    version_key: VersionKey,
    *,
    chroms: list[str] | None = None,
    workers: int = 4,
    samtools_threads: int = 2,
    work_dir: Path | None = None,
    bam_glob: str = "unique_reads.*.bam",
) -> Path:
    """Build the L0b Parquet store from all prefix-partition BAMs.

    Parameters
    ----------
    partition_dir:
        Root of the ``global_partitioned`` directory (contains per-prefix subdirs).
    out_dir:
        Where to write L0b shards.  Created if absent.
    version_key:
        The 0.4 version key; written to ``out_dir/_meta.json``.
    chroms:
        Subset of chromosomes to build (default: all in BAM header).
    workers:
        Number of parallel chrom workers.
    samtools_threads:
        Threads per ``samtools`` call within each worker.
    work_dir:
        Scratch directory for per-chrom sorted BAMs.  Defaults to a tmpdir
        inside ``out_dir``.
    bam_glob:
        Glob pattern relative to each partition subdir to find the BAM.

    Returns
    -------
    Path to the L0b output directory.
    """
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Discover BAMs — exactly one per partition subdir.
    bam_paths: list[Path] = []
    ambiguous: list[tuple[str, list[str]]] = []
    for subdir in sorted(p for p in partition_dir.iterdir() if p.is_dir()):
        matches = sorted(subdir.glob(bam_glob))
        if len(matches) > 1:
            ambiguous.append((subdir.name, [m.name for m in matches]))
        bam_paths.extend(matches)
    if ambiguous:
        sample = ambiguous[0]
        raise ValueError(
            f"bam_glob {bam_glob!r} matches multiple BAMs in {len(ambiguous)} "
            f"partition dir(s) — this would double-count reads. "
            f"e.g. {sample[0]}/: {sample[1]}. "
            f"Pass a more specific --bam-glob (e.g. 'unique_reads.*.filtered.bam')."
        )
    if not bam_paths:
        raise FileNotFoundError(f"No BAMs matching {bam_glob!r} under {partition_dir}")
    log.info("Found %d partition BAMs (one per partition)", len(bam_paths))

    # samtools merge -R (per-chrom region fetch) needs every input indexed.
    _ensure_indexed(bam_paths, threads=samtools_threads)

    # Chromosome list
    all_chroms = _get_chroms(bam_paths)
    target_chroms: list[str] = chroms if chroms is not None else all_chroms
    log.info("Building L0b for %d chromosomes", len(target_chroms))

    # Work dir for per-chrom sorted BAMs
    _tmp_ctx = None
    if work_dir is None:
        _tmp_ctx = tempfile.TemporaryDirectory(dir=out_dir, prefix="l0b_work_")
        work_dir = Path(_tmp_ctx.name)
    else:
        work_dir = work_dir.resolve()
        work_dir.mkdir(parents=True, exist_ok=True)

    # Write meta before we start (so a crash mid-build is detectable)
    write_meta(out_dir, version_key, {"status": "building", "n_chroms": len(target_chroms)})

    # Launch parallel workers
    bam_strs = [str(p) for p in bam_paths]
    out_str = str(out_dir)
    work_str = str(work_dir)

    total_rows = 0
    errors: list[str] = []

    # Recycle each worker after one chromosome so per-process RSS (pyarrow /
    # malloc arenas that don't return memory to the OS) cannot accumulate across
    # the ~hundreds of chroms a long-lived worker would otherwise handle.
    # max_tasks_per_child is Python 3.11+; degrade gracefully on older runtimes.
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
                bam_strs,
                chrom,
                work_str,
                out_str,
                samtools_threads,
            ): chrom
            for chrom in target_chroms
        }
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                _, n = fut.result()
                total_rows += n
                log.info("✓ %s  (%d rows)", chrom, n)
            except Exception as exc:
                log.error("✗ %s  %s", chrom, exc)
                errors.append(f"{chrom}: {exc}")

    if _tmp_ctx is not None:
        _tmp_ctx.cleanup()

    # Update meta to reflect final status
    write_meta(
        out_dir,
        version_key,
        {
            "status": "complete" if not errors else "partial",
            "n_chroms_built": len(target_chroms) - len(errors),
            "n_chroms_failed": len(errors),
            "total_rows": total_rows,
            "errors": errors,
        },
    )

    if errors:
        raise RuntimeError(
            f"L0b build finished with {len(errors)} failed chromosomes:\n" + "\n".join(errors)
        )

    log.info("L0b build complete: %d total rows → %s", total_rows, out_dir)
    return out_dir

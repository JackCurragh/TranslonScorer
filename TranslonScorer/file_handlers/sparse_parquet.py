"""Sparse Parquet unique-read matrix support.

The sparse matrix stores read counts by sample/run; the companion BAM stores
the genomic location for each unique read. This module joins those two pieces
without materializing the full matrix.
"""

from __future__ import annotations

import json
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import polars as pl
import pysam

READ_NAME_RE = re.compile(r"(?:^|[^\w])read_(\d+)(?:$|[^\d])")


@dataclass(frozen=True)
class Region:
    """0-based, half-open genomic interval."""

    chrom: str
    start: int
    end: int


def _load_manifest(path: str | Path) -> tuple[Path, dict]:
    path = Path(path)
    if path.is_dir():
        # Partition directories use a prefixed name: global.<PREFIX>_matrix_manifest.json
        candidates = sorted(path.glob("global.*_matrix_manifest.json"))
        if not candidates:
            # Fallback: unprefixed name used by older builds
            candidates = sorted(path.glob("global_matrix_manifest.json"))
        if not candidates:
            raise FileNotFoundError(f"No matrix manifest found in {path}")
        manifest_path = candidates[0]
    else:
        manifest_path = path
    with open(manifest_path, "r") as handle:
        manifest = json.load(handle)
    if manifest.get("matrix_format") != "sparse-parquet":
        raise ValueError(f"{manifest_path} is not a sparse-parquet manifest")
    return manifest_path, manifest


def _generation_paths(manifest_path: Path, manifest: dict) -> list[Path]:
    generations = manifest.get("generations") or [
        {
            "generation": manifest.get("generation", 1),
            "path": manifest.get("path", "global_counts/generation=000001"),
        }
    ]
    return [manifest_path.parent / item["path"] for item in generations]


def _parse_read_id(query_name: Optional[str]) -> Optional[int]:
    if query_name is None:
        return None
    match = READ_NAME_RE.search(query_name)
    return int(match.group(1)) if match else None


def _load_tombstone_filters(
    manifest_path: Path, manifest: dict
) -> tuple[set[int], set[str], set[str]]:
    raw = manifest.get("tombstones_path", "global_tombstones.parquet")
    tombstones_path = manifest_path.parent / raw
    if not tombstones_path.exists():
        # Partition builds use a prefixed tombstones file; discover it
        candidates = sorted(manifest_path.parent.glob("global.*_tombstones.parquet"))
        tombstones_path = candidates[0] if candidates else tombstones_path
    if not tombstones_path.exists():
        return set(), set(), set()
    frame = pl.read_parquet(tombstones_path)
    if frame.is_empty():
        return set(), set(), set()
    frame = frame.filter(pl.col("active") == True)
    read_ids: set[int] = set()
    sample_names: set[str] = set()
    study_ids: set[str] = set()
    for row in frame.iter_rows(named=True):
        if row.get("read_id") is not None:
            read_ids.add(int(row["read_id"]))
        if row.get("sample_id") is not None:
            sample_names.add(str(row["sample_id"]))
        if row.get("study_id") is not None:
            study_ids.add(str(row["study_id"]))
    return read_ids, sample_names, study_ids


def _load_samples(manifest_path: Path, manifest: dict) -> pl.LazyFrame:
    samples_path = manifest_path.parent / manifest.get("lookup_tables", {}).get(
        "samples", "global_samples.parquet"
    )
    if not samples_path.exists():
        raise FileNotFoundError(f"Sample lookup table does not exist: {samples_path}")
    return pl.scan_parquet(str(samples_path)).select(
        "sample_id",
        "sample_name",
        "study_id",
        "study_id_int",
    )


def _count_bucket_paths(generation_roots: Iterable[Path], read_bucket: int) -> list[Path]:
    bucket_name = f"read_bucket={read_bucket:06d}"
    paths: list[Path] = []
    for root in generation_roots:
        paths.extend(sorted((root / bucket_name).glob("*.parquet")))
    return paths


def _scan_event_bucket(path: Path) -> pl.LazyFrame:
    schema = {
        "read_bucket": pl.UInt32,
        "read_id": pl.UInt64,
        "chr": pl.Utf8,
        "start": pl.Int64,
        "stop": pl.Int64,
        "strand": pl.Utf8,
        "length": pl.Int64,
        "event_count": pl.UInt32,
    }
    # `dtypes=` was the pre-1.0 spelling of `schema_overrides=`; polars is
    # pinned ==1.36.1 so the old fallback was unreachable.
    return pl.scan_csv(path, separator="\t", schema_overrides=schema)


def _regions_from_exons(exon_df: pl.DataFrame) -> list[Region]:
    regions: list[Region] = []
    for row in exon_df.select(["chr", "start", "stop"]).iter_rows(named=True):
        starts = row["start"] if isinstance(row["start"], list) else [row["start"]]
        stops = row["stop"] if isinstance(row["stop"], list) else [row["stop"]]
        for start, stop in zip(starts, stops):
            regions.append(Region(str(row["chr"]), int(start), int(stop)))
    return regions


def _spool_bam_events(
    bam_path: str | Path,
    read_bucket_size: int,
    out_dir: Path,
    *,
    exclude_read_ids: set[int],
    regions: Optional[list[Region]],
) -> dict[int, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    handles: dict = {}
    paths: dict[int, Path] = {}
    seen_region_events = set()
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            if regions:
                if not bam.has_index():
                    raise ValueError(
                        f"Regional sparse matrix profiling requires an index for {bam_path}. "
                        "Create one with samtools index."
                    )
                # Normalise chromosome names: if BAM uses chr-prefix and regions don't (or vice versa)
                bam_chroms = set(bam.references)

                def _normalise_chrom(chrom: str) -> Optional[str]:
                    if chrom in bam_chroms:
                        return chrom
                    alt = f"chr{chrom}" if not chrom.startswith("chr") else chrom[3:]
                    return alt if alt in bam_chroms else None

                records = (
                    record
                    for region in regions
                    for _chrom in [_normalise_chrom(region.chrom)]
                    if _chrom is not None
                    for record in bam.fetch(_chrom, region.start, region.end)
                )
            else:
                records = bam.fetch(until_eof=True)  # type: ignore[assignment]

            for record in records:
                if record.is_unmapped or record.reference_name is None:
                    continue
                read_id = _parse_read_id(record.query_name)
                if read_id is None or read_id in exclude_read_ids:
                    continue
                strand = "-" if record.is_reverse else "+"
                start = int(record.reference_start)
                stop = int(record.reference_end or record.reference_start)
                length = int(record.query_length or max(0, stop - start))
                event_key = (
                    read_id,
                    record.reference_name,
                    start,
                    stop,
                    strand,
                    record.cigarstring,
                )
                if regions and event_key in seen_region_events:
                    continue
                seen_region_events.add(event_key)
                read_bucket = read_id // read_bucket_size
                handle = handles.get(read_bucket)
                if handle is None:
                    path = out_dir / f"events.read_bucket={read_bucket:06d}.tsv"
                    handle = open(path, "w")
                    handle.write(
                        "read_bucket\tread_id\tchr\tstart\tstop\tstrand\tlength\tevent_count\n"
                    )
                    handles[read_bucket] = handle
                    paths[read_bucket] = path
                handle.write(
                    f"{read_bucket}\t{read_id}\t{record.reference_name}\t{start}\t{stop}\t{strand}\t{length}\t1\n"
                )
    finally:
        for handle in handles.values():
            handle.close()
    return paths


def sparse_matrix_genomic_counts(
    *,
    bam_path: str | Path,
    manifest_path: str | Path,
    exon_df: Optional[pl.DataFrame] = None,
    regions: Optional[list[Region]] = None,
    sample_names: Optional[list[str]] = None,
) -> pl.DataFrame:
    """Return sample-resolved genomic starts from a sparse Parquet matrix.

    Output columns are ``sample_id`` (run/sample name), ``sample_index``,
    ``study_id``, ``study_id_int``, ``chr``, ``start``, ``stop``, ``strand``,
    ``length``, and ``count``.
    """
    manifest_file, manifest = _load_manifest(manifest_path)
    read_bucket_size = int(manifest["read_bucket_size"])
    generation_roots = _generation_paths(manifest_file, manifest)
    samples = _load_samples(manifest_file, manifest)
    read_tombstones, sample_tombstones, study_tombstones = _load_tombstone_filters(
        manifest_file, manifest
    )
    if sample_names:
        samples = samples.filter(pl.col("sample_name").is_in(sample_names))
    if sample_tombstones:
        samples = samples.filter(~pl.col("sample_name").is_in(sorted(sample_tombstones)))
    if study_tombstones:
        samples = samples.filter(~pl.col("study_id").is_in(sorted(study_tombstones)))

    fetch_regions = (
        regions
        if regions is not None
        else (_regions_from_exons(exon_df) if exon_df is not None else None)
    )
    with tempfile.TemporaryDirectory(prefix="translonscorer_sparse_matrix_") as tmp:
        event_paths = _spool_bam_events(
            bam_path,
            read_bucket_size,
            Path(tmp),
            exclude_read_ids=read_tombstones,
            regions=fetch_regions,
        )
        parts: list[pl.DataFrame] = []
        for read_bucket, event_path in sorted(event_paths.items()):
            count_paths = _count_bucket_paths(generation_roots, read_bucket)
            if not count_paths:
                continue
            events = _scan_event_bucket(event_path)
            counts = pl.scan_parquet([str(path) for path in count_paths])
            joined = (
                events.join(counts, on="read_id", how="inner")
                .join(samples, on=["sample_id", "study_id_int"], how="inner")
                .with_columns(
                    (
                        pl.col("event_count").cast(pl.Float64) * pl.col("count").cast(pl.Float64)
                    ).alias("_count")
                )
                .group_by(
                    "sample_name",
                    "sample_id",
                    "study_id",
                    "study_id_int",
                    "chr",
                    "start",
                    "stop",
                    "strand",
                    "length",
                )
                .agg(pl.col("_count").sum().alias("count"))
                .select(
                    pl.col("sample_name").alias("sample_id"),
                    pl.col("sample_id").cast(pl.UInt32).alias("sample_index"),
                    "study_id",
                    pl.col("study_id_int").cast(pl.UInt32),
                    "chr",
                    pl.col("start").cast(pl.Int64),
                    pl.col("stop").cast(pl.Int64),
                    "strand",
                    pl.col("length").cast(pl.Int64),
                    "count",
                )
            ).collect()
            if not joined.is_empty():
                parts.append(joined)

    if not parts:
        return pl.DataFrame(
            schema={
                "sample_id": pl.Utf8,
                "sample_index": pl.UInt32,
                "study_id": pl.Utf8,
                "study_id_int": pl.UInt32,
                "chr": pl.Utf8,
                "start": pl.Int64,
                "stop": pl.Int64,
                "strand": pl.Utf8,
                "length": pl.Int64,
                "count": pl.Float64,
            }
        )
    return (
        pl.concat(parts)
        .group_by(
            "sample_id",
            "sample_index",
            "study_id",
            "study_id_int",
            "chr",
            "start",
            "stop",
            "strand",
            "length",
        )
        .agg(pl.col("count").sum())
        .sort(["sample_id", "chr", "start", "strand", "length"])
    )

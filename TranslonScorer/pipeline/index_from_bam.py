from __future__ import annotations

"""
Build a Zarr read alignment index (read_id -> genomic alignment) from a BAM file.

Supports three mapping strategies to align Zarr rows to BAM reads:
  - qname-based (preferred) using metadata/parquet or FASTA headers
  - sequence-based using stable hashes
  - row-order-based (brittle) if row order matches BAM iteration

Outputs a Parquet file with columns: read_id, chr, start, stop, strand, length
Optionally can produce a splits index (per-read junctions) in a follow-up.
"""

from typing import Dict, Optional, Tuple, Iterable
import os
import io
import hashlib

import polars as pl


def _open_counts(zarr_root: str):
    # Lazy import to avoid hard dependency at CLI import time
    from ..file_handlers.zarr import _open_counts  # type: ignore
    return _open_counts(zarr_root)


def _counts_is_samples_first(arr) -> bool:
    from ..file_handlers.zarr import _counts_is_samples_first  # type: ignore
    return _counts_is_samples_first(arr)


def _detect_row_count(zarr_root: str) -> int:
    arr = _open_counts(zarr_root)
    # rows correspond to reads dimension
    return arr.shape[1] if _counts_is_samples_first(arr) else arr.shape[0]


def _iter_fasta(path: str) -> Iterable[Tuple[str, str]]:
    """Yield (header, sequence) from a FASTA file without external deps."""
    header = None
    seq_chunks = []
    with open(path, 'r') as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line.strip()[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(seq_chunks)


def _hash_seq(seq: str, alg: str = 'sha1') -> str:
    if alg == 'md5':
        return hashlib.md5(seq.encode('utf-8')).hexdigest()
    if alg == 'xxh64':
        try:
            import xxhash  # type: ignore
            return xxhash.xxh64(seq).hexdigest()
        except Exception:
            # Fallback to sha1 if xxhash missing
            return hashlib.sha1(seq.encode('utf-8')).hexdigest()
    return hashlib.sha1(seq.encode('utf-8')).hexdigest()


def _read_zarr_metadata(
    zarr_root: str,
    zarr_metadata: Optional[str],
    zarr_reads_fasta: Optional[str],
    zarr_read_key: str,
    hash_alg: str,
) -> Tuple[str, Dict[str, int]]:
    """
    Build a mapping from a join key -> row_id, and return which key was used.

    Returns (mode, mapping) where mode in {'qname','sequence','row_id'}.
    """
    n_rows = _detect_row_count(zarr_root)

    # Try metadata parquet for qname or sequence with explicit row_id
    if zarr_metadata and os.path.isfile(zarr_metadata):
        try:
            meta = pl.read_parquet(zarr_metadata)
            cols = set(meta.columns)
            # Identify a row id column
            row_col = 'row_id'
            for candidate in ['row_id', 'read_id', 'index', 'row']:
                if candidate in cols:
                    row_col = candidate
                    break
            if 'qname' in cols and (zarr_read_key in ('auto','qname')):
                mapping = {str(q): int(i) for i, q in meta.select([row_col, 'qname']).iter_rows()}
                return 'qname', mapping
            if 'sequence' in cols and (zarr_read_key in ('auto','sequence')):
                mapping = {str(_hash_seq(seq, hash_alg)): int(i) for i, seq in meta.select([row_col, 'sequence']).iter_rows()}
                return 'sequence', mapping
        except Exception:
            pass

    # Try FASTA: headers for qname, or sequences for hashing
    if zarr_reads_fasta and os.path.isfile(zarr_reads_fasta):
        mapping_q: Dict[str,int] = {}
        mapping_h: Dict[str,int] = {}
        for idx, (hdr, seq) in enumerate(_iter_fasta(zarr_reads_fasta)):
            if zarr_read_key in ('auto','qname'):
                mapping_q[str(hdr)] = idx
            if zarr_read_key in ('auto','sequence'):
                mapping_h[_hash_seq(seq, hash_alg)] = idx
        if zarr_read_key in ('auto','qname') and mapping_q:
            return 'qname', mapping_q
        if zarr_read_key in ('auto','sequence') and mapping_h:
            return 'sequence', mapping_h

    # Fallback: row-order mapping only if explicitly requested or auto with no better info
    if zarr_read_key in ('auto','row_id'):
        # identity mapping by position (0..n_rows-1)
        return 'row_id', {str(i): i for i in range(n_rows)}

    raise ValueError(
        "Unable to infer Zarr read mapping. Provide --zarr-metadata-parquet with qname/sequence, "
        "or --zarr-reads-fasta, or set --zarr-read-key=row_id if row order matches BAM."
    )


def build_read_index_from_bam(
    *,
    bam_path: str,
    zarr_root: str,
    out_index: str,
    zarr_metadata: Optional[str] = None,
    zarr_reads_fasta: Optional[str] = None,
    zarr_read_key: str = 'auto',
    bam_key: str = 'auto',
    hash_alg: str = 'sha1',
    chunk_size: int = 1_000_000,
) -> str:
    """Build read_index.parquet aligned to Zarr rows.

    Returns the path to the written Parquet file.
    """
    import pysam
    from pyarrow import Table as _ArrowTable
    from pyarrow import parquet as pq

    # Determine Zarr mapping key -> row_id
    mode, zmap = _read_zarr_metadata(
        zarr_root=zarr_root,
        zarr_metadata=zarr_metadata,
        zarr_reads_fasta=zarr_reads_fasta,
        zarr_read_key=zarr_read_key,
        hash_alg=hash_alg,
    )

    # Decide BAM key usage
    if bam_key == 'auto':
        bam_key = 'qname' if mode == 'qname' else ('sequence' if mode == 'sequence' else 'row_id')

    # Prepare writer
    os.makedirs(os.path.dirname(out_index) or '.', exist_ok=True)
    writer: Optional[pq.ParquetWriter] = None

    def _write_chunk(df: pl.DataFrame):
        nonlocal writer
        if df.is_empty():
            return
        tbl = df.to_arrow()
        if writer is None:
            writer = pq.ParquetWriter(out_index, schema=tbl.schema)
        writer.write_table(tbl)  # type: ignore

    # Row-order fast path (brittle; rely on Zarr row count)
    if mode == 'row_id' and bam_key == 'row_id':
        n_rows = _detect_row_count(zarr_root)
        curr = 0
        rows = []
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            for aln in bam.fetch(until_eof=True):
                if aln.is_unmapped:
                    continue
                if curr >= n_rows:
                    break
                strand = '-' if aln.is_reverse else '+'
                rows.append({
                    'read_id': curr,
                    'chr': bam.get_reference_name(aln.reference_id),
                    'start': int(aln.reference_start),
                    'stop': int(aln.reference_end),
                    'strand': strand,
                    'length': int(aln.query_length or 0),
                })
                if len(rows) >= chunk_size:
                    _write_chunk(pl.from_dicts(rows))
                    rows = []
                curr += 1
        if rows:
            _write_chunk(pl.from_dicts(rows))
        if writer:
            writer.close()
        return out_index

    # qname / sequence mapping
    rows = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            if bam_key == 'qname':
                key = str(aln.query_name)
            elif bam_key == 'sequence':
                key = _hash_seq(aln.query_sequence or '', hash_alg)
            else:
                raise ValueError("Unsupported bam_key for mapping: %s" % bam_key)
            if key not in zmap:
                continue  # unmatched read; acceptable if counts are unique reads
            read_id = int(zmap[key])
            strand = '-' if aln.is_reverse else '+'
            rows.append({
                'read_id': read_id,
                'chr': bam.get_reference_name(aln.reference_id),
                'start': int(aln.reference_start),
                'stop': int(aln.reference_end),
                'strand': strand,
                'length': int(aln.query_length or 0),
            })
            if len(rows) >= chunk_size:
                _write_chunk(pl.from_dicts(rows))
                rows = []
    if rows:
        _write_chunk(pl.from_dicts(rows))
    if writer:
        writer.close()
    return out_index


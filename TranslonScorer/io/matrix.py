"""Sparse-Parquet matrix I/O helpers: manifest, count-parquet, samples, BAM discovery.

Pure I/O adapters — open → read → close, return basic Python / Polars objects.
No computation, no frame assignment.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import List, Optional, Tuple

import polars as pl

# ---------------------------------------------------------------------------
# Manifest + lookup helpers
# ---------------------------------------------------------------------------


def _manifest(partition_dir: "str | Path") -> Tuple[Path, dict]:
    d = Path(partition_dir)
    candidates = sorted(d.glob("global.*_matrix_manifest.json"))
    if not candidates:
        candidates = sorted(d.glob("global_matrix_manifest.json"))
    if not candidates:
        raise FileNotFoundError(f"No matrix manifest in {d}")
    mp = candidates[0]
    return mp, json.loads(mp.read_text())


def _count_parquets(mp: Path, manifest: dict) -> List[str]:
    generations = manifest.get("generations") or [
        {"path": manifest.get("path", "global.AAAA_counts/generation=000001")}
    ]
    paths = []
    for gen in generations:
        gen_dir = mp.parent / gen["path"]
        paths.extend(str(p) for p in sorted(gen_dir.rglob("*.parquet")))
    return paths


def _samples_df(mp: Path, manifest: dict) -> pl.DataFrame:
    key = manifest.get("lookup_tables", {}).get("samples", "")
    sp = mp.parent / key if key else None
    if sp is None or not sp.exists():
        candidates = sorted(mp.parent.glob("global.*_samples.parquet"))
        sp = candidates[0] if candidates else None
    if sp is None:
        raise FileNotFoundError(f"samples.parquet not found in {mp.parent}")
    return pl.read_parquet(str(sp))


def _reads_parquet_path(mp: Path, manifest: dict) -> str:
    key = manifest.get("global_reads", "")
    rp = mp.parent / key if key else None
    if rp is None or not rp.exists():
        candidates = sorted(mp.parent.glob("global.*_reads.parquet"))
        rp = candidates[0] if candidates else None
    if rp is None:
        raise FileNotFoundError(f"reads.parquet not found in {mp.parent}")
    return str(rp)


def _discover_bam(partition_dir: "str | Path") -> Optional[Path]:
    """Find the unique_reads BAM for a partition directory."""
    d = Path(partition_dir)
    candidates = sorted(d.glob("unique_reads.*.bam"))
    return candidates[0] if candidates else None


# ---------------------------------------------------------------------------
# Read-ID parsing
# ---------------------------------------------------------------------------

_READ_ID_RE = re.compile(r"read_(\d+)")


def _parse_read_id(qname: str) -> Optional[int]:
    m = _READ_ID_RE.search(qname)
    return int(m.group(1)) if m else None

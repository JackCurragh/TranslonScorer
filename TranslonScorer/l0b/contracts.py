"""Phase 0 contracts: version key and L0b Parquet schema.

These two artefacts are the foundation everything else is built on.  They are
defined *once* here and imported everywhere so there is a single source of truth.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pyarrow as pa

# ---------------------------------------------------------------------------
# 0.4  Version key
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VersionKey:
    """Immutable identity of an L0b (and downstream) artefact.

    Two artefacts with the same key are bit-for-bit interchangeable; any field
    change forces a rebuild.  ``multimap_policy`` is stored here (not only in L1)
    because the policy annotation *could* change what we record in L0b in future
    (e.g. EM weights); storing it now avoids ambiguity.
    """

    genome_id: str  # e.g. "GRCh38.p14"
    aligner_cfg_hash: str  # sha256[:12] of the aligner config / index params
    junction_set_id: str  # e.g. "GENCODE_v44"
    l0a_version: str  # e.g. "2026-06"
    multimap_policy: str  # e.g. "unique_only"

    def as_dict(self) -> dict[str, str]:
        return asdict(self)

    def fingerprint(self) -> str:
        """Stable 12-char hex digest of the key, usable in directory names."""
        blob = json.dumps(self.as_dict(), sort_keys=True).encode()
        return hashlib.sha256(blob).hexdigest()[:12]

    def dir_fragment(self) -> str:
        """Human-readable path component for directory naming."""
        return (
            f"genome={self.genome_id}"
            f"/jset={self.junction_set_id}"
            f"/l0a={self.l0a_version}"
            f"/policy={self.multimap_policy}"
        )


def write_meta(root: Path, key: VersionKey, extra: dict[str, Any] | None = None) -> None:
    """Write ``_meta.json`` to *root*."""
    root.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {"version_key": key.as_dict(), "fingerprint": key.fingerprint()}
    if extra:
        payload.update(extra)
    (root / "_meta.json").write_text(json.dumps(payload, indent=2))


def read_meta(root: Path) -> dict[str, Any]:
    """Read ``_meta.json`` from *root*; raise if missing."""
    meta_path = root / "_meta.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"No _meta.json in {root}")
    return json.loads(meta_path.read_text())


# ---------------------------------------------------------------------------
# 1.1  L0b Parquet schema
# ---------------------------------------------------------------------------

# Junctions are stored as a list of (donor, acceptor) int64 pairs.
# We encode them as a list<struct<donor: int64, acceptor: int64>>.
_JUNCTION_TYPE = pa.list_(pa.struct([("donor", pa.int64()), ("acceptor", pa.int64())]))

L0B_SCHEMA = pa.schema(
    [
        # Read identity
        pa.field("read_id", pa.int64(), nullable=False),
        # Genomic locus
        pa.field("chrom", pa.large_utf8(), nullable=False),
        pa.field("pos5", pa.int64(), nullable=False),  # 0-based 5′ end
        pa.field("end", pa.int64(), nullable=False),  # 0-based exclusive
        pa.field("strand", pa.int8(), nullable=False),  # +1 / -1
        pa.field("length", pa.int16(), nullable=False),  # read length (bp)
        # Alignment metadata
        pa.field("cigar", pa.large_utf8(), nullable=True),  # None for simple M-only
        pa.field("mapq", pa.uint8(), nullable=False),
        pa.field("nh", pa.int32(), nullable=False),  # from NH tag
        pa.field("is_secondary", pa.bool_(), nullable=False),
        pa.field("aln_score", pa.int32(), nullable=True),  # AS tag
        pa.field("mismatches", pa.int32(), nullable=True),  # NM tag
        # Splice junctions crossed (null / empty for unspliced reads)
        pa.field("junctions_crossed", _JUNCTION_TYPE, nullable=True),
        # Placement weight — null at build time; filled by a policy pass
        pa.field("weight", pa.float32(), nullable=True),
    ]
)

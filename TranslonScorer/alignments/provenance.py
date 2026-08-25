"""Artefact identity: the version key that makes a build reproducible.

Defined *once* here and imported everywhere so there is a single source of
truth for what identifies a built artefact.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Version key
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VersionKey:
    """Immutable identity of an alignment table (and everything built on it).

    Two artefacts with the same key are bit-for-bit interchangeable; any field
    change forces a rebuild.  ``multimap_policy`` is stored here (not only on
    the count tables) because the policy *could* change what we record in the
    alignment table in future (e.g. EM weights); storing it now avoids
    ambiguity.
    """

    genome_id: str  # e.g. "GRCh38.p14"
    aligner_cfg_hash: str  # sha256[:12] of the aligner config / index params
    junction_set_id: str  # e.g. "GENCODE_v44"
    source_data_version: str  # upstream read-data release, e.g. "2026-06"
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
            f"/source={self.source_data_version}"
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

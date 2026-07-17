"""L0b: multimapping-aware genomic alignment loci table.

One row per (unique_read × alignment). Built once from the partitioned prefix
BAMs; never rebuilt unless the junction set or aligner config changes.
"""

from .builder import build_l0b
from .contracts import L0B_SCHEMA, VersionKey, read_meta, write_meta

__all__ = ["L0B_SCHEMA", "VersionKey", "write_meta", "read_meta", "build_l0b"]

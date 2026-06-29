"""L0b: multimapping-aware genomic alignment loci table.

One row per (unique_read × alignment). Built once from the partitioned prefix
BAMs; never rebuilt unless the junction set or aligner config changes.
"""
from .contracts import L0B_SCHEMA, VersionKey, write_meta, read_meta
from .builder import build_l0b

__all__ = ["L0B_SCHEMA", "VersionKey", "write_meta", "read_meta", "build_l0b"]

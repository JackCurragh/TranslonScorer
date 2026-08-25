"""Multimapping-aware genomic alignment table.

One row per (unique_read × alignment). Built once from the partitioned prefix
BAMs; never rebuilt unless the junction set or aligner config changes.
"""

from .builder import build_alignments
from .provenance import VersionKey, read_meta, write_meta
from .schema import ALIGNMENT_SCHEMA

__all__ = ["ALIGNMENT_SCHEMA", "VersionKey", "write_meta", "read_meta", "build_alignments"]

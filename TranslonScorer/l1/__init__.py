"""L1: raw 5′ positional cache — unique mappers only, pre-offset.

Two streams per chromosome:
  {chrom}_pos.parquet   — (pos5, strand, length, sample_id) → count
  {chrom}_junc.parquet  — (donor, acceptor, strand, length, sample_id) → count
"""

from .builder import build_l1
from .schema import L1_JUNC_SCHEMA, L1_POS_SCHEMA

__all__ = ["L1_POS_SCHEMA", "L1_JUNC_SCHEMA", "build_l1"]

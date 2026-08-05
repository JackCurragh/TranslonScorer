"""Raw 5′-end count tables — unique mappers only, pre-offset.

Two streams per chromosome:
  {chrom}_pos.parquet   — (pos5, strand, length, sample_id) → count
  {chrom}_junc.parquet  — (donor, acceptor, strand, length, sample_id) → count
"""

from .builder import build_counts
from .schema import FIVE_PRIME_COUNT_SCHEMA, JUNCTION_COUNT_SCHEMA

__all__ = ["FIVE_PRIME_COUNT_SCHEMA", "JUNCTION_COUNT_SCHEMA", "build_counts"]

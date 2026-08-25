"""Parquet schema for the alignment table: one row per aligned read record.

Note this is per *alignment*, not per read — ``is_secondary`` and ``nh`` mean a
multimapping read contributes several rows.
"""

from __future__ import annotations

import pyarrow as pa

# Junctions are stored as a list of (donor, acceptor) int64 pairs.
# We encode them as a list<struct<donor: int64, acceptor: int64>>.
_JUNCTION_TYPE = pa.list_(pa.struct([("donor", pa.int64()), ("acceptor", pa.int64())]))

ALIGNMENT_SCHEMA = pa.schema(
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

"""L1 Parquet schemas — positional stream and junction stream."""

import pyarrow as pa

L1_POS_SCHEMA = pa.schema(
    [
        pa.field("pos5", pa.int64(), nullable=False),
        pa.field("strand", pa.int8(), nullable=False),  # +1 / -1
        pa.field("length", pa.int16(), nullable=False),
        pa.field("sample_id", pa.uint32(), nullable=False),
        pa.field("count", pa.uint32(), nullable=False),
    ]
)

L1_JUNC_SCHEMA = pa.schema(
    [
        pa.field("donor", pa.int64(), nullable=False),
        pa.field("acceptor", pa.int64(), nullable=False),
        pa.field("strand", pa.int8(), nullable=False),
        pa.field("length", pa.int16(), nullable=False),
        pa.field("sample_id", pa.uint32(), nullable=False),
        pa.field("count", pa.uint32(), nullable=False),
    ]
)

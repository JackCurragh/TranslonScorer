from pathlib import Path

import polars as pl

from TranslonScorer.io.store import read_event_overlap


def test_read_event_overlap_normalises_mixed_integer_widths(tmp_path: Path) -> None:
    root = tmp_path / "event_overlap"
    root.mkdir()
    pl.DataFrame(
        {"event_id": [1], "other_event_id": [2], "overlap_start": [10], "overlap_end": [20]}
    ).with_columns(
        pl.col("event_id").cast(pl.Int128), pl.col("other_event_id").cast(pl.Int128)
    ).write_parquet(
        root / "chr1.parquet"
    )
    pl.DataFrame(
        {"event_id": [3], "other_event_id": [4], "overlap_start": [30], "overlap_end": [40]}
    ).with_columns(
        pl.col("event_id").cast(pl.UInt64), pl.col("other_event_id").cast(pl.UInt64)
    ).write_parquet(
        root / "chrM.parquet"
    )

    out = read_event_overlap(str(root))
    assert out.height == 2
    assert out["event_id"].to_list() == [1, 3]

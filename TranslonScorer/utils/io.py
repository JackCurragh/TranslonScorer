from __future__ import annotations

"""
I/O helpers for safe writes.

Utilities here avoid common errors where parent directories do not exist
before writing Parquet/CSV outputs.
"""

import os
from typing import Any


def ensure_dir_for_file(path: str) -> str:
    """Ensure parent directory exists for a file path and return the path."""
    parent = os.path.dirname(path) or "."
    os.makedirs(parent, exist_ok=True)
    return path


def write_parquet_safe(df: Any, path: str) -> str:
    """Create parent dir then write a Polars DataFrame to Parquet."""
    ensure_dir_for_file(path)
    df.write_parquet(path)
    return path


def write_csv_safe(df: Any, path: str, **kwargs) -> str:
    """Create parent dir then write a Polars/Pandas DataFrame to CSV."""
    ensure_dir_for_file(path)
    df.write_csv(path, **kwargs)
    return path


"""Consequentiality policy — pure transform.

Given per-translon score reports and a ConsequentialityPolicy, produce a
consequentiality label per translon.  No file I/O.

Public API
----------
apply_policy  — reports + policy → DataFrame with consequentiality column
"""
from __future__ import annotations

import polars as pl

from TranslonScorer.model import ConsequentialityPolicy


def apply_policy(
    report: pl.DataFrame,
    policy: ConsequentialityPolicy = ConsequentialityPolicy(),
) -> pl.DataFrame:
    """Apply a consequentiality policy to a per-translon report DataFrame.

    Adds a boolean `consequential` column.  Currently a pass-through stub;
    full implementation lands in Phase 4 (T13).
    """
    _ = policy
    return report.with_columns(pl.lit(False).alias("consequential"))

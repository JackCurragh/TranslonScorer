"""Coverage-provider capability convention (duck-typed, no Protocol hierarchy).

The three provider implementations (BamSetProvider, BigwigSetProvider,
MatrixProvider) have genuinely different capabilities. Rather than a set of
runtime-checkable Protocols, capability is expressed by which methods a provider
exposes and honoured by duck typing at the single call site that needs it:

- ``coverage(regions, *, site="A", by_sample=False) -> DataFrame(pos, count[, ...])``
    every provider — the minimal contract (``site`` is ignored by bigwig).
- ``size_factors() -> {sample_id: float}``
    every provider.
- ``junction_support(junctions, *, by_sample=False) -> DataFrame(junction_id, count)``
    BAM/matrix only. BigwigSetProvider defines it to raise NotImplementedError
    (no CIGAR), so callers gate on ``hasattr(provider, "junction_support")`` AND
    catch NotImplementedError (see workflows._junction_support_for_chrom).
- ``mappability_ledger(events) -> DataFrame(MAPPABILITY_LEDGER_SCHEMA)``
    BAM/matrix only.

This module keeps the shared MAPPABILITY_LEDGER_SCHEMA that ledger callers need.
"""

from __future__ import annotations

from typing import Dict

import polars as pl

# ---------------------------------------------------------------------------
# Mappability ledger schema
# ---------------------------------------------------------------------------

MAPPABILITY_LEDGER_SCHEMA: Dict[str, type] = {
    "event_id": pl.UInt64,
    "unique_reads": pl.Float64,
    "multimapper_reads": pl.Float64,
    "mean_loci_per_multiread": pl.Float64,
    "competing_event_ids": pl.Utf8,
}

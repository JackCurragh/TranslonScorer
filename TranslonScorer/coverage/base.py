"""Capability Protocols for coverage providers.

Three provider implementations (BamSetProvider, BigwigSetProvider,
MatrixProvider) have genuinely different capabilities.  A single ABC would
force junction_support to raise NotImplementedError on bigwig and site to be
silently ignored.  Instead, runtime-checkable Protocols express capability in
the type system so a bigwig provider passed where SupportsJunctions is
expected is a *static type error*, not a runtime crash.

Usage
-----
from TranslonScorer.coverage.base import (
    CoverageProvider, SupportsSites, SupportsJunctions, SupportsMappability,
    MAPPABILITY_LEDGER_SCHEMA,
)

isinstance(provider, SupportsJunctions)   # runtime capability check
"""
from __future__ import annotations

from typing import Dict, List, Tuple
from typing_extensions import Protocol, runtime_checkable

import polars as pl

from TranslonScorer.model import Region


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


# ---------------------------------------------------------------------------
# Capability Protocols
# ---------------------------------------------------------------------------

@runtime_checkable
class CoverageProvider(Protocol):
    """Minimal contract satisfied by all coverage sources.

    coverage() returns a tidy DataFrame with at least columns:
        pos   (int)   genomic or transcript position
        count (float) total reads at that position (aggregate across samples)
    and optionally sample_id (str) when by_sample=True.

    size_factors() returns {sample_id: float} depth-normalisation factors
    (median-ratio or library-size; 1.0 for single-sample sources).
    """

    def coverage(
        self,
        regions: List[Region],
        *,
        by_sample: bool = False,
    ) -> pl.DataFrame: ...

    def size_factors(self) -> Dict[str, float]: ...


@runtime_checkable
class SupportsSites(Protocol):
    """Provider that distinguishes P-site (initiation) from A-site (elong/term).

    BAM and matrix providers implement this; bigwig providers do not (they carry
    pre-computed coverage with no read-length information).

    site ∈ {"P", "A"}.  A-site default is P-site offset + 3 nt unless the
    provider's offset table specifies a distinct A-site offset.
    """

    def coverage(
        self,
        regions: List[Region],
        *,
        site: str = "A",
        by_sample: bool = False,
    ) -> pl.DataFrame: ...

    def size_factors(self) -> Dict[str, float]: ...


@runtime_checkable
class SupportsJunctions(Protocol):
    """Provider that can count reads spanning splice junctions (requires CIGAR).

    BAM and matrix providers satisfy this; bigwig providers do not.

    junctions: list of (chrom, donor_pos, acceptor_pos, strand_int, junction_id)
    Returns a DataFrame with columns: junction_id, count[, sample_id].
    """

    def junction_support(
        self,
        junctions: List[Tuple[str, int, int, int, int]],
        *,
        by_sample: bool = False,
    ) -> pl.DataFrame: ...


@runtime_checkable
class SupportsMappability(Protocol):
    """Provider that can emit per-event unique/multimapper accounting.

    Returns a DataFrame with MAPPABILITY_LEDGER_SCHEMA columns.
    """

    def mappability_ledger(self, events: pl.DataFrame) -> pl.DataFrame: ...

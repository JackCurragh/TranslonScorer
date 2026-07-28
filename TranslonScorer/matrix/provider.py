"""MatrixProvider — coverage from the sparse annotation-scale matrix.

Coverage comes from the P-site index built by `build-psite-index`, which
stores per-(sample, length) offset-corrected genomic P-site positions.
`site` selection shifts by one codon (A = P + 3).

This is the matrix half of the two coverage strategies: `MatrixProvider`
scans a pre-built index over thousands of samples, `BamSetProvider` /
`BigwigSetProvider` read a handful of files directly. They diverge only in
how per-locus coverage is obtained — the events, the scorer, the report and
the store below them are identical (see workflows._score_events_over_provider).

Capabilities (duck-typed; see coverage/base.py)
-----------------------------------------------
coverage(…, site="P"|"A", …)  — per-position P/A-site coverage
size_factors()                — per-sample depth factors
junction_support()            — spanning-read counts per junction
mappability_ledger()          — per-event unique/multimapper accounting (stub)
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import polars as pl

from TranslonScorer.coverage.base import MAPPABILITY_LEDGER_SCHEMA
from TranslonScorer.model import Region

PathLike = Union[str, Path]


class MatrixProvider:
    """CoverageProvider backed by the sparse annotation-scale BAM matrix.

    Parameters
    ----------
    partition_dirs  : path(s) to partition directories produced by the matrix
                      pipeline (each directory contains a unique-read BAM and
                      count/manifest parquets). Required — used for
                      junction_support()/mappability_ledger(), which read the
                      partition BAMs rather than the index.
    psite_index_dir : directory produced by `build-psite-index` (Phase 2).
                      Required. coverage() reads per-(sample, length)
                      offset-corrected genomic P-site positions from it. A
                      single flat offset across read lengths smears the P-site
                      across frames and pins elong_in_frame to the ~0.33 random
                      floor, so no flat-offset alternative is offered.
                      See psite_index.query_genomic_coverage.
    sample_names    : optional list of sample names to include (None = all).
    n_workers       : worker processes for parallel partition scanning
                      (junction_support only — coverage reads pre-built
                      Parquet shards directly, no scanning).
    """

    def __init__(
        self,
        partition_dirs: Union[PathLike, List[PathLike]],
        *,
        psite_index_dir: PathLike,
        sample_names: Optional[List[str]] = None,
        n_workers: Optional[int] = None,
    ) -> None:
        if isinstance(partition_dirs, (str, Path)):
            self._dirs: List[Path] = [Path(partition_dirs)]
        else:
            self._dirs = [Path(d) for d in partition_dirs]
        if not psite_index_dir:
            raise ValueError(
                "MatrixProvider requires psite_index_dir (build it with "
                "`translonscorer build-psite-index`). Matrix coverage must use "
                "per-(sample, length) P-site offsets; a flat offset collapses "
                "elong_in_frame to the ~0.33 random floor."
            )
        self._psite_index_dir = Path(psite_index_dir)
        self._sample_names = sample_names
        self._n_workers = n_workers

    # ------------------------------------------------------------------
    # CoverageProvider
    # ------------------------------------------------------------------

    def coverage(
        self,
        regions: List[Region],
        *,
        site: str = "A",
        by_sample: bool = False,
    ) -> pl.DataFrame:
        """Return per-position coverage over genomic regions.

        Positions come from the P-site index already offset-corrected per
        (sample, length), so site="P" needs no shift and site="A" applies the
        +3 nt (one codon) shift.

        Parameters
        ----------
        regions    : list of Region(chrom, start, end) namedtuples.
        site       : "A" (elongation/termination) or "P" (initiation).
        by_sample  : False → aggregate across samples; True → per-sample rows.

        Returns
        -------
        DataFrame with columns pos, count, strand[, sample_name].
        """
        if site not in {"P", "A"}:
            raise ValueError(f"site must be 'P' or 'A', got {site!r}")

        from TranslonScorer.matrix.psite_index import query_genomic_coverage

        return query_genomic_coverage(
            self._psite_index_dir,
            [(r.chrom, r.start, r.end) for r in regions],
            site=site,
            sample_names=self._sample_names,
            group_level="sample" if by_sample else "aggregate",
        )

    def size_factors(self) -> Dict[str, float]:
        """Per-sample library-size factors (stub: all 1.0 until ledger is built)."""
        return {}

    # ------------------------------------------------------------------
    # SupportsJunctions
    # ------------------------------------------------------------------

    def junction_support(
        self,
        junctions: List[Tuple[str, int, int, int, int]],
        *,
        by_sample: bool = False,
    ) -> pl.DataFrame:
        """Count reads spanning each junction event.

        Parameters
        ----------
        junctions : list of (chrom, donor_pos, acceptor_pos, strand_int, junction_id).
        by_sample : True → per-sample rows; False → aggregate.

        Returns
        -------
        DataFrame with columns junction_id, kind, count[, group].
        """
        from TranslonScorer.matrix.rollup import tabulate_junctions

        group_level = "sample" if by_sample else "aggregate"
        return tabulate_junctions(
            self._dirs,
            junctions,
            group_level=group_level,
            sample_names=self._sample_names,
            n_workers=self._n_workers,
        )

    # ------------------------------------------------------------------
    # SupportsMappability
    # ------------------------------------------------------------------

    def mappability_ledger(self, events: pl.DataFrame) -> pl.DataFrame:
        """Return empty mappability ledger (full implementation deferred to T12)."""
        return pl.DataFrame(schema=MAPPABILITY_LEDGER_SCHEMA)

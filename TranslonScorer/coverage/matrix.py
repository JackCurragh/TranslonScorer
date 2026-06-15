"""MatrixProvider — coverage from the sparse annotation-scale matrix.

Wraps `matrix_rollup.region_coverage` and `tabulate_junctions` behind
the capability-Protocol interface defined in `coverage/base.py`.

The matrix stores pre-computed A-site positions at a configurable offset;
`site` selection is supported by shifting the query region when asking for
P-site coverage (P = A − 3, so the query is shifted by +3 to align the
stored A-site positions).

Implements
----------
CoverageProvider    — coverage(), size_factors()
SupportsSites       — coverage(…, site="P"|"A", …)
SupportsJunctions   — junction_support()
SupportsMappability — mappability_ledger() (stub; returns empty ledger)
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import polars as pl

from TranslonScorer.model import Region
from TranslonScorer.coverage.base import MAPPABILITY_LEDGER_SCHEMA


PathLike = Union[str, Path]


class MatrixProvider:
    """CoverageProvider backed by the sparse annotation-scale BAM matrix.

    Parameters
    ----------
    partition_dirs : path(s) to partition directories produced by the matrix
                     pipeline (each directory contains a unique-read BAM and
                     count/manifest parquets).
    ref_offset     : default P-site read offset applied when building the
                     matrix (used to interpret stored positions).
    sample_names   : optional list of sample names to include (None = all).
    n_workers      : worker processes for parallel partition scanning.
    """

    def __init__(
        self,
        partition_dirs: Union[PathLike, List[PathLike]],
        *,
        ref_offset: int = 15,
        sample_names: Optional[List[str]] = None,
        n_workers: Optional[int] = None,
    ) -> None:
        if isinstance(partition_dirs, (str, Path)):
            self._dirs: List[Path] = [Path(partition_dirs)]
        else:
            self._dirs = [Path(d) for d in partition_dirs]
        self._ref_offset = int(ref_offset)
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

        The matrix stores positions at the P-site offset used during matrix
        construction (`ref_offset`).  When site="A" (default) positions are
        shifted forward by 3 nt to convert stored P-site to A-site coordinates.

        Parameters
        ----------
        regions    : list of Region(chrom, start, end) namedtuples.
        site       : "A" (elongation/termination) or "P" (initiation).
        by_sample  : False → aggregate across samples; True → per-sample rows.

        Returns
        -------
        DataFrame with columns pos, count[, strand][, sample_name].
        """
        if site not in {"P", "A"}:
            raise ValueError(f"site must be 'P' or 'A', got {site!r}")

        from TranslonScorer.matrix_rollup import region_coverage

        a_shift = 3 if site == "A" else 0
        query_regions = [Region(r.chrom, r.start - a_shift, r.end - a_shift) for r in regions]

        group_level = "sample" if by_sample else "aggregate"
        df = region_coverage(
            self._dirs,
            [(r.chrom, r.start, r.end) for r in query_regions],
            ref_offset=self._ref_offset,
            group_level=group_level,
            sample_names=self._sample_names,
            n_workers=self._n_workers,
        )
        if df.is_empty():
            schema: Dict[str, type] = {"pos": pl.Int64, "count": pl.Float64}
            if by_sample:
                schema["sample_name"] = pl.Utf8
            return pl.DataFrame(schema=schema)

        if a_shift and "pos" in df.columns:
            df = df.with_columns((pl.col("pos") + a_shift).alias("pos"))
        return df

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
        from TranslonScorer.matrix_rollup import tabulate_junctions

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


# ---------------------------------------------------------------------------
# coverage() adapter: dict → DataFrame (for testing / lightweight paths)
# ---------------------------------------------------------------------------


class DictCoverageProvider:
    """Minimal CoverageProvider wrapping a {pos: count} dict.

    Useful for unit tests and offline scoring (no matrix files required).
    Satisfies the CoverageProvider protocol; does NOT satisfy SupportsSites,
    SupportsJunctions, or SupportsMappability.
    """

    def __init__(self, cov_dict: Dict[int, float]) -> None:
        self._dict = cov_dict

    def coverage(
        self,
        regions: List[Region],
        *,
        site: str = "A",
        by_sample: bool = False,
    ) -> pl.DataFrame:
        """Return coverage for positions in the requested regions."""
        positions = []
        counts = []
        for r in regions:
            for pos in range(r.start, r.end):
                c = self._dict.get(pos)
                if c is not None:
                    positions.append(pos)
                    counts.append(float(c))
        return pl.DataFrame(
            {"pos": positions, "count": counts},
            schema={"pos": pl.Int64, "count": pl.Float64},
        )

    def size_factors(self) -> Dict[str, float]:
        return {"": 1.0}

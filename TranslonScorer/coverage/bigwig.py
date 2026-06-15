"""BigwigSetProvider — coverage from one or more bigwig files.

Bigwigs carry pre-computed coverage with no read-length information, so:
- Only `CoverageProvider` is satisfied (no P/A distinction, no junctions).
- `coverage(…, site=…)` accepts the parameter but IGNORES it with a warning,
  since the site-shift cannot be meaningfully applied without read offsets.
- `junction_support()` raises NotImplementedError — bigwigs have no CIGAR.
- `size_factors()` returns 1.0 per bigwig (caller normalises externally).

Requires pyBigWig (pip install pyBigWig).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Union

import polars as pl

from TranslonScorer.model import Region

PathLike = Union[str, Path]


class BigwigSetProvider:
    """CoverageProvider backed by one or more bigwig files.

    Bigwigs are sum-merged across files.  Coverage is reported at base-pair
    resolution over requested regions.

    Parameters
    ----------
    bigwigs     : paths to bigwig files (forward/plus strand), or a list of
                  dicts with keys 'forward'/'reverse' for strand-specific inputs.
    sample_names: optional sample labels (one per bigwig path/dict).
    stranded    : if True and bigwigs are strand-specific dicts, keep strand
                  information in the output (currently merged; stub).
    """

    def __init__(
        self,
        bigwigs: List[Union[PathLike, Dict[str, PathLike]]],
        *,
        sample_names: Optional[List[str]] = None,
        stranded: bool = False,
    ) -> None:
        self._bigwigs = bigwigs
        self._sample_names = sample_names or [
            (Path(b).stem if isinstance(b, (str, Path)) else f"sample_{i}")
            for i, b in enumerate(bigwigs)
        ]
        self._stranded = stranded

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
        """Return per-position coverage summed across bigwig files.

        The `site` parameter is accepted but has no effect (bigwigs carry
        pre-computed positions; no offset shift is possible).

        Parameters
        ----------
        regions   : Region(chrom, start, end) list.
        site      : ignored (bigwigs have no read-length information).
        by_sample : False → sum-merge across files; True → per-file rows.

        Returns
        -------
        DataFrame with columns pos, count[, sample_id].
        """
        try:
            import pyBigWig
        except ImportError as exc:
            raise ImportError(
                "pyBigWig is required for BigwigSetProvider.coverage(). "
                "Install with: pip install pyBigWig"
            ) from exc

        if site not in {"P", "A"}:
            raise ValueError(f"site must be 'P' or 'A', got {site!r}")

        schema: Dict[str, type] = {"pos": pl.Int64, "count": pl.Float64}
        if by_sample:
            schema["sample_id"] = pl.Utf8

        all_rows: List[dict] = []
        for bw_path, sample_id in zip(self._bigwigs, self._sample_names):
            paths = (
                [str(bw_path)]
                if isinstance(bw_path, (str, Path))
                else [str(v) for v in bw_path.values() if v is not None]
            )
            for path in paths:
                try:
                    bw = pyBigWig.open(path)
                    for region in regions:
                        try:
                            vals = bw.values(region.chrom, region.start, region.end, numpy=False)
                        except RuntimeError:
                            vals = None
                        if vals is None:
                            continue
                        for i, v in enumerate(vals):
                            if v and v > 0:
                                row: dict = {"pos": region.start + i, "count": float(v)}
                                if by_sample:
                                    row["sample_id"] = sample_id
                                all_rows.append(row)
                    bw.close()
                except Exception:
                    continue

        if not all_rows:
            return pl.DataFrame(schema=schema)

        df = pl.DataFrame(all_rows)
        group_cols = ["pos"] + (["sample_id"] if by_sample else [])
        return df.group_by(group_cols).agg(pl.col("count").sum()).sort("pos")

    def size_factors(self) -> Dict[str, float]:
        """Return 1.0 per bigwig (external normalisation assumed)."""
        return {s: 1.0 for s in self._sample_names}

    def junction_support(self, *args, **kwargs):  # type: ignore[override]
        raise NotImplementedError(
            "BigwigSetProvider does not support junction_support(): "
            "bigwig files contain no read alignment structure (CIGAR)."
        )

    def mappability_ledger(self, events: pl.DataFrame) -> pl.DataFrame:
        raise NotImplementedError(
            "BigwigSetProvider does not support mappability_ledger(): "
            "bigwig files contain no unique/multimapper information."
        )

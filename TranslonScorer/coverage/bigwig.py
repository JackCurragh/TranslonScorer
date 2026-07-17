"""BigwigSetProvider — coverage from one or more bigwig files.

Bigwigs carry pre-computed coverage with no read-length information, so:
- Only `CoverageProvider` is satisfied (no P/A distinction, no junctions).
- `coverage(…, site=…)` accepts the parameter but IGNORES it, since the
  site-shift cannot be meaningfully applied without read offsets.
- `junction_support()` raises NotImplementedError — bigwigs have no CIGAR.
- `size_factors()` returns 1.0 per bigwig (caller normalises externally).

Strand convention: a stranded entry is a dict with 'forward'/'reverse' keys
pointing at two separate bigwig files, both holding POSITIVE values (the
convention used elsewhere in this codebase, e.g. STAR's
Signal.Unique.str1/str2.out.bw) — strand is inferred from which file a value
came from, not from its sign. If your bigwigs instead encode strand via
negative values in a single file, negate the reverse-strand file before
passing it in (or pre-split it) — this provider does not sign-flip.

Requires pyBigWig (pip install pyBigWig).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

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
    stranded    : if True, every entry in `bigwigs` MUST be a
                  {'forward','reverse'} dict, and coverage() emits a `strand`
                  column (1 for forward, -1 for reverse) so
                  workflows._score_events_over_provider scores each strand
                  against its own coverage — required for correct Ribo-seq
                  scoring (see module docstring for the sign convention). If
                  False, entries are sum-merged regardless of shape (a
                  stranded dict's forward+reverse values are combined into
                  one unstranded pileup — only useful for coverage QC, not
                  strand-sensitive event scoring).
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
        DataFrame with columns pos, count[, strand (if stranded), sample_id
        (if by_sample)].
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
        if self._stranded:
            schema["strand"] = pl.Int64
        if by_sample:
            schema["sample_id"] = pl.Utf8

        all_rows: List[dict] = []
        for bw_entry, sample_id in zip(self._bigwigs, self._sample_names):
            # (path, strand_value) pairs to read for this entry. strand_value
            # is None when unstranded (no strand column emitted).
            if self._stranded:
                if not isinstance(bw_entry, dict) or not {"forward", "reverse"} <= bw_entry.keys():
                    raise ValueError(
                        "stranded=True requires every bigwig entry to be a dict with "
                        f"'forward'/'reverse' keys; got {bw_entry!r} for sample {sample_id!r}"
                    )
                sources: List[Tuple[Union[str, Path], Optional[int]]] = [
                    (bw_entry["forward"], 1),
                    (bw_entry["reverse"], -1),
                ]
            elif isinstance(bw_entry, (str, Path)):
                sources = [(bw_entry, None)]
            else:
                sources = [(v, None) for v in bw_entry.values() if v is not None]

            for path, strand_val in sources:
                if path is None:
                    continue
                try:
                    bw = pyBigWig.open(str(path))
                    chrom_lens = bw.chroms()
                    for region in regions:
                        # Clamp the (padded) window to the chromosome bounds: an
                        # event within _EVENT_FLANK_PAD of a contig end pushes
                        # region.end past the chrom length, and pyBigWig raises on
                        # out-of-bounds intervals — which would otherwise drop all
                        # coverage for the whole (merged) region. Skip chroms the
                        # bigwig doesn't carry.
                        clen = chrom_lens.get(region.chrom)
                        if clen is None:
                            continue
                        rstart = max(0, region.start)
                        rend = min(region.end, int(clen))
                        if rstart >= rend:
                            continue
                        try:
                            vals = bw.values(region.chrom, rstart, rend, numpy=False)
                        except RuntimeError:
                            vals = None
                        if vals is None:
                            continue
                        for i, v in enumerate(vals):
                            if v and v > 0:
                                row: dict = {"pos": rstart + i, "count": float(v)}
                                if self._stranded:
                                    row["strand"] = strand_val
                                if by_sample:
                                    row["sample_id"] = sample_id
                                all_rows.append(row)
                    bw.close()
                except Exception:
                    continue

        if not all_rows:
            return pl.DataFrame(schema=schema)

        df = pl.DataFrame(all_rows)
        group_cols = (
            ["pos"] + (["strand"] if self._stranded else []) + (["sample_id"] if by_sample else [])
        )
        return df.group_by(group_cols).agg(pl.col("count").sum()).sort(group_cols)

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

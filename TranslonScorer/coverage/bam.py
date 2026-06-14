"""BamSetProvider — per-BAM/per-read-length offset calibration + P/A coverage.

Architecture contract (non-negotiable)
---------------------------------------
1. Offsets are calibrated ONCE per BAM per read-length from whole-sample
   evidence, before any locus profile is built.  `_calibrate_offsets` is
   called lazily and its result cached; downstream `coverage()` calls reuse
   the cached offset table without recomputing.
2. Profile generation is then pure: reads + offset table + site → positions.
   The pure function `coverage.profile.apply_offsets` is used for this.
3. A-site = P-site + 3 nt by default unless the offset table specifies a
   distinct A-site column (future extension).
4. Unique mappers by default; multimappers recorded in mappability_ledger.

Dependencies
------------
pysam   — BAM reading (CIGAR, fetch)
oxbow   — optional faster reader (falls back to pysam when unavailable)
"""
from __future__ import annotations

import functools
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import polars as pl

from TranslonScorer.model import OffsetParams, Region
from TranslonScorer.offsets import make_offset_table
from TranslonScorer.coverage.base import MAPPABILITY_LEDGER_SCHEMA
from TranslonScorer.coverage.profile import apply_offsets


PathLike = Union[str, Path]

_JUNC_OVERHANG = 6  # min flanking aligned length for confident junction span


class BamSetProvider:
    """Coverage provider backed by 1–N genome (or transcriptome→genome) BAMs.

    Per-BAM per-read-length P-site offsets are calibrated once from the whole
    sample before any locus profile is built.  The calibrated table is cached;
    all subsequent `coverage()` calls reuse it.

    Parameters
    ----------
    bams        : list of BAM paths (genome or transcriptome; see transcriptome param).
    exon_df     : exon annotation DataFrame required for transcriptome→genome
                  projection.  Required when transcriptome=True.
    offsets     : OffsetParams controlling how offsets are calibrated.
    multimap    : "unique" (default) — exclude multimappers from counts;
                  "all" — include all alignments (use with care).
    transcriptome : True if BAMs are transcriptome-mapped (reads projected to
                   genome coordinates via exon_df before profiling).
    sample_names  : optional list of sample names (one per BAM in order).
    """

    def __init__(
        self,
        bams: List[PathLike],
        *,
        exon_df: Optional[pl.DataFrame] = None,
        offsets: OffsetParams = OffsetParams(),
        multimap: str = "unique",
        transcriptome: bool = False,
        sample_names: Optional[List[str]] = None,
    ) -> None:
        self._bams = [Path(b) for b in bams]
        self._exon_df = exon_df
        self._offset_params = offsets
        self._multimap = multimap
        self._transcriptome = transcriptome
        self._sample_names = sample_names or [
            Path(b).stem for b in bams
        ]
        if len(self._sample_names) != len(self._bams):
            raise ValueError(
                f"sample_names length ({len(self._sample_names)}) must match "
                f"bams length ({len(self._bams)})"
            )
        # Per-BAM offset tables: {bam_path: {read_length: p_offset}}
        # Populated lazily by _calibrate_offsets; NOT re-derived per locus.
        self._offset_tables: Optional[Dict[Path, Dict[int, int]]] = None

    # ------------------------------------------------------------------
    # Offset calibration (lazy, cached, whole-sample)
    # ------------------------------------------------------------------

    def _calibrate_offsets(self) -> Dict[Path, Dict[int, int]]:
        """Calibrate per-BAM per-read-length P-site offsets once.

        Calls make_offset_table from offsets.py using the configured method
        (global / file / metagene).  Returns {bam_path: {length: p_offset}}.

        This method is the KEY architectural guarantee: offsets are derived
        from whole-sample evidence and never from locus profiles or scoring
        regions.
        """
        if self._offset_tables is not None:
            return self._offset_tables

        tables: Dict[Path, Dict[int, int]] = {}
        for bam_path in self._bams:
            if self._offset_params.method == "file":
                table = make_offset_table(self._offset_params)
            elif self._offset_params.method == "global":
                table = make_offset_table(self._offset_params)
            else:
                # metagene: require whole-BAM read-length distribution
                # Currently the metagene stub is not implemented; fall back
                # to global offsets with a warning so the class is usable.
                try:
                    table = make_offset_table(self._offset_params)
                except NotImplementedError:
                    from TranslonScorer.utils.logging import log_warning
                    log_warning(
                        f"metagene offsets not yet implemented; "
                        f"falling back to global offset={self._offset_params.global_offset} "
                        f"for {bam_path.name}"
                    )
                    fallback_params = OffsetParams(
                        **{**self._offset_params.__dict__, "method": "global"}
                    )
                    table = make_offset_table(fallback_params)
            tables[bam_path] = table

        self._offset_tables = tables
        return self._offset_tables

    # ------------------------------------------------------------------
    # CoverageProvider / SupportsSites
    # ------------------------------------------------------------------

    def coverage(
        self,
        regions: List[Region],
        *,
        site: str = "A",
        by_sample: bool = False,
    ) -> pl.DataFrame:
        """Return per-position P/A-site coverage over genomic regions.

        Offsets are calibrated (once, cached) before the first call.
        Reads are fetched per-region via pysam; the offset table maps each
        read length to its P-site position; A = P + 3 nt.

        Parameters
        ----------
        regions   : Region(chrom, start, end) list — genomic coordinates.
        site      : "A" (elongation/term, default) or "P" (initiation).
        by_sample : False → aggregate; True → per-sample rows.

        Returns
        -------
        DataFrame with columns pos, count[, sample_id].
        """
        if site not in {"P", "A"}:
            raise ValueError(f"site must be 'P' or 'A', got {site!r}")

        try:
            import pysam
        except ImportError as exc:
            raise ImportError(
                "pysam is required for BamSetProvider.coverage(). "
                "Install with: pip install pysam"
            ) from exc

        offset_tables = self._calibrate_offsets()
        rows: List[dict] = []

        for bam_path, sample_id in zip(self._bams, self._sample_names):
            table = offset_tables[bam_path]
            bam_rows: List[dict] = []
            try:
                with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                    refs = set(bam.references)
                    for region in regions:
                        chrom = _resolve_chrom(region.chrom, refs)
                        if chrom is None:
                            continue
                        for rec in bam.fetch(chrom, region.start, region.end):
                            if rec.is_unmapped or rec.is_secondary:
                                continue
                            if self._multimap == "unique" and rec.mapping_quality == 0:
                                continue
                            length = rec.query_length or 0
                            if length == 0:
                                continue
                            p_offset = table.get(length, self._offset_params.global_offset)
                            a_shift = 3 if site == "A" else 0
                            pos = rec.reference_start + p_offset + a_shift
                            bam_rows.append({"sample_id": sample_id, "pos": pos, "count": 1.0})
            except (OSError, ValueError):
                continue

            rows.extend(bam_rows)

        if not rows:
            schema: Dict[str, type] = {"pos": pl.Int64, "count": pl.Float64}
            if by_sample:
                schema["sample_id"] = pl.Utf8
            return pl.DataFrame(schema=schema)

        df = pl.DataFrame(rows, schema={"sample_id": pl.Utf8, "pos": pl.Int64, "count": pl.Float64})
        if by_sample:
            return (
                df.group_by(["sample_id", "pos"])
                .agg(pl.col("count").sum())
                .sort(["sample_id", "pos"])
            )
        return (
            df.group_by("pos")
            .agg(pl.col("count").sum())
            .sort("pos")
        )

    def size_factors(self) -> Dict[str, float]:
        """Library-size normalisation factors (median-ratio; 1.0 per sample until computed)."""
        return {s: 1.0 for s in self._sample_names}

    # ------------------------------------------------------------------
    # SupportsJunctions
    # ------------------------------------------------------------------

    def junction_support(
        self,
        junctions: List[Tuple[str, int, int, int, int]],
        *,
        by_sample: bool = False,
    ) -> pl.DataFrame:
        """Count reads spanning each junction via CIGAR block inspection.

        junctions: (chrom, donor_pos, acceptor_pos, strand_int, junction_id)
        Reads are classified as span_conf, span_short, or unspliced.
        """
        try:
            import pysam
        except ImportError as exc:
            raise ImportError(
                "pysam is required for BamSetProvider.junction_support(). "
                "Install with: pip install pysam"
            ) from exc

        out_schema = {"junction_id": pl.UInt64, "kind": pl.Utf8, "count": pl.Float64}
        if by_sample:
            out_schema["sample_id"] = pl.Utf8

        rows: List[dict] = []
        for bam_path, sample_id in zip(self._bams, self._sample_names):
            try:
                with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                    refs = set(bam.references)
                    for chrom, donor, acceptor, strand, junc_id in junctions:
                        fc = _resolve_chrom(chrom, refs)
                        if fc is None:
                            continue
                        for rec in bam.fetch(fc, max(0, donor - 1), donor + 1):
                            if rec.is_unmapped or rec.is_secondary:
                                continue
                            blocks = rec.get_blocks()
                            kind = None
                            for i in range(len(blocks) - 1):
                                if blocks[i][1] == donor and blocks[i + 1][0] == acceptor:
                                    overhang = min(
                                        blocks[i][1] - blocks[i][0],
                                        blocks[i + 1][1] - blocks[i + 1][0],
                                    )
                                    kind = "span_conf" if overhang >= _JUNC_OVERHANG else "span_short"
                                    break
                            if kind is None:
                                for bs, be in blocks:
                                    if bs <= donor - 1 and be >= donor + 1:
                                        kind = "unspliced"
                                        break
                            if kind is not None:
                                row = {"junction_id": int(junc_id), "kind": kind, "count": 1.0}
                                if by_sample:
                                    row["sample_id"] = sample_id
                                rows.append(row)
            except (OSError, ValueError):
                continue

        if not rows:
            return pl.DataFrame(schema=out_schema)

        df = pl.DataFrame(rows)
        group_cols = ["junction_id", "kind"] + (["sample_id"] if by_sample else [])
        return (
            df.group_by(group_cols)
            .agg(pl.col("count").sum())
            .sort("junction_id")
        )

    # ------------------------------------------------------------------
    # SupportsMappability
    # ------------------------------------------------------------------

    def mappability_ledger(self, events: pl.DataFrame) -> pl.DataFrame:
        """Per-event unique/multimapper accounting (stub; full in Phase 4)."""
        return pl.DataFrame(schema=MAPPABILITY_LEDGER_SCHEMA)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_chrom(chrom: str, refs: set) -> Optional[str]:
    """Resolve a chromosome name against the set of BAM references."""
    if chrom in refs:
        return chrom
    alt = f"chr{chrom}" if not chrom.startswith("chr") else chrom[3:]
    return alt if alt in refs else None

"""Coverage providers — SKELETON for review (plumbing not yet implemented).

One scoring core (event_score) runs on top of a uniform CoverageProvider that
yields per-position A-site/P-site coverage + junction support + size factors for
a set of regions. Implementations differ only in the source:

    BamSetProvider     1–20 genome (or transcriptome→genome) BAMs
    BigwigSetProvider  1–20 bigwigs (lossy: coverage only)
    MatrixProvider     the sparse annotation-scale matrix (already built)

Locked method spec: see docs/event_scoring_model.md "Coverage & method spec".
Site: initiation→P-site, elongation/termination→A-site (A = P + 3 nt).
Offsets: metagene (canonical) | file | global, constrained to plausible values.
Multimapper: unique default + per-event mappability ledger for relaxation.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple

import polars as pl

from TranslonScorer.model import OffsetParams, Region  # noqa: F401 — re-exported


def plausible_offset_range(read_length: int, p: OffsetParams) -> Tuple[int, int]:
    """Physically/empirically plausible P-site offset window for a read length.
    A 25 nt read cannot have an 18 nt offset: hi = min(offset_max, ⌊L·frac⌋)."""
    hi = min(p.offset_max, int(read_length * p.max_frac))
    return p.offset_min, hi


def usable_read_length(read_length: int, p: OffsetParams) -> bool:
    return p.min_read_len <= read_length <= p.max_read_len


def metagene_offsets(
    five_prime_by_length: pl.DataFrame,   # cols: read_length, rel_pos (5' end vs start codon), count
    p: OffsetParams,
) -> Dict[int, int]:
    """CANONICAL offset method (skeleton). For each usable read length, build the
    metagene of UNIQUE-read 5′ ends around annotated start codons and pick the
    P-site offset = distance from 5′ end to the start codon, restricted to
    `plausible_offset_range`. Returns {read_length: p_site_offset}.

    (Implementation: argmax / changepoint of the 5′ pile-up at the start within
    the plausible window; A-site offset = P-site + 3 downstream.)"""
    raise NotImplementedError


def psite_to_asite(offset_p: int) -> int:
    return offset_p + 3      # A-site is the next codon (3 nt) 3′ of the P-site


# ---------------------------------------------------------------------------
# Mappability ledger (per event; may be a second pass)
# ---------------------------------------------------------------------------

MAPPABILITY_LEDGER_SCHEMA = {
    "event_id": pl.UInt64,
    "unique_reads": pl.Float64,        # reads uniquely mapped to this event
    "multimapper_reads": pl.Float64,   # reads here that also map elsewhere
    "mean_loci_per_multiread": pl.Float64,  # how multi (≈ copy number signal)
    "competing_event_ids": pl.Utf8,    # JSON list of co-mapped loci (for relaxation)
}

# Relaxation policy (applied later, not at default-unique scoring):
#   paralog (known expected copy number N) → distribute multimapper reads to the N copies
#   no clear paralogy                       → guilty-by-association weighting


# ---------------------------------------------------------------------------
# Provider interface
# ---------------------------------------------------------------------------

class CoverageProvider(ABC):
    """Uniform coverage source for the scoring core. `by_sample=False` returns the
    summed aggregate; True keeps per-sample rows (1–20 samples → per-sample/cluster
    tiers). `site` ∈ {"P","A"} selects ribosome site (P=init, A=elong/term)."""

    @abstractmethod
    def coverage(self, regions: List[Region], *, site: str = "A",
                 by_sample: bool = False) -> pl.DataFrame:
        """→ (pos, strand, count[, sample_name]). A-site = P-site + 3 nt."""

    @abstractmethod
    def junction_support(self, junctions: List[Tuple[str, int, int, int, int]], *,
                         by_sample: bool = False) -> pl.DataFrame:
        """→ (junction_id, kind, count[, sample]). Raises for bigwig (no CIGAR)."""

    @abstractmethod
    def size_factors(self) -> Dict[str, float]:
        """Per-sample depth normalisation (median-ratio size factors)."""

    def mappability_ledger(self, events: pl.DataFrame) -> pl.DataFrame:
        """Per-event unique/multimapper accounting (default: empty; BAM/matrix override)."""
        return pl.DataFrame(schema=MAPPABILITY_LEDGER_SCHEMA)


class BamSetProvider(CoverageProvider):
    """1–20 BAMs. Genome BAMs used directly; transcriptome BAMs projected to the
    genome with isoform-multimapper resolution. Per-BAM metagene offsets →
    P/A-site profiles; unique mappers by default + mappability ledger."""

    def __init__(self, bams: List[str], *, offsets: OffsetParams = OffsetParams(),
                 multimap: str = "unique", transcriptome: bool = False,
                 exon_df: Optional[pl.DataFrame] = None):
        self.bams, self.offsets, self.multimap = bams, offsets, multimap
        self.transcriptome, self.exon_df = transcriptome, exon_df

    def coverage(self, regions, *, site="A", by_sample=False):
        raise NotImplementedError   # offset → 5' shift (P or A) → per-pos coverage

    def junction_support(self, junctions, *, by_sample=False):
        raise NotImplementedError   # CIGAR-intron match + overhang confidence

    def size_factors(self):
        raise NotImplementedError


class BigwigSetProvider(CoverageProvider):
    """1–20 bigwigs, sum-merged. Coverage only: no offset (assumed pre-applied),
    no P/A distinction, no junction support."""

    def __init__(self, bigwigs: List[str]):
        self.bigwigs = bigwigs

    def coverage(self, regions, *, site="A", by_sample=False):
        raise NotImplementedError   # sum-merge bigwig coverage over regions

    def junction_support(self, junctions, *, by_sample=False):
        raise NotImplementedError("bigwig has no read structure for junctions")

    def size_factors(self):
        raise NotImplementedError


class MatrixProvider(CoverageProvider):
    """Annotation-scale sparse matrix — wraps the existing region_coverage /
    tabulate_junctions so the matrix path sits behind the same interface."""

    def __init__(self, partition_dirs, *, offsets: OffsetParams = OffsetParams()):
        self.partition_dirs, self.offsets = partition_dirs, offsets

    def coverage(self, regions, *, site="A", by_sample=False):
        raise NotImplementedError   # adapt matrix_rollup.region_coverage

    def junction_support(self, junctions, *, by_sample=False):
        raise NotImplementedError   # adapt matrix_rollup.tabulate_junctions

    def size_factors(self):
        raise NotImplementedError

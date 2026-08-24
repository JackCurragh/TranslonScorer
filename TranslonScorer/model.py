"""Central dataclasses: zero runtime dependencies, importable from anywhere.

All config and immutable records live here so the rest of the codebase has a
single, dep-free place to import shared types.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

# ---------------------------------------------------------------------------
# Coverage / offset config
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int


@dataclass(frozen=True)
class OffsetParams:
    offset_min: int = 8
    offset_max: int = 20
    max_frac: float = 0.6667  # offset ≤ ⌊read_len × max_frac⌋
    min_read_len: int = 25
    max_read_len: int = 35
    method: str = "metagene"  # "metagene" | "file" | "global"
    global_offset: int = 12  # used when method == "global" (P-site)
    offsets_file: Optional[str] = None


# ---------------------------------------------------------------------------
# Scoring config
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ScoreThresholds:
    version: str = "v1"
    min_reads: float = 20.0  # eligibility floor (all aspects)
    init_rise: float = 1.00  # log2 fold-change (body vs leader) for SUPPORTED start
    term_drop: float = 1.00  # log2 fold-change (body vs UTR) for SUPPORTED stop
    step_pseudocount: float = 1.0  # density floor α in log2FC — de-saturates
    flank_peakiness: float = 6.0  # evidence-only review flag
    instability: float = 0.40  # evidence-only review flag
    elong_in_frame: float = 0.50  # in-frame fraction for SUPPORTED elongation
    elong_breadth: float = 0.20  # fraction of in-frame codons covered
    elong_identifiability: float = 0.50  # contended attribution below this → AMBIGUOUS
    junc_min_spanning: float = 20.0  # spanning reads for SUPPORTED junction
    mappability_low: float = 0.50  # map_track_mean below this -> map_track_low=True (evidence-only)
    periodicity_min_codons: int = 5  # per-side codon floor for the boundary significance test
    periodicity_significance_alpha: float = 0.05  # periodicity_p below this counts as "significant"
    periodicity_min_agree_frac: float = 0.5  # fraction of flank lengths that must agree to
    # resolve an AMBIGUOUS init/term call to SUPPORTED via periodicity (see _decide_step)


# ---------------------------------------------------------------------------
# Frame support config
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FrameSupportParams:
    frame_method: str = "none"  # "none"|"linear"|"linear+hmm"|"deblur+linear+hmm"|"latent"
    frame_by_length: bool = False
    frame_background: str = "uniform"
    frame_hmm_lambda: float = 1.0


# ---------------------------------------------------------------------------
# Consequentiality policy (weights for tier × confidence × expression × context)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConsequentialityPolicy:
    min_tier_confidence: float = 0.0
    min_expression_percentile: float = 0.0
    context_weight: float = 1.0
    # consequentiality = tier_confidence
    #                    * (expression_floor + expression_weight * expression_pct)
    #                    * context_weight
    # floor=0.5/weight=0.5 means expression can at most double a translon's rank
    # and never drives it to zero on its own.
    expression_floor: float = 0.5
    expression_weight: float = 0.5
    # Which aspects count towards tier_confidence.  None = every aspect present
    # in the report.  Defaults to the linear translation chain, so a translon is
    # not penalised in confidence for an unsupported junction unless asked.
    chain_aspects: Optional[Tuple[str, ...]] = ("init", "elongation", "term")

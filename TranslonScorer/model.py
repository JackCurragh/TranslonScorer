"""Central dataclasses: zero runtime dependencies, importable from anywhere.

All config and immutable records live here so the rest of the codebase has a
single, dep-free place to import shared types.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


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
    max_frac: float = 0.6667        # offset ≤ ⌊read_len × max_frac⌋
    min_read_len: int = 25
    max_read_len: int = 35
    method: str = "metagene"        # "metagene" | "file" | "global"
    global_offset: int = 12         # used when method == "global" (P-site)
    offsets_file: Optional[str] = None


# ---------------------------------------------------------------------------
# Scoring config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ScoreThresholds:
    version: str = "v1"
    min_reads: float = 20.0          # eligibility floor (all aspects)
    init_rise: float = 1.00          # log2 fold-change (body vs leader) for SUPPORTED start
    term_drop: float = 1.00          # log2 fold-change (body vs UTR) for SUPPORTED stop
    step_pseudocount: float = 1.0    # density floor α in log2FC — de-saturates
    flank_peakiness: float = 6.0     # evidence-only review flag
    instability: float = 0.40        # evidence-only review flag
    elong_in_frame: float = 0.50     # in-frame fraction for SUPPORTED elongation
    elong_breadth: float = 0.20      # fraction of in-frame codons covered
    elong_identifiability: float = 0.50  # contended attribution below this → AMBIGUOUS
    junc_min_spanning: float = 20.0  # spanning reads for SUPPORTED junction


# ---------------------------------------------------------------------------
# Score record (typed mirror of _RECORD_SCHEMA in scoring/run.py)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ScoreRecord:
    event_id: int
    aspect: str
    group: str
    tier: str
    n_reads: float
    metric: Optional[float]
    metric_name: Optional[str]
    eligibility: str
    call: Optional[str]
    evidence: str
    thresholds_version: str


# ---------------------------------------------------------------------------
# Frame support config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FrameSupportParams:
    frame_method: str = "none"        # "none"|"linear"|"linear+hmm"|"deblur+linear+hmm"|"latent"
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

"""
Core functionality for TranslonScorer.

This module contains the core algorithms and functionality including:
- ORF prediction and classification
- Scoring algorithms and metrics
- Coordinate transformations
"""

from .coordinates import change_point_analysis, classify_orf, orfrelativeposition
from .orf_scoring import (
    assigningscore,
    calculate_scores,
    existingscore,
    globalscores,
    newscoring,
    oldscoring,
    sru_score,
)
from .orffinder import build_codon_automaton, find_all_positions, find_orfs, preporfs

__all__ = [
    # Scoring functions
    "sru_score",
    "calculate_scores",
    "oldscoring",
    "newscoring",
    "globalscores",
    "existingscore",
    "assigningscore",
    # Coordinate functions
    "classify_orf",
    "orfrelativeposition",
    "change_point_analysis",
    # ORF finding functions
    "find_orfs",
    "preporfs",
    "find_all_positions",
    "build_codon_automaton",
]

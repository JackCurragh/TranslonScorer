"""
Core functionality for TranslonScorer.

This module contains the core algorithms and functionality including:
- ORF prediction and classification
- Scoring algorithms and metrics
- Coordinate transformations
"""

from .scoring import (
    sru_score,
    calculate_scores,
    oldscoring,
    newscoring,
    globalscores,
    existingscore,
    assigningscore
)

from .coordinates import (
    classify_orf,
    orfrelativeposition
)

from .orffinder import (
    find_orfs,
    preporfs,
    find_all_positions,
    build_codon_automaton
)

__all__ = [
    # Scoring functions
    'sru_score',
    'calculate_scores',
    'oldscoring',
    'newscoring',
    'globalscores',
    'existingscore',
    'assigningscore',
    
    # Coordinate functions
    'classify_orf',
    'orfrelativeposition',
    
    # ORF finding functions
    'find_orfs',
    'preporfs',
    'find_all_positions',
    'build_codon_automaton'
] 
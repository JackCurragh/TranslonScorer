
"""
File handling functionality for TranslonScorer.

This package contains modules for handling different file formats:
- BAM file processing
- BED file processing
- BigWig file processing
"""

from .workflow import (
    find_orfs_workflow,
    score_orfs_workflow,
    plot_workflow,
    all_workflow,
)
from .config import Config
from .validator import validate_config

__all__ = [
    'find_orfs_workflow',
    'score_orfs_workflow',
    'plot_workflow',
    'all_workflow',
    'Config',
    'validate_config',
] 
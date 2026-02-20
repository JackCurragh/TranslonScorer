"""
File handling functionality for TranslonScorer.

Avoid importing heavy optional dependencies at package import time to
prevent unnecessary ImportError/ABI issues (e.g., pyBigWig vs NumPy).

Import submodules explicitly where used, e.g.:

    from TranslonScorer.file_handlers import bam as bam_handlers
    from TranslonScorer.file_handlers import bigwig as bigwig_handlers
    from TranslonScorer.file_handlers import bed as bed_handlers

This module intentionally does not import submodules eagerly.
"""

__all__ = []

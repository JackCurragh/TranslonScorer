"""
File handling functionality for TranslonScorer.

This package contains modules for handling different file formats:
- BAM file processing
- BED file processing
- BigWig file processing
"""

from .bam import (
    get_bam_tran,
    bamtranscript,
    process_transcriptomic_bam,
    calculate_differences,
    detect_bam_type,
)
from .bed import (
    asitecalc,
    bedtobigwig,
)
from .bigwig import (
    transcriptreads,
    scoring,
)

__all__ = [
    'get_bam_tran',
    'bamtranscript',
    'process_transcriptomic_bam',
    'calculate_differences',
    'detect_bam_type',
    'asitecalc',
    'bedtobigwig',
    'transcriptreads',
    'scoring',
] 
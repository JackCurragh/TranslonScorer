"""Thin re-export shim — implementation moved to TranslonScorer/scoring/.

All symbols previously defined here are now imported from scoring/aspects.py,
scoring/evidence.py, and scoring/run.py.  This shim keeps existing callers
working until the final cleanup task (T14) removes pipeline/.
"""
from __future__ import annotations

# Re-export model types (already re-exported from T1)
from TranslonScorer.model import ScoreThresholds  # noqa: F401

# Re-export io helper (already re-exported from T5)
from TranslonScorer.io.store import persist_scores  # noqa: F401

# Re-export all scoring symbols
from TranslonScorer.scoring.evidence import (  # noqa: F401
    _RECORD_SCHEMA,
    _decide_step,
    _elong_evidence,
    _jsonable,
    event_record,
)
from TranslonScorer.scoring.aspects import (  # noqa: F401
    abs_frame,
    score_elongation_event,
    score_initiation_event,
    score_junction_event,
    score_termination_event,
)
from TranslonScorer.scoring.run import (  # noqa: F401
    DEFAULT_THRESHOLDS,
    _prefix_sums,
    _range_sums,
    score_elongation_batch,
    score_events,
    score_events_vectorised,
)

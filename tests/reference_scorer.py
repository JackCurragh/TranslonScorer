"""Scalar reference scorer — test-only independent implementation.

This is the straightforward "loop over events, score each one" scorer. It used
to ship in ``TranslonScorer/scoring/run.py`` alongside the vectorised scorer,
with a docstring asking the two to be kept in sync by hand. Two implementations
of the same decision logic in product code is a standing drift risk, so the
scalar one lives here instead: it exists purely to be compared against the
shipped ``scoring.run.score_events``.

Keeping it (rather than deleting it) is deliberate. Elongation in the product
scorer is a prefix-sum kernel whose arithmetic is not obvious by inspection;
an independent, obviously-correct implementation to diff against is worth more
than the lines it costs. It must never be imported by product code.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import polars as pl

from TranslonScorer.events import SpliceContext
from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.aspects import (
    score_elongation_event,
    score_initiation_event,
    score_junction_event,
    score_termination_event,
)
from TranslonScorer.scoring.evidence import _RECORD_SCHEMA, event_record
from TranslonScorer.scoring.run import DEFAULT_THRESHOLDS


def score_events_scalar(
    events: pl.DataFrame,
    coverage: Dict[int, float],
    *,
    overlaps: Optional[Dict[int, list]] = None,
    junction_support: Optional[Dict[int, dict]] = None,
    group: str = "aggregate",
    tier: str = "aggregate",
    thr: ScoreThresholds = DEFAULT_THRESHOLDS,
    splice_context: Optional[SpliceContext] = None,
    map_track: Optional[Dict[int, dict]] = None,
) -> pl.DataFrame:
    """Score every event one at a time → the same long-form table as the product.

    Takes coverage as a plain ``{pos: count}`` dict rather than the product's
    ``cov_df``; otherwise the signature and output mirror
    ``scoring.run.score_events`` exactly.
    """
    overlaps = overlaps or {}
    junction_support = junction_support or {}
    map_track = map_track or {}
    rows: List[dict] = []
    for r in events.iter_rows(named=True):
        t, eid = r["type"], r["event_id"]
        if t == "init":
            raw = score_initiation_event(
                r["start"],
                r["strand"],
                coverage,
                thr=thr,
                chrom=r.get("chrom"),
                splice_context=splice_context,
            )
        elif t == "term":
            raw = score_termination_event(
                r["start"],
                r["strand"],
                coverage,
                thr=thr,
                chrom=r.get("chrom"),
                splice_context=splice_context,
            )
        elif t == "elongation":
            raw = score_elongation_event(
                r["start"],
                r["end"],
                r["phase"],
                r["strand"],
                coverage,
                overlaps.get(eid),
                thr=thr,
            )
        elif t == "junction":
            raw = score_junction_event(junction_support.get(eid, {}), thr=thr)
        else:
            continue
        if eid in map_track:
            raw = {**raw, **map_track[eid]}
        rows.append(event_record(eid, t, group, tier, raw, thr.version))
    return (
        pl.from_dicts(rows, schema=_RECORD_SCHEMA) if rows else pl.DataFrame(schema=_RECORD_SCHEMA)
    )

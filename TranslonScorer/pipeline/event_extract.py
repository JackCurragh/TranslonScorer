"""Extract deduplicated genomic *events* from the translon annotation.

Translons (features) decompose into shared genomic events; scoring each event
once (rather than per translon) is ~50x less work and makes shared evidence
consistent. This module turns ``translon_blocks`` + ``translons`` into:

  events        — unique genomic events (init / term / elongation / junction)
  feature_event — translon → event membership (role, rank, phase)
  event_overlap — elongation events that overlap another frame (contention)

Event identity is a deterministic hash of the genomic key, so a rebuilt
annotation reproduces the same ``event_id`` and prior scores still join.

Elongation phase is matrix-engine-consistent:
    + strand: phase = (ts - start) % 3
    - strand: phase = (ts + end - 1) % 3
where ``ts`` is the CDS-relative offset of the block's 5' end (cumulative length
of prior blocks in translation order). Two translons using a genomic segment in
the same reading-frame register share the elongation event.
"""
from __future__ import annotations

# Re-export shims — implementations live in TranslonScorer.events
from TranslonScorer.events import (  # noqa: F401
    _contention,
    _eid,
    _mod3,
    extract_events,
    run_extract,
)

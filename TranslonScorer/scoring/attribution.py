"""Per-aspect lens batteries and neighbor attribution
(docs/significance_testing_plan.md §0, §6, "Sequencing" step 5).

Two separate jobs live here, deliberately kept apart:

1. **Battery / confidence** (`*_battery` functions): given one event's already-
   computed evidence dict, how many of that aspect's *distinct, convergent*
   lenses agree? This is the `confidence` basis for §6 (see
   ``scoring/evidence.py``'s ``_RECORD_SCHEMA`` and ``compose_confidence``),
   not a new eligibility/call gate -- batteries never override
   ``_decide_step``/``_elong_evidence``'s existing decisions.

2. **Neighbor attribution** (`classify_neighbor_outcome` + the
   ``pair_*`` functions): given two candidate events in different frames at
   the same locus, EACH scored independently and on its OWN terms (never
   head-to-head), which of the three outcomes from §0 applies? This never
   invents a winner -- both events keep their own call/confidence; this
   only adds a label describing the relationship between them.

Every threshold used below lives on ``ScoreThresholds``, tagged first-pass/
unvalidated in ``model.py``, and is deliberately open to retuning once a
real cohort exists to validate against (see
docs/significance_testing_results.md).
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple

from TranslonScorer.model import ScoreThresholds

# ---------------------------------------------------------------------------
# Battery / confidence
# ---------------------------------------------------------------------------


def _agree(flags: Iterable[Optional[bool]]) -> Tuple[int, int]:
    """(n_agree, n_total) over a battery's lens flags, where None means the
    lens could not be evaluated (e.g. below its own codon/read floor) and is
    excluded from both counts -- an untested lens is neither agreement nor
    disagreement."""
    present = [f for f in flags if f is not None]
    return sum(1 for f in present if f), len(present)


def init_battery(evidence: dict, thr: ScoreThresholds) -> dict:
    """Initiation: Level (did density rise at all), Phase (periodicity_p
    significant at most flanks -- reuses `_decide_step`'s own majority-of-
    flanks logic), and two independent Uniformity readings that answer
    different questions (see significance_testing_plan.md §1's
    implementation note): `uniformity_gini` (magnitude concentration, low
    gini_body_inframe at most flanks) and `uniformity_share`
    (body_share_consistency_p NOT significant at most flanks -- the
    near/far halves of the body flank agree on frame-0 share)."""
    axes = list((evidence.get("boundary_axes_by_flank") or {}).values())
    level = (
        None if evidence.get("consensus_rise") is None else evidence["consensus_rise"] > 0
    )

    def _flank_majority(field: str, is_pass) -> Optional[bool]:
        vals = [a[field] for a in axes if a.get(field) is not None]
        if not vals:
            return None
        return (sum(1 for v in vals if is_pass(v)) / len(vals)) >= thr.periodicity_min_agree_frac

    phase = _flank_majority("periodicity_p", lambda p: p < thr.periodicity_significance_alpha)
    uniformity_gini = _flank_majority("gini_body_inframe", lambda g: g < thr.init_uniformity_gini_max)
    uniformity_share = _flank_majority(
        "body_share_consistency_p", lambda p: p >= thr.periodicity_significance_alpha
    )

    n_agree, n_total = _agree([level, phase, uniformity_gini, uniformity_share])
    return {
        "lenses": {
            "level": level,
            "phase": phase,
            "uniformity_gini": uniformity_gini,
            "uniformity_share": uniformity_share,
        },
        "n_agree": n_agree,
        "n_total": n_total,
        "confidence": (n_agree / n_total) if n_total else None,
        "battery_passed": (n_total > 0) and (n_agree / n_total) >= thr.lens_agree_frac,
    }


def term_battery(evidence: dict, thr: ScoreThresholds) -> dict:
    """Termination: Level (dropoff_significance_p significant), Continuity
    (dropoff_continuity_p NOT significant -- large p means the before-stop
    window looks like ordinary upstream elongation rather than an isolated
    pileup). Downstream/readthrough attribution is handled separately by
    `classify_termination_downstream`, not folded into this battery -- it
    answers a different question (what is the after-signal, not how
    confident is the drop itself)."""
    sig_p = evidence.get("dropoff_significance_p")
    level = None if sig_p is None else sig_p < thr.periodicity_significance_alpha

    cont_p = evidence.get("dropoff_continuity_p")
    continuity = None if cont_p is None else cont_p >= thr.periodicity_significance_alpha

    n_agree, n_total = _agree([level, continuity])
    return {
        "lenses": {"level": level, "continuity": continuity},
        "n_agree": n_agree,
        "n_total": n_total,
        "confidence": (n_agree / n_total) if n_total else None,
        "battery_passed": (n_total > 0) and (n_agree / n_total) >= thr.lens_agree_frac,
    }


def elong_battery(evidence: dict, thr: ScoreThresholds) -> dict:
    """Elongation: Level (frame_chisq_p significant), Uniformity
    (body_uniformity_p NOT significant -- consistent across both halves).
    Mappability is a separate downweight (see `compose_confidence`), not a
    battery lens: it never "agrees" or "disagrees" with a claim, it discounts
    trust in the whole measurement."""
    chisq_p = evidence.get("frame_chisq_p")
    level = None if chisq_p is None else chisq_p < thr.periodicity_significance_alpha

    uni_p = evidence.get("body_uniformity_p")
    uniformity = None if uni_p is None else uni_p >= thr.periodicity_significance_alpha

    n_agree, n_total = _agree([level, uniformity])
    return {
        "lenses": {"level": level, "uniformity": uniformity},
        "n_agree": n_agree,
        "n_total": n_total,
        "confidence": (n_agree / n_total) if n_total else None,
        "battery_passed": (n_total > 0) and (n_agree / n_total) >= thr.lens_agree_frac,
    }


def cif_battery(evidence: dict, thr: ScoreThresholds) -> dict:
    """CIF: Level (frac_significant over half), Uniformity (significant
    codons scattered, not clustered into one run), Cross-aspect agreement
    (CIF's own significance call matches the elongation chi-square call on
    the same event -- two independently-computed measures over the same
    underlying per-codon tallies agreeing is stronger evidence than either
    alone)."""
    sig = evidence.get("cif_significance") or {}
    frac = sig.get("frac_significant")
    level = None if frac is None else frac >= thr.cif_level_min_frac

    contig = evidence.get("cif_contiguity")
    uniformity = None if contig is None else contig["max_run_frac"] < thr.cif_contiguity_max_run_frac

    chisq_p = evidence.get("frame_chisq_p")
    cross_agree = None
    if level is not None and chisq_p is not None:
        elong_sig = chisq_p < thr.periodicity_significance_alpha
        cross_agree = level == elong_sig

    n_agree, n_total = _agree([level, uniformity, cross_agree])
    return {
        "lenses": {"level": level, "uniformity": uniformity, "cross_agree": cross_agree},
        "n_agree": n_agree,
        "n_total": n_total,
        "confidence": (n_agree / n_total) if n_total else None,
        "battery_passed": (n_total > 0) and (n_agree / n_total) >= thr.lens_agree_frac,
    }


def compose_confidence(
    aspect: str, evidence: dict, thr: ScoreThresholds
) -> Optional[float]:
    """`confidence` for §6: fraction of this aspect's lenses that agree, with
    elongation additionally downweighted (never gated) when
    `map_track_low` is set -- docs/significance_testing_plan.md §2 "Cross-
    check against existing region flags" / item 8. Multiplicative, not a
    veto: a low-mappability region with strong signal keeps most of its
    confidence rather than losing it outright.

    For "elongation", CIF's own battery (§3) is pooled in alongside
    elongation's: CIF has no separate row in `_RECORD_SCHEMA` (its evidence
    already rides on the elongation event, see `cif`/`n_codons`), so its
    lens agreement is folded into the one `confidence` value the elongation
    row actually gets. `elong_battery`/`cif_battery` stay separately
    callable for diagnostics (e.g. the results report breaks them out), but
    the schema-level number is their combined n_agree/n_total.
    """
    if aspect == "init":
        conf = init_battery(evidence, thr)["confidence"]
    elif aspect == "term":
        conf = term_battery(evidence, thr)["confidence"]
    elif aspect == "elongation":
        e, c = elong_battery(evidence, thr), cif_battery(evidence, thr)
        n_agree, n_total = e["n_agree"] + c["n_agree"], e["n_total"] + c["n_total"]
        conf = (n_agree / n_total) if n_total else None
    else:
        return None
    if conf is None:
        return None
    if aspect == "elongation" and evidence.get("map_track_low"):
        conf *= thr.mappability_confidence_penalty
    return conf


# ---------------------------------------------------------------------------
# Neighbor attribution (§0 "attribution not exclusion")
# ---------------------------------------------------------------------------


def classify_neighbor_outcome(
    own_passed: Optional[bool],
    neighbor_passed: Optional[bool],
    *,
    leakage: Optional[bool] = None,
) -> str:
    """The three outcomes from §0/§1, computed from batteries already run
    independently per frame -- never a head-to-head comparison of raw
    metrics (a head-to-head framing structurally can't return "both pass").

    Returns one of:
      "no_neighbor"        no competing frame at this locus.
      "only_this_frame"    this frame passes, neighbor doesn't -- ordinary
                            single-frame event; the neighboring signal was
                            noise or unrelated.
      "both_independent"   both frames independently pass their own battery
                            -- report both as genuine overlapping events.
      "leakage"             this frame's own battery fails, the neighbor's
                            passes, AND `leakage` is True (a caller-supplied
                            signal, e.g. elongation's `identifiability`
                            falling below `elong_identifiability` -- the
                            contended reads mostly belong to the neighbor's
                            frame) -- this frame's apparent signal is really
                            the neighbor's tail, not a real event here.
      "neither_supported"  neither frame's battery passes, or this frame
                            fails with no leakage signal to explain why --
                            attribution has nothing further to add.
    """
    if neighbor_passed is None:
        return "no_neighbor"
    if own_passed and neighbor_passed:
        return "both_independent"
    if own_passed and not neighbor_passed:
        return "only_this_frame"
    if not own_passed and neighbor_passed and leakage:
        return "leakage"
    return "neither_supported"


def classify_termination_downstream(
    evidence: dict,
    thr: ScoreThresholds,
    *,
    neighbor_passed: Optional[bool] = None,
) -> str:
    """The three outcomes from §4 "Downstream attribution": clean drop,
    genuine readthrough, or a distinct nearby event, using
    `dropoff_after_share` (magnitude of in-frame signal past the stop).

    `neighbor_passed`: the battery result (any aspect -- typically a
    downstream init event's own battery) of a distinct candidate event
    found overlapping the after-window, if the caller has one to offer
    (e.g. via `pair_boundary_neighbors` scoped to the downstream side).
    Without one, elevated after-signal is reported as "readthrough" rather
    than silently downgraded to "clean_drop" -- sustained in-frame signal
    past a stop is real, documented biology (§4), not noise to explain away
    by default.
    """
    share = evidence.get("dropoff_after_share")
    if share is None or share < thr.elong_in_frame:
        return "clean_drop"
    if neighbor_passed:
        return "distinct_downstream"
    return "readthrough"


def pair_elongation_neighbors(
    scored: Dict[int, dict],
    thr: ScoreThresholds,
) -> Dict[int, dict]:
    """Attribution outcome per elongation event against its competitors,
    using the identifiability/competitor_share fields `_elong_evidence`
    already computes from `overlaps_df` -- no new neighbor-finding needed,
    since overlapping elongation candidates are already separate,
    independently-scored events (`score_elongation_batch`) by construction.

    `evidence["competitor_share"]` is `{competitor_event_id: frame_share}`.
    Two-way overlaps get a real classification; three or more distinct
    competitor frames at one event get flagged `multi_way_overlap=True`
    with each event's own (already-independent) result reported as-is --
    deliberately NOT resolved, per docs/significance_testing_plan.md's own
    scope (no N>2 policy invented here; see the significance_testing_results
    report for whether this fired on real data).
    """
    out: Dict[int, dict] = {}
    for eid, ev in scored.items():
        battery = elong_battery(ev, thr)
        own_passed = battery["battery_passed"]
        competitors = ev.get("competitor_share") or {}
        identifiability = ev.get("identifiability")
        leakage = identifiability is not None and identifiability < thr.elong_identifiability
        if not competitors:
            out[eid] = {"outcome": "no_neighbor", "battery": battery}
            continue
        if len(competitors) > 1:
            out[eid] = {
                "outcome": "multi_way_overlap",
                "battery": battery,
                "n_competitors": len(competitors),
            }
            continue
        (comp_id,) = competitors.keys()
        comp_ev = scored.get(int(comp_id))
        comp_passed = elong_battery(comp_ev, thr)["battery_passed"] if comp_ev else None
        outcome = classify_neighbor_outcome(own_passed, comp_passed, leakage=leakage)
        out[eid] = {"outcome": outcome, "battery": battery, "competitor_event_id": int(comp_id)}
    return out


def pair_cif_neighbors(scored: Dict[int, dict], thr: ScoreThresholds) -> Dict[int, dict]:
    """Same shape as `pair_elongation_neighbors`, using `cif_battery` instead
    of `elong_battery` -- docs/significance_testing_plan.md §3's own
    attribution note: a stretch of codons significantly dominant in a
    DIFFERENT frame is a candidate for a real local overlapping element,
    checked on that frame's own per-codon battery, not treated as
    contamination of this ORF's CIF. Reuses the same overlaps/competitor
    data as elongation (CIF is computed over the same event span)."""
    out: Dict[int, dict] = {}
    for eid, ev in scored.items():
        battery = cif_battery(ev, thr)
        own_passed = battery["battery_passed"]
        competitors = ev.get("competitor_share") or {}
        if not competitors:
            out[eid] = {"outcome": "no_neighbor", "battery": battery}
            continue
        if len(competitors) > 1:
            out[eid] = {
                "outcome": "multi_way_overlap",
                "battery": battery,
                "n_competitors": len(competitors),
            }
            continue
        (comp_id,) = competitors.keys()
        comp_ev = scored.get(int(comp_id))
        comp_passed = cif_battery(comp_ev, thr)["battery_passed"] if comp_ev else None
        outcome = classify_neighbor_outcome(own_passed, comp_passed)
        out[eid] = {"outcome": outcome, "battery": battery, "competitor_event_id": int(comp_id)}
    return out


def pair_boundary_neighbors(
    events: List[dict],
    scored: Dict[int, dict],
    battery_fn,
    thr: ScoreThresholds,
    *,
    window_nt: int = 60,
) -> Dict[int, dict]:
    """Attribution outcome for init/term events: find other events of the
    SAME type on the same (chrom, strand) whose position falls within
    `window_nt` of this one in a DIFFERENT frame, and classify the pair.

    No upstream "candidate neighbor" concept exists for init/term the way
    overlaps_df already gives elongation (see events.py) -- this is a
    pragmatic proximity heuristic scoped to this session
    (docs/significance_testing_plan.md leaves "how neighbors are found" for
    init/term unspecified), using `window_nt` matching the largest step-score
    flank (60 nt) as the locality assumption. `events` rows need
    event_id/type/chrom/strand/start/phase; events missing `chrom` or
    `phase` are skipped (evidence-only degrade, not a crash).
    """
    from TranslonScorer.scoring.aspects import abs_frame

    by_key: Dict[Tuple[str, int, str], List[dict]] = {}
    for r in events:
        if r.get("chrom") is None or r.get("phase") is None or r["event_id"] not in scored:
            continue
        by_key.setdefault((r["chrom"], r["strand"], r["type"]), []).append(r)

    out: Dict[int, dict] = {}
    for key, rows in by_key.items():
        rows = sorted(rows, key=lambda r: r["start"])
        for i, r in enumerate(rows):
            eid = r["event_id"]
            ev = scored[eid]
            battery = battery_fn(ev, thr)
            own_passed = battery["battery_passed"]
            r_frame = abs_frame(r["phase"], r["strand"])
            neighbors = [
                o
                for o in rows
                if o["event_id"] != eid
                and abs(o["start"] - r["start"]) <= window_nt
                and abs_frame(o["phase"], o["strand"]) != r_frame
            ]
            if not neighbors:
                out[eid] = {"outcome": "no_neighbor", "battery": battery}
                continue
            if len(neighbors) > 1:
                out[eid] = {
                    "outcome": "multi_way_overlap",
                    "battery": battery,
                    "n_competitors": len(neighbors),
                }
                continue
            other = neighbors[0]
            other_ev = scored.get(other["event_id"])
            other_passed = battery_fn(other_ev, thr)["battery_passed"] if other_ev else None
            outcome = classify_neighbor_outcome(own_passed, other_passed)
            out[eid] = {
                "outcome": outcome,
                "battery": battery,
                "neighbor_event_id": other["event_id"],
            }
    return out

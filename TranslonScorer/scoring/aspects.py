"""Per-aspect event scorers: initiation, elongation, termination, junction.

Each scorer receives pre-extracted coverage and returns a raw evidence dict.
Decision logic (eligibility/call) is applied inside each scorer via helpers
in evidence.py so that the measurement and decision layers are clearly
separated: changing thresholds re-derives calls without re-measuring coverage.

Public API
----------
abs_frame               — phase × strand → genomic frame carrying in-frame signal
score_initiation_event  — P-site coverage step UP at start codon
score_termination_event — A-site coverage step DOWN at stop codon
score_elongation_event  — frame-resolved in-frame fraction + identifiability
score_junction_event    — splice-junction confident-spanning count + psi
"""

from __future__ import annotations

import math
import statistics
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as _np

from TranslonScorer.events import SpliceContext
from TranslonScorer.model import ScoreThresholds
from TranslonScorer.scoring.evidence import _decide_step, _elong_evidence
from TranslonScorer.scoring.signature import (
    DROPOFF_DOWNSTREAM_NT,
    DROPOFF_UPSTREAM_NT,
    DROPOFF_WINDOW_NT,
    cif,
    cif_codon_contiguity,
    cif_codon_significance,
    dropoff,
    gini,
    orf_signal_vector,
)

# ---------------------------------------------------------------------------
# Frame utility
# ---------------------------------------------------------------------------


def abs_frame(phase: int, strand: int) -> int:
    """Genomic frame (p % 3) carrying this event's in-frame (codon-0) signal."""
    return (-phase) % 3 if strand > 0 else phase % 3


# ---------------------------------------------------------------------------
# Initiation / termination: smoothed, flank-robust step scoring
# ---------------------------------------------------------------------------
# Each flank is binned into codons (3 nt, periodicity-aware) and the MEDIAN
# bin is taken as the robust level.  The step is measured at several flank
# lengths; the consensus (median of log2FCs) is the headline metric.


def _flank_positions(pos: int, strand: int, length: int, *, body_side: bool) -> List[int]:
    """Genomic positions for a flank of `length` nt on the body or outer side of
    `pos`, in transcript orientation. body_side=True → into the ORF.

    Flat genomic window — wrong whenever an intron falls inside the flank
    (the leader/UTR side is spliced). Used directly only when no splice
    context is available; see _project_flank for the splice-aware version.

    Minus-strand results are returned in DESCENDING genomic order, because
    that is transcript order there. This used to return ascending genomic on
    both strands, which made the fallback the exact reverse of the
    splice-aware path for every minus-strand flank. Nothing noticed: the only
    consumer was _codon_levels, which bins in 3s and takes median/max/total,
    and every flank length in use is a multiple of 3 — so the bin set, and
    hence every score, was identical either way. It became visible only with
    the first order-sensitive consumer (_dropoff_at, whose window is
    frame-locked to a specific index).
    """
    if strand > 0:
        rng = range(pos, pos + length) if body_side else range(pos - length, pos)
        return list(rng)
    rng = range(pos - length + 1, pos + 1) if body_side else range(pos + 1, pos + length + 1)
    return list(reversed(rng))


def _project_flank(
    pos: int,
    strand: int,
    length: int,
    *,
    body_side: bool,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
) -> Tuple[List[int], bool]:
    """`length` positions spanning a transcript-distance flank on the body or
    outer side of `pos`, jumping across any intron the flank would otherwise
    read straight through.

    Without this, a leader/UTR flank near a splice site reads raw genomic
    bases across the intron — mostly zero-coverage sequence a spliced read
    never touches — which spuriously inflates init_rise/term_drop (the
    leader/UTR looks emptier than it is). This walks in transcript order
    instead, using `splice_context[(chrom, strand)]` (the same genome-wide,
    strand-scoped intron list — sorted (donor, acceptor) genomic intervals —
    that `events._context_junction_events` already uses; it is not scoped to
    a single gene, so a nearby unrelated intron on the same strand can in
    principle be picked up, and where isoforms disagree on the immediate
    upstream/downstream exon this follows exactly one path — the design
    already defers isoform-of-origin attribution, see
    docs/reannotation_engine_design.md §6.3).

    Falls back to the flat genomic window (identical positions, in the same
    order) when no splice_context/chrom is given or no intron is nearby, so
    behaviour is unchanged unless a caller opts in.

    Returns (positions, spliced) — spliced=True iff at least one intron was
    crossed while building this flank.
    """
    if not splice_context or chrom is None:
        return list(_flank_positions(pos, strand, length, body_side=body_side)), False

    introns = splice_context.get((chrom, "+" if strand > 0 else "-"))
    if not introns:
        return list(_flank_positions(pos, strand, length, body_side=body_side)), False

    # donor = first intronic base, acceptor = first exonic base after the
    # intron (both ascending-genomic, per events.build_splice_context). Where
    # two introns share a donor (alternative acceptor usage) the later one
    # (by sorted order) wins — deterministic, not a resolved ambiguity.
    donor_to_acceptor = dict(introns)
    acceptor_to_donor = {a: d for d, a in introns}

    def step_up(p: int) -> int:
        a = donor_to_acceptor.get(p + 1)
        return a if a is not None else p + 1

    def step_down(p: int) -> int:
        d = acceptor_to_donor.get(p)
        return (d - 1) if d is not None else p - 1

    forward = step_up if strand > 0 else step_down  # transcript 5'->3'
    backward = step_down if strand > 0 else step_up

    spliced = False
    if body_side:
        out = [pos]
        p = pos
        for _ in range(length - 1):
            nxt = forward(p)
            spliced = spliced or abs(nxt - p) != 1
            p = nxt
            out.append(p)
        return out, spliced

    rev: List[int] = []
    p = pos
    for _ in range(length):
        nxt = backward(p)
        spliced = spliced or abs(nxt - p) != 1
        p = nxt
        rev.append(p)
    return list(reversed(rev)), spliced


def _codon_levels(
    coverage: Dict[int, float],
    positions: Sequence[int],
) -> Tuple[float, float, float]:
    """(median codon-bin sum, max codon-bin sum, total) over a window.

    Median is robust to a single spike; codon bins fold in 3-nt periodicity.
    """
    pos = list(positions)
    bins = [
        sum(coverage.get(p, 0.0) for p in pos[i : i + 3]) for i in range(0, max(len(pos) - 2, 0), 3)
    ]
    total = sum(coverage.get(p, 0.0) for p in pos)
    if not bins:
        return total, total, total
    return statistics.median(bins), max(bins), total


def _codon_bins(
    coverage: Dict[int, float],
    positions: Sequence[int],
) -> Tuple[_np.ndarray, _np.ndarray]:
    """(per-codon-bin total, per-codon-bin frame-0-only) over a
    transcript-ordered position list whose index 0 is a codon's first base.

    Every flank list this is called on is built by ``_project_flank``/
    ``_flank_positions`` anchored at ``pos`` — the start codon's first base
    for initiation, the terminal nucleotide for termination — both of which
    are codon boundaries by construction, and every flank length in use
    (9/18/30/60) is a multiple of 3. So index 0, 3, 6, ... is always the
    frame-0 position of its bin on BOTH sides: the body flank because it
    starts at ``pos`` itself, and the outer flank because it is contiguous
    ascending-transcript-order nucleotides immediately adjacent to ``pos``
    with a length divisible by 3 — no separate frame bookkeeping needed, and
    no genomic ``pos % 3`` math, matching how ``orf_signal_vector`` reindexes
    to transcript order before treating index % 3 as frame.

    Trailing 1-2 positions that don't complete a codon are dropped, matching
    ``_codon_levels``. Returns arrays of length ``len(positions) // 3``.
    """
    pos = list(positions)
    n_bins = len(pos) // 3
    total_bins = _np.zeros(n_bins, dtype=float)
    frame0_bins = _np.zeros(n_bins, dtype=float)
    for i in range(n_bins):
        triple = pos[i * 3 : i * 3 + 3]
        total_bins[i] = sum(coverage.get(p, 0.0) for p in triple)
        frame0_bins[i] = coverage.get(triple[0], 0.0)
    return total_bins, frame0_bins


def _peakiness(total_bins: _np.ndarray) -> float:
    """max/median codon-bin ratio -- the existing flank_peakiness formula,
    generalised to any bin array instead of one fixed reference flank."""
    if len(total_bins) == 0:
        return 0.0
    med = float(_np.median(total_bins))
    mx = float(_np.max(total_bins))
    return (mx / med) if med > 0 else (float("inf") if mx > 0 else 0.0)


def _breadth(total_bins: _np.ndarray) -> float:
    """Fraction of codon bins with any signal at all."""
    if len(total_bins) == 0:
        return 0.0
    return float(_np.count_nonzero(total_bins > 0) / len(total_bins))


def _periodicity(total_bins: _np.ndarray, frame0_bins: _np.ndarray) -> Optional[float]:
    """Frame-0 share of this window's total signal -- PIF's definition,
    applied to a boundary flank instead of a whole ORF span. None when the
    window has no signal at all (undefined, not 0)."""
    total = float(_np.sum(total_bins))
    if total <= 0:
        return None
    return float(_np.sum(frame0_bins) / total)


def _periodicity_significance(
    body_total: _np.ndarray,
    body_frame0: _np.ndarray,
    out_total: _np.ndarray,
    out_frame0: _np.ndarray,
    *,
    min_codons: int = 5,
) -> Optional[float]:
    """One-sided Mann-Whitney U p-value: is per-codon frame-0 share higher
    in the body flank's codons than the outer flank's?

    Inspired by RiboCode's use of a nonparametric periodicity test, not a
    reproduction of its method -- this compares per-codon frame-0 SHARE
    (frame0/total) between the two flanks directly, rather than testing one
    side against a fixed null, which is what "is there a change at the
    boundary" actually asks.

    Codons with zero total signal have an undefined share and are dropped
    from both samples. None (not a p-value, and not 0/1) when either side
    has fewer than `min_codons` codons left after dropping, or when scipy
    is unavailable -- a persuasive p-value needs enough codons to rank, and
    fabricating one from a handful of points, or from a hard scipy
    dependency this codebase otherwise keeps optional, is worse than an
    honest "don't know."
    """
    with _np.errstate(divide="ignore", invalid="ignore"):
        body_share = body_frame0 / body_total
        out_share = out_frame0 / out_total
    body_share = body_share[_np.isfinite(body_share)]
    out_share = out_share[_np.isfinite(out_share)]
    if len(body_share) < min_codons or len(out_share) < min_codons:
        return None
    try:
        from scipy.stats import mannwhitneyu
    except ImportError:
        return None
    try:
        _, p = mannwhitneyu(body_share, out_share, alternative="greater")
    except ValueError:
        return None
    if not math.isfinite(p):
        return None
    return float(p)


def _share_consistency(
    a_total: _np.ndarray,
    a_frame0: _np.ndarray,
    b_total: _np.ndarray,
    b_frame0: _np.ndarray,
    *,
    min_codons: int = 5,
) -> Optional[float]:
    """Two-sided sibling of `_periodicity_significance`: are two per-codon
    frame-0-share samples statistically indistinguishable?

    Same share computation and codon-floor/scipy-availability honesty
    convention as `_periodicity_significance` -- deliberately NOT
    implemented by calling it, because the alternative hypothesis differs
    (two-sided "are these different" vs one-sided "is A greater than B") and
    scipy's `mannwhitneyu` needs that as an argument, not a post-hoc
    reinterpretation of a one-sided p-value.

    A LARGE p-value is the "consistent" outcome here (fail to reject equal
    distributions) -- used where the question is uniformity/continuity, not
    a directional step (docs/significance_testing_plan.md §2, §4).
    """
    with _np.errstate(divide="ignore", invalid="ignore"):
        a_share = a_frame0 / a_total
        b_share = b_frame0 / b_total
    a_share = a_share[_np.isfinite(a_share)]
    b_share = b_share[_np.isfinite(b_share)]
    if len(a_share) < min_codons or len(b_share) < min_codons:
        return None
    try:
        from scipy.stats import mannwhitneyu
    except ImportError:
        return None
    try:
        _, p = mannwhitneyu(a_share, b_share, alternative="two-sided")
    except ValueError:
        return None
    if not math.isfinite(p):
        return None
    return float(p)


def _split_half_consistency(
    total: _np.ndarray, frame0: _np.ndarray, *, min_codons: int
) -> Optional[float]:
    """Split a codon-bin array in half and test frame-0-SHARE consistency
    between the two halves with `_share_consistency` -- the same spatial
    "does this hold up across the window, not just in aggregate" question
    `_elong_body_uniformity` asks of the elongation body and
    `_dropoff_continuity` asks of the termination before-window, applied
    here to a boundary flank's body side (docs/significance_testing_plan.md
    §1).

    Distinct from `gini_body_inframe` on the same flank: Gini measures
    magnitude concentration (is the mass piled into one bin), this measures
    spatial consistency (does the SHARE look the same in the near half vs
    the far half). A flank can score well on one and poorly on the other --
    e.g. broadly spread signal that is nonetheless much more in-frame near
    the boundary than further out would have unremarkable Gini but a small
    `_share_consistency` p-value.
    """
    n = len(total)
    if n < min_codons * 2:
        return None
    mid = n // 2
    return _share_consistency(
        total[:mid], frame0[:mid], total[mid:], frame0[mid:], min_codons=min_codons
    )


def _boundary_axes_for_flank(
    coverage: Dict[int, float],
    body_pos: Sequence[int],
    out_pos: Sequence[int],
    *,
    min_codons: int,
) -> dict:
    """All new (periodicity/uniformity/breadth/significance) axes for one
    flank length, both sides. Evidence only -- see _step_score."""
    body_total, body_frame0 = _codon_bins(coverage, body_pos)
    out_total, out_frame0 = _codon_bins(coverage, out_pos)
    peri_body = _periodicity(body_total, body_frame0)
    peri_outer = _periodicity(out_total, out_frame0)
    return {
        "periodicity_body": peri_body,
        "periodicity_outer": peri_outer,
        "periodicity_delta": (
            (peri_body - peri_outer) if peri_body is not None and peri_outer is not None else None
        ),
        "peakiness_body": _peakiness(body_total),
        "peakiness_outer": _peakiness(out_total),
        "gini_body": gini(body_total),
        "gini_outer": gini(out_total),
        # docs/significance_testing_plan.md §1 "Uniformity" lens -- Gini of
        # the FRAME-0-ONLY codon series, not the total-signal series above.
        # A genuine start settles into steady, roughly uniform in-frame
        # elongation downstream (low gini_body_inframe); a stray pileup
        # tends to stay peaky/patchy even if its total-signal Gini looks
        # unremarkable, because gini_body mixes in off-frame noise.
        "gini_body_inframe": gini(body_frame0),
        "gini_outer_inframe": gini(out_frame0),
        "breadth_body": _breadth(body_total),
        "breadth_outer": _breadth(out_total),
        "periodicity_p": _periodicity_significance(
            body_total, body_frame0, out_total, out_frame0, min_codons=min_codons
        ),
        # Second, distinct "uniformity" reading -- see _split_half_consistency
        # docstring for why this is not redundant with gini_body_inframe.
        "body_share_consistency_p": _split_half_consistency(
            body_total, body_frame0, min_codons=min_codons
        ),
    }


def _step_score(
    pos: int,
    strand: int,
    coverage: Dict[int, float],
    flanks: Sequence[int],
    min_reads: float,
    alpha: float,
    *,
    forward_is_body: bool = True,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
    min_codons: int = 5,
) -> dict:
    """Robust step (body_level - outer_level) across `pos`, over several flank
    lengths. Metric = log2 fold-change (body vs outer) with pseudocount α.

    `forward_is_body` picks which side of `pos` plays the "body" role:
    True (initiation) — downstream/into-ORF is body, upstream/leader is outer.
    False (termination) — upstream/into-ORF is body, downstream/UTR is outer.

    When `splice_context`/`chrom` are given, both flanks are built via
    _project_flank so a nearby intron is jumped rather than read straight
    through (see _project_flank docstring); otherwise this is unchanged flat
    genomic windowing.

    Alongside the depth-only rise (unchanged from before), each flank length
    also gets the frame-resolved axes from `_boundary_axes_for_flank` —
    periodicity, uniformity (peakiness + Gini), breadth, and a periodicity
    significance test, each computed both sides of `pos`. These are
    evidence, exactly like `stability`/`flank_peakiness` already are — they
    do not feed `consensus_rise` and do not gate eligibility/call here (see
    `_decide_step`).
    """
    rises: Dict[int, float] = {}
    flank_spliced = False
    axes_by_flank: Dict[int, dict] = {}
    for L in flanks:
        body_pos, sp1 = _project_flank(
            pos, strand, L, body_side=forward_is_body, chrom=chrom, splice_context=splice_context
        )
        out_pos, sp2 = _project_flank(
            pos,
            strand,
            L,
            body_side=not forward_is_body,
            chrom=chrom,
            splice_context=splice_context,
        )
        flank_spliced = flank_spliced or sp1 or sp2
        body_med, _, _ = _codon_levels(coverage, body_pos)
        out_med, _, _ = _codon_levels(coverage, out_pos)
        rises[L] = math.log2((body_med + alpha) / (out_med + alpha))
        axes_by_flank[L] = _boundary_axes_for_flank(
            coverage, body_pos, out_pos, min_codons=min_codons
        )
    vals = list(rises.values())
    consensus = statistics.median(vals)
    stability = (max(vals) - min(vals)) if len(vals) > 1 else 0.0
    refL = sorted(flanks)[len(flanks) // 2]
    out_pos_ref, _ = _project_flank(
        pos, strand, refL, body_side=not forward_is_body, chrom=chrom, splice_context=splice_context
    )
    out_med, out_max, _ = _codon_levels(coverage, out_pos_ref)
    flank_peakiness = (out_max / out_med) if out_med > 0 else (float("inf") if out_max > 0 else 0.0)
    body_pos_min, _ = _project_flank(
        pos,
        strand,
        min(flanks),
        body_side=forward_is_body,
        chrom=chrom,
        splice_context=splice_context,
    )
    n_reads = _codon_levels(coverage, body_pos_min)[2]
    return {
        "rise_by_flank": rises,
        "consensus_rise": consensus,
        "stability": stability,
        "flank_peakiness": flank_peakiness,
        "n_reads": n_reads,
        "flank_spliced": flank_spliced,
        "boundary_axes_by_flank": axes_by_flank,
    }


_DEFAULT_THR = ScoreThresholds()


def score_initiation_event(
    start_pos: int,
    strand: int,
    coverage: Dict[int, float],
    *,
    flanks: Sequence[int] = (9, 18, 30, 60),
    thr: ScoreThresholds = _DEFAULT_THR,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
) -> dict:
    """Initiation = smoothed, flank-robust step UP at the start codon.

    Headline metric: consensus log2 fold-change (body / outer flank).

    `chrom`/`splice_context`: when given, the leader flank is projected
    across a nearby intron instead of read as raw flanking genomic bases —
    see _project_flank. Omit for the legacy flat-genomic behaviour.
    """
    s = _step_score(
        start_pos,
        strand,
        coverage,
        flanks,
        thr.min_reads,
        thr.step_pseudocount,
        forward_is_body=True,
        chrom=chrom,
        splice_context=splice_context,
        min_codons=thr.periodicity_min_codons,
    )
    s["metric"] = s["consensus_rise"]
    s["metric_name"] = "init_rise"
    s["eligibility"], s["call"] = _decide_step(s, rise_thr=thr.init_rise, thr=thr)
    return s


def score_termination_event(
    term_pos: int,
    strand: int,
    coverage: Dict[int, float],
    *,
    flanks: Sequence[int] = (9, 18, 30, 60),
    thr: ScoreThresholds = _DEFAULT_THR,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
) -> dict:
    """Termination = step DOWN after the stop codon (body before > UTR after).

    Headline metric: consensus log2 fold-change (body / UTR) — positive = drop.

    `chrom`/`splice_context`: when given, the UTR flank is projected across a
    nearby intron instead of read as raw flanking genomic bases — see
    _project_flank. Omit for the legacy flat-genomic behaviour.
    """
    s = _step_score(
        term_pos,
        strand,
        coverage,
        flanks,
        thr.min_reads,
        thr.step_pseudocount,
        forward_is_body=False,
        chrom=chrom,
        splice_context=splice_context,
        min_codons=thr.periodicity_min_codons,
    )
    s["metric"] = s["consensus_rise"]
    s["metric_name"] = "term_drop"
    s["eligibility"], s["call"] = _decide_step(s, rise_thr=thr.term_drop, thr=thr)
    s["dropoff"] = _dropoff_at(
        term_pos, strand, coverage, chrom=chrom, splice_context=splice_context
    )
    s["dropoff_significance_p"] = _dropoff_significance(
        term_pos,
        strand,
        coverage,
        chrom=chrom,
        splice_context=splice_context,
        min_codons=thr.periodicity_min_codons,
    )
    s["dropoff_continuity_p"] = _dropoff_continuity(
        term_pos,
        strand,
        coverage,
        chrom=chrom,
        splice_context=splice_context,
        min_codons=thr.periodicity_min_codons,
    )
    s["dropoff_after_share"] = _dropoff_after_share(
        term_pos, strand, coverage, chrom=chrom, splice_context=splice_context
    )
    return s


def _dropoff_at(
    term_pos: int,
    strand: int,
    coverage: Dict[int, float],
    *,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
) -> Optional[float]:
    """Chothani et al. (2022) ribosome drop-off across the stop, as evidence.

    Rides on the termination aspect rather than elongation because it needs
    positions PAST the ORF end, and this is where the splice-aware flank
    machinery already lives — a flat genomic window would read straight
    through an intron for any stop near a 3' exon boundary.

    Distinct from the headline ``term_drop``, which is a multi-flank consensus
    log2 fold-change over all positions with flanks out to 60 nt. This is a
    bounded [0,1] ratio over a fixed 33-nt window using frame-0 positions only.
    Both answer "does signal fall after the stop?"; they are not
    interconvertible and neither supersedes the other.

    ASSUMPTION, and it is load-bearing: ``term_pos`` is the translon's terminal
    nucleotide in transcript orientation (events.py: ``bed_end - 1`` on +,
    ``bed_start`` on -). The reference's window is frame-locked so that index
    17 is the terminal STOP nucleotide, so the two coincide only if the
    translon span includes its stop codon. If an upstream annotation excludes
    it, every drop-off here is 3 nt out of register — so check it per source.
    Verified on translon_db/translons.sqlite (8,852,481 rows): 99.9%
    terminal_codon_class='stop', 100.0% length_mod3=0.

    None when the window runs off the contig — the reference flags that case
    rather than truncating, and a truncated window would silently change which
    positions are in frame.
    """
    before, _ = _project_flank(
        term_pos,
        strand,
        DROPOFF_UPSTREAM_NT,
        body_side=False,
        chrom=chrom,
        splice_context=splice_context,
    )
    after, _ = _project_flank(
        term_pos,
        strand,
        DROPOFF_DOWNSTREAM_NT + 1,  # includes term_pos itself at index 0
        body_side=True,
        chrom=chrom,
        splice_context=splice_context,
    )
    window = before + after
    if len(window) != DROPOFF_WINDOW_NT or min(window) < 0:
        return None
    return dropoff(_np.array([coverage.get(p, 0.0) for p in window], dtype=float))


def _dropoff_window_positions(
    term_pos: int,
    strand: int,
    *,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
) -> Optional[List[int]]:
    """The same 33-position, frame-locked window `_dropoff_at` builds, as
    positions rather than a coverage-sampled array -- shared by
    `_dropoff_significance`/`_dropoff_continuity` so both stay pinned to the
    reference's exact frame alignment (index 17 == the terminal stop
    nucleotide) instead of re-deriving it from a truncated flank, which is
    not itself codon-boundary-aligned (17 is not a multiple of 3)."""
    before, _ = _project_flank(
        term_pos,
        strand,
        DROPOFF_UPSTREAM_NT,
        body_side=False,
        chrom=chrom,
        splice_context=splice_context,
    )
    after, _ = _project_flank(
        term_pos,
        strand,
        DROPOFF_DOWNSTREAM_NT + 1,
        body_side=True,
        chrom=chrom,
        splice_context=splice_context,
    )
    window = before + after
    if len(window) != DROPOFF_WINDOW_NT or min(window) < 0:
        return None
    return window


def _dropoff_significance(
    term_pos: int,
    strand: int,
    coverage: Dict[int, float],
    *,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
    min_codons: int = 5,
) -> Optional[float]:
    """One-sided significance mirror of `dropoff()`: is per-codon frame-0 share
    higher in the before-stop window than the after-stop window?

    Reuses `_periodicity_significance` itself (docs/significance_testing_plan.md
    §4, "Level" lens) rather than reimplementing the comparison, over the same
    frame-locked window `dropoff()` uses: 6 codons up to and including the
    stop vs. the 5 codons after it. `min_codons` applies per side, so with the
    default of 5 this is already close to running at the floor -- expect
    `None` more often than the boundary-flank version of this test, which has
    more codons to work with at larger flank lengths.
    """
    window = _dropoff_window_positions(term_pos, strand, chrom=chrom, splice_context=splice_context)
    if window is None:
        return None
    before_pos, after_pos = window[:18], window[18:]
    before_total, before_frame0 = _codon_bins(coverage, before_pos)
    after_total, after_frame0 = _codon_bins(coverage, after_pos)
    return _periodicity_significance(
        before_total, before_frame0, after_total, after_frame0, min_codons=min_codons
    )


def _dropoff_continuity(
    term_pos: int,
    strand: int,
    coverage: Dict[int, float],
    *,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
    min_codons: int = 5,
    extend_nt: int = 18,
) -> Optional[float]:
    """Does the before-stop window look like ordinary upstream elongation, or
    an isolated pileup right at the stop (docs/significance_testing_plan.md
    §4, "Continuity" lens)?

    Builds a window immediately upstream of the drop-off's 6-codon before-stop
    window (`extend_nt` nt, default 18 == 6 codons, so it is itself a whole
    number of codons away and stays frame-aligned) and compares per-codon
    frame-0 share between the two with `_periodicity_significance`, "near"
    (closer to the stop) vs "far" (further upstream). A small p-value means
    the near window is significantly MORE frame-0-dominant than ordinary
    upstream elongation looks -- consistent with a stop-proximal pileup
    rather than a genuine continuation of the ORF's own signature. A large
    p-value is the "looks like normal elongation" outcome.
    """
    window = _dropoff_window_positions(term_pos, strand, chrom=chrom, splice_context=splice_context)
    if window is None:
        return None
    near_pos = window[:18]
    far_pos, _ = _project_flank(
        near_pos[0], strand, extend_nt, body_side=False, chrom=chrom, splice_context=splice_context
    )
    if len(far_pos) != extend_nt:
        return None
    near_total, near_frame0 = _codon_bins(coverage, near_pos)
    far_total, far_frame0 = _codon_bins(coverage, far_pos)
    return _periodicity_significance(
        near_total, near_frame0, far_total, far_frame0, min_codons=min_codons
    )


def _dropoff_after_share(
    term_pos: int,
    strand: int,
    coverage: Dict[int, float],
    *,
    chrom: Optional[str] = None,
    splice_context: Optional[SpliceContext] = None,
) -> Optional[float]:
    """In-frame share of the 5-codon after-stop window (docs/
    significance_testing_plan.md §4 "Downstream attribution"). Feeds
    `attribution.classify_termination_downstream`'s clean-drop /
    readthrough / distinct-downstream call: unlike `dropoff_significance_p`
    (which only says whether before beats after), this is a magnitude
    reading of the after side alone, needed because "no signal after the
    stop" and "clear signal after the stop that just isn't statistically
    beaten by before" are different situations for that classification.
    """
    window = _dropoff_window_positions(term_pos, strand, chrom=chrom, splice_context=splice_context)
    if window is None:
        return None
    after_total, after_frame0 = _codon_bins(coverage, window[18:])
    return _periodicity(after_total, after_frame0)


# ---------------------------------------------------------------------------
# Elongation
# ---------------------------------------------------------------------------


def _elong_frame_chisq(full_by_frame: Sequence[float]) -> Optional[float]:
    """Chi-square goodness-of-fit of the 3-way frame tally against uniform
    1/3 (docs/significance_testing_plan.md §2, "Level" lens).

    `full_by_frame` is every covered position's read count binned by
    genomic ``p % 3`` across the whole span (contended and clean together --
    unlike `cont_by_frame`, which is contended-only). Preferred over a plain
    frame-0-vs-rest binomial test because it also catches signal split
    unevenly between the two non-frame-0 positions, which a binomial test
    against frame-0 alone cannot distinguish from a real uniform null.

    None when there is no signal, when any expected cell would be zero
    (chisquare's degenerate case), or when scipy is unavailable -- same
    honesty convention as `_periodicity_significance`.
    """
    total = float(sum(full_by_frame))
    if total <= 0:
        return None
    try:
        from scipy.stats import chisquare
    except ImportError:
        return None
    try:
        _, p = chisquare(list(full_by_frame))
    except ValueError:
        return None
    if not math.isfinite(p):
        return None
    return float(p)


def _elong_body_uniformity(vec: _np.ndarray, *, min_codons: int = 5) -> Optional[float]:
    """Split the ORF body into two halves and test whether frame-0
    dominance is consistent along its length (docs/significance_testing_plan.md
    §2, "Uniformity" lens).

    `vec` is the same per-nucleotide, frame-0-at-index-0 signal vector CIF is
    computed from (`orf_signal_vector`) -- reused rather than re-summed. A
    confound that only overlaps part of the ORF shows up as a local patch
    (small p, halves differ); genuine elongation should look similar
    throughout (large p). None when either half has fewer than `min_codons`
    codons after dropping zero-signal ones, or on an odd number of trailing
    nt (dropped, matching `_codon_bins`).

    Thin wrapper over `_split_half_consistency` (the same split+test
    initiation's `body_share_consistency_p` uses) after binning `vec` into
    per-codon total/frame0 arrays -- the two were computing the same thing
    from different starting representations (raw per-nt vector here, an
    already-binned array there).
    """
    n = len(vec) - (len(vec) % 3)
    if n <= 0:
        return None
    codons = vec[:n].reshape(-1, 3)
    total, frame0 = codons.sum(axis=1), codons[:, 0]
    return _split_half_consistency(total, frame0, min_codons=min_codons)


def score_elongation_event(
    start: int,
    end: int,
    phase: int,
    strand: int,
    coverage: Dict[int, float],
    overlaps: Optional[List[Tuple[int, int, int, int]]] = None,
    *,
    thr: ScoreThresholds = _DEFAULT_THR,
) -> dict:
    """Score one elongation event.

    overlaps: [(overlap_start, overlap_end, competitor_phase, competitor_event_id)]
    for every competing event sharing genomic span in a different frame.

    Returns clean vs contended in-frame fractions, per-overlap identifiability,
    and per-competitor signal shares via _elong_evidence.
    """
    overlaps = overlaps or []
    a_e = abs_frame(phase, strand)

    contended: Dict[int, None] = {}
    comp_frame: Dict[int, int] = {}
    for os, oe, cph, cid in overlaps:
        for p in range(max(os, start), min(oe, end)):
            contended[p] = None
        comp_frame[cid] = abs_frame(cph, strand)

    n = clean_tot = clean_inf = cont_tot = 0.0
    cont_by_frame = [0.0, 0.0, 0.0]
    full_by_frame = [0.0, 0.0, 0.0]
    covered = 0
    for p in range(start, end):
        c = coverage.get(p, 0.0)
        if c <= 0:
            continue
        covered += 1
        n += c
        full_by_frame[p % 3] += c
        if p in contended:
            cont_tot += c
            cont_by_frame[p % 3] += c
        else:
            clean_tot += c
            if p % 3 == a_e:
                clean_inf += c

    # CIF needs the per-codon vector, which the frame sums above cannot give.
    # Built from the same coverage map, in transcript order.
    _cov_pos = _np.array([p for p in range(start, end) if coverage.get(p)], dtype=_np.int64)
    _cov_cnt = _np.array([float(coverage[int(p)]) for p in _cov_pos], dtype=float)
    _vec = orf_signal_vector(_cov_pos, _cov_cnt, start, end, a_e, strand)
    _cif_sig = cif_codon_significance(
        _vec, min_reads=thr.cif_codon_min_reads, alpha=thr.periodicity_significance_alpha
    )
    _cif_contig = cif_codon_contiguity(_cif_sig["sig_mask"]) if _cif_sig else None

    return _elong_evidence(
        n,
        covered,
        clean_tot,
        clean_inf,
        cont_tot,
        cont_by_frame,
        a_e,
        comp_frame,
        len(contended),
        end - start,
        thr,
        cif_value=cif(_vec),
        n_codons=(len(_vec) // 3) if len(_vec) else None,
        frame_chisq_p=_elong_frame_chisq(full_by_frame),
        cif_significance=_cif_sig,
        cif_contiguity=_cif_contig,
        body_uniformity_p=_elong_body_uniformity(_vec, min_codons=thr.elong_uniformity_min_codons),
    )


# ---------------------------------------------------------------------------
# Junction
# ---------------------------------------------------------------------------


def score_junction_event(
    support: dict,
    *,
    thr: ScoreThresholds = _DEFAULT_THR,
) -> dict:
    """Splice support with an overhang-confidence (identifiability) layer.

    support keys: span_conf, span_short, unspliced (weighted read counts).
      span_conf  — intron matches with ≥6 nt aligned both sides (confident)
      span_short — matches with short overhang (~ambiguous)
      unspliced  — aligned straight through the donor (intron retention)

    Headline metric: confident spanning count. psi = spliced-in ratio.
    """
    conf = float(support.get("span_conf", 0.0))
    short = float(support.get("span_short", 0.0))
    unspl = float(support.get("unspliced", 0.0))
    spanning = conf + short
    psi = (spanning / (spanning + unspl)) if (spanning + unspl) > 0 else None
    if conf >= thr.junc_min_spanning:
        elig, call = "ELIGIBLE", "SUPPORTED"
    elif spanning >= thr.junc_min_spanning:
        elig, call = "ELIGIBLE", "AMBIGUOUS"
    else:
        elig, call = "INSUFFICIENT", None
    return {
        "n_reads": spanning,
        "metric": conf,
        "metric_name": "junc_confident_spanning",
        "confident_spanning": conf,
        "short_spanning": short,
        "unspliced": unspl,
        "psi": psi,
        "eligibility": elig,
        "call": call,
    }

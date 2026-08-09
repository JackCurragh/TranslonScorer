"""Translation signature scores (Chothani et al., Molecular Cell 2022).

Pure transforms over an ORF's per-nucleotide P-site signal vector.  No I/O, no
coverage access, no thresholds applied — the caller supplies the vector and gets
numbers back, so these can be tested directly against the R reference in
``Sonia_translation signature scores.txt``.

The three scores:

  PIF      frame-0 share of total signal
  CIF      fraction of codons whose first position carries >33.33% of the codon
  dropoff  in-frame signal before the stop, over before + after

Faithful-to-R details that are easy to get wrong, each covered by a test:

* **CIF's denominator is every codon**, including codons with no signal at all.
  In R, ``0/0`` is ``NaN`` and ``which(NaN > 33.33)`` drops it from the
  numerator, while the denominator stays ``nrow(psites)/3``.  Uncovered codons
  therefore *penalise* CIF rather than being excluded from it.

* **The CIF threshold is the literal 33.33, not 100/3.**  A codon with exactly
  equal signal in all three frames scores 33.333… which is ``> 33.33``, so
  perfectly uniform codons count as frame-0 dominant.  Almost certainly not the
  original intent, but it is what the reference does and what published numbers
  reflect, so it is preserved deliberately.

* **Drop-off window indices are frame-locked to the stop codon.**  In the
  33-position window, index 17 (0-based) is the terminal stop nucleotide and
  index 15 is the stop codon's first base, which is why the numerator runs
  0,3,…,15 and the denominator continues 18,21,…,30.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

# Sonia's literal threshold. NOT 100/3 — see the module docstring; the
# difference decides how perfectly uniform codons are counted.
CIF_FRAME_DOMINANCE_PCT = 33.33

# Drop-off window: 17 nt upstream of the terminal stop nucleotide, the stop
# nucleotide itself, and 15 nt downstream.
DROPOFF_UPSTREAM_NT = 17
DROPOFF_DOWNSTREAM_NT = 15
DROPOFF_WINDOW_NT = DROPOFF_UPSTREAM_NT + 1 + DROPOFF_DOWNSTREAM_NT  # 33

# Frame-0 offsets within that window, before and after the stop codon.
_DROPOFF_BEFORE = tuple(range(0, 18, 3))  # 0,3,6,9,12,15 — last 6 codons incl. stop
_DROPOFF_AFTER = tuple(range(18, 33, 3))  # 18,21,24,27,30 — 5 codons past the stop


def pif(signal: np.ndarray) -> Optional[float]:
    """Periodicity In Frame: frame-0 share of the ORF's total signal.

    ``sum(x[0::3]) / sum(x)``.  None when there is no signal at all, which is
    a different statement from 0.0 (no *in-frame* signal).
    """
    total = float(np.sum(signal))
    if total <= 0:
        return None
    return float(np.sum(signal[::3]) / total)


def cif(signal: np.ndarray) -> Optional[float]:
    """Codons In Frame: fraction of codons whose first position dominates.

    Denominator is *every* codon in the ORF, including zero-signal ones — see
    the module docstring.  None when the length is not a whole number of
    codons, since the reference silently recycles vectors in that case rather
    than erroring, and a wrong number is worse than a missing one.
    """
    n = len(signal)
    if n == 0 or n % 3:
        return None
    codons = np.asarray(signal, dtype=float).reshape(-1, 3)
    codon_total = codons.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        share = codons[:, 0] * 100.0 / codon_total
    # NaN (an all-zero codon) must fail the comparison, matching R's which().
    dominant = np.nan_to_num(share, nan=-np.inf, posinf=np.inf, neginf=-np.inf)
    return float(np.count_nonzero(dominant > CIF_FRAME_DOMINANCE_PCT) / codons.shape[0])


def dropoff(window: np.ndarray) -> Optional[float]:
    """Ribosome drop-off across the stop codon, from a 33-nt window.

    ``window`` must be the signal over ``DROPOFF_WINDOW_NT`` transcript
    positions positioned so index 17 is the terminal stop nucleotide.  Returns
    in-frame signal up to and including the stop codon, divided by that plus
    the in-frame signal of the 5 codons past it.  None when the denominator is
    zero — no signal either side is not a drop-off of 0.
    """
    if len(window) != DROPOFF_WINDOW_NT:
        raise ValueError(f"drop-off window must be {DROPOFF_WINDOW_NT} nt, got {len(window)}")
    w = np.asarray(window, dtype=float)
    before = float(w[list(_DROPOFF_BEFORE)].sum())
    after = float(w[list(_DROPOFF_AFTER)].sum())
    denom = before + after
    if denom <= 0:
        return None
    return before / denom


def dropoff_window_indices(stop_index: int, n_transcript_positions: int) -> Optional[range]:
    """Transcript-coordinate slice for the drop-off window, or None if it does not fit.

    ``stop_index`` is the 0-based index of the terminal stop nucleotide within
    the strand-ordered transcript position vector.  Needs 17 positions before
    and 15 after; the reference flags both cases rather than truncating.
    """
    start = stop_index - DROPOFF_UPSTREAM_NT
    end = stop_index + DROPOFF_DOWNSTREAM_NT + 1
    if start < 0 or end > n_transcript_positions:
        return None
    return range(start, end)


def orf_signal_vector(
    pos_sorted: np.ndarray,
    cnt_sorted: np.ndarray,
    start: int,
    end: int,
    expected_frame: int,
    strand: int,
) -> np.ndarray:
    """Per-nucleotide signal over ``[start, end)`` in TRANSCRIPT order.

    Trimmed at the 5' end so index 0 is a codon's first base, and at the 3' end
    so the length is a whole number of codons.  That is what lets ``pif`` and
    ``cif`` stay strand-agnostic: they can assume ``0, 3, 6, …`` is frame 0.

    The strand handling is the part that is easy to get wrong.  A codon is
    three consecutive *transcript* positions, so on the minus strand its first
    base sits at the HIGHEST genomic coordinate of the three.  Grouping raw
    genomic triples the same way on both strands would put the codon boundary
    one or two bases off for every minus-strand ORF, and CIF would be measuring
    the wrong position's share.  Reversing first sidesteps that entirely.

    ``expected_frame`` is the genomic ``pos % 3`` residue of frame-0 bases —
    i.e. ``aspects.abs_frame(phase, strand)``.

    Coverage is looked up from the sorted sparse arrays the batch scorer
    already holds; absent positions are 0.0.
    """
    n = end - start
    if n <= 0:
        return np.zeros(0, dtype=float)

    lo = int(np.searchsorted(pos_sorted, start, side="left"))
    hi = int(np.searchsorted(pos_sorted, end, side="left"))
    v = np.zeros(n, dtype=float)
    if hi > lo:
        v[pos_sorted[lo:hi] - start] = cnt_sorted[lo:hi]

    if strand < 0:
        v = v[::-1]
        # v[i] is genomic position end-1-i; frame 0 where (end-1-i) % 3 == e.
        offset = (end - 1 - expected_frame) % 3
    else:
        # v[i] is genomic position start+i; frame 0 where (start+i) % 3 == e.
        offset = (expected_frame - start) % 3

    v = v[offset:]
    return v[: len(v) - (len(v) % 3)]


def signal_from_intervals(
    positions: Sequence[int],
    starts: np.ndarray,
    ends: np.ndarray,
    values: np.ndarray,
    *,
    strict_r_compat: bool = False,
) -> np.ndarray:
    """Look up per-position signal from half-open intervals (``start < p <= end``).

    Uncovered positions are 0.0, which is what a sparse bedGraph means.

    ``strict_r_compat=True`` instead *drops* uncovered positions, reproducing
    the reference's ``unlist(lapply(...))``: R returns ``numeric(0)`` for a
    position with no covering interval and ``unlist`` silently discards it, so
    the vector comes back short.  The returned length is the caller's signal
    that it fired.

    What that costs depends on the gap, and the distinction matters:

    * **Gap length divisible by 3** — frame is preserved, so PIF is unaffected.
      CIF still changes, because the dropped codons leave the denominator
      instead of penalising the score (see ``cif``), which inflates it.
    * **Gap length not divisible by 3** — every downstream position changes
      frame, ``seq(1,n,3)`` no longer selects frame 0, and both PIF and CIF
      become meaningless rather than merely biased.

    So this is a bug, but not a uniformly catastrophic one, and on a fully
    dense track it never fires at all — which is the likeliest reason it went
    unnoticed.  Kept reproducible because it is what produced the published
    numbers, so the gap between the two readings can be quantified per dataset
    rather than argued about.
    """
    pos = np.asarray(positions, dtype=np.int64)
    out = np.zeros(len(pos), dtype=float)
    covered = np.zeros(len(pos), dtype=bool)
    order = np.argsort(starts, kind="stable")
    s, e, v = starts[order], ends[order], values[order]
    # Half-open (start, end]: the last interval starting strictly before p.
    idx = np.searchsorted(s, pos, side="left") - 1
    valid = idx >= 0
    if np.any(valid):
        hit = valid.copy()
        hit[valid] = pos[valid] <= e[idx[valid]]
        out[hit] = v[idx[hit]]
        covered[hit] = True
    if strict_r_compat:
        return out[covered]
    return out

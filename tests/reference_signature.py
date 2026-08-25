"""Literal transliteration of the R signature scores, as an independent oracle.

Deliberately NOT refactored: this mirrors ``Sonia_translation signature
scores.txt`` line by line, including its indexing and its treatment of NaN, so
that ``TranslonScorer.scoring.signature`` has something to disagree with.  The
vectorised implementation is not obviously correct by inspection; a second,
dumber implementation is worth more than the lines it costs.

Must never be imported by product code — same rule as tests/reference_scorer.py.

R is 1-based and its ``seq(a, b, by)`` is inclusive of both ends.  Every index
expression below carries the R it came from.
"""

from __future__ import annotations

from typing import List, Optional, Sequence


def _r_seq(start: int, stop: int, by: int) -> List[int]:
    """R's seq(start, stop, by) -> 0-based Python indices."""
    out, i = [], start
    while i <= stop:
        out.append(i - 1)  # 1-based -> 0-based
        i += by
    return out


def r_pif(psites: Sequence[float]) -> Optional[float]:
    """pif = sum(psites[seq(1,nrow(psites),3),2]) / sum(psites[,2])"""
    n = len(psites)
    total = sum(psites)
    if total == 0:
        return None  # R yields NaN here; we surface it as missing
    return sum(psites[i] for i in _r_seq(1, n, 3)) / total


def r_cif(psites: Sequence[float]) -> Optional[float]:
    """tot = p[seq(2,n,3)] + p[seq(1,n,3)] + p[seq(3,n,3)]
    cif = length(which(p[seq(1,n,3)]*100/tot > 33.33)) / (nrow(psites)/3)
    """
    n = len(psites)
    if n == 0 or n % 3:
        return None  # R would recycle mismatched seq() lengths; refuse instead
    i1 = _r_seq(1, n, 3)
    i2 = _r_seq(2, n, 3)
    i3 = _r_seq(3, n, 3)
    tot = [psites[a] + psites[b] + psites[c] for a, b, c in zip(i1, i2, i3)]
    hits = 0
    for a, t in zip(i1, tot):
        if t == 0:
            continue  # 0/0 -> NaN, and which() drops NaN
        if psites[a] * 100 / t > 33.33:
            hits += 1
    return hits / (n / 3)


def r_dropoff(window: Sequence[float]) -> Optional[float]:
    """dropoff = sum(w[seq(1,18,3)]) / sum(w[c(seq(1,18,3), seq(19,33,3))])"""
    num_idx = _r_seq(1, 18, 3)
    den_idx = _r_seq(1, 18, 3) + _r_seq(19, 33, 3)
    denom = sum(window[i] for i in den_idx)
    if denom == 0:
        return None
    return sum(window[i] for i in num_idx) / denom

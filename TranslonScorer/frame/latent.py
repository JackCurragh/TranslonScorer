from __future__ import annotations

"""
Probabilistic latent frame model for Ribo-seq.

Method 3 from the frame-assignment framework:
  "Model observed signal as a mixture of latent processes:
     O = M T + B
   where T = true frame signal, M = frame leakage model, B = background."

Each codon's observed frame-count vector O^(c) = [o0, o1, o2] arises from:

    E[o_i^(c)] = sum_j  M[i,j] * t_j^(c)  +  b_i

where:
  - M[i,j] = P(observe frame i | true frame j)   (shared across codons)
  - T^(c)  = [t0, t1, t2]                         (latent per-codon support)
  - B      = [b0, b1, b2]                          (global background rate)

EM algorithm
------------
Starting from init_M we alternate over expected true-frame counts T:

E-step (Richardson-Lucy update for T, per codon):
    Lambda^(c)   = M @ T^(c) + B
    T_new^(c)[j] = T^(c)[j] * sum_i  M[i,j] * O^(c)[i] / Lambda^(c)[i]
    In matrix form: T_new = T * (ratio @ M)   where ratio = O / Lambda

Optional M-step (update M from all codons):
    M_new[i,j] proportional to  M[i,j] * sum_c  ratio[c,i] * T[c,j]
             = M * (ratio.T @ T)
    Columns of M_new normalised to sum 1.

By default the leakage matrix is held fixed. That makes this model a
count-weighted probabilistic refinement of the CDS-learned linear correction.
Updating M from all codons is available, but should be regularised toward
init_M; otherwise sparse or non-canonical regions can move the leakage model
away from the calibration evidence.

After convergence T is row-normalised to posterior probabilities.
"""

from typing import Tuple
import numpy as np

_EPS = 1e-12


def _normalise_columns(M: np.ndarray) -> np.ndarray:
    out = M.astype(np.float64).copy()
    col_s = out.sum(axis=0, keepdims=True)
    col_s[col_s < _EPS] = 1.0
    return out / col_s


def _initialise_true_counts(O: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Initialise T by regularised pseudo-inverse and preserve row totals."""
    try:
        M_pinv = np.linalg.pinv(M + 1e-3 * np.eye(3))
        T = O @ M_pinv.T
    except np.linalg.LinAlgError:
        T = O.copy()

    T[T < 0] = _EPS
    observed_total = O.sum(axis=1, keepdims=True)
    corrected_total = T.sum(axis=1, keepdims=True)
    nonzero = corrected_total.squeeze() > _EPS
    if np.any(nonzero):
        T[nonzero] *= observed_total[nonzero] / corrected_total[nonzero]
    return T


def fit_latent(
    O: np.ndarray,
    init_M: np.ndarray,
    background: str = "flat",
    max_iter: int = 50,
    tol: float = 1e-5,
    background_frac: float = 0.005,
    update_M: bool = False,
    m_prior_strength: float = 1000.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit latent frame model via EM (Poisson mixture).

    Parameters
    ----------
    O             : (N, 3) observed frame-count vectors per codon.
    init_M        : (3, 3) initial confusion matrix (rows=observed, cols=true).
    background    : 'flat' (small uniform B) or 'zero' (no background).
    max_iter      : maximum EM iterations.
    tol           : relative log-likelihood change for early stopping.
    background_frac: fraction of mean signal used as background rate.
    update_M      : if true, re-estimate the leakage matrix from all codons.
                    Default false keeps the CDS-calibrated leakage fixed.
    m_prior_strength: pseudo-count strength anchoring M to init_M when update_M
                    is true.

    Returns
    -------
    P     : (N, 3) normalised latent frame posteriors.
    M_fit : (3, 3) fitted confusion matrix.
    B     : (3,)   fitted background rates.
    """
    if O.size == 0:
        return O.copy(), init_M.copy(), np.zeros(3, dtype=np.float64)

    O_counts = O.astype(np.float64).copy()
    O_counts[O_counts < 0] = 0.0

    # Initialise M (column-normalised so each col sums to 1)
    M_prior = _normalise_columns(init_M)
    M = M_prior.copy()

    # Background
    if background == "zero":
        B = np.zeros(3, dtype=np.float64)
    else:
        B = np.full(3, background_frac * float(O_counts.mean()), dtype=np.float64)

    # T is in count space, not per-codon proportions. This keeps high-depth
    # codons more informative than sparse/noisy codons.
    T = _initialise_true_counts(O_counts, M)

    prev_ll = -np.inf

    for iteration in range(max_iter):
        # E-step: Lambda[c,i] = sum_j M[i,j]*T[c,j] + B[i]
        Lambda = T @ M.T + B[np.newaxis, :]   # (N, 3)
        Lambda[Lambda < _EPS] = _EPS

        ratio  = O_counts / Lambda             # (N, 3)

        # T_new[c,j] = T[c,j] * sum_i M[i,j] * ratio[c,i]
        #            = (T * (ratio @ M))
        T_new  = T * (ratio @ M)               # (N, 3)
        T_new[T_new < 0] = _EPS

        if update_M:
            # M-step: expected allocations, regularised by CDS-learned M.
            expected = M * (ratio.T @ T)
            if m_prior_strength > 0:
                expected = expected + (float(m_prior_strength) * M_prior)
            M_new = _normalise_columns(expected)
        else:
            M_new = M

        # Convergence check on log-likelihood
        Lambda_new = T_new @ M_new.T + B[np.newaxis, :]
        Lambda_new[Lambda_new < _EPS] = _EPS
        ll = float(np.sum(O_counts * np.log(Lambda_new) - Lambda_new))

        T = T_new
        M = M_new

        if iteration > 0 and abs(ll - prev_ll) < tol * (abs(prev_ll) + _EPS):
            break
        prev_ll = ll

    # Row-normalise T to posterior probabilities
    P   = T.copy()
    p_s = P.sum(axis=1, keepdims=True)
    p_s[p_s == 0] = 1.0
    P  /= p_s

    return P, M, B

from __future__ import annotations

from typing import Tuple
import numpy as np


def _transition_matrix(lambda_switch: float) -> np.ndarray:
    """Build a 3x3 transition matrix with strong stay probability.

    Map lambda to stay prob s via s = 1 - 0.5*exp(-lambda). For lambda=2 → s≈0.932.
    Off-diagonals share remaining mass equally.
    """
    s = float(1.0 - 0.5 * np.exp(-max(0.0, lambda_switch)))
    s = min(max(s, 0.5), 0.999)
    off = (1.0 - s) / 2.0
    A = np.array(
        [
            [s, off, off],
            [off, s, off],
            [off, off, s],
        ],
        dtype=float,
    )
    return A


def smooth_posteriors(S: np.ndarray, lambda_switch: float) -> Tuple[np.ndarray, np.ndarray]:
    """Forward-backward smoothing of 3-state posteriors with simple transition prior.

    S: shape (N,3) emissions (unnormalized ok; will be normalized per row)
    Returns (P, path): smoothed marginals and Viterbi path over {0,1,2}.
    """
    if S.size == 0:
        return S, np.zeros((0,), dtype=int)

    N = S.shape[0]
    E = S.copy()
    # Normalize emissions per row to avoid underflow
    row_sum = E.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    E = E / row_sum

    A = _transition_matrix(lambda_switch)
    pi = np.array([1 / 3, 1 / 3, 1 / 3], dtype=float)

    # Forward pass with scaling
    alpha = np.zeros_like(E)
    c = np.zeros((N,), dtype=float)
    alpha[0] = pi * E[0]
    c[0] = alpha[0].sum() or 1.0
    alpha[0] /= c[0]
    for t in range(1, N):
        alpha[t] = (alpha[t - 1] @ A) * E[t]
        c[t] = alpha[t].sum() or 1.0
        alpha[t] /= c[t]

    # Backward pass
    beta = np.zeros_like(E)
    beta[-1] = 1.0
    for t in range(N - 2, -1, -1):
        beta[t] = A @ (E[t + 1] * beta[t + 1])
        s = beta[t].sum() or 1.0
        beta[t] /= s

    P = alpha * beta
    P_sum = P.sum(axis=1, keepdims=True)
    P_sum[P_sum == 0] = 1.0
    P /= P_sum

    # Viterbi path
    delta = np.zeros_like(E)
    psi = np.zeros((N, 3), dtype=int)
    delta[0] = np.log(pi + 1e-12) + np.log(E[0] + 1e-12)
    logA = np.log(A + 1e-12)
    for t in range(1, N):
        vals = delta[t - 1][:, None] + logA
        psi[t] = np.argmax(vals, axis=0)
        delta[t] = np.max(vals, axis=0) + np.log(E[t] + 1e-12)
    path = np.zeros((N,), dtype=int)
    path[-1] = int(np.argmax(delta[-1]))
    for t in range(N - 2, -1, -1):
        path[t] = int(psi[t + 1, path[t + 1]])

    return P, path

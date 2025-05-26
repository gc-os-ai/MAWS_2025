"""
routines.py – numerical helpers used by MAWS.

Most functions operate on **iterables of energies** (floats).  NumPy is used
for speed; set `use_mp=True` if you really need >64-bit precision.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from mpmath import mp

# --------------------------------------------------------------------------- #
# Global high-precision control
# --------------------------------------------------------------------------- #
mp.dps = 60  # Decimal places for mpmath when requested


# --------------------------------------------------------------------------- #
# Core mathematics
# --------------------------------------------------------------------------- #
def kullback_leibler_divergence(
    sample: Sequence[float],
    reference: Sequence[float],
    *,
    use_mp: bool = False,
) -> float:
    """
    KL(P‖Q) = ∑ Pᵢ log(Pᵢ / Qᵢ).

    Parameters
    ----------
    sample, reference
        Probability distributions **P** and **Q**.  Must be the same length and
        sum to 1.
    use_mp
        If *True*, compute with mpmath `mp.mpf` high precision.

    Returns
    -------
    float
        KL divergence (≥ 0); `np.inf` if any reference probability is 0 while
        sample is non-zero.

    Examples
    --------
    >>> kullback_leibler_divergence([0.5, 0.5], [0.4, 0.6])  # doctest: +ELLIPSIS
    0.020...
    """
    if len(sample) != len(reference):
        raise ValueError("sample and reference must have the same length")

    if use_mp:
        p = list(map(mp.mpf, sample))
        q = list(map(mp.mpf, reference))
        return float(-mp.fsum(pi * mp.log(pi / qi) for pi, qi in zip(p, q)))

    p_arr = np.asarray(sample, dtype=float)
    q_arr = np.asarray(reference, dtype=float)
    if np.any(q_arr == 0.0) and np.any((p_arr > 0) & (q_arr == 0)):
        return float("inf")

    with np.errstate(divide="ignore", invalid="ignore"):
        log_term = np.log(p_arr / q_arr, where=(p_arr > 0))
    return float(-np.sum(p_arr * log_term, where=(p_arr > 0)))


def entropy(
    probabilities: Sequence[float],
    *,
    use_mp: bool = False,
) -> float:
    r"""
    Shannon entropy: S = –∑ pᵢ log(pᵢ).

    Notes
    -----
    If the input is *not* a probability distribution yet, you can still call
    this function – the values are first normalised to sum 1.
    """
    if use_mp:
        probs = list(map(mp.mpf, probabilities))
        total = mp.fsum(probs)
        probs = [p / total for p in probs]
        return float(-mp.fsum(p * mp.log(p) for p in probs))

    p_arr = np.asarray(probabilities, dtype=float)
    p_arr = p_arr / p_arr.sum()
    with np.errstate(divide="ignore"):
        log_p = np.log(p_arr, where=(p_arr > 0))
    return float(-np.sum(p_arr * log_p, where=(p_arr > 0)))


# --------------------------------------------------------------------------- #
# Energy-dependent helpers
# --------------------------------------------------------------------------- #
def good_energy(energy: float, threshold: float) -> float:
    """Return +1 if *energy* ≤ *threshold*, else –1."""
    return 1.0 if energy <= threshold else -1.0


def zps(
    energies: Sequence[float],
    *,
    beta: float = 0.01,
) -> Tuple[float, NDArray[np.floating], float]:
    """
    Partition function **Z**, Boltzmann probabilities **P**, and entropy **S**.

    Z = ∑ exp(–βEᵢ),   Pᵢ = exp(–βEᵢ) / Z,   S = –∑ Pᵢ log Pᵢ
    """
    e_arr = np.asarray(energies, dtype=float)
    boltz = np.exp(-beta * e_arr)
    Z = boltz.sum()
    P = boltz / Z
    S = entropy(P)
    return float(Z), P, S


def s(energies: Sequence[float], *, beta: float = 0.01) -> float:
    """Wrapper that returns only the entropy part of :func:`zps`."""
    return zps(energies, beta=beta)[2]


def score(energies: Sequence[float], *, beta: float = 0.01) -> float:
    """
    Optimisation score balancing low energy *and* high entropy.

    Defined as ``–S * min(E)``.
    """
    return -s(energies, beta=beta) * min(energies)


def score_threshold(
    energies: Sequence[float],
    *,
    beta: float = 0.01,
    threshold: float = 0.0,
) -> float:
    """
    Variant that penalises energies above *threshold*.

    Returns ``S * good_energy(min(E), threshold)``.
    """
    return s(energies, beta=beta) * good_energy(min(energies), threshold)


# --------------------------------------------------------------------------- #
# Candidate selection helpers
# --------------------------------------------------------------------------- #
def best_position(
    positions: NDArray[np.floating] | List[NDArray[np.floating]],
    free_energies: Sequence[float],
) -> NDArray[np.floating]:
    """
    Return the position with the **lowest** free energy.

    Parameters
    ----------
    positions
        Array-like list of 3-D coordinates (or higher-dimensional features).
    free_energies
        Energies aligned with *positions*.

    Returns
    -------
    numpy.ndarray
        The position with minimal energy.
    """
    idx = int(np.argmin(np.asarray(free_energies)))
    return np.asarray(positions)[idx]


def choose_candidates(
    entropies: Sequence[float],
    sequences: Sequence[str],
    *,
    threshold: float = 0.0,
) -> List[str]:
    """
    Pick sequences whose entropy lies within *threshold* of the best one.

    Returned list is **sorted** by increasing entropy.
    """
    ent_arr = np.asarray(entropies, dtype=float)
    order = np.argsort(ent_arr)
    best_ent = ent_arr[order[0]]
    good_mask = ent_arr <= best_ent + threshold
    return list(np.asarray(sequences)[order][good_mask])


# --------------------------------------------------------------------------- #
# Backward-compatibility camelCase aliases
# --------------------------------------------------------------------------- #
kullbackLeiblerDivergence = kullback_leibler_divergence
Entropy = entropy
GoodEnergy = good_energy
ZPS = zps
S = s
Score = score
ScoreThreshold = score_threshold
best_position_k = best_position  # to avoid name clash; adjust as needed
choose_candidates_k = choose_candidates

__all__ = [
    # snake_case
    "kullback_leibler_divergence",
    "entropy",
    "good_energy",
    "zps",
    "s",
    "score",
    "score_threshold",
    "best_position",
    "choose_candidates",
    # original camelCase
    "kullbackLeiblerDivergence",
    "Entropy",
    "GoodEnergy",
    "ZPS",
    "S",
    "Score",
    "ScoreThreshold",
    "best_position_k",
    "choose_candidates_k",
]
########################################################
# if __name__ == "__main__":
#     energies = [1.2, 0.9, 1.5, 0.7]
#     print("Entropy:", s(energies))
#     print("Score  :", score(energies))
#     print("Best pos idx:", best_position(np.arange(4), energies))

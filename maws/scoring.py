"""
maws.scoring
============

Turning many energies into one number that says how promising a strand is.

MAWS tries thousands of shapes for a candidate strand and gets an energy for
each. It then has to reduce that pile of numbers to a single score, so that one
candidate can be compared with another. That is what :func:`entropy_score`
does, and it is the heart of how MAWS chooses.

Why not simply take the lowest energy? Because a strand that has one good shape
and thousands of bad ones is not a good binder. Binding means the strand
settles into a small family of shapes and stays there. So the score measures
how *concentrated* the energies are, not how low they get.

Concretely, the energies are turned into probabilities — each shape is weighted
by ``exp(-beta * energy)``, so lower-energy shapes count for more — and the
score is the spread of that distribution. A strand whose weight piles onto a
few shapes scores low. One whose weight is spread evenly over everything scores
high. Lower is better.

Examples
--------
Energies that are all equal give no preference to any shape at all, the worst
possible outcome, and score zero:

>>> abs(entropy_score([100.0, 100.0, 100.0])) < 1e-12
True

One shape far better than the rest concentrates the weight, and scores well
below zero:

>>> entropy_score([0.0] + [1000.0] * 9) < entropy_score([0.0] + [10.0] * 9)
True
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol, runtime_checkable

import numpy as np
from mpmath import mp

from maws.errors import ConfigurationError

__all__ = ["Scorer", "boltzmann_weights", "entropy_score"]

mp.dps = 60
"""Digits of precision used for the sums below.

Energies of a whole molecule run to thousands of kJ/mol, and the weights are
exponentials of those numbers, which ordinary floating point cannot hold: the
largest weights overflow and the smallest vanish to zero. Arbitrary-precision
arithmetic sidesteps the problem outright.
"""


@runtime_checkable
class Scorer(Protocol):
    """Anything that reduces a list of energies to one number, lower is better.

    :func:`entropy_score` satisfies this. Writing it as an interface means an
    alternative — one that also takes into account how well the strand folds on
    its own, say — can be dropped in without touching the search.
    """

    def __call__(self, energies: Sequence[float], *, beta: float) -> float:
        """Return the score for one candidate.

        Parameters
        ----------
        energies : sequence of float
            Energies of the shapes tried for this candidate, in kJ/mol.
        beta : float
            How sharply lower energies are favoured.

        Returns
        -------
        float
            The score. Lower means a more promising candidate.
        """
        ...


def entropy_score(energies: Sequence[float], *, beta: float = 0.01) -> float:
    """Return how concentrated a candidate's energies are. Lower is better.

    Parameters
    ----------
    energies : sequence of float
        Energies of the shapes tried for one candidate strand, in kJ/mol. At
        least one is required.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ. Raising it makes
        the score depend mostly on the few best shapes; lowering it lets every
        shape count roughly equally, and at zero the score is always zero
        because every shape then weighs the same.

    Returns
    -------
    float
        Zero when every shape is equally likely, and increasingly negative as
        the weight concentrates onto fewer shapes.

    Raises
    ------
    maws.errors.ConfigurationError
        If `energies` is empty or `beta` is negative.

    Notes
    -----
    The calculation is:

    1. weight each shape by ``exp(-beta * energy)``;
    2. divide by the total, giving a probability ``p`` per shape;
    3. return ``-sum(p * log(p * N))``, where ``N`` is the number of shapes.

    The ``* N`` inside the logarithm is what puts the zero point at "every
    shape equally likely" regardless of how many shapes were tried, so scores
    from runs that sampled different numbers of shapes can still be compared.

    Examples
    --------
    Ten shapes, one of them far better than the rest, versus ten shapes that
    are all much of a muchness. The first scores lower, meaning more promising:

    >>> one_clear_winner = entropy_score([0.0] + [1000.0] * 9)
    >>> nothing_to_choose = entropy_score([0.0] + [10.0] * 9)
    >>> one_clear_winner < nothing_to_choose
    True
    """
    values = list(energies)
    if not values:
        raise ConfigurationError(
            "entropy_score needs at least one energy; it was given none, "
            "which usually means the sampling loop never ran"
        )
    if beta < 0:
        raise ConfigurationError(f"beta must not be negative, got {beta}")

    count = len(values)
    weights = [mp.exp(mp.fmul(-beta, value)) for value in values]
    total = mp.fsum(weights)
    probabilities = [mp.fdiv(weight, total) for weight in weights]
    entropy = -mp.fsum(mp.fmul(p, mp.log(mp.fmul(p, count))) for p in probabilities)
    return float(entropy)


def boltzmann_weights(energies: Sequence[float], *, beta: float = 0.01) -> np.ndarray:
    """Return the probability of each shape, as :func:`entropy_score` sees it.

    Exposed for inspecting a run: it answers "how much of the score came from
    which shapes", which the score alone cannot say.

    Parameters
    ----------
    energies : sequence of float
        Energies of the shapes tried, in kJ/mol.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ.

    Returns
    -------
    numpy.ndarray
        Shape ``(len(energies),)``, adding up to one.

    Raises
    ------
    maws.errors.ConfigurationError
        If `energies` is empty or `beta` is negative.

    Examples
    --------
    >>> weights = boltzmann_weights([0.0, 0.0])
    >>> weights
    array([0.5, 0.5])
    """
    values = list(energies)
    if not values:
        raise ConfigurationError("boltzmann_weights needs at least one energy")
    if beta < 0:
        raise ConfigurationError(f"beta must not be negative, got {beta}")
    weights = [mp.exp(mp.fmul(-beta, value)) for value in values]
    total = mp.fsum(weights)
    return np.array([float(mp.fdiv(weight, total)) for weight in weights])

r"""
maws.scoring
============

Turning many energies into one number that says how promising a strand is.

MAWS tries thousands of shapes for a candidate strand and gets an energy for
each. It then has to reduce that pile of numbers to a single score, so that one
candidate can be compared with another. Lower always means more promising.

Two ways of doing that live here, and the difference between them matters.

:func:`free_energy_score` is what the search uses. It asks how much of the
strand's freedom to move ends up in shapes that are low in energy. A strand
scores well when its energies are low *and* many of the shapes tried reach
them. It is the ordinary free energy of the shapes that were tried, so it is
measured in kJ/mol and can be compared against an energy directly.

:func:`entropy_score` is the original MAWS score and is kept because runs are
compared against published ones. It asks a narrower question: how *concentrated*
the weight is across the shapes tried, ignoring how low the energies are. That
is a real difference, not a matter of scale — see the warning on that function
for the case it gets wrong.

Examples
--------
A strand whose shapes are all around -1000 kJ/mol is more promising than one
whose shapes are all around -100, and the free energy says so:

>>> free_energy_score([-1000.0] * 4) < free_energy_score([-100.0] * 4)
True

Between two strands whose best shapes are equally good, the one that reaches
that energy more often scores lower:

>>> free_energy_score([-500.0] * 4) < free_energy_score([-500.0] + [0.0] * 3)
True
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol, runtime_checkable

import numpy as np
from mpmath import mp

from maws.errors import ConfigurationError

__all__ = ["Scorer", "boltzmann_weights", "entropy_score", "free_energy_score"]

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


def free_energy_score(energies: Sequence[float], *, beta: float = 0.01) -> float:
    r"""Return the free energy of the shapes tried, in kJ/mol. Lower is better.

    Parameters
    ----------
    energies : sequence of float
        Energies of the shapes tried for one candidate strand, in kJ/mol. At
        least one is required.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ. Raising it makes
        the answer depend mostly on the few best shapes, approaching the
        lowest energy of all. Lowering it lets every shape count roughly
        equally, approaching their average.

    Returns
    -------
    float
        A value in kJ/mol, always between the lowest energy given and their
        average.

    Raises
    ------
    maws.errors.ConfigurationError
        If `energies` is empty or `beta` is not positive. Zero would mean
        dividing by it.

    See Also
    --------
    entropy_score : The original MAWS score, which asks a different question.
    boltzmann_weights : How much each shape counted towards the answer.

    Notes
    -----
    The free energy of a set of sampled shapes is

    .. math::
        F = -\frac{1}{\beta} \ln\!\left(
            \frac{1}{N} \sum_{i=1}^{N} e^{-\beta E_i}
        \right).

    Two different things push it down, which is what makes it a fair
    comparison between candidates. Low energies push it down, because each
    term grows as its energy falls. And *many* low energies push it down
    further, because the sum has more large terms in it. A strand that reaches
    -1000 kJ/mol in one shape out of a thousand and clashes in the other 999
    is therefore beaten by one that reaches -1000 kJ/mol in most of them,
    which is the behaviour wanted: the first has found a coincidence, the
    second has found a fit.

    Dividing by :math:`N` before taking the logarithm is what lets runs that
    tried different numbers of shapes be compared. Without it, simply trying
    more shapes would lower the score.

    Examples
    --------
    Shapes that are all equally good give back that energy exactly, whatever
    `beta` is:

    >>> free_energy_score([-250.0] * 8)
    -250.0

    The answer never goes below the best shape, nor above the average:

    >>> energies = [-900.0, -100.0, -100.0, -100.0]
    >>> min(energies) <= free_energy_score(energies) <= sum(energies) / 4
    True

    Raising `beta` moves it towards the best shape:

    >>> sharp = free_energy_score(energies, beta=1.0)
    >>> gentle = free_energy_score(energies, beta=0.001)
    >>> sharp < gentle
    True
    """
    values = list(energies)
    if not values:
        raise ConfigurationError(
            "free_energy_score needs at least one energy; it was given none, "
            "which usually means the sampling loop never ran"
        )
    if beta <= 0:
        raise ConfigurationError(f"beta must be greater than zero, got {beta}")

    weights = [mp.exp(mp.fmul(-beta, value)) for value in values]
    average = mp.fdiv(mp.fsum(weights), len(values))
    return float(mp.fdiv(-mp.log(average), beta))


def entropy_score(energies: Sequence[float], *, beta: float = 0.01) -> float:
    """Return how concentrated a candidate's energies are. Lower is better.

    .. warning::
        This is not what the search uses, and it should not be chosen without
        reading the rest of this docstring. It answers "how concentrated",
        never "how low", and those come apart badly. A candidate whose shapes
        almost all crash atoms into each other scores *better* here than one
        whose shapes are all reasonable, because a few survivors among many
        disasters is a very concentrated distribution. See the second example
        below. :func:`free_energy_score` is the default for that reason.

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

    See Also
    --------
    free_energy_score : What the search uses instead.

    Notes
    -----
    The calculation is:

    1. weight each shape by ``exp(-beta * energy)``;
    2. divide by the total, giving a probability ``p`` per shape;
    3. return ``-sum(p * log(p * N))``, where ``N`` is the number of shapes.

    The ``* N`` inside the logarithm is what puts the zero point at "every
    shape equally likely" regardless of how many shapes were tried, so scores
    from runs that sampled different numbers of shapes can still be compared.

    Step 2 is also where the trouble comes from. Dividing by the total throws
    away the scale of the energies: adding the same amount to every energy
    leaves the answer untouched, so the score cannot tell a strand that
    reaches -1000 kJ/mol from one that never gets below zero.

    Examples
    --------
    Ten shapes, one of them far better than the rest, versus ten shapes that
    are all much of a muchness. The first scores lower, meaning more promising:

    >>> one_clear_winner = entropy_score([0.0] + [1000.0] * 9)
    >>> nothing_to_choose = entropy_score([0.0] + [10.0] * 9)
    >>> one_clear_winner < nothing_to_choose
    True

    The same arithmetic run on a candidate that mostly crashes atoms into each
    other. Eight of its ten shapes are impossible and two are ordinary, and it
    still scores far better than a candidate whose ten shapes are all good:

    >>> mostly_impossible = entropy_score([-100.0, -100.0] + [1e6] * 8)
    >>> all_reasonable = entropy_score([-100.0] * 10)
    >>> mostly_impossible < all_reasonable
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

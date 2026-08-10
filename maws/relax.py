"""
maws.relax
==========

Shaking a structure loose from an awkward arrangement.

Joining a new residue onto a strand leaves the atoms near the join in positions
that are geometrically correct but physically strained: bonds at slightly wrong
lengths, atoms slightly too close. Settling the structure fixes the worst of
that by rolling downhill in energy, but downhill from a strained start can lead
into a poor arrangement that is merely the nearest one, not a good one.

:func:`perturb_and_minimize` gets past that by nudging every atom a little at
random and settling again, several times over. The nudges are small enough not
to destroy the shape and large enough to escape a shallow dip.

Examples
--------
>>> import numpy as np
>>> from maws.energy import StubEnergy
>>> from maws.pose import Pose
>>> from maws.values import AtomRange
>>> pose = Pose(np.array([[0.0, 0, 0], [3.0, 0, 0]]))
>>> model = StubEnergy(AtomRange(0, 1), AtomRange(1, 2))
>>> settled = perturb_and_minimize(
...     pose, model, iterations=3, rng=np.random.default_rng(0)
... )
>>> settled.pose.n_atoms
2
"""

from __future__ import annotations

import numpy as np

from maws.energy import EnergyModel, Relaxed
from maws.errors import ConfigurationError
from maws.pose import Pose

__all__ = ["perturb_and_minimize"]


def perturb_and_minimize(
    pose: Pose,
    energy: EnergyModel,
    *,
    size: float = 0.1,
    iterations: int = 50,
    max_iterations: int = 100,
    rng: np.random.Generator | None = None,
) -> Relaxed:
    """perturb_and_minimize(pose, energy, *, size=0.1, iterations=50) -> Relaxed

    Nudge every atom at random and settle the structure, repeatedly.

    Parameters
    ----------
    pose : maws.pose.Pose
        The atom positions to start from.
    energy : maws.energy.EnergyModel
        What scores and settles the structure.
    size : float, default=0.1
        How far each atom may be nudged along each axis, in ångström. Drawn
        evenly from ``-size`` to ``+size``. Larger values escape a poor
        arrangement more readily but risk pulling the structure apart; 0.5 Å is
        a reasonable choice straight after a new residue has been joined on.
    iterations : int, default=50
        How many nudge-and-settle rounds to run. More rounds explore further at
        a proportional cost.
    max_iterations : int, default=100
        How many adjustment steps each settling may take.
    rng : numpy.random.Generator, optional
        Source of randomness. Pass one built with a fixed seed to make a run
        repeatable. Defaults to a fresh generator, so runs differ.

    Returns
    -------
    maws.energy.Relaxed
        The positions after the final round, and their energy.

    Raises
    ------
    maws.errors.ConfigurationError
        If `size` is negative, or `iterations` is negative.

    See Also
    --------
    maws.energy.EnergyModel.minimize : One settling step on its own.

    Notes
    -----
    .. warning::
        The result is the last round's, not the best round's. Each round starts
        from where the previous one finished, so the sequence is a walk rather
        than a search, and its final energy can be higher than one it passed
        through. This matches how the structure is meant to be used afterwards:
        as a starting point for sampling, not as an answer.

    Examples
    --------
    >>> import numpy as np
    >>> from maws.energy import StubEnergy
    >>> from maws.pose import Pose
    >>> from maws.values import AtomRange
    >>> pose = Pose(np.array([[0.0, 0, 0], [3.0, 0, 0]]))
    >>> model = StubEnergy(AtomRange(0, 1), AtomRange(1, 2))
    >>> result = perturb_and_minimize(
    ...     pose, model, size=0.5, iterations=2, rng=np.random.default_rng(0)
    ... )
    >>> result.energy > 0
    True
    """
    if size < 0:
        raise ConfigurationError(f"size must not be negative, got {size}")
    if iterations < 0:
        raise ConfigurationError(f"iterations must not be negative, got {iterations}")

    generator = rng if rng is not None else np.random.default_rng()
    result = Relaxed(energy.evaluate(pose), pose)
    for _ in range(iterations):
        nudged = result.pose.jittered(
            generator.uniform(-size, size, result.pose.xyz.shape)
        )
        result = energy.minimize(nudged, max_iterations=max_iterations)
    return result

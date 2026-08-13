"""
maws.relax
==========

Shaking a structure loose from an awkward arrangement.

Joining a new *residue* — a single nucleotide of the strand — onto an aptamer
leaves the atoms near the join in positions that are geometrically correct but
physically strained: bonds at slightly wrong lengths, atoms slightly too close.
Settling the structure fixes the worst of that by rolling downhill in energy,
but downhill from a strained start can lead into a poor arrangement that is
merely the nearest one, not a good one.

:func:`perturb_and_minimize` gets past that by nudging every atom a little at
random and settling again, several times over. The nudges are small enough not
to destroy the shape and large enough to escape a shallow dip. It is the only
thing in this module.

Distances here are in ångström and energies in kJ/mol.

See Also
--------
maws.energy.EnergyModel : Supplies both the score and the settling.
maws.pose.Pose.jittered : Applies one round of nudges.

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
from maws.pose import ChainView, Pose
from maws.values import AtomRange

__all__ = ["perturb_and_minimize"]


def perturb_and_minimize(
    pose: Pose,
    energy: EnergyModel,
    *,
    size: float = 0.1,
    iterations: int = 50,
    max_iterations: int = 100,
    moving: AtomRange | ChainView | None = None,
    rng: np.random.Generator | None = None,
) -> Relaxed:
    """perturb_and_minimize(pose, energy, *, size=0.1, iterations=50,
                         max_iterations=100, moving=None, rng=None) -> Relaxed

    Nudge every atom at random and settle the structure, repeatedly.

    One round adds an independent random displacement to every atom and then
    lets the energy model roll the whole structure downhill again. The rounds
    run one after another, each starting from where the last one finished.

    Parameters
    ----------
    pose : maws.pose.Pose
        The atom positions to work from.
    energy : maws.energy.EnergyModel
        What scores and settles the structure. Its energies are in kJ/mol.
    size : float, default=0.1
        How far each atom may be nudged along each axis, in ångström. Drawn
        evenly from ``-size`` to ``+size``. Larger values escape a poor
        arrangement more readily but risk pulling the structure apart; 0.5 Å is
        a reasonable choice straight after a new residue has been joined on.
    iterations : int, default=50
        How many nudge-and-settle rounds to run. More rounds explore further at
        a proportional cost. Zero does no nudging at all and simply reports the
        energy of `pose` as it stands.
    max_iterations : int, default=100
        How many adjustment steps to allow before stopping, within each single
        settling. Raising it settles each round more thoroughly and costs more.
    moving : maws.values.AtomRange or maws.pose.ChainView, optional
        Which atoms to nudge. Normally the strand, so that the target it is
        being fitted against is not shaken about as well. Every atom is nudged
        by default.
    rng : numpy.random.Generator, optional
        Source of randomness. Pass one built with a fixed seed to make a run
        repeatable. Defaults to a fresh generator, so runs differ.

    Returns
    -------
    maws.energy.Relaxed
        The positions after the final round, and their energy in kJ/mol. A new
        pose; `pose` is unchanged.

    Raises
    ------
    maws.errors.ConfigurationError
        If `size` is negative, or `iterations` is negative.

    See Also
    --------
    maws.energy.EnergyModel.minimize : One settling step on its own.
    maws.pose.Pose.jittered : Applies the random displacements of one round.

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
        offsets = np.zeros(result.pose.xyz.shape)
        rows = slice(None) if moving is None else _span_of(moving).as_slice()
        offsets[rows] = generator.uniform(-size, size, offsets[rows].shape)
        result = energy.minimize(
            result.pose.jittered(offsets), max_iterations=max_iterations
        )
    return result


def _span_of(where: AtomRange | ChainView) -> AtomRange:
    """Return the run of atoms named by a span or by a chain.

    Parameters
    ----------
    where : maws.values.AtomRange or maws.pose.ChainView
        Either a run of atom numbers, or a chain that has one.

    Returns
    -------
    maws.values.AtomRange
        The run of atom numbers.
    """
    return where if isinstance(where, AtomRange) else where.span

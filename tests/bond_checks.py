"""
Checking that a move of the atoms left the molecule in one piece.

Turning a bond is supposed to change a molecule's *shape* and nothing else. No
two bonded atoms may end up further apart, or closer together, than they
started. A rotation that moves one atom of a bond but not the other stretches
that bond, and a structure with a stretched bond is not a molecule any more —
it will still produce an energy, a score, and a ranking, all of them
meaningless, and nothing else in the package notices.

These helpers are the check for that. They are used by
``tests/unit/test_geometry_invariants.py`` on shapes written out by hand, and
by ``tests/integration/test_real_geometry.py`` on structures built by
AmberTools.

Which atoms are bonded has to be said, not guessed, so
:func:`assert_bonds_survive` takes the list of pairs. There are two ways to get
one. :func:`bonds_from_topology` reads the real list out of what AmberTools
worked out, and is what the integration tier uses. :func:`bonded_pairs` guesses
it from how close the atoms are, which is only safe on a shape written out by
hand: a structure straight out of ``tleap`` has not been settled yet and can
hold two unbonded atoms 1.1 Å apart, which the guess would call a bond and then
report as broken the moment anything moved.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from maws.pose import Pose

if TYPE_CHECKING:  # pragma: no cover - needed by type checkers only
    from openmm.app import Topology

__all__ = [
    "BOND_CUTOFF",
    "TOLERANCE",
    "assert_bonds_survive",
    "bond_lengths",
    "bonded_pairs",
    "bonds_from_topology",
]

BOND_CUTOFF = 1.8
"""Atoms closer than this many ångström are treated as bonded.

Comfortably above a carbon-carbon single bond (1.54 Å), and used only on shapes
written out by hand where nothing unbonded is placed that close.
"""

TOLERANCE = 1e-6
"""How far a bond length may move before it counts as broken, in ångström.

A rotation is an exact rigid motion of the atoms it moves, so a bond with both
its atoms inside one moving group should not change length at all. Anything
above rounding error means the two atoms were moved differently.
"""


def bonded_pairs(xyz: np.ndarray, cutoff: float = BOND_CUTOFF) -> np.ndarray:
    """Guess which atoms are bonded from how close together they are.

    Only for shapes written out by hand. See the module docstring for why a
    built structure needs :func:`bonds_from_topology` instead.

    Parameters
    ----------
    xyz : numpy.ndarray
        Shape ``(N, 3)`` positions in ångström.
    cutoff : float, default=1.8
        How close two atoms must be, in ångström, to count as bonded.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_pairs, 2)`` of atom index pairs, each listed once with the
        lower index first.

    Examples
    --------
    Three atoms in a line, 1 Å apart: the two neighbouring pairs are bonded and
    the end-to-end pair is not.

    >>> import numpy as np
    >>> bonded_pairs(np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]]))
    array([[0, 1],
           [1, 2]])
    """
    gaps = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    close = (gaps < cutoff) & (gaps > 0.0)
    return np.argwhere(np.triu(close))


def bonds_from_topology(topology: Topology) -> np.ndarray:
    """Return the bonds AmberTools worked out, as pairs of atom numbers.

    Parameters
    ----------
    topology : openmm.app.Topology
        The atom-and-bond listing that came with a built structure, reached as
        ``system.amber.topology``.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_bonds, 2)`` of atom index pairs, lower index first. The
        numbering is the same one a :class:`~maws.pose.Pose` uses, because both
        come from the same build.
    """
    pairs = [sorted((bond.atom1.index, bond.atom2.index)) for bond in topology.bonds()]
    return np.array(sorted(pairs), dtype=int)


def bond_lengths(xyz: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """Return the distance across each listed pair of atoms.

    Parameters
    ----------
    xyz : numpy.ndarray
        Shape ``(N, 3)`` positions in ångström.
    pairs : numpy.ndarray
        Shape ``(n_pairs, 2)`` of atom index pairs.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_pairs,)`` distances in ångström.

    Examples
    --------
    >>> import numpy as np
    >>> xyz = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]])
    >>> bond_lengths(xyz, np.array([[0, 2]]))
    array([2.])
    """
    return np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=-1)


def assert_bonds_survive(
    before: Pose, after: Pose, what: str, *, pairs: np.ndarray
) -> None:
    """Fail if any listed bond has a different length in `after` than in `before`.

    Parameters
    ----------
    before : maws.pose.Pose
        Positions before the move.
    after : maws.pose.Pose
        Positions after it. Must describe the same atoms in the same order.
    what : str
        What was done, quoted in the failure message.
    pairs : numpy.ndarray
        Shape ``(n_pairs, 2)`` of bonded atom index pairs, from
        :func:`bonds_from_topology` or :func:`bonded_pairs`.

    Raises
    ------
    AssertionError
        If any bond length changed by more than :data:`TOLERANCE`. The message
        names the worst pair and by how much it moved, so the move that caused
        it can be found without re-running anything.

    Examples
    --------
    Turning both atoms of a bond together leaves it alone:

    >>> import numpy as np
    >>> from maws.pose import Pose
    >>> from maws.values import AtomRange, Torsion
    >>> pose = Pose(np.array([[0.0, 0, 0], [1.0, 0, 0], [1.0, 1, 0]]))
    >>> turned = pose.rotate(Torsion(0, 1, AtomRange(1, 3)), 1.0)
    >>> assert_bonds_survive(
    ...     pose, turned, "a quarter turn", pairs=bonded_pairs(pose.xyz)
    ... )
    """
    assert len(pairs) > 0, "no bonds were given, so nothing would be checked"

    lengths_before = bond_lengths(before.xyz, pairs)
    lengths_after = bond_lengths(after.xyz, pairs)
    drift = np.abs(lengths_after - lengths_before)
    worst = int(np.argmax(drift))
    assert drift[worst] < TOLERANCE, (
        f"{what} changed a bond length. Atoms "
        f"{pairs[worst, 0]}-{pairs[worst, 1]} went from "
        f"{lengths_before[worst]:.3f} to {lengths_after[worst]:.3f} Å "
        f"({drift[worst]:.3f} Å); {int((drift > TOLERANCE).sum())} of "
        f"{len(pairs)} bonds moved."
    )

"""
Invariants that any move of the atoms must respect.

Every other test in this suite checks that some operation does what its
docstring says. These check something cruder: that the molecule survives it.

Three things are checked here. That a move which is supposed to be rigid
really is one, so that sliding or turning a whole chain leaves every distance
inside it alone. That no turnable bond offered to the search has nothing to
move. And that the bond-length check itself both passes on a turn that keeps a
molecule together and fails on one that does not.

The shapes used are written out by hand or produced by the stand-in builder,
which lays atoms out on a grid. Grid positions are not chemistry, so what is
checked here is only what does not depend on the positions being chemically
sensible. Everything that does — above all, that turning any real bond of a
real strand leaves every other bond alone — is in
``tests/integration/test_real_geometry.py``, which builds with AmberTools.

See Also
--------
bond_checks : The bond-length check itself, shared with the integration tier.
"""

from __future__ import annotations

import numpy as np
import pytest
from bond_checks import assert_bonds_survive

from maws.build import FakeBuilder
from maws.forcefield import ForceField
from maws.libraries import rna
from maws.pose import Pose
from maws.topology import Assembly
from maws.values import AtomRange, Torsion


@pytest.fixture
def strand():
    """A four-nucleotide strand, built without AmberTools.

    The grid positions are not chemistry, but the residue spans, the atom
    counts and the turnable bonds are the real ones, which is what the checks
    below are about.
    """
    return FakeBuilder(spacing=1.5).build(
        Assembly().with_aptamer(rna(), "G A U C"),
        ForceField.for_target("RNA", "protein"),
    )


def all_distances(pose: Pose) -> np.ndarray:
    """Return the distance between every pair of atoms.

    Parameters
    ----------
    pose : maws.pose.Pose
        The positions to measure.

    Returns
    -------
    numpy.ndarray
        Shape ``(N, N)``, in ångström, symmetric with a zero diagonal.

    Notes
    -----
    Used instead of a list of bonds where every atom moves together. Nothing
    then needs to be known about which atoms are bonded: if the move was rigid,
    *no* distance changed, which is a stronger statement and one that holds
    whatever the positions happen to be.
    """
    return np.linalg.norm(pose.xyz[:, None, :] - pose.xyz[None, :, :], axis=-1)


class TestMovingAWholeChainChangesNoDistance:
    """Sliding or turning an entire chain cannot change its internal shape."""

    def test_sliding_a_chain_keeps_every_distance(self, strand):
        """Moving every atom of a chain by the same amount changes no distance."""
        span = strand.chain("aptamer").span
        moved = strand.pose.translate(span, np.array([12.0, -3.0, 7.5]))
        np.testing.assert_allclose(
            all_distances(moved), all_distances(strand.pose), atol=1e-9
        )

    def test_turning_a_whole_chain_keeps_every_distance(self, strand):
        """Turning a chain as one rigid body changes where it is, not what it is."""
        span = strand.chain("aptamer").span
        turned = strand.pose.rotate_about(
            span,
            axis=np.array([0.3, 1.0, -0.4]),
            angle=1.1,
            origin=strand.pose.centroid(span),
        )
        np.testing.assert_allclose(
            all_distances(turned), all_distances(strand.pose), atol=1e-9
        )

    def test_turning_by_a_full_circle_puts_every_atom_back(self, strand):
        """A whole turn is the identity, which is the cheapest sanity check."""
        torsion = strand.chain("aptamer").residue(0).torsion(0)
        turned = strand.pose.rotate(torsion, 2.0 * np.pi)
        np.testing.assert_allclose(turned.xyz, strand.pose.xyz, atol=1e-9)


class TestTorsionsMoveSomething:
    """A bond that moves nothing is a wasted degree of freedom.

    The search spends one of its random angles on each bond it is handed. One
    that cannot move anything silently reduces how much of the strand's shape
    is explored — and at the ends of a strand it does so in one direction only,
    which would quietly make growing at one end explore less than growing at
    the other.

    What is checked here is the bookkeeping: that every bond offered has some
    atom to move that is not one of the two the axis runs through. Whether that
    atom is really off the axis depends on where it is, so the matching check
    on real positions is in ``tests/integration/test_real_geometry.py``.
    """

    @pytest.mark.parametrize("direction", ["3prime", "5prime"])
    def test_every_bond_offered_has_an_atom_to_move(self, strand, direction):
        """A bond whose only moving atoms are its own two ends cannot turn anything.

        The axis of a turn runs through the pivot atom and the bond atom, so
        neither of them moves however far the bond is turned.
        """
        chain = strand.chain("aptamer")
        for residue_index in range(chain.n_residues):
            for torsion in chain.residue(residue_index).torsions(direction):
                movable = set(torsion.moving) - {torsion.pivot, torsion.bond}
                assert movable, (
                    f"residue {residue_index} bond {torsion.pivot}-{torsion.bond} "
                    f"read towards the {direction} end has nothing to move but "
                    f"the two atoms of the bond itself"
                )


HOOK = Pose(np.array([[0.0, 0, 0], [1.0, 0, 0], [1.0, 1, 0], [1.0, 2, 0]]))
"""Four atoms in a hook: two along x, then two more turning off along y.

Atoms 0 and 1 define the axis of the turn used below. Atoms 2 and 3 are off
that axis and one ångström apart, so the bond between them is the one a turn
can stretch — and it is stretched exactly when one of the two moves and the
other does not.
"""

HOOK_BONDS = np.array([[0, 1], [1, 2], [2, 3]])
"""Which atoms of :data:`HOOK` are bonded, written out rather than guessed.

A guess from the distances would also call atoms 0 and 2 bonded, since they
are 1.41 Å apart. Saying it explicitly keeps the two tests below about the
turn being checked rather than about the guess.
"""


class TestHandBuiltCases:
    """Small shapes whose answer can be worked out on paper."""

    def test_a_turn_that_moves_both_ends_of_a_bond_keeps_it(self):
        """Both atoms of a bond moving together leaves the bond alone."""
        turned = HOOK.rotate(Torsion(0, 1, AtomRange(2, 4)), np.pi / 2)
        assert_bonds_survive(
            HOOK, turned, "turning both atoms of a bond", pairs=HOOK_BONDS
        )

    def test_a_turn_that_moves_one_end_of_a_bond_is_caught(self):
        """The check has to actually fail when a bond is broken.

        Without this, a bug in the check itself would let every test that uses
        it pass on a broken structure. Here only atom 3 is moved, so it swings
        away from atom 2 and the bond between them goes from 1 Å to 2.24 Å.
        """
        broken = HOOK.rotate(Torsion(0, 1, AtomRange(3, 4)), np.pi / 2)
        with pytest.raises(AssertionError, match="changed a bond length"):
            assert_bonds_survive(
                HOOK, broken, "a deliberately broken turn", pairs=HOOK_BONDS
            )

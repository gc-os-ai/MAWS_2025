"""
Tests for the superposition in :mod:`maws.regrow`.

Growing a strand means building it again from scratch, which produces correct
chemistry in an arbitrary place. :func:`~maws.regrow.kabsch` finds the single
rotation and shift that carries the new build onto the shape already found,
and :func:`~maws.regrow.splice` applies it and then puts the atoms that
already existed back exactly where they were.

The point sets here are four corners of a small box, so each expected answer
is a quarter turn and a whole-number shift that can be checked by eye.
"""

from __future__ import annotations

import numpy as np
import pytest

from maws.errors import ConfigurationError
from maws.pose import Pose
from maws.regrow import kabsch, splice
from maws.values import AtomRange

QUARTER_TURN_ABOUT_Z = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
"""Rotation sending the x direction to y, the y direction to -x, z unchanged.

Written out rather than built with :func:`~maws.pose.rodrigues`, so that a test
recovering it is not comparing a function against itself.
"""

BOX_CORNERS = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 2, 0], [0.0, 0, 3]])
"""Four points spanning all three directions, with a different spread in each.

Three points would be the bare minimum; the fourth keeps them off a plane, and
the unequal spread keeps the covariance matrix free of repeated singular
values, which would leave the fitted rotation only partly determined.
"""


class TestKabsch:
    """Finding the rotation and shift that best match two sets of points."""

    def test_a_pure_shift_is_recovered_with_no_rotation(self):
        """Points moved sideways come back as the identity and that shift."""
        rotation, shift = kabsch(BOX_CORNERS, BOX_CORNERS + [5.0, -2.0, 0.0])
        np.testing.assert_allclose(rotation, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(shift, [5.0, -2.0, 0.0], atol=1e-12)

    def test_a_pure_rotation_is_recovered_with_no_shift(self):
        """Points turned about the origin come back as that turn and no shift."""
        turned = BOX_CORNERS @ QUARTER_TURN_ABOUT_Z.T
        rotation, shift = kabsch(BOX_CORNERS, turned)
        np.testing.assert_allclose(rotation, QUARTER_TURN_ABOUT_Z, atol=1e-12)
        np.testing.assert_allclose(shift, [0.0, 0.0, 0.0], atol=1e-12)

    def test_a_turn_and_a_shift_together_are_recovered(self):
        """Applying what comes back reproduces the target points."""
        target = BOX_CORNERS @ QUARTER_TURN_ABOUT_Z.T + [10.0, 20.0, 30.0]
        rotation, shift = kabsch(BOX_CORNERS, target)
        np.testing.assert_allclose(BOX_CORNERS @ rotation.T + shift, target, atol=1e-12)

    def test_the_result_is_always_a_rotation_and_never_a_reflection(self):
        """Points and their mirror image still give a turn, not a flip.

        Reflecting the z coordinate would match the two sets perfectly, and a
        plain decomposition returns exactly that: a matrix of determinant -1,
        which turns a molecule into its mirror image and makes nonsense of its
        chemistry. The fit corrects for it, so the determinant is +1 and the
        mirrored points are matched as well as a genuine rotation can.
        """
        mirrored = BOX_CORNERS * [1.0, 1.0, -1.0]
        rotation, _ = kabsch(BOX_CORNERS, mirrored)
        assert np.linalg.det(rotation) == pytest.approx(1.0)

    def test_two_point_sets_of_different_sizes_cannot_be_matched(self):
        """Row *i* of one must describe the same atom as row *i* of the other."""
        with pytest.raises(ConfigurationError, match="same shape"):
            kabsch(BOX_CORNERS, BOX_CORNERS[:3])

    @pytest.mark.parametrize("n_points", [0, 1, 2])
    def test_fewer_than_three_points_leave_the_orientation_undecided(self, n_points):
        """Two points can be spun freely about the line joining them.

        With one or two points there is no unique turn that matches them, so
        this is refused rather than returning an arbitrary one.
        """
        few = BOX_CORNERS[:n_points]
        with pytest.raises(ConfigurationError, match="at least three points"):
            kabsch(few, few)


class TestSplice:
    """Fitting a freshly built strand onto the shape an older one had."""

    def _old(self) -> Pose:
        """Return the structure from before the strand grew.

        Five atoms: a three-atom strand at the origin, then a two-atom target
        far away from it.

        Returns
        -------
        maws.pose.Pose
            The shape worth keeping.
        """
        return Pose(
            np.array(
                [[0.0, 0, 0], [1.0, 0, 0], [0.0, 1, 0], [20.0, 20, 20], [21.0, 20, 20]]
            )
        )

    def _fresh(self) -> Pose:
        """Return the structure the builder produced after the strand grew.

        Eight atoms: the same three strand atoms as :meth:`_old` but turned a
        quarter turn about z and moved 100 along x, then three new atoms for
        the residue just added, then the two target atoms wherever the builder
        happened to leave them.

        Returns
        -------
        maws.pose.Pose
            Correct internal geometry, arbitrary overall placement.
        """
        return Pose(
            np.array(
                [
                    [100.0, 0, 0],
                    [100.0, 1, 0],
                    [99.0, 0, 0],
                    [101.0, 0, 0],
                    [102.0, 0, 0],
                    [101.0, 1, 0],
                    [0.0, 0, 0],
                    [0.0, 0, 0],
                ]
            )
        )

    def _splice(self, **overrides) -> Pose:
        """Fit :meth:`_fresh` onto :meth:`_old`, with any argument replaced.

        The strand here is one residue that grew into two. Its one old residue
        is the three atoms at the front of both structures, so there is a
        single entry in `matches` and the fit is computed from it.

        Parameters
        ----------
        **overrides
            Argument name to replacement value.

        Returns
        -------
        maws.pose.Pose
            Positions for the grown structure.
        """
        arguments = {
            "fresh_chain": AtomRange(0, 6),
            "matches": [(AtomRange(0, 3), AtomRange(0, 3))],
            "anchor": 0,
            "others": [(AtomRange(6, 8), AtomRange(3, 5))],
        }
        arguments.update(overrides)
        return splice(self._fresh(), self._old(), **arguments)

    def test_atoms_that_already_existed_land_on_their_old_positions_exactly(self):
        """Those rows are copied over, not fitted, so they match to the last bit."""
        spliced = self._splice()
        np.testing.assert_array_equal(
            spliced.xyz[AtomRange(0, 3).as_slice()], self._old().xyz[:3]
        )

    def test_the_new_residue_is_carried_along_by_the_same_fit(self):
        """The added atoms are turned back a quarter turn and moved 100 along -x.

        Atom 3 sat one step along x from the strand's first atom in the fresh
        build; undoing the quarter turn puts it one step along -y instead.
        """
        np.testing.assert_allclose(
            self._splice().xyz[3:6],
            [[0.0, -1, 0], [0.0, -2, 0], [1.0, -1, 0]],
            atol=1e-12,
        )

    def test_the_shape_of_the_new_residue_survives_the_fit(self):
        """Fitting turns and shifts the strand, so distances within it do not change."""
        fresh = self._fresh()
        spliced = self._splice()
        assert np.linalg.norm(spliced.xyz[3] - spliced.xyz[5]) == pytest.approx(
            float(np.linalg.norm(fresh.xyz[3] - fresh.xyz[5]))
        )

    def test_a_chain_that_did_not_grow_is_copied_into_its_new_span(self):
        """The target keeps its old positions, at whatever indices it now has."""
        np.testing.assert_array_equal(self._splice().xyz[6:8], self._old().xyz[3:5])

    def test_the_result_covers_the_grown_structure(self):
        """There is one row per atom of the new build, not of the old one."""
        assert self._splice().n_atoms == 8

    @pytest.mark.parametrize("n_kept", [0, 1, 2])
    def test_too_few_kept_atoms_leaves_the_fresh_positions_alone(self, n_kept):
        """With no orientation to fit, the build is used exactly as it came."""
        spliced = self._splice(matches=[(AtomRange(0, n_kept), AtomRange(0, n_kept))])
        np.testing.assert_array_equal(spliced.xyz[0:6], self._fresh().xyz[0:6])

    def test_a_chain_that_did_not_grow_is_copied_even_with_nothing_kept(self):
        """The other chains are placed whether or not a fit happens."""
        spliced = self._splice(matches=[(AtomRange(0, 0), AtomRange(0, 0))])
        np.testing.assert_array_equal(spliced.xyz[6:8], self._old().xyz[3:5])

    def test_an_empty_strand_has_nothing_to_fit_onto(self):
        """The first residue of a run is placed by the sampler, not by a fit.

        There is no earlier shape to keep, which the caller signals with an
        anchor of -1, so the fresh build is left exactly as it came.
        """
        spliced = self._splice(matches=[], anchor=-1)
        np.testing.assert_array_equal(spliced.xyz[0:6], self._fresh().xyz[0:6])

    def test_only_the_anchor_residue_decides_where_the_new_one_goes(self):
        """Residues further from the join are copied back, but do not steer the fit.

        A fresh build comes out fully extended while the shape being kept has
        been folded, so no single turn and shift carries the whole strand onto
        it. Fitting from the neighbouring residue alone keeps the new residue
        correctly joined; fitting from all of them would smear that error over
        the one atom group being positioned.

        The second entry here is deliberately mismatched — its old positions
        are nowhere near its fresh ones — and the answer must not change.
        """
        with_distraction = self._splice(
            matches=[
                (AtomRange(0, 3), AtomRange(0, 3)),
                (AtomRange(3, 5), AtomRange(3, 5)),
            ],
            anchor=0,
        )
        np.testing.assert_allclose(
            with_distraction.xyz[5], self._splice().xyz[5], atol=1e-12
        )

    def test_the_kept_part_must_be_the_same_size_in_both_structures(self):
        """Different counts mean the two spans are not the same atoms."""
        with pytest.raises(ConfigurationError, match="must name the same atoms"):
            self._splice(matches=[(AtomRange(0, 3), AtomRange(0, 2))])

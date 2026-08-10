"""
Tests for :mod:`maws.pose`.

A pose is a set of atom positions that never changes: every method returns a
new pose and leaves the one it was called on exactly as it was. Most of what
follows checks that promise one method at a time, because everything else in
MAWS is written assuming it.

The poses here are three or four atoms on the axes, so each expected position
is a quarter turn worked out by hand rather than a number copied out of a
previous run.
"""

from __future__ import annotations

import numpy as np
import pytest

from maws.errors import ConfigurationError
from maws.libraries import rna
from maws.pose import ChainView, Pose, rodrigues
from maws.sampling import Placement
from maws.topology import BuiltSystem
from maws.values import AtomRange, NucleotideSequence, Torsion

ELBOW = Torsion(0, 1, AtomRange(2, 3))
"""The bond from atom 0 to atom 1 of :func:`corner`, swinging atom 2.

Its axis runs along z, so a quarter turn about it moves atom 2 from the x
direction to the y direction.
"""


@pytest.fixture
def corner() -> Pose:
    """Three atoms in an L: the corner, one atom up z, one atom out along x."""
    return Pose(np.array([[0.0, 0, 0], [0, 0, 1], [1, 0, 1]]))


@pytest.fixture
def two_chains() -> Pose:
    """Six atoms in two chains of three, one at the origin and one far off."""
    return Pose(
        np.array(
            [
                [0.0, 0, 0],
                [1.0, 0, 0],
                [0.0, 1, 0],
                [50.0, 0, 0],
                [51.0, 0, 0],
                [50.0, 1, 0],
            ]
        )
    )


FIRST_CHAIN = AtomRange(0, 3)
"""The atoms of the first chain of :func:`two_chains`."""

SECOND_CHAIN = AtomRange(3, 6)
"""The atoms of the second chain of :func:`two_chains`."""


class TestPoseConstruction:
    """Making a pose out of an array of positions."""

    def test_a_pose_knows_how_many_atoms_it_holds(self):
        """One row per atom, so the atom count is the number of rows."""
        assert Pose(np.zeros((4, 3))).n_atoms == 4

    def test_a_pose_is_as_long_as_its_atom_count(self):
        """``len(pose)`` and ``pose.n_atoms`` are the same number."""
        assert len(Pose(np.zeros((4, 3)))) == 4

    @pytest.mark.parametrize("shape", [(4,), (4, 2), (2, 3, 3)])
    def test_positions_must_be_one_row_of_three_numbers_per_atom(self, shape):
        """Anything that is not ``(N, 3)`` is a mistake, caught where it is made."""
        with pytest.raises(ConfigurationError, match="must be shaped"):
            Pose(np.zeros(shape))

    def test_a_design_with_no_atoms_is_allowed(self):
        """An empty pose is what a search starts from before anything is built."""
        assert Pose(np.zeros((0, 3))).n_atoms == 0

    def test_the_positions_cannot_be_written_to(self):
        """Writing to the array is refused, which is what makes a pose a value."""
        assert Pose(np.zeros((4, 3))).xyz.flags.writeable is False

    def test_the_array_handed_in_is_left_alone(self):
        """Building a pose does not take ownership of the caller's array."""
        source = np.zeros((2, 3))
        Pose(source)
        source[0, 0] = 99.0
        assert source[0, 0] == 99.0

    def test_changing_the_array_afterwards_does_not_change_the_pose(self):
        """A pose holds its own copy, so nothing can edit it from outside."""
        source = np.zeros((2, 3))
        pose = Pose(source)
        source[0, 0] = 99.0
        assert pose.xyz[0, 0] == 0.0

    def test_a_pose_built_from_a_slice_is_still_its_own(self):
        """A view shares memory with what it looks at, so it has to be copied.

        Without a copy the pose would change whenever the larger array did,
        which would quietly break the promise that a pose never changes.
        """
        base = np.zeros((4, 3))
        pose = Pose(base[:2])
        base[0, 0] = 99.0
        assert pose.xyz[0, 0] == 0.0

    def test_positions_can_be_given_as_plain_lists(self):
        """Anything shaped like a list of triples is accepted."""
        assert Pose([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]).n_atoms == 2

    def test_whole_numbers_are_stored_as_fractions(self):
        """An array of integers does not make later moves round to integers."""
        assert Pose(np.array([[1, 2, 3]])).xyz.dtype == np.float64

    def test_a_pose_says_how_big_it_is_when_printed(self):
        """The atom count is the one thing worth seeing at a prompt."""
        assert repr(Pose(np.zeros((7, 3)))) == "<Pose 7 atoms>"


class TestRodrigues:
    """The matrix that turns points about a line through the origin."""

    def test_a_quarter_turn_about_z_sends_x_to_y_by_the_right_hand_rule(self):
        """A positive angle turns x towards y, not away from it."""
        matrix = rodrigues(np.array([0.0, 0, 1]), np.pi / 2)
        np.testing.assert_allclose(
            np.array([1.0, 0, 0]) @ matrix.T, [0.0, 1.0, 0.0], atol=1e-12
        )

    def test_a_full_turn_brings_a_point_back_to_itself(self):
        """Turning by 2π is the same as not turning, to floating-point accuracy."""
        matrix = rodrigues(np.array([1.0, 2, 3]), 2 * np.pi)
        assert np.allclose(matrix, np.eye(3))

    def test_turning_by_nothing_is_exactly_the_identity(self):
        """At angle zero the two terms vanish outright, so this is exact."""
        np.testing.assert_array_equal(rodrigues(np.array([1.0, 2, 3]), 0.0), np.eye(3))

    def test_the_length_of_the_axis_is_ignored(self):
        """The axis can be the difference between two atom positions as it comes."""
        short = rodrigues(np.array([0.0, 0, 1]), 0.7)
        long = rodrigues(np.array([0.0, 0, 40]), 0.7)
        assert np.allclose(short, long)

    def test_turning_the_other_way_undoes_the_turn(self):
        """The matrix for -θ is the inverse of the matrix for θ."""
        axis = np.array([1.0, 1, 0])
        assert np.allclose(rodrigues(axis, 0.4) @ rodrigues(axis, -0.4), np.eye(3))

    def test_an_axis_of_zero_length_leaves_no_line_to_turn_about(self):
        """Two atoms at the same position define no bond to turn."""
        with pytest.raises(ConfigurationError, match="zero-length axis"):
            rodrigues(np.zeros(3), np.pi / 2)


class TestPoseRotate:
    """Turning one bond, swinging the atoms attached to one side of it."""

    def test_a_quarter_turn_swings_the_moving_atom_round_the_axis(self, corner):
        """Atom 2 starts out along x and a quarter turn about z puts it along y."""
        turned = corner.rotate(ELBOW, np.pi / 2)
        np.testing.assert_allclose(turned.xyz[2], [0.0, 1.0, 1.0], atol=1e-12)

    def test_a_full_turn_returns_every_atom_to_where_it_was(self, corner):
        """Turning by 2π changes nothing, to floating-point accuracy."""
        assert np.allclose(corner.rotate(ELBOW, 2 * np.pi).xyz, corner.xyz)

    def test_turning_one_way_then_back_returns_every_atom_to_where_it_was(self, corner):
        """+θ followed by -θ is the same as doing nothing."""
        there_and_back = corner.rotate(ELBOW, 0.9).rotate(ELBOW, -0.9)
        assert np.allclose(there_and_back.xyz, corner.xyz)

    def test_atoms_outside_the_moving_range_never_move(self, corner):
        """Only the atoms the bond names swing; the rest of the design sits still."""
        turned = corner.rotate(ELBOW, 1.3)
        np.testing.assert_array_equal(turned.xyz[:2], corner.xyz[:2])

    def test_the_pivot_atom_stays_where_it_is(self, corner):
        """The turn is centred on the pivot, so that atom is the one fixed point."""
        turned = corner.rotate(ELBOW, 1.3)
        np.testing.assert_array_equal(turned.xyz[ELBOW.pivot], corner.xyz[ELBOW.pivot])

    def test_a_bond_between_two_atoms_at_the_same_position_is_refused(self):
        """There is no axis to turn about, so this is reported rather than nan."""
        flat = Pose(np.array([[0.0, 0, 0], [0.0, 0, 0], [1.0, 0, 0]]))
        with pytest.raises(ConfigurationError, match="zero-length axis"):
            flat.rotate(ELBOW, np.pi / 2)


class TestPoseRotateAll:
    """Turning several bonds in one call."""

    def _elbows(self) -> tuple[Torsion, Torsion]:
        """Return two bonds of :func:`four_atoms`, the first moving the second.

        The first bond swings atom 2, which is one of the two atoms defining
        the second bond's axis. Turning the first therefore changes what the
        second turns about.

        Returns
        -------
        tuple of maws.values.Torsion
            The bond from atom 0 to atom 1, then the one from atom 1 to atom 2.
        """
        return (
            Torsion(0, 1, AtomRange(2, 3)),
            Torsion(1, 2, AtomRange(3, 4)),
        )

    @pytest.fixture
    def four_atoms(self) -> Pose:
        """Four atoms in a staircase, so that two bonds in a row can be turned."""
        return Pose(np.array([[0.0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1]]))

    def test_the_turns_are_applied_one_after_another(self, four_atoms):
        """The result is the same as calling ``rotate`` twice in the same order."""
        first, second = self._elbows()
        together = four_atoms.rotate_all([first, second], [0.4, 1.1])
        separately = four_atoms.rotate(first, 0.4).rotate(second, 1.1)
        np.testing.assert_array_equal(together.xyz, separately.xyz)

    def test_each_turn_works_from_where_the_turn_before_it_left_the_atoms(
        self, four_atoms
    ):
        """A bond turns about wherever an earlier bond has already put its atoms.

        The same two turns in the other order therefore give a different shape,
        because the second bond's axis is only where it is once the first has
        been made.
        """
        first, second = self._elbows()
        forwards = four_atoms.rotate_all([first, second], [0.4, 1.1])
        backwards = four_atoms.rotate_all([second, first], [1.1, 0.4])
        assert not np.allclose(forwards.xyz, backwards.xyz)

    def test_turning_no_bonds_at_all_changes_nothing(self, four_atoms):
        """An empty list of bonds is allowed and leaves every atom alone."""
        np.testing.assert_array_equal(four_atoms.rotate_all([], []).xyz, four_atoms.xyz)

    def test_there_must_be_exactly_one_angle_per_bond(self, four_atoms):
        """A short list of angles means a mistake, not bonds to leave unturned."""
        first, second = self._elbows()
        with pytest.raises(ConfigurationError, match="one angle per bond"):
            four_atoms.rotate_all([first, second], [0.4])


class TestPoseTranslate:
    """Sliding a run of atoms without turning them."""

    def test_the_named_atoms_all_move_by_the_same_amount(self, two_chains):
        """Every atom of the run shifts along each axis by what it was given."""
        moved = two_chains.translate(FIRST_CHAIN, np.array([0.0, 0, 10]))
        np.testing.assert_allclose(
            moved.atoms(FIRST_CHAIN), two_chains.atoms(FIRST_CHAIN) + [0, 0, 10]
        )

    def test_atoms_outside_the_run_are_left_exactly_where_they_were(self, two_chains):
        """The other chain is copied across untouched, not recomputed."""
        moved = two_chains.translate(FIRST_CHAIN, np.array([0.0, 0, 10]))
        np.testing.assert_array_equal(
            moved.atoms(SECOND_CHAIN), two_chains.atoms(SECOND_CHAIN)
        )


class TestPosePlace:
    """Putting a whole chain somewhere a sampler proposed."""

    def test_the_chain_is_slid_first_and_then_turned_about_its_first_atom(
        self, two_chains
    ):
        """A shift of 10 along x, then a quarter turn about z through the new corner.

        Atom 0 lands at ``[10, 0, 0]`` and stays there, because the turn is
        centred on it. Atom 1 sat one step along x from it and ends one step
        along y; atom 2 sat one step along y and ends one step back along x.
        """
        placement = Placement(
            position=np.array([10.0, 0, 0]),
            axis=np.array([0.0, 0, 1]),
            angle=np.pi / 2,
        )
        placed = two_chains.place(FIRST_CHAIN, placement)
        np.testing.assert_allclose(
            placed.atoms(FIRST_CHAIN),
            [[10.0, 0, 0], [10.0, 1, 0], [9.0, 0, 0]],
            atol=1e-12,
        )

    def test_every_other_chain_stays_exactly_where_it_was(self, two_chains):
        """Placing one chain is not allowed to disturb the rest of the design."""
        placement = Placement(
            position=np.array([10.0, 0, 0]),
            axis=np.array([0.0, 0, 1]),
            angle=np.pi / 2,
        )
        placed = two_chains.place(FIRST_CHAIN, placement)
        np.testing.assert_array_equal(
            placed.atoms(SECOND_CHAIN), two_chains.atoms(SECOND_CHAIN)
        )

    def test_a_chain_can_be_named_by_its_view_instead_of_its_span(
        self, two_residue_system: BuiltSystem
    ):
        """Passing a chain view moves that chain and leaves the target alone."""
        pose = two_residue_system.pose
        ligand = two_residue_system.chain("ligand")
        placement = Placement(
            position=np.array([0.0, 0, 7]),
            axis=np.array([0.0, 0, 1]),
            angle=0.0,
        )
        placed = pose.place(two_residue_system.chain("aptamer"), placement)
        np.testing.assert_array_equal(placed.atoms(ligand), pose.atoms(ligand))


class TestPoseJittered:
    """Nudging every atom by its own small displacement."""

    def test_each_atom_moves_by_its_own_row(self, corner):
        """The offsets are added position by position, not applied as one shift."""
        offsets = np.array([[1.0, 0, 0], [0.0, 2, 0], [0.0, 0, 3]])
        np.testing.assert_allclose(corner.jittered(offsets).xyz, corner.xyz + offsets)

    def test_a_random_nudge_moves_every_atom_a_little(self, corner, rng):
        """A jitter drawn at random leaves no atom exactly where it was."""
        offsets = rng.normal(scale=0.1, size=corner.xyz.shape)
        assert not np.any(corner.jittered(offsets).xyz == corner.xyz)

    def test_there_must_be_one_offset_per_atom(self, corner):
        """An offset array of the wrong size cannot be lined up with the atoms."""
        with pytest.raises(ConfigurationError, match="offsets must be shaped"):
            corner.jittered(np.zeros((2, 3)))


class TestPoseWithSpan:
    """Giving one run of atoms entirely new positions."""

    def test_the_named_atoms_take_the_positions_given(self, two_chains):
        """The replacement rows are written in as they are."""
        replacement = np.array([[7.0, 7, 7], [8.0, 8, 8], [9.0, 9, 9]])
        updated = two_chains.with_span(FIRST_CHAIN, replacement)
        np.testing.assert_array_equal(updated.atoms(FIRST_CHAIN), replacement)

    def test_atoms_outside_the_run_are_left_exactly_where_they_were(self, two_chains):
        """Replacing one chain does not disturb any other."""
        replacement = np.zeros((3, 3))
        updated = two_chains.with_span(FIRST_CHAIN, replacement)
        np.testing.assert_array_equal(
            updated.atoms(SECOND_CHAIN), two_chains.atoms(SECOND_CHAIN)
        )

    def test_there_must_be_one_row_per_atom_of_the_run(self, two_chains):
        """Two rows cannot fill a three-atom run, and are refused rather than padded."""
        with pytest.raises(ConfigurationError, match="rows of positions"):
            two_chains.with_span(FIRST_CHAIN, np.zeros((2, 3)))


class TestEveryMoveReturnsANewPose:
    """The promise the whole module rests on: nothing changes a pose in place."""

    def test_turning_a_bond_leaves_the_pose_it_was_called_on_alone(self, corner):
        """``rotate`` gives back a new pose."""
        before = corner.xyz.copy()
        turned = corner.rotate(ELBOW, 1.3)
        assert turned is not corner
        np.testing.assert_array_equal(corner.xyz, before)

    def test_turning_several_bonds_leaves_the_pose_it_was_called_on_alone(self, corner):
        """``rotate_all`` gives back a new pose."""
        before = corner.xyz.copy()
        turned = corner.rotate_all([ELBOW], [1.3])
        assert turned is not corner
        np.testing.assert_array_equal(corner.xyz, before)

    def test_sliding_atoms_leaves_the_pose_it_was_called_on_alone(self, two_chains):
        """``translate`` gives back a new pose."""
        before = two_chains.xyz.copy()
        moved = two_chains.translate(FIRST_CHAIN, np.array([1.0, 2, 3]))
        assert moved is not two_chains
        np.testing.assert_array_equal(two_chains.xyz, before)

    def test_placing_a_chain_leaves_the_pose_it_was_called_on_alone(self, two_chains):
        """``place`` gives back a new pose."""
        before = two_chains.xyz.copy()
        placed = two_chains.place(
            FIRST_CHAIN,
            Placement(
                position=np.array([1.0, 2, 3]),
                axis=np.array([0.0, 0, 1]),
                angle=0.8,
            ),
        )
        assert placed is not two_chains
        np.testing.assert_array_equal(two_chains.xyz, before)

    def test_nudging_atoms_leaves_the_pose_it_was_called_on_alone(self, corner):
        """``jittered`` gives back a new pose."""
        before = corner.xyz.copy()
        nudged = corner.jittered(np.full((3, 3), 0.05))
        assert nudged is not corner
        np.testing.assert_array_equal(corner.xyz, before)

    def test_replacing_a_run_leaves_the_pose_it_was_called_on_alone(self, two_chains):
        """``with_span`` gives back a new pose."""
        before = two_chains.xyz.copy()
        updated = two_chains.with_span(FIRST_CHAIN, np.zeros((3, 3)))
        assert updated is not two_chains
        np.testing.assert_array_equal(two_chains.xyz, before)


class TestChainView:
    """A named window onto one chain, and where its residues sit inside it."""

    def _chain(self, start: int = 0) -> ChainView:
        """Return a two-nucleotide RNA chain whose atoms begin at `start`.

        Parameters
        ----------
        start : int, default=0
            Index of the chain's first atom within the whole design.

        Returns
        -------
        ChainView
            A chain of G followed by A, 66 atoms in all.
        """
        return ChainView.build(
            "aptamer",
            AtomRange(start, start + 66),
            NucleotideSequence.parse("G A"),
            rna(),
        )

    def test_a_two_nucleotide_strand_has_two_residues(self):
        """One residue per nucleotide written."""
        assert self._chain().n_residues == 2

    def test_the_ends_of_a_strand_take_their_end_specific_forms(self):
        """G at the 5' end is G5 and A at the 3' end is A3."""
        assert self._chain().canonical == ("G5", "A3")

    def test_each_residue_starts_where_the_one_before_it_ended(self):
        """G5 has 32 atoms, so the second residue starts at offset 32."""
        assert self._chain().residue_offsets == (0, 32)

    def test_the_first_residue_covers_its_own_atom_count(self):
        """G5 is 32 atoms, counted from the chain's first atom."""
        assert self._chain().residue(0).span == AtomRange(0, 32)

    def test_the_last_residue_runs_to_the_end_of_the_chain(self):
        """A3 is 34 atoms, and 32 plus 34 is the chain's 66."""
        assert self._chain().residue(1).span == AtomRange(32, 66)

    def test_a_residue_span_is_counted_from_the_start_of_the_design(self):
        """A chain that does not start at atom 0 shifts all of its residues."""
        assert self._chain(start=100).residue(1).span == AtomRange(132, 166)

    def test_counting_back_from_the_three_prime_end(self):
        """``-1`` names the last residue, as with a Python list."""
        assert self._chain().residue(-1).name == "A3"

    def test_every_residue_can_be_asked_for_at_once(self):
        """``residues`` lists them 5' end first."""
        assert [r.name for r in self._chain().residues()] == ["G5", "A3"]

    @pytest.mark.parametrize("index", [2, -3])
    def test_a_residue_index_off_either_end_is_rejected(self, index):
        """Asking for a residue a two-residue chain does not have is a mistake."""
        with pytest.raises(ConfigurationError, match="does not name one"):
            self._chain().residue(index)


class TestResidueViewTorsion:
    """A residue's turnable bonds, renumbered for the whole design."""

    def _chain(self, start: int = 100) -> ChainView:
        """Return a two-nucleotide RNA chain beginning at atom `start`.

        Parameters
        ----------
        start : int, default=100
            Index of the chain's first atom within the whole design. Non-zero
            by default, so that a forgotten offset shows up.

        Returns
        -------
        ChainView
            A chain of G followed by A, 66 atoms in all.
        """
        return ChainView.build(
            "aptamer",
            AtomRange(start, start + 66),
            NucleotideSequence.parse("G A"),
            rna(),
        )

    def test_a_bond_is_numbered_from_the_start_of_the_design(self):
        """The residue tables count from a residue; a bond counts from atom 0."""
        assert self._chain().residue(0).torsion(0) == Torsion(
            pivot=100, bond=101, moving=AtomRange(101, 166)
        )

    def test_the_same_bond_of_an_unshifted_chain_starts_at_zero(self):
        """With the chain at the start of the design, the offset is nothing."""
        assert self._chain(start=0).residue(0).torsion(0) == Torsion(
            pivot=0, bond=1, moving=AtomRange(1, 66)
        )

    def test_a_later_residue_is_offset_by_the_residues_before_it(self):
        """The second residue's atoms begin 32 further along than the chain does."""
        assert self._chain().residue(1).torsion(0) == Torsion(
            pivot=132, bond=135, moving=AtomRange(135, 166)
        )

    def test_a_backbone_bond_swings_everything_to_the_three_prime_end(self):
        """The moving run reaches the last atom of the chain, not of the residue."""
        assert self._chain().residue(0).torsion(0).moving.stop == 166

    def test_the_five_prime_form_swings_the_other_part_of_the_strand(self):
        """The two directions turn the same bond but move opposite parts."""
        residue = self._chain().residue(1)
        assert residue.torsion(0).moving == AtomRange(135, 166)
        assert residue.torsion(0, "5prime").moving == AtomRange(100, 135)

    def test_the_five_prime_form_runs_its_axis_the_other_way(self):
        """Swapping the pivot and bond atoms is what reverses the turn."""
        residue = self._chain().residue(1)
        forwards = residue.torsion(0)
        backwards = residue.torsion(0, "5prime")
        assert (forwards.pivot, forwards.bond) == (backwards.bond, backwards.pivot)

    def test_a_bond_index_the_residue_does_not_have_is_rejected(self):
        """G5 has four turnable bonds, so there is no fifth to ask for."""
        with pytest.raises(ConfigurationError, match="does not name one"):
            self._chain().residue(0).torsion(4)

    def test_a_direction_that_is_neither_end_is_rejected(self):
        """Only the two ends of a strand can be named, and a typo is caught here."""
        with pytest.raises(ConfigurationError, match="must be '3prime' or '5prime'"):
            self._chain().residue(0).torsion(0, "sideways")


class TestResidueViewTorsions:
    """Asking for several of a residue's bonds at once."""

    def _residue(self):
        """Return the first residue of a lone guanine chain.

        Returns
        -------
        maws.pose.ResidueView
            A window onto GN, which has four turnable bonds.
        """
        return ChainView.build(
            "aptamer", AtomRange(0, 33), NucleotideSequence.parse("G"), rna()
        ).residue(0)

    def test_all_the_bonds_come_back_when_no_limit_is_given(self):
        """A lone guanine declares four turnable bonds."""
        assert len(self._residue().torsions()) == 4

    def test_a_limit_returns_that_many_bonds(self):
        """The first three are the ones nearest the residue's own numbering start."""
        assert len(self._residue().torsions(limit=3)) == 3

    def test_a_limit_larger_than_the_residue_has_returns_what_there_is(self):
        """Asking for ten bonds from a residue with four gives four."""
        assert len(self._residue().torsions(limit=10)) == 4

    def test_the_direction_is_passed_on_to_every_bond(self):
        """One direction argument covers the whole set."""
        forwards = self._residue().torsions()
        backwards = self._residue().torsions("5prime")
        assert all(f.pivot == b.bond for f, b in zip(forwards, backwards, strict=True))

"""
Tests for :mod:`maws.relax`.

Shaking a structure loose means nudging every atom at random and settling the
result, over and over. The scorer used throughout is
:class:`~maws.energy.StubEnergy`, whose settling step hands the positions back
untouched, so what is left to check is the nudging: how far it moves atoms, how
many times, and whether the same generator gives the same answer twice.
"""

from __future__ import annotations

import numpy as np
import pytest

from maws.energy import Relaxed, StubEnergy
from maws.errors import ConfigurationError
from maws.relax import perturb_and_minimize


@pytest.fixture
def scorer(one_residue_system) -> StubEnergy:
    """A scorer that pushes the strand and the target apart."""
    return StubEnergy(
        one_residue_system.chain("aptamer").span,
        one_residue_system.chain("ligand").span,
    )


@pytest.fixture
def start(one_residue_system):
    """The positions a shake starts from."""
    return one_residue_system.pose


class TestPerturbAndMinimize:
    """Nudging every atom and settling, several rounds over."""

    def test_it_hands_back_positions_and_the_energy_of_those_positions(
        self, start, scorer, rng
    ):
        """The result is a settled structure together with its score."""
        result = perturb_and_minimize(start, scorer, iterations=3, rng=rng)
        assert isinstance(result, Relaxed)
        assert result.energy == pytest.approx(scorer.evaluate(result.pose))

    def test_no_atom_is_lost_along_the_way(self, start, scorer, rng):
        """Every atom is nudged, and none is added or dropped."""
        result = perturb_and_minimize(start, scorer, iterations=5, rng=rng)
        assert result.pose.n_atoms == start.n_atoms

    def test_the_same_seed_gives_the_same_shape(self, start, scorer):
        """A run can be repeated exactly by fixing the seed."""
        first = perturb_and_minimize(
            start, scorer, iterations=4, rng=np.random.default_rng(11)
        )
        second = perturb_and_minimize(
            start, scorer, iterations=4, rng=np.random.default_rng(11)
        )
        assert first.pose.xyz == pytest.approx(second.pose.xyz)
        assert first.energy == pytest.approx(second.energy)

    def test_different_seeds_give_different_shapes(self, start, scorer):
        """The nudges really are random, not a fixed pattern."""
        first = perturb_and_minimize(
            start, scorer, iterations=4, rng=np.random.default_rng(11)
        )
        second = perturb_and_minimize(
            start, scorer, iterations=4, rng=np.random.default_rng(12)
        )
        assert not np.allclose(first.pose.xyz, second.pose.xyz)

    def test_asking_for_no_rounds_gives_back_the_starting_shape(
        self, start, scorer, rng
    ):
        """With nothing to do it scores what it was given and stops there."""
        result = perturb_and_minimize(start, scorer, iterations=0, rng=rng)
        assert result.pose is start
        assert result.energy == pytest.approx(scorer.evaluate(start))

    def test_a_nudge_of_zero_leaves_every_atom_where_it_was(self, start, scorer, rng):
        """The rounds still run; they simply move nothing."""
        result = perturb_and_minimize(start, scorer, size=0.0, iterations=5, rng=rng)
        assert result.pose.xyz == pytest.approx(start.xyz)

    def test_no_atom_moves_further_than_the_size_allows_in_one_round(
        self, start, scorer, rng
    ):
        """Each nudge is drawn evenly from ``-size`` to ``+size`` per axis, so
        one round of 0.25 Å cannot shift an atom further than that on any
        axis."""
        result = perturb_and_minimize(start, scorer, size=0.25, iterations=1, rng=rng)
        assert np.abs(result.pose.xyz - start.xyz).max() <= 0.25

    def test_a_larger_nudge_moves_the_structure_further(self, start, scorer):
        """Size is the lever for escaping a poor arrangement."""
        small = perturb_and_minimize(
            start, scorer, size=0.01, iterations=1, rng=np.random.default_rng(5)
        )
        large = perturb_and_minimize(
            start, scorer, size=1.0, iterations=1, rng=np.random.default_rng(5)
        )
        moved_a_little = np.abs(small.pose.xyz - start.xyz).max()
        moved_a_lot = np.abs(large.pose.xyz - start.xyz).max()
        assert moved_a_lot > moved_a_little

    def test_the_starting_positions_are_left_as_they_were(self, start, scorer, rng):
        """Every round produces new positions rather than editing the ones it
        was handed, so the shape a search started from survives the shake."""
        before = start.xyz.copy()
        perturb_and_minimize(start, scorer, size=0.5, iterations=5, rng=rng)
        assert start.xyz == pytest.approx(before)

    def test_a_negative_nudge_is_rejected(self, start, scorer):
        """A distance below zero has no meaning."""
        with pytest.raises(ConfigurationError, match="size must not be negative"):
            perturb_and_minimize(start, scorer, size=-0.1)

    def test_a_negative_number_of_rounds_is_rejected(self, start, scorer):
        """Minus three rounds is a mistake, not an instruction to do nothing."""
        with pytest.raises(ConfigurationError, match="iterations must not be negative"):
            perturb_and_minimize(start, scorer, iterations=-3)

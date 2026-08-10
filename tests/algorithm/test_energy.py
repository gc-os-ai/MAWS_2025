"""
Tests for :mod:`maws.energy`.

Only the scorers that are pure arithmetic are covered here: the stand-in, the
one that adds several scorers together, and the pair of numbers they hand back.
Scoring with real physics needs OpenMM installed and belongs with the
integration tests.

The stand-in's energy is a sum of ``scale / distance²`` over every pair of atoms
taken one from each of two groups, so every number in this file is worked out by
hand in the test that asserts it.
"""

from __future__ import annotations

import numpy as np
import pytest

from maws.energy import CompositeEnergy, EnergyModel, Relaxed, StubEnergy
from maws.errors import ConfigurationError
from maws.pose import Pose
from maws.values import AtomRange

FIRST = AtomRange(0, 1)
"""The single atom making up the first group in most of these tests."""

SECOND = AtomRange(1, 2)
"""The single atom making up the second group in most of these tests."""


def pair(separation: float, *, offset: float = 0.0) -> Pose:
    """Return two atoms `separation` ångström apart along the x axis.

    Parameters
    ----------
    separation : float
        How far apart to put them, in ångström.
    offset : float, default=0.0
        How far along x to shift the whole pair, in ångström. Moving the pair
        without changing the gap must not change the score.

    Returns
    -------
    maws.pose.Pose
        A two-atom pose.
    """
    return Pose(np.array([[offset, 0.0, 0.0], [offset + separation, 0.0, 0.0]]))


class TestStubEnergyDistance:
    """How the stand-in score depends on how far apart the two groups are."""

    def test_doubling_the_distance_quarters_the_score(self):
        """Two atoms 5 Å apart score 1/25 = 0.04; at 10 Å they score
        1/100 = 0.01."""
        model = StubEnergy(FIRST, SECOND)
        assert model.evaluate(pair(5.0)) == pytest.approx(0.04)
        assert model.evaluate(pair(10.0)) == pytest.approx(0.01)

    def test_moving_the_groups_apart_lowers_the_score(self):
        """A lower score is a better arrangement, so separation is rewarded."""
        model = StubEnergy(FIRST, SECOND)
        assert model.evaluate(pair(12.0)) < model.evaluate(pair(6.0))

    def test_moving_the_groups_together_raises_the_score(self):
        """Atoms pushed into each other are the arrangement to avoid."""
        model = StubEnergy(FIRST, SECOND)
        assert model.evaluate(pair(2.0)) > model.evaluate(pair(6.0))

    def test_two_atoms_in_the_same_place_give_a_large_number(self):
        """Squared distances are floored at 1e-6, so coincident atoms score
        1e6 rather than dividing by zero."""
        model = StubEnergy(FIRST, SECOND)
        score = model.evaluate(pair(0.0))
        assert np.isfinite(score)
        assert score == pytest.approx(1e6)

    def test_only_the_gap_matters_not_where_the_pair_sits(self):
        """The same two atoms moved 500 Å along x score exactly the same."""
        model = StubEnergy(FIRST, SECOND)
        assert model.evaluate(pair(4.0)) == pytest.approx(
            model.evaluate(pair(4.0, offset=500.0))
        )

    def test_every_pair_of_atoms_across_the_groups_contributes(self):
        """Two atoms placed either side of a third, each 5 Å from it, score
        1/25 twice over: 0.08."""
        model = StubEnergy(AtomRange(0, 2), AtomRange(2, 3))
        pose = Pose(np.array([[-5.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 0.0, 0.0]]))
        assert model.evaluate(pose) == pytest.approx(0.08)


class TestStubEnergyScale:
    """The multiplier that says how strongly the two groups repel."""

    @pytest.mark.parametrize("scale", [0.5, 1.0, 3.0])
    def test_the_scale_multiplies_the_whole_score(self, scale):
        """Two atoms 5 Å apart score ``scale`` times 0.04."""
        model = StubEnergy(FIRST, SECOND, scale=scale)
        assert model.evaluate(pair(5.0)) == pytest.approx(0.04 * scale)


class TestStubEnergyMinimize:
    """What settling a structure does when the scorer has no forces in it."""

    def test_the_positions_come_back_as_they_went_in(self):
        """Nothing is adjusted, so this is the pose that was handed over."""
        model = StubEnergy(FIRST, SECOND)
        pose = pair(5.0)
        assert model.minimize(pose).pose is pose

    def test_the_energy_reported_is_the_energy_of_that_pose(self):
        """Settling and scoring agree, so the two can be compared freely."""
        model = StubEnergy(FIRST, SECOND)
        pose = pair(5.0)
        assert model.minimize(pose).energy == pytest.approx(model.evaluate(pose))


class TestRelaxed:
    """The pair of values a settling step hands back."""

    def test_it_unpacks_as_energy_then_pose(self):
        """``energy, pose = model.minimize(pose)`` reads left to right."""
        model = StubEnergy(FIRST, SECOND)
        started_from = pair(5.0)
        energy, pose = model.minimize(started_from)
        assert energy == pytest.approx(0.04)
        assert pose is started_from

    def test_the_same_two_values_have_names(self):
        """Keeping the result whole and reading its parts works as well."""
        model = StubEnergy(FIRST, SECOND)
        result = model.minimize(pair(5.0))
        assert isinstance(result, Relaxed)
        assert result.energy == pytest.approx(0.04)
        assert result.pose.n_atoms == 2


class TestCompositeEnergy:
    """Adding several scorers together, each with its own weight."""

    def test_each_part_is_multiplied_by_its_weight_and_the_results_added(self):
        """One part at weight 2 and the same part at weight 3 come to 5 times
        the part on its own."""
        part = StubEnergy(FIRST, SECOND)
        combined = CompositeEnergy((part, 2.0), (part, 3.0))
        pose = pair(5.0)
        assert combined.evaluate(pose) == pytest.approx(5 * part.evaluate(pose))

    def test_a_weight_of_zero_removes_a_part_from_the_total(self):
        """The score is then whatever the remaining parts say."""
        part = StubEnergy(FIRST, SECOND)
        combined = CompositeEnergy((part, 1.0), (part, 0.0))
        pose = pair(5.0)
        assert combined.evaluate(pose) == pytest.approx(part.evaluate(pose))

    def test_settling_reports_the_weighted_total_of_every_part(self):
        """The energy handed back can be compared with a plain score."""
        part = StubEnergy(FIRST, SECOND)
        combined = CompositeEnergy((part, 2.0))
        pose = pair(5.0)
        assert combined.minimize(pose).energy == pytest.approx(combined.evaluate(pose))

    def test_a_composite_with_no_parts_is_rejected(self):
        """It would score every arrangement zero, so nothing could be chosen."""
        with pytest.raises(ConfigurationError, match="at least one part"):
            CompositeEnergy()


class TestEnergyModelProtocol:
    """What a class has to do to be usable as a scorer."""

    def test_the_stand_in_scorer_counts_as_one(self):
        """It has both methods, and there is nothing else to register."""
        assert isinstance(StubEnergy(FIRST, SECOND), EnergyModel)

    def test_a_composite_counts_as_one_too(self):
        """A composite can therefore be a part of another composite."""
        assert isinstance(
            CompositeEnergy((StubEnergy(FIRST, SECOND), 1.0)), EnergyModel
        )

    def test_something_with_neither_method_does_not(self):
        """A plain object cannot be handed to a search as its scorer."""
        assert not isinstance(object(), EnergyModel)

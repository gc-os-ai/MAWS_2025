"""
Tests for :mod:`maws.scoring`.

The score reduces the energies of every shape tried for one candidate strand to
a single number, lower being more promising. What the tests here pin down is
where its zero point sits, which direction it moves in, and which changes to
the energies it deliberately ignores — because those are the properties that
make two candidates, scored in two different runs, comparable at all.

Energies are in kJ/mol and ``beta`` is in mol/kJ throughout.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from maws.errors import ConfigurationError
from maws.scoring import boltzmann_weights, entropy_score


class TestTheZeroPoint:
    """Where the score sits when no shape is preferred to any other."""

    def test_equal_energies_score_exactly_zero(self):
        """Zero is defined as "every shape equally likely", so it is exact.

        The comparison is against a tolerance rather than ``== 0`` only because
        the sums are floating point at the end; there is no rounding in the
        answer itself, which is why the tolerance can be this tight.
        """
        assert abs(entropy_score([100.0, 100.0, 100.0])) < 1e-12

    @pytest.mark.parametrize("count", [1, 3, 300])
    def test_the_zero_point_does_not_move_with_the_number_of_shapes(self, count):
        """Runs that tried different numbers of shapes stay comparable."""
        assert abs(entropy_score([100.0] * count)) < 1e-12

    def test_a_single_shape_scores_zero(self):
        """One shape is trivially the only shape, so nothing is preferred."""
        assert abs(entropy_score([-4321.0])) < 1e-12

    def test_no_preference_at_all_scores_zero_whatever_the_energies(self):
        """At ``beta`` of zero every shape weighs the same, however bad it is."""
        assert entropy_score([0.0, 500.0, -900.0, 12.5], beta=0.0) == 0.0


class TestWhichWayTheScoreMoves:
    """Lower means the weight has piled onto fewer shapes."""

    def test_one_dominant_shape_scores_well_below_zero(self):
        """A strand that settles into one shape is what a binder looks like."""
        assert entropy_score([0.0] + [1000.0] * 9) < -2.0

    def test_shapes_that_are_much_of_a_muchness_score_near_zero(self):
        """Energies within a thousandth of each other give almost no preference."""
        assert abs(entropy_score([0.0, 0.001, 0.002])) < 1e-9

    def test_concentrating_the_weight_scores_lower_than_spreading_it(self):
        """The whole point of the score: concentration beats a flat spread."""
        one_clear_winner = entropy_score([0.0] + [1000.0] * 9)
        nothing_to_choose = entropy_score([0.0] + [10.0] * 9)
        assert one_clear_winner < nothing_to_choose

    @pytest.mark.parametrize(
        "energies",
        [
            [0.0, 1.0],
            [0.0, 1000.0],
            [5.0, 5.0, 5.0],
            [-200.0, -199.0, 340.0, 12.0],
            [0.0] * 50 + [-3000.0],
        ],
    )
    def test_the_score_never_rises_above_zero(self, energies):
        """Nothing is worse than no preference, so zero is the ceiling."""
        assert entropy_score(energies) <= 0.0

    def test_raising_beta_makes_the_same_energies_score_lower(self):
        """A larger ``beta`` lets the best shapes count for more."""
        energies = [0.0, 50.0, 100.0]
        assert entropy_score(energies, beta=0.1) < entropy_score(energies, beta=0.01)

    @pytest.mark.parametrize("beta", [0.001, 0.01, 0.1, 1.0])
    def test_the_score_stays_at_or_below_zero_for_any_beta(self, beta):
        """Sharpening the preference cannot turn the score positive."""
        assert entropy_score([0.0, 50.0, 100.0], beta=beta) <= 0.0


class TestWhatTheScoreIgnores:
    """Changes to the energies that must leave the answer alone."""

    @pytest.mark.parametrize("shift", [1000.0, -1000.0, 1e5])
    def test_adding_the_same_amount_to_every_energy_changes_nothing(self, shift):
        """Only the gaps between energies matter, never where they sit.

        Two runs can report energies offset from each other — a different
        reference point for the target on its own, say — and still be compared,
        because the shift cancels when the weights are divided by their total.
        """
        energies = [0.0, 50.0, 120.0, 7.5]
        shifted = [energy + shift for energy in energies]
        assert entropy_score(shifted) == pytest.approx(
            entropy_score(energies), abs=1e-12
        )

    def test_reordering_the_energies_changes_nothing(self):
        """The shapes are a set; the order they were tried in carries nothing."""
        assert entropy_score([0.0, 50.0, 120.0]) == pytest.approx(
            entropy_score([120.0, 0.0, 50.0]), abs=1e-12
        )


class TestVeryLargeEnergies:
    """Whole-molecule energies run to thousands of kJ/mol, and must still work."""

    def test_a_spread_of_a_hundred_thousand_still_gives_a_real_number(self):
        """Ordinary floating point loses these weights; this arithmetic does not."""
        assert math.isfinite(entropy_score([0.0, 1e5]))

    def test_a_hopeless_shape_takes_the_score_to_its_floor(self):
        """One shape holding all the weight out of two scores minus ln 2."""
        assert entropy_score([0.0, 1e5]) == pytest.approx(-math.log(2.0), abs=1e-9)

    @pytest.mark.parametrize(
        "energies",
        [[0.0, 1e5], [-1e5, 1e5], [0.0, 1e5, 2e5, 5e4], [1e5] * 4],
    )
    def test_no_spread_of_energies_produces_a_broken_number(self, energies):
        """Never a nan and never an infinity, whatever the numbers were."""
        assert math.isfinite(entropy_score(energies))


class TestScoreArguments:
    """What the score refuses to be given."""

    def test_no_energies_at_all_is_rejected(self):
        """An empty list usually means the sampling loop never ran."""
        with pytest.raises(ConfigurationError, match="at least one energy"):
            entropy_score([])

    def test_a_negative_beta_is_rejected(self):
        """A negative ``beta`` would prefer the worst shapes over the best."""
        with pytest.raises(ConfigurationError, match="must not be negative"):
            entropy_score([0.0, 1.0], beta=-0.01)


class TestBoltzmannWeights:
    """The per-shape probabilities the score is worked out from."""

    def test_the_weights_add_up_to_one(self):
        """They are probabilities, so they account for the whole of the weight."""
        weights = boltzmann_weights([0.0, 50.0, 120.0, -30.0])
        assert float(weights.sum()) == pytest.approx(1.0, abs=1e-12)

    def test_there_is_one_weight_per_shape(self):
        """Each energy given gets its own share of the total back."""
        assert boltzmann_weights([0.0, 50.0, 120.0]).shape == (3,)

    def test_equal_energies_share_the_weight_evenly(self):
        """Four shapes nobody prefers take a quarter each."""
        assert boltzmann_weights([7.0] * 4) == pytest.approx(np.full(4, 0.25))

    def test_the_lower_energy_of_a_pair_takes_the_larger_share(self):
        """Lower energy means a more settled shape, and so more weight."""
        weights = boltzmann_weights([0.0, 100.0])
        assert weights[0] > weights[1]

    def test_a_pair_a_hundred_apart_splits_roughly_three_to_one(self):
        """At ``beta`` of 0.01 the worse shape weighs ``exp(-1)`` as much."""
        weights = boltzmann_weights([0.0, 100.0])
        assert weights[1] / weights[0] == pytest.approx(math.exp(-1.0))

    def test_no_energies_at_all_is_rejected(self):
        """Nothing to weigh means the sampling loop produced nothing."""
        with pytest.raises(ConfigurationError, match="at least one energy"):
            boltzmann_weights([])

    def test_a_negative_beta_is_rejected(self):
        """The same rule as for the score, for the same reason."""
        with pytest.raises(ConfigurationError, match="must not be negative"):
            boltzmann_weights([0.0, 1.0], beta=-1.0)

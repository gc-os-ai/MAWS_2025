"""
Tests for :mod:`maws.scoring`.

A score reduces the energies of every shape tried for one candidate strand to a
single number, lower being more promising. Two of them live in that module and
they do not agree, so the tests come in three parts.

:func:`~maws.scoring.free_energy_score` is what the search uses. What is pinned
down for it is that it tracks the energies themselves, that it rewards a
candidate for reaching a low energy often rather than once, and that it stays
between the best energy and the average.

:func:`~maws.scoring.entropy_score` is the original MAWS score, kept so that
runs can be compared against published ones. What is pinned down for it is
where its zero point sits, which direction it moves in, and what it ignores.

The last part puts the two side by side on the case that separates them.

Energies are in kJ/mol and ``beta`` is in mol/kJ throughout.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from maws.errors import ConfigurationError
from maws.scoring import boltzmann_weights, entropy_score, free_energy_score


class TestTheFreeEnergyTracksTheEnergies:
    """The score the search uses moves with the energies it is given."""

    def test_shapes_that_are_all_the_same_score_exactly_that_energy(self):
        """With nothing to choose between the shapes, the answer is their energy."""
        assert free_energy_score([-250.0] * 8) == pytest.approx(-250.0)

    @pytest.mark.parametrize("beta", [0.001, 0.01, 1.0])
    def test_that_holds_whatever_beta_is(self, beta):
        """``beta`` decides how the shapes are weighed against each other.

        When they are all the same there is nothing to weigh, so it drops out.
        """
        assert free_energy_score([-250.0] * 8, beta=beta) == pytest.approx(-250.0)

    @pytest.mark.parametrize("shift", [1000.0, -1000.0, 1e5])
    def test_moving_every_energy_moves_the_score_with_them(self, shift):
        """A strand that binds 1000 kJ/mol more strongly must score 1000 lower.

        This is the property the concentration score does not have, and it is
        why the search does not use that one: without it, a candidate that
        never gets below zero cannot be told apart from one that reaches
        -1000 kJ/mol.
        """
        energies = [-900.0, -100.0, -100.0, -100.0]
        shifted = [energy + shift for energy in energies]
        assert free_energy_score(shifted) == pytest.approx(
            free_energy_score(energies) + shift
        )

    def test_lower_energies_score_lower(self):
        """The whole point: a strand that binds harder is more promising."""
        assert free_energy_score([-1000.0] * 4) < free_energy_score([-100.0] * 4)

    def test_reaching_a_good_energy_more_often_scores_lower(self):
        """Two candidates with the same best shape are separated by the rest.

        Binding means settling into a family of good shapes, not stumbling on
        one. The candidate that reaches -500 kJ/mol in every shape it tried is
        the one that has found a fit.
        """
        every_time = free_energy_score([-500.0] * 4)
        just_once = free_energy_score([-500.0, 0.0, 0.0, 0.0])
        assert every_time < just_once


class TestWhereTheFreeEnergyCanLie:
    """The answer is bounded by the numbers it was given, at either end."""

    ENERGIES = [-900.0, -100.0, -100.0, -100.0]
    """One good shape among three ordinary ones, so the bounds are far apart."""

    @pytest.mark.parametrize("beta", [0.001, 0.01, 0.1, 1.0])
    def test_it_never_beats_the_best_shape(self, beta):
        """No weighting of the shapes can be better than the best of them."""
        assert free_energy_score(self.ENERGIES, beta=beta) >= min(self.ENERGIES)

    @pytest.mark.parametrize("beta", [0.001, 0.01, 0.1, 1.0])
    def test_it_is_never_worse_than_the_average(self, beta):
        """Weighting by energy can only help, since low energies weigh more."""
        average = sum(self.ENERGIES) / len(self.ENERGIES)
        assert free_energy_score(self.ENERGIES, beta=beta) <= average + 1e-9

    def test_a_large_beta_pulls_it_towards_the_best_shape(self):
        """At a high enough ``beta`` only the best shape is counted at all."""
        assert free_energy_score(self.ENERGIES, beta=10.0) == pytest.approx(
            min(self.ENERGIES), abs=0.5
        )

    def test_a_small_beta_pulls_it_towards_the_average(self):
        """At a low enough ``beta`` every shape counts the same."""
        average = sum(self.ENERGIES) / len(self.ENERGIES)
        assert free_energy_score(self.ENERGIES, beta=1e-6) == pytest.approx(
            average, abs=0.5
        )


class TestFreeEnergyArguments:
    """What the free energy refuses to be given."""

    def test_no_energies_at_all_is_rejected(self):
        """An empty list usually means the sampling loop never ran."""
        with pytest.raises(ConfigurationError, match="at least one energy"):
            free_energy_score([])

    @pytest.mark.parametrize("beta", [0.0, -0.01])
    def test_a_beta_of_zero_or_less_is_rejected(self, beta):
        """The answer is divided by ``beta``, and a negative one inverts it."""
        with pytest.raises(ConfigurationError, match="greater than zero"):
            free_energy_score([0.0, 1.0], beta=beta)

    @pytest.mark.parametrize(
        "energies",
        [[0.0, 1e5], [-1e5, 1e5], [0.0, 1e5, 2e5, 5e4], [1e5] * 4],
    )
    def test_no_spread_of_energies_produces_a_broken_number(self, energies):
        """Never a nan and never an infinity, whatever the numbers were.

        A shape where two atoms end up on top of each other produces an energy
        in the millions, and the weights are exponentials of these numbers.
        """
        assert math.isfinite(free_energy_score(energies))


class TestTheTwoScoresDisagree:
    """The case that decides which of the two the search should use.

    One candidate crashes atoms into each other in almost every shape tried
    and is ordinary in the rest. The other is ordinary in all of them, and its
    best shape is no worse. The first is a strand that has found nothing; the
    second is a strand that fits.
    """

    MOSTLY_IMPOSSIBLE = [-100.0, -100.0] + [1e6] * 8
    """Eight shapes with atoms on top of each other, and two ordinary ones."""

    ALL_REASONABLE = [-100.0] * 10
    """Ten shapes, none of them impossible, reaching the same best energy."""

    def test_the_concentration_score_prefers_the_one_that_mostly_clashes(self):
        """Two survivors among eight disasters is a very concentrated spread.

        That is what the concentration score measures, so it scores it well.
        This is not a bug in the arithmetic; it is what the question "how
        concentrated" answers when asked about a candidate like this.
        """
        assert entropy_score(self.MOSTLY_IMPOSSIBLE) < entropy_score(
            self.ALL_REASONABLE
        )

    def test_the_free_energy_prefers_the_one_that_never_clashes(self):
        """The impossible shapes contribute nothing, and their count counts.

        Averaging over ten shapes when only two of them are reachable is what
        pushes the mostly-impossible candidate up, so the search picks the
        other one.
        """
        assert free_energy_score(self.ALL_REASONABLE) < free_energy_score(
            self.MOSTLY_IMPOSSIBLE
        )


class TestTheZeroPoint:
    """Where the concentration score sits when no shape is preferred."""

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


class TestWhichWayTheConcentrationScoreMoves:
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


class TestWhatTheConcentrationScoreIgnores:
    """Changes to the energies that leave the concentration score alone.

    The first of these is the reason it is not the default. Ignoring where the
    energies sit is what lets it rate a candidate that never binds as highly
    as one that does.
    """

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


class TestTheConcentrationScoreOnVeryLargeEnergies:
    """Whole-molecule energies run to thousands of kJ/mol, and must still work."""

    def test_a_spread_of_a_hundred_thousand_still_gives_a_real_number(self):
        """Ordinary floating point loses these weights; this arithmetic does not."""
        assert math.isfinite(entropy_score([0.0, 1e5]))

    def test_a_hopeless_shape_takes_the_score_to_its_floor(self):
        """One shape holding all the weight out of two scores minus ln 2.

        Written down because it is the arithmetic, not because it is wanted:
        the second energy here is an impossible shape, and adding it to the
        list *improves* the concentration score from zero to -0.69. See
        ``TestTheTwoScoresDisagree`` for what that costs.
        """
        assert entropy_score([0.0, 1e5]) == pytest.approx(-math.log(2.0), abs=1e-9)

    @pytest.mark.parametrize(
        "energies",
        [[0.0, 1e5], [-1e5, 1e5], [0.0, 1e5, 2e5, 5e4], [1e5] * 4],
    )
    def test_no_spread_of_energies_produces_a_broken_number(self, energies):
        """Never a nan and never an infinity, whatever the numbers were."""
        assert math.isfinite(entropy_score(energies))


class TestConcentrationScoreArguments:
    """What the concentration score refuses to be given."""

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

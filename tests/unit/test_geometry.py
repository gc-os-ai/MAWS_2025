"""
Tests for :mod:`maws.geometry`.

Three small vector calculations with no state and no chemistry, so every
answer here is worked out by hand: right-angled triangles, quarter turns, and
a two-atom balance point.
"""

from __future__ import annotations

import numpy as np
import pytest

from maws.errors import ConfigurationError
from maws.geometry import angle_between, centre_of_mass, unit_vector


class TestUnitVector:
    """Shrinking a vector to length one without changing where it points."""

    def test_the_result_has_length_one(self):
        """Whatever goes in, what comes out measures one."""
        assert np.linalg.norm(unit_vector(np.array([1.0, 2.0, 3.0]))) == pytest.approx(
            1.0
        )

    def test_a_three_four_five_triangle_divides_by_five(self):
        """A vector of length 5 is divided by 5, component by component."""
        np.testing.assert_allclose(
            unit_vector(np.array([3.0, 4.0, 0.0])), [0.6, 0.8, 0.0]
        )

    def test_lengthening_a_vector_does_not_change_its_direction(self):
        """Only direction survives, so a scaled copy gives the same answer."""
        vector = np.array([1.0, 2.0, 3.0])
        np.testing.assert_allclose(unit_vector(vector), unit_vector(100.0 * vector))

    def test_a_vector_of_integers_comes_back_as_fractions(self):
        """Whole numbers in do not truncate the division to whole numbers."""
        np.testing.assert_allclose(unit_vector(np.array([0, 3, 0])), [0.0, 1.0, 0.0])

    def test_a_vector_of_zero_length_points_nowhere(self):
        """There is no direction to return, so this is refused rather than nan."""
        with pytest.raises(ConfigurationError, match="no direction"):
            unit_vector(np.zeros(3))


class TestAngleBetween:
    """The angle between two directions, in radians."""

    def test_two_axes_are_a_quarter_turn_apart(self):
        """The x and y directions meet at a right angle."""
        angle = angle_between(np.array([1.0, 0, 0]), np.array([0.0, 1, 0]))
        assert angle == pytest.approx(np.pi / 2)

    def test_opposite_directions_are_half_a_turn_apart(self):
        """The largest angle two directions can make is π."""
        angle = angle_between(np.array([1.0, 0, 0]), np.array([-1.0, 0, 0]))
        assert angle == pytest.approx(np.pi)

    def test_length_is_ignored(self):
        """A long vector and a short one meet at the same angle as two short ones."""
        short = angle_between(np.array([1.0, 0, 0]), np.array([0.0, 1, 0]))
        mixed = angle_between(np.array([1000.0, 0, 0]), np.array([0.0, 0.001, 0]))
        assert mixed == pytest.approx(short)

    def test_a_vector_makes_no_angle_with_a_scaled_copy_of_itself(self):
        """Two vectors pointing the same way are zero radians apart.

        ``[1, 1, 1]`` and ``[0.7, 0.7, 0.7]`` are the interesting case: their
        cosine works out to 1.0000000000000002 in floating-point arithmetic,
        and an arc-cosine of that is not-a-number. Clipping the cosine to
        ``[-1, 1]`` first is what keeps this a number.
        """
        angle = angle_between(np.array([1.0, 1.0, 1.0]), np.array([0.7, 0.7, 0.7]))
        assert not np.isnan(angle)
        assert angle == 0.0

    @pytest.mark.parametrize("position", [0, 1])
    def test_a_vector_of_zero_length_is_rejected(self, position):
        """Either argument being zero-length leaves no angle to measure."""
        arguments = [np.array([1.0, 0, 0]), np.array([1.0, 0, 0])]
        arguments[position] = np.zeros(3)
        with pytest.raises(ConfigurationError, match="no direction"):
            angle_between(*arguments)


class TestCentreOfMass:
    """The balance point of a group of atoms."""

    def test_equal_masses_give_the_plain_average_position(self):
        """When every atom weighs the same, this is the midpoint."""
        centre = centre_of_mass(
            np.array([[0.0, 0, 0], [2.0, 0, 0]]), np.array([5.0, 5.0])
        )
        np.testing.assert_allclose(centre, [1.0, 0.0, 0.0])

    def test_the_balance_point_sits_nearer_the_heavier_atom(self):
        """A nine-to-one split lands nine tenths of the way along."""
        centre = centre_of_mass(
            np.array([[0.0, 0, 0], [10.0, 0, 0]]), np.array([1.0, 9.0])
        )
        np.testing.assert_allclose(centre, [9.0, 0.0, 0.0])

    def test_every_axis_is_averaged_separately(self):
        """The three coordinates do not interact."""
        centre = centre_of_mass(
            np.array([[0.0, 4.0, 8.0], [2.0, 0.0, 0.0]]), np.array([1.0, 1.0])
        )
        np.testing.assert_allclose(centre, [1.0, 2.0, 4.0])

    def test_one_mass_per_atom_is_required(self):
        """A missing mass is a mistake worth naming, not one to average over."""
        with pytest.raises(ConfigurationError, match="positions but"):
            centre_of_mass(np.zeros((3, 3)), np.array([1.0, 1.0]))

    def test_massless_atoms_have_no_balance_point(self):
        """Dividing by a total mass of zero is refused rather than returning inf."""
        with pytest.raises(ConfigurationError, match="total mass is zero"):
            centre_of_mass(np.zeros((2, 3)), np.zeros(2))

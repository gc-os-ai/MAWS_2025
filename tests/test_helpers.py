"""
Tests for maws.helpers module.

This module tests the unit conversion and vector utility functions:
- angstrom/nostrom: Attach/strip Ångström units
- kJ/noJ: Attach/strip energy units
- angle: Compute angle between vectors
- directed_angle: Signed angle about an axis
- center_of_mass: Compute centroid of coordinates
"""

import math

import numpy as np
import pytest
from openmm import unit

from maws.helpers import (
    angle,
    angstrom,
    center_of_mass,
    directed_angle,
    kJ,
    noJ,
    nostrom,
)


class TestUnitConversions:
    """Tests for unit attachment/stripping functions."""

    def test_angstrom_attaches_units(self):
        """angstrom() attaches Å units to numeric array."""
        vec = angstrom([1.0, 2.0, 3.0])
        assert vec.unit == unit.angstrom
        # Check values preserved
        values = vec.value_in_unit(unit.angstrom)
        np.testing.assert_array_equal(values, [1.0, 2.0, 3.0])

    def test_nostrom_strips_units(self):
        """nostrom() strips Å units and returns raw numbers."""
        vec = angstrom([5.0, 10.0, 15.0])
        raw = nostrom(vec)
        np.testing.assert_array_equal(raw, [5.0, 10.0, 15.0])

    def test_nostrom_converts_nanometers(self):
        """nostrom() converts other length units to Ångströms."""
        # 1 nm = 10 Å
        vec = 1.0 * unit.nanometers
        raw = nostrom(vec)
        assert raw == pytest.approx(10.0)

    def test_kJ_attaches_energy_units(self):  # noqa: N802
        """kJ() attaches kJ/mol units to numeric values."""
        energy = kJ(100.5)
        assert energy.unit == unit.kilojoules_per_mole
        assert energy.value_in_unit(unit.kilojoules_per_mole) == pytest.approx(100.5)

    def test_noJ_strips_energy_units(self):  # noqa: N802
        """noJ() strips kJ/mol units and returns raw numbers."""
        energy = kJ(250.0)
        raw = noJ(energy)
        assert float(raw) == pytest.approx(250.0)


class TestAngleCalculations:
    """Tests for angle computation functions."""

    def test_angle_perpendicular_vectors(self):
        """angle() returns π/2 for perpendicular vectors."""
        result = angle([1, 0, 0], [0, 1, 0])
        assert result == pytest.approx(math.pi / 2)

    def test_angle_parallel_vectors(self):
        """angle() returns 0 for parallel vectors."""
        result = angle([1, 0, 0], [1, 0, 0])
        assert result == pytest.approx(0.0)

    def test_angle_antiparallel_vectors(self):
        """angle() returns π for antiparallel vectors."""
        result = angle([1, 0, 0], [-1, 0, 0])
        assert result == pytest.approx(math.pi)

    def test_directed_angle_positive(self):
        """directed_angle() returns positive angle for CCW rotation."""
        # Rotate [1,0,0] to [0,1,0] about z-axis = +90 degrees
        result = directed_angle([1, 0, 0], [0, 1, 0], [0, 0, 1])
        assert result == pytest.approx(math.pi / 2)

    def test_directed_angle_negative(self):
        """directed_angle() returns negative angle for CW rotation."""
        # Rotate [0,1,0] to [1,0,0] about z-axis = -90 degrees
        result = directed_angle([0, 1, 0], [1, 0, 0], [0, 0, 1])
        assert result == pytest.approx(-math.pi / 2)


class TestCenterOfMass:
    """Tests for center_of_mass function."""

    def test_center_of_mass_simple(self):
        """center_of_mass() returns correct centroid."""
        coords = np.array([[0, 0, 0], [2, 0, 0], [1, 3, 0]])
        result = center_of_mass(coords)
        expected = np.array([1.0, 1.0, 0.0])
        np.testing.assert_array_almost_equal(result, expected)

    def test_center_of_mass_single_point(self):
        """center_of_mass() of single point is that point."""
        coords = np.array([[5.0, 3.0, 1.0]])
        result = center_of_mass(coords)
        np.testing.assert_array_equal(result, [5.0, 3.0, 1.0])

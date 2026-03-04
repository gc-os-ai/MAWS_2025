import numpy as np
import pytest

from maws import helpers
from openmm import unit


def test_angstrom_and_nostrom():
    arr = [1, 2, 3]
    q = helpers.angstrom(arr)
    assert q.unit == unit.angstrom
    np.testing.assert_allclose(q.value_in_unit(unit.angstrom), np.array([1.0, 2.0, 3.0]))
    # round-trip
    assert np.allclose(helpers.nostrom(q), np.array(arr))


def test_kJ_and_noJ():
    val = 5.0
    q = helpers.kJ(val)
    assert q.unit == unit.kilojoules_per_mole
    assert helpers.noJ(q) == pytest.approx(val)


def test_angle_directed_angle():
    a = [1, 0, 0]
    b = [0, 1, 0]
    assert helpers.angle(a, b) == pytest.approx(np.pi / 2)
    # same direction should be zero
    assert helpers.angle(a, a) == pytest.approx(0.0)
    # directed angle around z axis
    axis = [0, 0, 1]
    assert helpers.directed_angle(a, b, axis) == pytest.approx(np.pi / 2)
    assert helpers.directed_angle(b, a, axis) == pytest.approx(-np.pi / 2)

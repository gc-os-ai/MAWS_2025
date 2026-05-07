"""
Tests for maws.space module.

Surface-aware sampling primitives:
- Cube / Sphere / Shell : envelope dataclasses
- NAngles : torsion angle sampling (U(1)^N)
- Excluder : SAS-style rejection
- SurfaceSampler : envelope + Excluder composer
- make_sampler : factory
"""

import math

from maws.space import NAngles


class TestNAngles:
    """Tests for NAngles (torsion angle) sampling space."""

    def test_nangles_generator_returns_n_elements(self):
        """NAngles.generator() returns array of N angles."""
        angles = NAngles(number=4)
        sample = angles.generator()
        assert len(sample) == 4

    def test_nangles_in_range(self):
        """NAngles samples are in [0, 2π)."""
        angles = NAngles(number=5)
        for _ in range(10):
            sample = angles.generator()
            for angle in sample:
                assert 0 <= angle < 2 * math.pi, f"Angle {angle} out of range"

    def test_nangles_different_counts(self):
        """NAngles works for various N values."""
        for n in [1, 3, 10]:
            angles = NAngles(number=n)
            sample = angles.generator()
            assert len(sample) == n


# ---------- New surface-aware sampler tests ----------

import numpy as np  # noqa: E402

from maws.space import _BONDI_VDW_RADII, _DEFAULT_VDW, Cube  # noqa: E402


class TestBondiTable:
    def test_has_common_biological_elements(self):
        for sym in ("H", "C", "N", "O", "P", "S"):
            assert sym in _BONDI_VDW_RADII

    def test_carbon_value(self):
        # Bondi 1964: C = 1.70 Å
        assert _BONDI_VDW_RADII["C"] == 1.70

    def test_default_vdw_is_carbon_like(self):
        assert _DEFAULT_VDW == 1.70


class TestCube:
    def test_generator_returns_7_elements(self):
        c = Cube(width=10.0, centre=np.array([0.0, 0.0, 0.0]))
        assert len(c.generator()) == 7

    def test_generator_position_within_bounds(self):
        c = Cube(width=10.0, centre=np.array([5.0, 5.0, 5.0]))
        for _ in range(50):
            s = c.generator()
            assert 0.0 <= s[0] <= 10.0
            assert 0.0 <= s[1] <= 10.0
            assert 0.0 <= s[2] <= 10.0

    def test_generator_rotation_axis_unit_length(self):
        c = Cube(width=10.0, centre=np.array([0.0, 0.0, 0.0]))
        for _ in range(20):
            s = c.generator()
            axis_norm = math.sqrt(s[3] ** 2 + s[4] ** 2 + s[5] ** 2)
            assert abs(axis_norm - 1.0) < 1e-9

    def test_generator_rotation_angle_in_range(self):
        c = Cube(width=10.0, centre=np.array([0.0, 0.0, 0.0]))
        for _ in range(20):
            s = c.generator()
            assert 0.0 <= s[6] <= math.pi

    def test_is_frozen_dataclass(self):
        c = Cube(width=10.0, centre=np.array([0.0, 0.0, 0.0]))
        try:
            c.width = 20.0  # should raise
        except Exception as e:
            assert "frozen" in str(e).lower() or "FrozenInstance" in type(e).__name__
            return
        raise AssertionError("expected FrozenInstanceError on attribute assignment")

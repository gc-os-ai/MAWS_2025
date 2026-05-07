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

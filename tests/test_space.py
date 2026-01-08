"""
Tests for maws.space module.

This module tests the sampling space classes used for molecular conformational
exploration in MAWS:
- Box/Cube: Rectangular sampling regions with random rotation
- Sphere/SphericalShell: Spherical sampling regions
- NAngles: Torsion angle sampling (U(1)^N)
"""

import math

import numpy as np

from maws.space import Box, Cube, NAngles, Sphere, SphericalShell


class TestBox:
    """Tests for Box sampling space."""

    def test_box_generator_returns_7_elements(self):
        """Box.generator() returns 7-element array [x,y,z,ax,ay,az,angle]."""
        box = Box(
            x_width=10.0,
            y_width=10.0,
            z_width=10.0,
            centre=[0, 0, 0],
        )
        sample = box.generator()
        assert len(sample) == 7

    def test_box_sample_within_bounds(self):
        """Box samples are within ±width/2 of center."""
        box = Box(
            x_width=10.0,
            y_width=20.0,
            z_width=30.0,
            centre=[0, 0, 0],
        )
        # Generate many samples and check bounds
        for _ in range(10):
            sample = box.generator()
            # Position is first 3 elements
            assert -5.0 <= sample[0] <= 5.0, "x out of bounds"
            assert -10.0 <= sample[1] <= 10.0, "y out of bounds"
            assert -15.0 <= sample[2] <= 15.0, "z out of bounds"

    def test_box_rotation_angle_in_range(self):
        """Box rotation angle is in [0, π]."""
        box = Box(x_width=10.0, y_width=10.0, z_width=10.0, centre=[0, 0, 0])
        for _ in range(10):
            sample = box.generator()
            # Rotation angle is last element
            assert 0 <= sample[6] <= math.pi


class TestCube:
    """Tests for Cube sampling space (equal width Box)."""

    def test_cube_generator_returns_7_elements(self):
        """Cube.generator() returns 7-element array."""
        cube = Cube(width=20.0, centre=[0, 0, 0])
        sample = cube.generator()
        assert len(sample) == 7

    def test_cube_sample_within_bounds(self):
        """Cube samples are within ±width/2 of center."""
        cube = Cube(width=20.0, centre=[5, 5, 5])
        for _ in range(10):
            sample = cube.generator()
            # Center is at [5,5,5], width=20 → range [-5, 15]
            assert -5.0 <= sample[0] <= 15.0
            assert -5.0 <= sample[1] <= 15.0
            assert -5.0 <= sample[2] <= 15.0


class TestSphere:
    """Tests for Sphere sampling space."""

    def test_sphere_generator_returns_7_elements(self):
        """Sphere.generator() returns 7-element array."""
        sphere = Sphere(radius=10.0, centre=[0, 0, 0])
        sample = sphere.generator()
        assert len(sample) == 7

    def test_sphere_sample_within_radius(self):
        """Sphere samples are within the specified radius."""
        sphere = Sphere(radius=10.0, centre=[0, 0, 0])
        for _ in range(20):
            sample = sphere.generator()
            # Distance from center
            pos = np.array(sample[:3])
            distance = np.linalg.norm(pos)
            assert distance <= 10.0, f"Sample {pos} outside sphere radius"


class TestSphericalShell:
    """Tests for SphericalShell sampling space."""

    def test_shell_generator_returns_7_elements(self):
        """SphericalShell.generator() returns 7-element array."""
        shell = SphericalShell(innerRadius=5.0, outerRadius=10.0, centre=[0, 0, 0])
        sample = shell.generator()
        assert len(sample) == 7

    def test_shell_sample_within_bounds(self):
        """SphericalShell samples are between inner and outer radii."""
        shell = SphericalShell(innerRadius=5.0, outerRadius=10.0, centre=[0, 0, 0])
        for _ in range(20):
            sample = shell.generator()
            pos = np.array(sample[:3])
            distance = np.linalg.norm(pos)
            assert 5.0 <= distance <= 10.0, f"Sample {pos} outside shell bounds"


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

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

    def test_sphere_sample_offset_by_centre(self):
        """Sphere samples are within radius of the (non-zero) centre, not the origin."""
        centre = np.array([50.0, -30.0, 12.0])
        sphere = Sphere(radius=10.0, centre=centre)
        for _ in range(50):
            sample = sphere.generator()
            pos = np.array(sample[:3])
            distance = np.linalg.norm(pos - centre)
            assert distance <= 10.0 + 1e-9, (
                f"Sample {pos} not within radius of centre {centre}"
            )

    def test_sphere_is_in_membership(self):
        """Sphere.is_in() correctly accepts inside points and rejects outside ones."""
        centre = np.array([10.0, 10.0, 10.0])
        sphere = Sphere(radius=5.0, centre=centre)
        assert bool(sphere.is_in(np.array([10.0, 10.0, 10.0])))
        assert bool(sphere.is_in(np.array([13.0, 10.0, 10.0])))
        assert not bool(sphere.is_in(np.array([100.0, 0.0, 0.0])))


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

    def test_shell_sample_offset_by_centre(self):
        """SphericalShell samples sit in the shell measured from centre, not origin."""
        centre = np.array([50.0, -30.0, 12.0])
        shell = SphericalShell(innerRadius=5.0, outerRadius=10.0, centre=centre)
        for _ in range(50):
            sample = shell.generator()
            pos = np.array(sample[:3])
            distance = np.linalg.norm(pos - centre)
            assert 5.0 - 1e-9 <= distance <= 10.0 + 1e-9, (
                f"Sample {pos} outside shell of centre {centre}"
            )

    def test_shell_radial_distribution_is_uniform_in_volume(self):
        """Shell sampling must be uniform per unit volume (∝ r²), not uniform in r.

        For R_in=5, R_out=10 with the volume-correct CDF
            r = (u·(R_out³ − R_in³) + R_in³)^(1/3),
        E[r] = (3/4)·(R_out⁴ − R_in⁴)/(R_out³ − R_in³) ≈ 8.036.
        A naive `r ~ Uniform(R_in, R_out)` gives E[r] = 7.5 — well outside
        the statistical band of the correct sampler at N=10_000.
        """
        rng = np.random.default_rng(0)
        np.random.seed(0)  # generator() uses np.random
        R_in, R_out = 5.0, 10.0
        shell = SphericalShell(innerRadius=R_in, outerRadius=R_out, centre=[0, 0, 0])
        radii = np.array([np.linalg.norm(shell.generator()[:3]) for _ in range(10_000)])
        expected = (3 / 4) * (R_out**4 - R_in**4) / (R_out**3 - R_in**3)
        # Tight band: correct sampler lands within ~0.05 of 8.036; biased
        # sampler is at 7.5, well outside this band.
        assert abs(radii.mean() - expected) < 0.1, (
            f"E[r] = {radii.mean():.3f}, expected ~{expected:.3f} (uniform-in-r bias?)"
        )
        _ = rng  # silence unused

    def test_shell_direction_is_uniform_on_sphere(self):
        """Shell directions must be uniform on the sphere (E[z²/r²] = 1/3).

        With `psi ~ Uniform(0, π)` the polar angle is biased toward the
        poles, giving E[cos²(psi)] = 1/2 instead of 1/3.
        """
        np.random.seed(1)
        R_in, R_out = 5.0, 10.0
        shell = SphericalShell(innerRadius=R_in, outerRadius=R_out, centre=[0, 0, 0])
        samples = np.array([shell.generator()[:3] for _ in range(10_000)])
        # Direction vectors (unit), then mean of z²
        radii = np.linalg.norm(samples, axis=1, keepdims=True)
        unit_z = (samples / radii)[:, 2]
        mean_z2 = float((unit_z**2).mean())
        # Uniform-on-sphere: 1/3 ≈ 0.333. Biased: 1/2 = 0.500. Tolerance
        # 0.04 cleanly separates the two at N=10_000.
        assert abs(mean_z2 - 1 / 3) < 0.04, (
            f"E[(z/r)²] = {mean_z2:.3f}, expected ~0.333 (polar-angle bias?)"
        )


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

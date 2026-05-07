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
        angles = NAngles(n=4)
        sample = angles.generator()
        assert len(sample) == 4

    def test_nangles_in_range(self):
        """NAngles samples are in [0, 2π)."""
        angles = NAngles(n=5)
        for _ in range(10):
            sample = angles.generator()
            for angle in sample:
                assert 0 <= angle < 2 * math.pi, f"Angle {angle} out of range"

    def test_nangles_different_counts(self):
        """NAngles works for various N values."""
        for n in [1, 3, 10]:
            angles = NAngles(n=n)
            sample = angles.generator()
            assert len(sample) == n

    def test_nangles_is_frozen_dataclass(self):
        a = NAngles(n=3)
        try:
            a.n = 5
        except Exception as e:
            assert "frozen" in str(e).lower() or "FrozenInstance" in type(e).__name__
            return
        raise AssertionError("expected FrozenInstanceError on attribute assignment")


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


class TestSphere:
    def test_generator_returns_7_elements(self):
        from maws.space import Sphere

        s = Sphere(radius=5.0, centre=np.array([0.0, 0.0, 0.0]))
        assert len(s.generator()) == 7

    def test_generator_within_radius_at_origin(self):
        from maws.space import Sphere

        s = Sphere(radius=5.0, centre=np.array([0.0, 0.0, 0.0]))
        for _ in range(50):
            sample = s.generator()
            assert np.linalg.norm(sample[:3]) <= 5.0 + 1e-9

    def test_generator_offset_by_centre(self):
        """Bug-fix from PR #38 carries over: samples must be offset by centre."""
        from maws.space import Sphere

        centre = np.array([50.0, -30.0, 12.0])
        s = Sphere(radius=10.0, centre=centre)
        for _ in range(50):
            sample = s.generator()
            assert np.linalg.norm(sample[:3] - centre) <= 10.0 + 1e-9

    def test_radial_distribution_volume_correct(self):
        """E[r] = 3R/4 for uniform-in-volume sampling (= 3.75 for R=5)."""
        from maws.space import Sphere

        np.random.seed(0)
        s = Sphere(radius=5.0, centre=np.array([0.0, 0.0, 0.0]))
        rs = np.array([np.linalg.norm(s.generator()[:3]) for _ in range(10_000)])
        assert abs(rs.mean() - 3.75) < 0.05

    def test_direction_uniform_on_sphere(self):
        """E[(z/r)^2] = 1/3 for uniform direction (vs 1/2 if biased to poles)."""
        from maws.space import Sphere

        np.random.seed(1)
        s = Sphere(radius=5.0, centre=np.array([0.0, 0.0, 0.0]))
        samples = np.array([s.generator()[:3] for _ in range(10_000)])
        rs = np.linalg.norm(samples, axis=1, keepdims=True)
        unit_z = (samples / rs)[:, 2]
        assert abs(float((unit_z**2).mean()) - 1 / 3) < 0.04


class TestShell:
    def test_generator_returns_7_elements(self):
        from maws.space import Shell

        sh = Shell(inner=5.0, outer=10.0, centre=np.array([0.0, 0.0, 0.0]))
        assert len(sh.generator()) == 7

    def test_generator_within_shell_at_origin(self):
        from maws.space import Shell

        sh = Shell(inner=5.0, outer=10.0, centre=np.array([0.0, 0.0, 0.0]))
        for _ in range(50):
            sample = sh.generator()
            r = np.linalg.norm(sample[:3])
            assert 5.0 - 1e-9 <= r <= 10.0 + 1e-9

    def test_generator_offset_by_centre(self):
        from maws.space import Shell

        centre = np.array([50.0, -30.0, 12.0])
        sh = Shell(inner=5.0, outer=10.0, centre=centre)
        for _ in range(50):
            sample = sh.generator()
            r = np.linalg.norm(sample[:3] - centre)
            assert 5.0 - 1e-9 <= r <= 10.0 + 1e-9

    def test_radial_distribution_uniform_in_volume(self):
        """E[r] = (3/4)·(R_out^4 - R_in^4)/(R_out^3 - R_in^3) ≈ 8.036 for [5, 10]."""
        from maws.space import Shell

        np.random.seed(0)
        sh = Shell(inner=5.0, outer=10.0, centre=np.array([0.0, 0.0, 0.0]))
        rs = np.array([np.linalg.norm(sh.generator()[:3]) for _ in range(10_000)])
        expected = (3 / 4) * (10**4 - 5**4) / (10**3 - 5**3)
        assert abs(rs.mean() - expected) < 0.1

    def test_direction_uniform_on_sphere(self):
        from maws.space import Shell

        np.random.seed(1)
        sh = Shell(inner=5.0, outer=10.0, centre=np.array([0.0, 0.0, 0.0]))
        samples = np.array([sh.generator()[:3] for _ in range(10_000)])
        rs = np.linalg.norm(samples, axis=1, keepdims=True)
        unit_z = (samples / rs)[:, 2]
        assert abs(float((unit_z**2).mean()) - 1 / 3) < 0.04


class TestExcluder:
    def test_clear_far_from_atoms(self, synthetic_two_carbon_complex):
        from maws.space import Excluder

        ex = Excluder(synthetic_two_carbon_complex, probe=1.4)
        assert ex.is_clear(np.array([100.0, 100.0, 100.0]))

    def test_blocked_at_atom_centre(self, synthetic_two_carbon_complex):
        from maws.space import Excluder

        ex = Excluder(synthetic_two_carbon_complex, probe=1.4)
        assert not ex.is_clear(np.array([0.0, 0.0, 0.0]))

    def test_boundary_just_inside_and_just_outside(self, synthetic_two_carbon_complex):
        """C vdW = 1.70, probe = 1.4 → blocked iff dist <= 3.10."""
        from maws.space import Excluder

        ex = Excluder(synthetic_two_carbon_complex, probe=1.4)
        # 3.0 from atom — inside the inflated sphere, blocked
        assert not ex.is_clear(np.array([3.0, 0.0, 0.0]))
        # 3.2 from atom — just outside, clear
        assert ex.is_clear(np.array([3.2, 0.0, 0.0]))

    def test_unknown_element_falls_back_with_warning(self):
        """Unknown element symbol uses _DEFAULT_VDW and emits a warning once."""
        import warnings

        from openmm import unit as ommunit

        from maws.space import _DEFAULT_VDW, Excluder

        positions = np.array([[0.0, 0.0, 0.0]]) * ommunit.angstrom
        atom = type("A", (), {"element": type("E", (), {"symbol": "Xx"})()})()
        topo = type("T", (), {"atoms": lambda self: iter([atom])})()
        cx = type("C", (), {"positions": positions, "topology": topo})()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ex = Excluder(cx, probe=0.0)
            # blocked iff dist < _DEFAULT_VDW (1.70)
            assert not ex.is_clear(np.array([1.0, 0.0, 0.0]))
            assert ex.is_clear(np.array([2.5, 0.0, 0.0]))
        assert any("Xx" in str(w.message) for w in caught), (
            f"expected warning mentioning 'Xx', got {[str(w.message) for w in caught]}"
        )
        # Touch the unused-but-meaningful default to silence linters
        assert _DEFAULT_VDW == 1.70


class TestComputeEnvelopeDims:
    def test_cube_dims_octahedron(self, synthetic_octahedron_complex):
        from maws.space import compute_envelope_dims

        d = compute_envelope_dims(synthetic_octahedron_complex, "cube", reach=10.0)
        # COM = origin; R_max = 5.0; width = 2*(5+10) = 30
        np.testing.assert_allclose(d["centre"], [0.0, 0.0, 0.0])
        assert d["width"] == 30.0

    def test_sphere_dims_octahedron(self, synthetic_octahedron_complex):
        from maws.space import compute_envelope_dims

        d = compute_envelope_dims(synthetic_octahedron_complex, "sphere", reach=10.0)
        # radius = R_max + reach = 15.0
        np.testing.assert_allclose(d["centre"], [0.0, 0.0, 0.0])
        assert d["radius"] == 15.0

    def test_shell_dims_octahedron(self, synthetic_octahedron_complex):
        from maws.space import compute_envelope_dims

        d = compute_envelope_dims(synthetic_octahedron_complex, "shell", reach=10.0)
        # R_min = R_max = 5.0 (octahedron); inner = max(0, 5 - 5) = 0; outer = 15
        np.testing.assert_allclose(d["centre"], [0.0, 0.0, 0.0])
        assert d["inner"] == 0.0
        assert d["outer"] == 15.0

    def test_unknown_shape_raises(self, synthetic_octahedron_complex):
        import pytest

        from maws.space import compute_envelope_dims

        with pytest.raises(ValueError, match="cube|sphere|shell"):
            compute_envelope_dims(synthetic_octahedron_complex, "blob", reach=10.0)


class TestSurfaceSampler:
    def test_all_samples_clear(self, synthetic_octahedron_complex):
        from maws.space import Excluder, Sphere, SurfaceSampler

        envelope = Sphere(radius=15.0, centre=np.array([0.0, 0.0, 0.0]))
        excluder = Excluder(synthetic_octahedron_complex, probe=1.4)
        sampler = SurfaceSampler(envelope=envelope, excluder=excluder)
        for _ in range(50):
            sample = sampler.generator()
            assert excluder.is_clear(sample[:3])

    def test_raises_when_envelope_buried(self, synthetic_two_carbon_complex):
        """A Cube of width 0.1 sitting on an atom is fully blocked."""
        import pytest

        from maws.space import Cube, Excluder, SamplingError, SurfaceSampler

        envelope = Cube(width=0.1, centre=np.array([0.0, 0.0, 0.0]))
        excluder = Excluder(synthetic_two_carbon_complex, probe=1.4)
        sampler = SurfaceSampler(
            envelope=envelope, excluder=excluder, max_rejections=20
        )
        with pytest.raises(SamplingError, match="buried|reach|probe"):
            sampler.generator()


class TestMakeSampler:
    def test_returns_surface_sampler(self, synthetic_octahedron_complex):
        from maws.space import SurfaceSampler, make_sampler

        s = make_sampler("shell", synthetic_octahedron_complex, reach=10.0, probe=1.4)
        assert isinstance(s, SurfaceSampler)

    def test_chooses_correct_envelope_per_shape(self, synthetic_octahedron_complex):
        from maws.space import Cube, Shell, Sphere, make_sampler

        assert isinstance(
            make_sampler("cube", synthetic_octahedron_complex).envelope, Cube
        )
        assert isinstance(
            make_sampler("sphere", synthetic_octahedron_complex).envelope, Sphere
        )
        assert isinstance(
            make_sampler("shell", synthetic_octahedron_complex).envelope, Shell
        )

    def test_rejects_unknown_shape(self, synthetic_octahedron_complex):
        import pytest

        from maws.space import make_sampler

        with pytest.raises(ValueError, match="cube|sphere|shell"):
            make_sampler("blob", synthetic_octahedron_complex)

    def test_passes_dims_through(self, synthetic_octahedron_complex):
        """For the octahedron with reach=10, sphere radius should be 15."""
        from maws.space import make_sampler

        s = make_sampler("sphere", synthetic_octahedron_complex, reach=10.0)
        assert s.envelope.radius == 15.0

    def test_returns_clear_samples(self, synthetic_octahedron_complex):
        from maws.space import make_sampler

        np.random.seed(0)
        s = make_sampler("shell", synthetic_octahedron_complex)
        for _ in range(20):
            sample = s.generator()
            assert s.excluder.is_clear(sample[:3])

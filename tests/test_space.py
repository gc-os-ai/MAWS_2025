"""
Tests for maws.space module.

Surface-aware sampling primitives:
- Sphere : envelope dataclass (only built-in)
- NAngles : torsion angle sampling (U(1)^N)
- Excluder : SAS-style rejection
- SurfaceSampler : envelope + Excluder composer
- make_sampler : factory
"""

import math

import numpy as np

from maws.space import _BONDI_VDW_RADII, _DEFAULT_VDW, NAngles


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


class TestBondiTable:
    def test_has_common_biological_elements(self):
        for sym in ("H", "C", "N", "O", "P", "S"):
            assert sym in _BONDI_VDW_RADII

    def test_carbon_value(self):
        # Bondi 1964: C = 1.70 Å
        assert _BONDI_VDW_RADII["C"] == 1.70

    def test_default_vdw_is_carbon_like(self):
        assert _DEFAULT_VDW == 1.70


class TestSphere:
    def test_generator_returns_sample(self):
        from maws.space import Sample, Sphere

        s = Sphere(radius=5.0, centre=np.array([0.0, 0.0, 0.0])).generator()
        assert isinstance(s, Sample)

    def test_generator_within_radius_at_origin(self):
        from maws.space import Sphere

        s = Sphere(radius=5.0, centre=np.array([0.0, 0.0, 0.0]))
        for _ in range(50):
            sample = s.generator()
            assert np.linalg.norm(sample.position) <= 5.0 + 1e-9

    def test_generator_offset_by_centre(self):
        """Bug-fix from PR #38 carries over: samples must be offset by centre."""
        from maws.space import Sphere

        centre = np.array([50.0, -30.0, 12.0])
        s = Sphere(radius=10.0, centre=centre)
        for _ in range(50):
            sample = s.generator()
            assert np.linalg.norm(sample.position - centre) <= 10.0 + 1e-9

    def test_radial_distribution_volume_correct(self):
        """E[r] = 3R/4 for uniform-in-volume sampling (= 3.75 for R=5)."""
        from maws.space import Sphere

        np.random.seed(0)
        s = Sphere(radius=5.0, centre=np.array([0.0, 0.0, 0.0]))
        rs = np.array([np.linalg.norm(s.generator().position) for _ in range(10_000)])
        assert abs(rs.mean() - 3.75) < 0.05

    def test_direction_uniform_on_sphere(self):
        """E[(z/r)^2] = 1/3 for uniform direction (vs 1/2 if biased to poles)."""
        from maws.space import Sphere

        np.random.seed(1)
        s = Sphere(radius=5.0, centre=np.array([0.0, 0.0, 0.0]))
        samples = np.array([s.generator().position for _ in range(10_000)])
        rs = np.linalg.norm(samples, axis=1, keepdims=True)
        unit_z = (samples / rs)[:, 2]
        assert abs(float((unit_z**2).mean()) - 1 / 3) < 0.04


class TestRandomUnitAxis:
    """Regression: ``_random_unit_axis`` must be isotropic on the unit sphere.

    A previous implementation drew uniformly from the cube
    ``[-1, 1]^3`` and normalised; that biased directions toward the
    cube's eight corners and away from its face-centres (the cube's
    corners are at distance √3 from the origin, faces at distance 1,
    so cube interior volume that normalises near a corner is larger
    than that which normalises near a face-centre).

    The fix is the standard Gaussian-then-normalise recipe (the
    multivariate standard normal is spherically symmetric).
    """

    def test_axes_are_unit_length(self):
        from maws.space import _random_unit_axis

        np.random.seed(0)
        for _ in range(50):
            axis = _random_unit_axis()
            assert abs(np.linalg.norm(axis) - 1.0) < 1e-9

    def test_isotropic_face_vs_corner_directions(self):
        """Statistic that cleanly separates uniform from cube-biased.

        Classify each axis by its absolute-value components:

        - ``face``: one component dominates, ``max|.|  > 0.85`` (axis
          lies near a ±x/±y/±z direction).
        - ``corner``: all three components similar size,
          ``min|.| > 0.45`` (axis lies near a (±1, ±1, ±1)/√3
          direction).

        At N = 50_000:
        - Uniform on sphere: face fraction ≈ 0.45, corner ≈ 0.07
        - Cube-then-normalise (buggy): face ≈ 0.30, corner ≈ 0.12

        A threshold of face ≥ 0.40 cleanly distinguishes the two.
        """
        from maws.space import _random_unit_axis

        np.random.seed(0)
        n = 50_000
        axes = np.abs(np.array([_random_unit_axis() for _ in range(n)]))
        max_per_axis = axes.max(axis=1)
        min_per_axis = axes.min(axis=1)
        face_fraction = float((max_per_axis > 0.85).mean())
        corner_fraction = float((min_per_axis > 0.45).mean())
        # Uniform target ≈ 0.45 / 0.07; the buggy cube-normalise gives
        # 0.30 / 0.12. The bands below cleanly separate them.
        assert face_fraction > 0.40, (
            f"face_fraction = {face_fraction:.3f}, expected > 0.40 "
            f"(cube-bias produces ~0.30)"
        )
        assert corner_fraction < 0.10, (
            f"corner_fraction = {corner_fraction:.3f}, expected < 0.10 "
            f"(cube-bias produces ~0.12)"
        )


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
    def test_sphere_dims_octahedron(self, synthetic_octahedron_complex):
        from maws.space import compute_envelope_dims

        d = compute_envelope_dims(synthetic_octahedron_complex, reach=10.0)
        # COM = origin; R_max = 5.0; radius = R_max + reach = 15.
        np.testing.assert_allclose(d["centre"], [0.0, 0.0, 0.0])
        assert d["radius"] == 15.0

    def test_radius_scales_with_reach(self, synthetic_octahedron_complex):
        from maws.space import compute_envelope_dims

        d_small = compute_envelope_dims(synthetic_octahedron_complex, reach=2.0)
        d_large = compute_envelope_dims(synthetic_octahedron_complex, reach=20.0)
        assert d_large["radius"] - d_small["radius"] == 18.0


class TestMakeSamplerValidation:
    def test_rejects_negative_reach(self, synthetic_octahedron_complex):
        import pytest

        from maws.space import make_sampler

        with pytest.raises(ValueError, match="reach must be >= 0"):
            make_sampler(synthetic_octahedron_complex, reach=-1.0)

    def test_rejects_negative_probe(self, synthetic_octahedron_complex):
        import pytest

        from maws.space import make_sampler

        with pytest.raises(ValueError, match="probe must be >= 0"):
            make_sampler(synthetic_octahedron_complex, probe=-1.0)


class TestSurfaceSampler:
    def test_all_samples_clear(self, synthetic_octahedron_complex):
        from maws.space import Excluder, Sphere, SurfaceSampler

        envelope = Sphere(radius=15.0, centre=np.array([0.0, 0.0, 0.0]))
        excluder = Excluder(synthetic_octahedron_complex, probe=1.4)
        sampler = SurfaceSampler(envelope=envelope, excluder=excluder)
        for _ in range(50):
            sample = sampler.generator()
            assert excluder.is_clear(sample.position)

    def test_raises_when_envelope_buried(self, synthetic_two_carbon_complex):
        """A tiny sphere sitting on an atom is fully blocked."""
        import pytest

        from maws.space import (
            Excluder,
            SamplingError,
            Sphere,
            SurfaceSampler,
        )

        envelope = Sphere(radius=0.05, centre=np.array([0.0, 0.0, 0.0]))
        excluder = Excluder(synthetic_two_carbon_complex, probe=1.4)
        sampler = SurfaceSampler(
            envelope=envelope, excluder=excluder, max_rejections=20
        )
        with pytest.raises(SamplingError, match="buried|reach|probe"):
            sampler.generator()


class TestMakeSampler:
    def test_returns_surface_sampler_with_sphere_envelope(
        self, synthetic_octahedron_complex
    ):
        from maws.space import Sphere, SurfaceSampler, make_sampler

        s = make_sampler(synthetic_octahedron_complex, reach=10.0, probe=1.4)
        assert isinstance(s, SurfaceSampler)
        assert isinstance(s.envelope, Sphere)

    def test_passes_dims_through(self, synthetic_octahedron_complex):
        """For the octahedron with reach=10, sphere radius should be 15."""
        from maws.space import make_sampler

        s = make_sampler(synthetic_octahedron_complex, reach=10.0)
        assert s.envelope.radius == 15.0

    def test_returns_clear_samples(self, synthetic_octahedron_complex):
        from maws.space import make_sampler

        np.random.seed(0)
        s = make_sampler(synthetic_octahedron_complex)
        for _ in range(20):
            sample = s.generator()
            assert s.excluder.is_clear(sample.position)


class TestSurfaceFollowingSampler:
    def test_generator_returns_sample(self, synthetic_octahedron_complex):
        from maws.space import Sample, SurfaceFollowingSampler

        np.random.seed(0)
        sf = SurfaceFollowingSampler(synthetic_octahedron_complex, d_max=6.0, probe=1.4)
        sample = sf.generator()
        assert isinstance(sample, Sample)
        assert sample.position.shape == (3,)
        assert sample.axis.shape == (3,)
        assert isinstance(sample.angle, float)

    def test_samples_within_d_max_of_some_atom(self, synthetic_octahedron_complex):
        from scipy.spatial import KDTree

        from maws.helpers import nostrom
        from maws.space import SurfaceFollowingSampler

        np.random.seed(0)
        d_max = 4.0
        positions = np.asarray(
            nostrom(synthetic_octahedron_complex.positions), dtype=float
        )
        tree = KDTree(positions)
        sf = SurfaceFollowingSampler(
            synthetic_octahedron_complex, d_max=d_max, probe=1.4
        )
        for _ in range(20):
            sample = sf.generator()
            nearest, _ = tree.query(sample.position, k=1)
            assert nearest <= d_max + 1e-9, (
                f"sample at {sample.position} is {nearest:.2f} Å "
                f"from nearest atom, exceeds d_max={d_max}"
            )

    def test_samples_are_sas_clear(self, synthetic_octahedron_complex):
        from maws.space import Excluder, SurfaceFollowingSampler

        np.random.seed(0)
        excluder = Excluder(synthetic_octahedron_complex, probe=1.4)
        sf = SurfaceFollowingSampler(synthetic_octahedron_complex, d_max=6.0, probe=1.4)
        for _ in range(20):
            sample = sf.generator()
            assert excluder.is_clear(sample.position)

    def test_rejects_nonpositive_d_max(self, synthetic_octahedron_complex):
        import pytest

        from maws.space import SurfaceFollowingSampler

        with pytest.raises(ValueError, match="d_max must be > 0"):
            SurfaceFollowingSampler(synthetic_octahedron_complex, d_max=0.0)
        with pytest.raises(ValueError, match="d_max must be > 0"):
            SurfaceFollowingSampler(synthetic_octahedron_complex, d_max=-1.0)


class TestMakeSamplerModes:
    def test_default_is_sphere(self, synthetic_octahedron_complex):
        from maws.space import SurfaceSampler, make_sampler

        s = make_sampler(synthetic_octahedron_complex)
        assert isinstance(s, SurfaceSampler)

    def test_explicit_sphere_mode(self, synthetic_octahedron_complex):
        from maws.space import SurfaceSampler, make_sampler

        s = make_sampler(synthetic_octahedron_complex, mode="sphere", reach=5.0)
        assert isinstance(s, SurfaceSampler)

    def test_surface_following_mode(self, synthetic_octahedron_complex):
        from maws.space import SurfaceFollowingSampler, make_sampler

        sf = make_sampler(
            synthetic_octahedron_complex, mode="surface-following", d_max=6.0
        )
        assert isinstance(sf, SurfaceFollowingSampler)

    def test_unknown_mode_raises(self, synthetic_octahedron_complex):
        import pytest

        from maws.space import make_sampler

        with pytest.raises(ValueError, match="sphere|surface-following"):
            make_sampler(synthetic_octahedron_complex, mode="bogus")

    def test_negative_probe_rejected_for_both_modes(self, synthetic_octahedron_complex):
        import pytest

        from maws.space import make_sampler

        with pytest.raises(ValueError, match="probe must be >= 0"):
            make_sampler(synthetic_octahedron_complex, mode="sphere", probe=-1.0)
        with pytest.raises(ValueError, match="probe must be >= 0"):
            make_sampler(
                synthetic_octahedron_complex,
                mode="surface-following",
                probe=-1.0,
            )

    def test_both_modes_yield_clear_samples(self, synthetic_octahedron_complex):
        from maws.space import Excluder, make_sampler

        excluder = Excluder(synthetic_octahedron_complex, probe=1.4)
        for mode in ("sphere", "surface-following"):
            np.random.seed(0)
            s = make_sampler(synthetic_octahedron_complex, mode=mode)
            for _ in range(10):
                assert excluder.is_clear(s.generator().position)

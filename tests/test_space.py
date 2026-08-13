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

        rng = np.random.default_rng(0)
        for _ in range(50):
            axis = _random_unit_axis(rng)
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

        rng = np.random.default_rng(0)
        n = 50_000
        axes = np.abs(np.array([_random_unit_axis(rng) for _ in range(n)]))
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


class TestClashFilter:
    """Tests for ClashFilter, the check on the atoms that actually get placed.

    The sampler's own filter asks whether a water-sized probe could sit at the
    sampled point. What gets placed there is a whole nucleotide, so a point
    can pass that check while the atoms put there overlap the target. This
    filter is the check on those atoms.

    Every atom in the fixture is carbon, Bondi radius 1.70 A. A moving atom
    and a rigid atom therefore count as clashing once they are closer than
    ``1.70 + 1.70 - tolerance``.

    See issue #48.
    """

    ELEMENT = [0, 1, 2]  # atoms 0-1 move, everything from index 2 is rigid

    @staticmethod
    def moved(x):
        """Return full positions with the moving pair starting at ``x`` on the x-axis.

        The two rigid atoms stay at x = 50 and x = 51, so the closest
        moving-rigid pair is atom 1 at ``x + 1`` against the rigid atom at 50.
        """
        return np.array(
            [[x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [50.0, 0.0, 0.0], [51.0, 0.0, 0.0]]
        )

    def test_clear_when_far_apart(self, synthetic_docking_complex):
        """A strand placed far from the target is accepted."""
        from maws.space import ClashFilter

        f = ClashFilter(synthetic_docking_complex, self.ELEMENT)
        assert f.is_clear(self.moved(0.0))

    def test_blocked_when_placed_on_top_of_a_rigid_atom(
        self, synthetic_docking_complex
    ):
        """A strand sharing coordinates with a target atom is rejected.

        This is the case the energy blow-ups come from: two atoms at the same
        place give a potential energy many orders of magnitude above a real
        pose, and that pose then wins on the selection score.
        """
        from maws.space import ClashFilter

        f = ClashFilter(synthetic_docking_complex, self.ELEMENT)
        assert not f.is_clear(self.moved(50.0))

    def test_the_cut_is_the_summed_radii_less_the_tolerance(
        self, synthetic_docking_complex
    ):
        """Two carbons are rejected below 3.00 A and accepted above it.

        Both radii are 1.70 A and the tolerance is 0.4 A, so the cut sits at
        1.70 + 1.70 - 0.4 = 3.00 A. Testing either side of it pins down that
        the filter compares against summed radii rather than a flat distance.
        """
        from maws.space import ClashFilter

        f = ClashFilter(synthetic_docking_complex, self.ELEMENT, tolerance=0.4)
        assert not f.is_clear(self.moved(46.1))  # atom 1 is 2.9 A away
        assert f.is_clear(self.moved(45.9))  # atom 1 is 3.1 A away

    def test_a_larger_tolerance_accepts_a_closer_contact(
        self, synthetic_docking_complex
    ):
        """Raising the tolerance lets through a contact a lower one rejects.

        The tolerance is how far two atoms may overlap before the pose is
        thrown away. Real complexes need some overlap allowed: a hydrogen bond
        puts a hydrogen about 0.8 A inside the summed radii.
        """
        from maws.space import ClashFilter

        strict = ClashFilter(synthetic_docking_complex, self.ELEMENT, tolerance=0.4)
        loose = ClashFilter(synthetic_docking_complex, self.ELEMENT, tolerance=1.5)
        assert not strict.is_clear(self.moved(46.1))
        assert loose.is_clear(self.moved(46.1))

    def test_two_rigid_atoms_overlapping_is_not_a_clash(
        self, synthetic_docking_complex
    ):
        """Only the moving atoms are tested, so the target's own geometry passes.

        The target comes from the input PDB and cannot be redrawn. Judging it
        would reject every pose in the run rather than the bad ones.
        """
        from maws.space import ClashFilter

        f = ClashFilter(synthetic_docking_complex, self.ELEMENT)
        positions = self.moved(0.0)
        positions[3] = positions[2]
        assert f.is_clear(positions)

    def test_rejects_negative_tolerance(self, synthetic_docking_complex):
        """A negative tolerance raises, rather than silently widening the cut.

        It would mean demanding a gap between the strand and the target that
        no bound pose can satisfy, so it is always a mistake.
        """
        import pytest

        from maws.space import ClashFilter

        with pytest.raises(ValueError, match="tolerance"):
            ClashFilter(synthetic_docking_complex, self.ELEMENT, tolerance=-1.0)


class _StubSampler:
    """Draws the same pose every time and counts how often it was asked."""

    def __init__(self):
        self.draws = 0

    def generator(self):
        from maws.space import Sample

        self.draws += 1
        return Sample(position=np.zeros(3), axis=np.array([0.0, 0.0, 1.0]), angle=0.0)


class _StubRotations:
    """Stands in for NAngles, handing back a fixed set of torsion angles."""

    def __init__(self, n=4):
        self.n = n

    def generator(self):
        return np.zeros(self.n)


class _StubChain:
    """Records the torsions applied to it. ``element`` covers atoms 0-1."""

    element = [0, 1, 2]

    def __init__(self):
        self.bends = 0

    def rotate_in_residue(self, residue, torsion, angle):
        self.bends += 1


class _StubComplex:
    """Moves atom 0 one A along x each time it is placed.

    ``offsets`` records where atom 0 sat at the start of every attempt, which
    is how a test sees whether the strand was reset between attempts.
    """

    def __init__(self):
        import openmm as mm
        from openmm import unit as ommunit

        self._mm = mm
        self._unit = ommunit
        self.positions = [mm.Vec3(0.0, 0.0, 0.0)] * 3 * ommunit.angstrom
        self.offsets = []

    def place_global(self, element, position, axis, angle):
        from maws.helpers import nostrom

        self.offsets.append(float(np.asarray(nostrom(self.positions))[0][0]))
        moved = self.positions[:]
        moved[0] += self._mm.Vec3(1.0, 0.0, 0.0) * self._unit.angstrom
        self.positions = moved


class _StubClash:
    """Rejects the first `rejections` conformations, then accepts everything.

    ``bends_at_check`` records how many torsions had been applied by the time
    each decision was made, which is how a test sees the order of operations.
    """

    def __init__(self, rejections, chain=None):
        self.rejections = rejections
        self.chain = chain
        self.bends_at_check = []

    def is_clear(self, positions):
        if self.chain is not None:
            self.bends_at_check.append(self.chain.bends)
        if self.rejections:
            self.rejections -= 1
            return False
        return True


class TestDrawClearConformation:
    """Tests for draw_clear_conformation, the reject-and-redraw loop.

    See issue #48.
    """

    def test_a_clear_conformation_is_returned_on_the_first_draw(self):
        """Nothing is redrawn when the first attempt already clears the target."""
        from maws.space import draw_clear_conformation

        sampler = _StubSampler()
        draw_clear_conformation(
            _StubComplex(), _StubChain(), sampler, _StubRotations(), _StubClash(0)
        )
        assert sampler.draws == 1

    def test_a_rejected_conformation_is_drawn_again(self):
        """Each rejection costs one more pose and one more set of torsions.

        A rejected attempt must not be scored, and must not be retried
        unchanged, or the loop would either keep the clash or never end.
        """
        from maws.space import draw_clear_conformation

        sampler, chain, rotations = _StubSampler(), _StubChain(), _StubRotations(4)
        draw_clear_conformation(
            _StubComplex(), chain, sampler, rotations, _StubClash(3)
        )
        assert sampler.draws == 4
        assert chain.bends == 4 * 4

    def test_the_strand_is_reset_before_each_attempt(self):
        """Every attempt starts from the coordinates the strand came in with.

        Without the reset each attempt would build on the rejected one before
        it, so the accepted conformation would be a composition of failures
        rather than the single draw the sampler reports.
        """
        from maws.space import draw_clear_conformation

        cx = _StubComplex()
        draw_clear_conformation(
            cx, _StubChain(), _StubSampler(), _StubRotations(), _StubClash(3)
        )
        assert cx.offsets == [0.0, 0.0, 0.0, 0.0]

    def test_the_torsions_are_applied_before_the_clash_check(self):
        """The check sees the bent strand, not the freshly placed one.

        Placing a strand clear of the target does not keep it clear: the
        torsions swing atoms far enough to reach back into it, so a check made
        before them would pass conformations that clash.
        """
        from maws.space import draw_clear_conformation

        chain = _StubChain()
        clash = _StubClash(0, chain=chain)
        draw_clear_conformation(
            _StubComplex(), chain, _StubSampler(), _StubRotations(4), clash
        )
        assert clash.bends_at_check == [4]

    def test_giving_up_raises(self):
        """A target nothing can clear fails loudly instead of looping forever."""
        import pytest

        from maws.space import SamplingError, draw_clear_conformation

        with pytest.raises(SamplingError, match="attempts"):
            draw_clear_conformation(
                _StubComplex(),
                _StubChain(),
                _StubSampler(),
                _StubRotations(),
                _StubClash(999),
                max_rejections=5,
            )


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

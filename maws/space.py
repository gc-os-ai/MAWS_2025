"""
maws.space
==========

Surface-aware sampling for the MAWS aptamer-design loop.

The single entry point is :func:`make_sampler`. It returns a sampler
whose ``.generator()`` yields candidate poses (:class:`Sample`) that
survive the configured filters. Two modes are available:

- ``mode="sphere"`` (default): bounding-sphere envelope + SAS rejection.
  Simple, fast, well-tested. Returns a :class:`SurfaceSampler`.
- ``mode="surface-following"`` (opt-in): also caps how far a candidate
  may sit from the nearest protein atom, so accepted samples concentrate
  near the molecular surface. Returns a :class:`SurfaceFollowingSampler`.

Public API
----------
Sphere : frozen dataclass
    Solid sphere envelope. Volume-correct uniform sampling.
NAngles : frozen dataclass
    N independent angles in [0, 2π) for in-residue rotations.
Excluder
    KDTree-backed SAS rejection: a candidate is "clear" iff its distance
    to every ligand atom exceeds (Bondi vdW + probe).
SurfaceSampler
    Composes an envelope + an Excluder into a rejection sampler.
SurfaceFollowingSampler
    Alternative sampler that concentrates poses near the molecular
    surface via an outer distance cap.
make_sampler(complex_obj, *, mode="sphere", reach=10.0, d_max=6.0, probe=1.4)
    Factory: builds the requested sampler for the given Complex.
compute_envelope_dims(complex_obj, reach)
    Returns the kwargs for the Sphere dataclass; useful in tests and
    notebooks.


Examples
--------
>>> from maws.space import NAngles
>>> sample = NAngles(n=4).generator()
>>> len(sample)
4
"""

import warnings
from dataclasses import dataclass
from typing import Literal, Protocol

import numpy as np
from openmm import unit
from scipy.spatial import KDTree

from maws.helpers import mass_weighted_center, nostrom

# Bondi (1964) "Van der Waals Volumes and Radii", J. Phys. Chem. 68(3):441-451.
# These are the values used by Chimera, PyMOL, FreeSASA. Hardcoded rather
# than pulled from a library to (a) avoid an extra dependency for a 25-line
# table and (b) make the surface a pure geometric concept, decoupled from
# whichever force field the run is using.
# Values in Å. Indexed by element symbol (case-sensitive).
_BONDI_VDW_RADII: dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
    "Si": 2.10,
    "B": 1.92,
    "Se": 1.90,
    "Na": 2.27,
    "Mg": 1.73,
    "K": 2.75,
    "Ca": 2.31,
    "Fe": 2.00,
    "Zn": 1.39,
    "Cu": 1.40,
    "Mn": 2.00,
    "Ni": 1.63,
}
# Carbon-equivalent fallback for unknown elements.
_DEFAULT_VDW = 1.70


@dataclass(frozen=True)
class Sample:
    """One pose drawn from an envelope: a position + a rotation."""

    position: np.ndarray  # shape (3,), Å
    axis: np.ndarray  # shape (3,), unit vector
    angle: float  # radians


def _random_unit_axis() -> np.ndarray:
    """Random unit vector uniformly distributed on the unit sphere.

    Uses the standard Gaussian-then-normalize recipe: the multivariate
    standard normal has spherically symmetric density, so normalizing
    yields a point uniform on the unit sphere. A naive uniform draw in
    the cube ``[-1, 1]^3`` followed by normalization would
    over-represent the cube's corner-aligned directions because the
    cube has greater radial extent along its diagonals than its faces.
    """
    v = np.random.standard_normal(3)
    return v / np.linalg.norm(v)


class Envelope(Protocol):
    """Anything with ``.generator() -> Sample`` is an envelope.

    :class:`Sphere` is the only built-in. Users may pass a custom envelope
    directly to :class:`SurfaceSampler` if they want a different
    geometric region.
    """

    def generator(self) -> Sample: ...


def _spherical_sample(centre: np.ndarray, r: float) -> np.ndarray:
    """Random point on a sphere of radius ``r`` around ``centre``
    (uniform direction)."""
    phi = np.random.uniform(0, 2 * np.pi)
    cos_psi = np.random.uniform(-1, 1)  # uniform on the sphere
    sin_psi = np.sqrt(1.0 - cos_psi * cos_psi)
    return np.array(
        [
            centre[0] + r * np.cos(phi) * sin_psi,
            centre[1] + r * np.sin(phi) * sin_psi,
            centre[2] + r * cos_psi,
        ]
    )


@dataclass(frozen=True)
class Sphere:
    """
    Solid sphere of given ``radius`` centered on ``centre`` (Å).

    Sampling
    --------
    :meth:`generator` returns a :class:`Sample` drawn uniformly **in
    volume** — every cubic ångström inside the sphere is equally likely.
    Implemented in two stages:

    1. Radial draw, ``r = radius · u^(1/3)`` where ``u ~ U(0, 1)``.
       The cube-root is the inverse of the volume CDF
       ``F(r) = (r/R)³``; without it the naive ``r ~ U(0, R)`` would
       concentrate samples toward the centre because outer shells have
       larger area.
    2. Direction draw on the unit sphere via :func:`_spherical_sample`:
       ``cos(ψ) ~ U(-1, 1)`` and ``φ ~ U(0, 2π)``, then
       ``(sin ψ cos φ, sin ψ sin φ, cos ψ)`` is the unit direction.
       This is the standard recipe for uniform points on a sphere; a
       naive ``ψ ~ U(0, π)`` would over-sample the poles.

    The final position is ``centre + r · direction``. ``axis`` is an
    independent random unit vector; ``angle`` is ``U(0, 2π)``.
    """

    radius: float
    centre: np.ndarray

    def generator(self) -> Sample:
        # Volume-correct: r = R · u^(1/3) for u ~ U(0,1).
        r = self.radius * np.random.uniform(0, 1) ** (1 / 3)
        return Sample(
            position=_spherical_sample(self.centre, r),
            axis=_random_unit_axis(),
            angle=float(np.random.uniform(0, 2 * np.pi)),
        )


@dataclass(frozen=True)
class NAngles:
    """N independent angles drawn from [0, 2π). Used for in-residue rotations."""

    n: int

    def generator(self) -> np.ndarray:
        return np.random.uniform(0, 2 * np.pi, self.n)


class Excluder:
    """
    SAS-style point-vs-protein membership test.

    What is the SAS?
    ----------------
    The **Solvent-Accessible Surface** (Lee & Richards, 1971) is defined by
    rolling a small spherical *probe* — historically a water molecule of
    radius ≈ 1.4 Å — over the union of the protein's atomic van der Waals
    spheres. The set of positions the probe's *centre* can occupy without
    overlapping any atom is the *accessible region*; its boundary is the
    SAS. Every common molecular-graphics tool (Chimera, PyMOL, FreeSASA)
    draws this surface with the same definition.

    Decision rule
    -------------
    For a candidate point ``p``, this class returns ``is_clear(p) = True``
    iff a probe-sized sphere could physically sit at ``p`` strictly
    outside every protein atom:

        ``dist(p, atom_i) > vdW(atom_i) + probe``  for every atom *i*.

    The strict inequality is a measure-zero convention for continuous
    rejection sampling: a point lying exactly on the SAS boundary is
    treated as blocked. With float64 positions this never matters in
    practice, but the docstring matches the code.

    The protein atomic vdW radii come from the Bondi (1964) table at the
    top of this module; unknown elements fall back to the carbon-equivalent
    ``_DEFAULT_VDW`` (1.70 Å) and emit a one-shot warning.


    Parameters
    ----------
    complex_obj
        Anything with ``.positions`` (openmm Quantity, Å-convertible,
        shape (N, 3)) and ``.topology.atoms()`` (callable yielding atoms
        whose ``.element.symbol`` is a string).
    probe : float, default 1.4
        Probe radius in Å. 1.4 Å is the water-equivalent SAS convention.
        Larger values are more conservative (only big pockets accessible);
        smaller values are more permissive (allows tighter fits).
    """

    # Process-wide: one warning per unknown element symbol.
    _warned: set[str] = set()

    def __init__(self, complex_obj, probe: float = 1.4):
        positions = np.asarray(nostrom(complex_obj.positions), dtype=float)
        radii = np.empty(len(positions), dtype=float)
        for i, atom in enumerate(complex_obj.topology.atoms()):
            sym = atom.element.symbol
            r = _BONDI_VDW_RADII.get(sym)
            if r is None:
                r = _DEFAULT_VDW
                if sym not in Excluder._warned:
                    Excluder._warned.add(sym)
                    warnings.warn(
                        f"Unknown element {sym!r} - using fallback "
                        f"vdW = {_DEFAULT_VDW} Å",
                        RuntimeWarning,
                        stacklevel=2,
                    )
            radii[i] = r
        self._positions = positions
        self._inflated = radii + probe
        self._max_inflated = float(self._inflated.max())
        self._tree = KDTree(positions)

    def is_clear(self, point: np.ndarray) -> bool:
        """True iff ``point`` lies outside every (vdW + probe) sphere."""
        idx = self._tree.query_ball_point(point, self._max_inflated)
        if not idx:
            return True
        diffs = self._positions[idx] - np.asarray(point)
        dists2 = (diffs**2).sum(axis=1)
        return bool((dists2 > self._inflated[idx] ** 2).all())


def _atom_mass_in_dalton(atom) -> float:
    m = atom.element.mass
    if hasattr(m, "value_in_unit"):
        return m.value_in_unit(unit.dalton)
    return float(m)


def compute_envelope_dims(complex_obj, reach: float) -> dict:
    """
    Compute auto-sized sphere envelope dimensions from the ligand atoms.

    Parameters
    ----------
    complex_obj
        Object with ``.positions`` (Quantity) and ``.topology.atoms()``
        (yielding atoms with ``.element.mass``).
    reach : float
        How far the envelope extends past the ligand's bounding radius (Å).

    Returns
    -------
    dict
        Kwargs for :class:`Sphere`: ``{"radius": R_max + reach, "centre": COM}``.
    """
    pos = np.asarray(nostrom(complex_obj.positions), dtype=float)
    masses = np.array(
        [_atom_mass_in_dalton(a) for a in complex_obj.topology.atoms()],
        dtype=float,
    )
    com = mass_weighted_center(pos, masses)
    dists = np.linalg.norm(pos - com, axis=1)
    return {"radius": float(dists.max()) + reach, "centre": com}


class SamplingError(RuntimeError):
    """Raised when SurfaceSampler cannot find a clear point in max_rejections tries."""


@dataclass
class SurfaceSampler:
    """
    Two-layer rejection sampler that combines an envelope with an
    :class:`Excluder`.

    Each call to :meth:`generator` applies two filters in order:

    1. **Envelope filter (cheap, geometric).** The envelope's own
       ``.generator()`` (built-in: :class:`Sphere`; custom envelopes are
       allowed as long as they satisfy :class:`Envelope`) draws a
       candidate ``Sample`` from its geometric region around the ligand.


    2. **SAS filter (precise, atom-aware).** The :class:`Excluder` checks
       the candidate's ``position`` against every protein atom's
       ``vdW + probe`` sphere. If the candidate would put an aptamer atom
       too close to any protein atom (steric clash by the SAS rule), it
       is rejected and we draw again. Cost: O(log N_atoms) per check
       (~50 µs for a ~3000-atom protein, dominated by the KDTree query).

    The loop continues until a candidate passes both filters or
    ``max_rejections`` attempts have failed; the latter raises
    :class:`SamplingError`. The cap exists so a fully-buried envelope
    (mis-sized, wrong centre, ``reach`` too small) fails fast with an
    actionable message instead of looping forever.

    Parameters
    ----------
    envelope
        Any object with ``.generator() -> Sample`` (built-in:
        :class:`Sphere`).
    excluder
        Configured :class:`Excluder` instance.
    max_rejections : int, default 1000
        Hard cap on consecutive rejected attempts before
        :class:`SamplingError` is raised.
    """

    envelope: Envelope
    excluder: Excluder
    max_rejections: int = 1000

    def generator(self) -> Sample:
        for _ in range(self.max_rejections):
            sample = self.envelope.generator()
            if self.excluder.is_clear(sample.position):
                return sample
        raise SamplingError(
            f"Could not draw a clear point in {self.max_rejections} attempts. "
            f"Envelope may be fully buried - increase --reach, decrease --probe, "
            f"or check ligand size."
        )


class SurfaceFollowingSampler:
    """
    Opt-in alternative to :class:`SurfaceSampler` that concentrates
    candidate poses near the molecular surface.

    Region
    ------
    Samples uniformly in the band

        ``{ p : p is SAS-clear AND dist(p, nearest_protein_atom) <= d_max }``

    i.e. a layer of thickness ``d_max`` wrapping the protein's atomic
    surface (including following into pockets, since "near any atom" is
    direction-agnostic), capped on the outside so we don't sample far
    solvent.

    Algorithm (per call to :meth:`generator`)
    -----------------------------------------
    1. Draw a uniform-in-volume point in the bounding sphere of radius
       ``R_max + d_max`` around the mass-weighted COM (same volume-correct
       cube-root radial law as :class:`Sphere`).
    2. Query a KDTree of protein atom positions for the nearest atom
       distance. Reject if greater than ``d_max`` (too far in solvent).
    3. Run an :class:`Excluder` SAS check. Reject if inside protein bulk.
    4. Repeat until a candidate passes both filters, or
       ``max_rejections`` attempts have failed (raises
       :class:`SamplingError`).

    Parameters
    ----------
    complex_obj
        Built ligand-only ``Complex`` (positions + topology).
    d_max : float, default 6.0
        Band thickness in Å. Smaller values concentrate samples
        tighter to the surface but raise rejection rate further.
    probe : float, default 1.4
        Probe radius for the SAS rejection (Å). Water-equivalent at
        the default.
    max_rejections : int, default 50000
        Hard cap on consecutive rejected attempts before
        :class:`SamplingError` is raised. The cap is higher than
        :class:`SurfaceSampler`'s because the band's rejection rate is
        higher.
    """

    def __init__(
        self,
        complex_obj,
        *,
        d_max: float = 6.0,
        probe: float = 1.4,
        max_rejections: int = 50_000,
    ):
        if d_max <= 0:
            raise ValueError(f"d_max must be > 0, got {d_max}")
        positions = np.asarray(nostrom(complex_obj.positions), dtype=float)
        masses = np.array(
            [_atom_mass_in_dalton(a) for a in complex_obj.topology.atoms()],
            dtype=float,
        )
        com = mass_weighted_center(positions, masses)
        R_max = float(np.linalg.norm(positions - com, axis=1).max())
        self._com = com
        self._R_bound = R_max + d_max
        self._tree = KDTree(positions)
        self._excluder = Excluder(complex_obj, probe=probe)
        self._d_max = d_max
        self._max_rejections = max_rejections

    def generator(self) -> Sample:
        for _ in range(self._max_rejections):
            # Uniform-in-volume draw inside the bounding sphere.
            r = self._R_bound * np.random.uniform(0, 1) ** (1 / 3)
            phi = np.random.uniform(0, 2 * np.pi)
            cos_psi = np.random.uniform(-1, 1)
            sin_psi = np.sqrt(1.0 - cos_psi * cos_psi)
            position = self._com + r * np.array(
                [np.cos(phi) * sin_psi, np.sin(phi) * sin_psi, cos_psi]
            )
            d_near, _ = self._tree.query(position, k=1)
            if d_near > self._d_max:
                continue  # too far from any atom
            if not self._excluder.is_clear(position):
                continue  # inside protein bulk
            return Sample(
                position=position,
                axis=_random_unit_axis(),
                angle=float(np.random.uniform(0, 2 * np.pi)),
            )
        raise SamplingError(
            f"Could not draw a surface-following sample in "
            f"{self._max_rejections} attempts. d_max may be too small or "
            f"the protein may be unusually dense."
        )


def make_sampler(
    complex_obj,
    *,
    mode: Literal["sphere", "surface-following"] = "sphere",
    reach: float = 10.0,
    d_max: float = 6.0,
    probe: float = 1.4,
):
    """
    Build a fully-configured surface-aware sampler for ``complex_obj``.

    Two sampling modes are available:

    - ``mode="sphere"`` (default): draws candidates uniformly in a
      bounding sphere around the ligand and rejects those inside the
      protein bulk via an SAS check. Returns a
      :class:`SurfaceSampler`. This is the simple, fast, well-tested
      default.

    - ``mode="surface-following"`` (opt-in): also rejects candidates
      that are more than ``d_max`` Å from any protein atom, so accepted
      samples concentrate near the molecular surface. Returns a
      :class:`SurfaceFollowingSampler`. Higher rejection rate; useful
      when sample-density-per-surface-area matters more than wall-clock
      speed.

    Parameters
    ----------
    complex_obj
        Built ligand-only ``Complex`` (positions + topology).
    mode : {"sphere", "surface-following"}, default "sphere"
        Which sampler to construct.
    reach : float, default 10.0
        ``mode="sphere"`` only. How far the bounding sphere extends
        past the ligand's bounding radius (Å). Must be ``>= 0``.
    d_max : float, default 6.0
        ``mode="surface-following"`` only. Band thickness from the
        molecular surface (Å). Must be ``> 0``.
    probe : float, default 1.4
        Probe radius for the SAS rejection (Å), water-equivalent by
        default. Must be ``>= 0``. Used by both modes.

    Returns
    -------
    SurfaceSampler or SurfaceFollowingSampler
        Concrete type depends on ``mode``. Both expose
        ``.generator() -> Sample``.

    Examples
    --------
    >>> # default (sphere) — most callers want this:
    >>> sampler = make_sampler(ligand_only)  # doctest: +SKIP
    >>> pose = sampler.generator()  # doctest: +SKIP

    >>> # opt in to surface-following with a 4 Å band:
    >>> sampler = make_sampler(  # doctest: +SKIP
    ...     ligand_only, mode="surface-following", d_max=4.0
    ... )
    """
    if probe < 0:
        raise ValueError(f"probe must be >= 0, got {probe}")
    if mode == "sphere":
        if reach < 0:
            raise ValueError(f"reach must be >= 0, got {reach}")
        dims = compute_envelope_dims(complex_obj, reach)
        envelope = Sphere(**dims)
        excluder = Excluder(complex_obj, probe=probe)
        return SurfaceSampler(envelope=envelope, excluder=excluder)
    if mode == "surface-following":
        return SurfaceFollowingSampler(complex_obj, d_max=d_max, probe=probe)
    raise ValueError(
        f"Unknown mode {mode!r}; expected one of 'sphere', 'surface-following'."
    )

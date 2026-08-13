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
ClashFilter
    Rejects a placed strand whose atoms overlap the docking target.
SurfaceSampler
    Composes an envelope + an Excluder into a rejection sampler.
SurfaceFollowingSampler
    Alternative sampler that concentrates poses near the molecular
    surface via an outer distance cap.
make_sampler(complex_obj, *, mode="sphere", reach=10.0, d_max=6.0, probe=1.4, rng=None)
    Factory: builds the requested sampler for the given Complex.
draw_clear_conformation(complex_obj, chain, sampler, rotations, clash)
    Draws poses and torsions until one sits clear of the target.
compute_envelope_dims(complex_obj, reach)
    Returns the kwargs for the Sphere dataclass; useful in tests and
    notebooks.

Reproducibility
---------------
Every random draw here comes from a :class:`numpy.random.Generator` given as
``rng``. Pass a seed to repeat a run. Leave it out and each sampler builds
its own generator, so runs differ.


Examples
--------
>>> from maws.space import NAngles
>>> sample = NAngles(n=4).generator()
>>> len(sample)
4
"""

import warnings
from dataclasses import InitVar, dataclass, field
from typing import Literal, Protocol

import numpy as np
from openmm import unit
from scipy.spatial import KDTree

from maws.helpers import atom_masses, mass_weighted_center, nostrom

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


def _random_unit_axis(rng: np.random.Generator) -> np.ndarray:
    """Random unit vector uniformly distributed on the unit sphere.

    Uses the standard Gaussian-then-normalize recipe: the multivariate
    standard normal has spherically symmetric density, so normalizing
    yields a point uniform on the unit sphere. A naive uniform draw in
    the cube ``[-1, 1]^3`` followed by normalization would
    over-represent the cube's corner-aligned directions because the
    cube has greater radial extent along its diagonals than its faces.
    """
    v = rng.standard_normal(3)
    return v / np.linalg.norm(v)


class Envelope(Protocol):
    """Anything with ``.generator() -> Sample`` is an envelope.

    :class:`Sphere` is the only built-in. Users may pass a custom envelope
    directly to :class:`SurfaceSampler` if they want a different
    geometric region.
    """

    def generator(self) -> Sample: ...


def _spherical_sample(
    centre: np.ndarray, r: float, rng: np.random.Generator
) -> np.ndarray:
    """Random point on a sphere of radius ``r`` around ``centre``
    (uniform direction)."""
    phi = rng.uniform(0, 2 * np.pi)
    cos_psi = rng.uniform(-1, 1)  # uniform on the sphere
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

    Parameters
    ----------
    radius : float
        How far the sphere reaches from `centre`, in ångström.
    centre : numpy.ndarray
        Shape ``(3,)`` point the sphere is built around, in ångström.
    rng : int or numpy.random.Generator, optional
        Source of randomness. Pass a seed to make a run repeatable.
        Defaults to a fresh generator, so runs differ.
    """

    radius: float
    centre: np.ndarray
    rng: InitVar[int | np.random.Generator | None] = None
    _rng: np.random.Generator = field(init=False, repr=False, compare=False)

    def __post_init__(self, rng):
        object.__setattr__(self, "_rng", np.random.default_rng(rng))

    def generator(self) -> Sample:
        # Volume-correct: r = R · u^(1/3) for u ~ U(0,1).
        r = self.radius * self._rng.uniform(0, 1) ** (1 / 3)
        return Sample(
            position=_spherical_sample(self.centre, r, self._rng),
            axis=_random_unit_axis(self._rng),
            angle=float(self._rng.uniform(0, 2 * np.pi)),
        )


@dataclass(frozen=True)
class NAngles:
    """N independent angles drawn from [0, 2π). Used for in-residue rotations.

    Parameters
    ----------
    n : int
        How many angles each draw returns, one per rotatable bond.
    rng : int or numpy.random.Generator, optional
        Source of randomness. Pass a seed to make a run repeatable.
        Defaults to a fresh generator, so runs differ.
    """

    n: int
    rng: InitVar[int | np.random.Generator | None] = None
    _rng: np.random.Generator = field(init=False, repr=False, compare=False)

    def __post_init__(self, rng):
        object.__setattr__(self, "_rng", np.random.default_rng(rng))

    def generator(self) -> np.ndarray:
        return self._rng.uniform(0, 2 * np.pi, self.n)


_WARNED_ELEMENTS: set[str] = set()
"""Element symbols already reported as missing from :data:`_BONDI_VDW_RADII`.

Process-wide, so a topology with thousands of atoms of an unlisted element
produces one warning rather than thousands.
"""


def _vdw_radii(topology) -> np.ndarray:
    """Return the van der Waals radius of every atom, in ångström.

    Parameters
    ----------
    topology : openmm.app.Topology
        Any object whose ``atoms()`` yields atoms carrying ``.element.symbol``.

    Returns
    -------
    numpy.ndarray
        Shape ``(N,)`` radii in ångström, ordered to match a positions array.

    Warns
    -----
    RuntimeWarning
        Once per element symbol missing from the Bondi table, which then falls
        back to :data:`_DEFAULT_VDW`.
    """
    radii = []
    for atom in topology.atoms():
        symbol = atom.element.symbol
        radius = _BONDI_VDW_RADII.get(symbol)
        if radius is None:
            radius = _DEFAULT_VDW
            if symbol not in _WARNED_ELEMENTS:
                _WARNED_ELEMENTS.add(symbol)
                warnings.warn(
                    f"Unknown element {symbol!r} - using fallback "
                    f"vdW = {_DEFAULT_VDW} Å",
                    RuntimeWarning,
                    stacklevel=2,
                )
        radii.append(radius)
    return np.array(radii, dtype=float)


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

    def __init__(self, complex_obj, probe: float = 1.4):
        positions = np.asarray(nostrom(complex_obj.positions), dtype=float)
        radii = _vdw_radii(complex_obj.topology)
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


class ClashFilter:
    """
    Steric test on the atoms a pose actually puts down.

    :class:`Excluder` answers a question about a *point*: could a water-sized
    probe sit here? What a pose puts at that point is a whole strand — tens of
    atoms spanning several ångström — so a point can pass that test while the
    atoms placed there sit inside the target. This class tests those atoms.

    Decision rule
    -------------
    The complex is split in two by `element`: the atoms of the strand being
    placed, which move, and everything else, which is the docking target and
    does not. A pose is rejected when any moving atom *i* and any rigid atom
    *j* are closer than their van der Waals spheres allow:

        ``dist(i, j) < vdW(i) + vdW(j) - tolerance``

    Only strand-target pairs are tested. A torsion swings a rigid group about
    a bond, so the strand keeps its own internal geometry, and the force field
    already resists what little self-overlap a torsion can create.

    Radii come from the Bondi (1964) table at the top of this module.

    Parameters
    ----------
    complex_obj
        Anything with ``.positions`` (openmm Quantity, Å-convertible, shape
        ``(N, 3)``) and ``.topology.atoms()``. The rigid atoms' coordinates
        are read once, here, so the target must not move afterwards.
    element : sequence of int
        ``[start, bond, end]`` atom indices of the strand being placed, with
        ``end`` exclusive. ``bond`` is ignored.
    tolerance : float, default=1.0
        How far two atoms may overlap their van der Waals spheres before the
        pose is rejected, in ångström. Some overlap has to be allowed: a
        hydrogen bond puts a hydrogen roughly 0.8 Å inside the summed radii,
        so a tolerance below that rejects poses that are bound, not clashing.
        Raising it further admits harder contacts.

    Raises
    ------
    ValueError
        If `tolerance` is negative.

    See Also
    --------
    Excluder : The point test this complements.
    maws.complex.Complex.place_global : Puts the strand down for this to judge.

    Examples
    --------
    Reject a pose, draw another, and keep the first that clears the target.

    >>> clash = ClashFilter(complex, chain.element)  # doctest: +SKIP
    >>> pose = sampler.generator()  # doctest: +SKIP
    >>> clash.is_clear(nostrom(complex.positions))  # doctest: +SKIP
    False
    """

    def __init__(self, complex_obj, element, *, tolerance: float = 1.0):
        if tolerance < 0:
            raise ValueError(f"tolerance must be >= 0, got {tolerance}")

        radii = _vdw_radii(complex_obj.topology)
        positions = np.asarray(nostrom(complex_obj.positions), dtype=float)
        self._start, self._end = element[0], element[2]

        rigid = np.ones(len(positions), dtype=bool)
        rigid[self._start : self._end] = False
        self._rigid_positions = positions[rigid]
        self._rigid_radii = radii[rigid]
        self._moving_radii = radii[self._start : self._end]
        self._tolerance = tolerance
        self._tree = KDTree(self._rigid_positions)
        self._reach = float(
            self._rigid_radii.max() + self._moving_radii.max() - tolerance
        )

    def is_clear(self, positions: np.ndarray) -> bool:
        """Return True when no atom of the strand overlaps the target.

        Parameters
        ----------
        positions : numpy.ndarray
            Shape ``(N, 3)`` coordinates of the whole complex in ångström,
            unitless. Only the rows `element` covers are read.

        Returns
        -------
        bool
            True when the pose is worth costing an energy evaluation.
        """
        moving = np.asarray(positions, dtype=float)[self._start : self._end]
        for k, near in enumerate(self._tree.query_ball_point(moving, self._reach)):
            if not near:
                continue
            gap = np.linalg.norm(self._rigid_positions[near] - moving[k], axis=1)
            limit = self._rigid_radii[near] + self._moving_radii[k] - self._tolerance
            if (gap < limit).any():
                return False
        return True


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
    com = mass_weighted_center(pos, atom_masses(complex_obj.topology))
    dists = np.linalg.norm(pos - com, axis=1)
    return {"radius": float(dists.max()) + reach, "centre": com}


class SamplingError(RuntimeError):
    """Raised when SurfaceSampler cannot find a clear point in max_rejections tries."""


def draw_clear_conformation(
    complex_obj,
    chain,
    sampler,
    rotations,
    clash: ClashFilter,
    *,
    residue: int = 0,
    max_rejections: int = 1000,
) -> Sample:
    """draw_clear_conformation(complex_obj, chain, sampler, rotations, clash, *,
    residue=0, max_rejections=1000)

    Draw whole conformations until one sits clear of the target.

    One attempt puts the strand at a drawn pose, bends one residue by a drawn
    set of torsion angles, and asks `clash` whether the result touches the
    target. A rejected attempt is neither scored nor counted, and the strand
    is returned to where it started before the next one, so only a single
    draw is accepted.

    Both steps have to happen before the test. Placing a strand clear of the
    target does not keep it clear: the torsions that follow swing atoms
    several ångström, far enough to reach back into it.

    Parameters
    ----------
    complex_obj : maws.complex.Complex
        Holds both the strand and the target. Left holding the accepted
        conformation.
    chain : maws.chain.Chain
        The strand being placed and bent.
    sampler
        Draws the pose. Any object with ``.generator() -> Sample``
        (built-in: :class:`SurfaceSampler`).
    rotations : NAngles
        Draws one torsion angle per rotatable backbone bond, in radians.
    clash : ClashFilter
        Built for the same `complex_obj` and ``chain.element``.
    residue : int, default=0
        Which residue of `chain` the torsions are applied to. The default is
        the only residue of a one-nucleotide strand.
    max_rejections : int, default=1000
        Hard cap on consecutive rejected attempts before giving up.

    Returns
    -------
    Sample
        The pose that was accepted. `complex_obj` already holds it.

    Raises
    ------
    SamplingError
        If nothing clears the target within `max_rejections` attempts.

    See Also
    --------
    ClashFilter : Makes the accept/reject decision.
    maws.complex.Complex.place_global : Puts the strand down.

    Examples
    --------
    >>> clash = ClashFilter(complex, chain.element)  # doctest: +SKIP
    >>> pose = draw_clear_conformation(  # doctest: +SKIP
    ...     complex, chain, sampler, NAngles(4), clash
    ... )
    """
    start = complex_obj.positions[:]
    for _ in range(max_rejections):
        complex_obj.positions = start[:]

        pose = sampler.generator()
        complex_obj.place_global(
            chain.element,
            pose.position * unit.angstrom,
            pose.axis * unit.angstrom,
            pose.angle,
        )
        for torsion, angle in enumerate(rotations.generator()):
            chain.rotate_in_residue(residue, torsion, angle)

        if clash.is_clear(nostrom(complex_obj.positions)):
            return pose

    complex_obj.positions = start[:]
    raise SamplingError(
        f"Could not fit the strand against the target in {max_rejections} "
        f"attempts. The sampling region may sit inside the target - increase "
        f"--reach, or raise the clash tolerance."
    )


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
            f"The region may lie inside the target - widen it with --reach, or "
            f"with --site-radius if a site centre was given, decrease --probe, "
            f"or check the ligand size."
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
    rng : int or numpy.random.Generator, optional
        Source of randomness. Pass a seed to make a run repeatable.
        Defaults to a fresh generator, so runs differ.
    """

    def __init__(
        self,
        complex_obj,
        *,
        d_max: float = 6.0,
        probe: float = 1.4,
        max_rejections: int = 50_000,
        rng: int | np.random.Generator | None = None,
    ):
        if d_max <= 0:
            raise ValueError(f"d_max must be > 0, got {d_max}")
        self._rng = np.random.default_rng(rng)
        positions = np.asarray(nostrom(complex_obj.positions), dtype=float)
        com = mass_weighted_center(positions, atom_masses(complex_obj.topology))
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
            r = self._R_bound * self._rng.uniform(0, 1) ** (1 / 3)
            phi = self._rng.uniform(0, 2 * np.pi)
            cos_psi = self._rng.uniform(-1, 1)
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
                axis=_random_unit_axis(self._rng),
                angle=float(self._rng.uniform(0, 2 * np.pi)),
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
    site_centre=None,
    site_radius: float | None = None,
    rng: int | np.random.Generator | None = None,
):
    """
    Build a fully-configured surface-aware sampler for ``complex_obj``.

    Both modes draw a point uniformly through the volume of a bounding
    sphere and keep it only if it lies outside the target's
    solvent-accessible surface, so both return points that sit in solvent.
    They differ in how far out into that solvent a point may sit.

    - ``mode="sphere"`` (default): the sphere is centred on the target's
      centre of mass with a radius of ``furthest atom + reach``, and every
      point outside the surface is kept. On a protein most of that volume
      is open solvent well clear of the target. Measured on a 2760-atom
      target with the shipped defaults: the sphere has a radius of 39 Å,
      half the accepted points sit more than 12 Å from the nearest target
      atom, and the furthest sits 25 Å out. A strand placed there floats
      free of the target. Returns a :class:`SurfaceSampler`.

    - ``mode="surface-following"``: a point is kept only while some target
      atom lies within ``d_max`` of it, so accepted points form a band over
      the surface that dips into pockets and wraps around protrusions. On
      the same target with ``d_max=6``, half the accepted points sit within
      4.4 Å of an atom and none beyond 6 Å. It spends more draws per
      accepted point. Returns a :class:`SurfaceFollowingSampler`.

    The two limits measure different things, which is what makes the second
    mode follow the shape: ``reach`` extends the target's overall bounding
    radius, one number for the whole target, while ``d_max`` is measured
    from whichever atom happens to be nearest the point.

    Sphere mode can be aimed at one site with `site_centre`, which replaces
    the auto-sized sphere with one of `site_radius` around a chosen point.

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
    site_centre : array-like, optional
        ``mode="sphere"`` only. Shape ``(3,)`` point in ångström to sample
        around, in the target's coordinate frame. Give it when the binding
        site is known: the whole sample budget is then spent on that site
        rather than spread over the target's entire surface. Omit it and the
        region is sized to hold the whole target.
    site_radius : float, default 15.0
        How far the region reaches from `site_centre`, in ångström. Applies
        only alongside `site_centre`. Must be ``> 0``.
    rng : int or numpy.random.Generator, optional
        Source of randomness. Pass a seed to make a run repeatable.
        Defaults to a fresh generator, so runs differ.

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
    if site_centre is None and site_radius is not None:
        raise ValueError("site_radius applies only alongside site_centre.")
    if site_centre is not None:
        if mode != "sphere":
            raise ValueError(
                f"site_centre applies only to mode='sphere', got mode={mode!r}. "
                f"Surface-following covers the whole surface."
            )
        site_centre = np.asarray(site_centre, dtype=float)
        if site_centre.shape != (3,):
            raise ValueError(
                f"site_centre must be a point of three coordinates, got shape "
                f"{site_centre.shape}"
            )
        site_radius = 15.0 if site_radius is None else float(site_radius)
        if site_radius <= 0:
            raise ValueError(f"site_radius must be > 0, got {site_radius}")

    if mode == "sphere":
        if reach < 0:
            raise ValueError(f"reach must be >= 0, got {reach}")
        if site_centre is None:
            dims = compute_envelope_dims(complex_obj, reach)
        else:
            dims = {"radius": site_radius, "centre": site_centre}
        envelope = Sphere(**dims, rng=rng)
        excluder = Excluder(complex_obj, probe=probe)
        return SurfaceSampler(envelope=envelope, excluder=excluder)
    if mode == "surface-following":
        return SurfaceFollowingSampler(complex_obj, d_max=d_max, probe=probe, rng=rng)
    raise ValueError(
        f"Unknown mode {mode!r}; expected one of 'sphere', 'surface-following'."
    )

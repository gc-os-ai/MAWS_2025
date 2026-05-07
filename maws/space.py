"""
maws.space
==========

Surface-aware sampling for the MAWS aptamer-design loop.

The intended public entry point is :func:`make_sampler`. It picks the
envelope shape, auto-sizes it from the ligand atoms (mass-weighted COM,
bounding radii), and wraps the result in a :class:`SurfaceSampler` whose
``.generator()`` rejects candidate poses inside the protein bulk via an
SAS-style Bondi-vdW + probe check.

Public API
----------
Cube, Sphere, Shell : frozen dataclasses
    Geometric envelopes; each ``.generator()`` returns a 7-element ndarray
    ``[x, y, z, axis_x, axis_y, axis_z, rotation_angle]``.
NAngles : frozen dataclass
    N independent angles in [0, 2π) for in-residue rotations.
Excluder
    KDTree-backed SAS rejection: a candidate is "clear" iff its distance
    to every ligand atom exceeds (Bondi vdW + probe).
SurfaceSampler
    Composes an envelope + an Excluder into a rejection sampler.
make_sampler(shape, complex_obj, *, reach=10.0, probe=1.4)
    Factory: builds envelope + Excluder for the given Complex.
compute_envelope_dims(complex_obj, shape, reach)
    Returns the kwargs for the matching envelope dataclass; useful
    in tests and notebooks.

Examples
--------
>>> from maws.space import NAngles
>>> sample = NAngles(n=4).generator()
>>> len(sample)
4
"""

import warnings
from dataclasses import dataclass
from typing import Protocol

import numpy as np
from openmm import unit
from scipy.spatial import KDTree

from maws.helpers import mass_weighted_center, nostrom

# ============================================================================
# Surface-aware sampling — public API
# ============================================================================
#
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
    v = np.random.uniform(-1, 1, 3)
    return v / np.linalg.norm(v)


class Envelope(Protocol):
    """Anything with `.generator() -> Sample` is an envelope."""

    def generator(self) -> Sample: ...


@dataclass(frozen=True)
class Cube:
    """Axis-aligned cube of side `width`, centered on `centre` (Å)."""

    width: float
    centre: np.ndarray

    def generator(self) -> Sample:
        half = self.width / 2.0
        position = np.array(
            [
                np.random.uniform(self.centre[0] - half, self.centre[0] + half),
                np.random.uniform(self.centre[1] - half, self.centre[1] + half),
                np.random.uniform(self.centre[2] - half, self.centre[2] + half),
            ]
        )
        return Sample(
            position=position,
            axis=_random_unit_axis(),
            angle=float(np.random.uniform(0, np.pi)),
        )


def _spherical_sample(centre: np.ndarray, r: float) -> np.ndarray:
    """Random point on a sphere of radius `r` around `centre` (uniform direction)."""
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
    """Sphere of given radius centered on `centre` (Å). Volume-correct sampling."""

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
class Shell:
    """Spherical shell between inner and outer radii, centered on `centre` (Å)."""

    inner: float
    outer: float
    centre: np.ndarray

    def generator(self) -> Sample:
        # Volume-correct shell: r = (u·(R_out^3 - R_in^3) + R_in^3)^(1/3).
        u = np.random.uniform(0, 1)
        r = (u * (self.outer**3 - self.inner**3) + self.inner**3) ** (1 / 3)
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
    Rejects candidate points lying inside (vdW + probe) of any ligand atom.

    Parameters
    ----------
    complex_obj
        Anything with `.positions` (openmm Quantity, Å-convertible, shape (N,3))
        and `.topology.atoms` (iterable yielding objects with `.element.symbol`).
    probe : float
        SAS probe radius in Å. Default 1.4 (water-like).
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
        """True iff `point` lies outside every (vdW + probe) sphere."""
        idx = self._tree.query_ball_point(point, self._max_inflated)
        if not idx:
            return True
        diffs = self._positions[idx] - np.asarray(point)
        dists2 = (diffs**2).sum(axis=1)
        return bool((dists2 > self._inflated[idx] ** 2).all())


# Per-shape formula for auto-sizing the envelope from the ligand bounding
# radii. Each entry maps shape name → (R_max, R_min, com, reach) → kwargs
# dict for the matching envelope dataclass. Single source of truth for
# both shape names AND sizing rules.
_DIM_FORMULAS: dict = {
    "cube": lambda R_max, R_min, com, reach: {
        "width": 2.0 * (R_max + reach),
        "centre": com,
    },
    "sphere": lambda R_max, R_min, com, reach: {
        "radius": R_max + reach,
        "centre": com,
    },
    "shell": lambda R_max, R_min, com, reach: {
        "inner": max(0.0, R_min - 5.0),
        "outer": R_max + reach,
        "centre": com,
    },
}


def _atom_mass_in_dalton(atom) -> float:
    m = atom.element.mass
    if hasattr(m, "value_in_unit"):
        return m.value_in_unit(unit.dalton)
    return float(m)


def compute_envelope_dims(complex_obj, shape: str, reach: float) -> dict:
    """
    Compute auto-sized envelope dimensions from the ligand atoms.

    Parameters
    ----------
    complex_obj
        Object with `.positions` (Quantity) and `.topology.atoms()` (yielding
        atoms with `.element.mass`).
    shape : {"cube", "sphere", "shell"}
        Envelope shape.
    reach : float
        Aptamer reach beyond the surface (Å).

    Returns
    -------
    dict
        Kwargs for the matching envelope dataclass.
    """
    formula = _DIM_FORMULAS.get(shape)
    if formula is None:
        raise ValueError(
            f"Unknown shape {shape!r}; expected one of {list(_DIM_FORMULAS)}."
        )
    pos = np.asarray(nostrom(complex_obj.positions), dtype=float)
    masses = np.array(
        [_atom_mass_in_dalton(a) for a in complex_obj.topology.atoms()],
        dtype=float,
    )
    com = mass_weighted_center(pos, masses)
    dists = np.linalg.norm(pos - com, axis=1)
    return formula(float(dists.max()), float(dists.min()), com, reach)


class SamplingError(RuntimeError):
    """Raised when SurfaceSampler cannot find a clear point in max_rejections tries."""


@dataclass
class SurfaceSampler:
    """
    Rejection sampler: draws from `envelope` until `excluder.is_clear` accepts.

    Parameters
    ----------
    envelope
        Any object with `.generator() -> 7-element ndarray` (Cube / Sphere / Shell).
    excluder
        Excluder instance.
    max_rejections : int
        Safety cap to fail fast on misconfigured envelopes.
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


# Single source of truth for supported shapes — used by both make_sampler
# (envelope construction) and compute_envelope_dims (validation).
_ENVELOPE_TYPES: dict[str, type] = {
    "cube": Cube,
    "sphere": Sphere,
    "shell": Shell,
}


def make_sampler(
    shape: str,
    complex_obj,
    *,
    reach: float = 10.0,
    probe: float = 1.4,
) -> SurfaceSampler:
    """
    Build a fully-configured surface-aware sampler for `complex_obj`.

    Parameters
    ----------
    shape : {"cube", "sphere", "shell"}
        Envelope shape (auto-sized from `complex_obj`).
    complex_obj
        Built ligand-only Complex (positions + topology).
    reach : float, default 10.0
        Distance the envelope extends beyond the ligand surface (Å).
    probe : float, default 1.4
        vdW probe radius for SAS rejection (Å). Water-like at the default.

    Returns
    -------
    SurfaceSampler
    """
    cls = _ENVELOPE_TYPES.get(shape)
    if cls is None:
        raise ValueError(
            f"Unknown shape {shape!r}; expected one of {list(_ENVELOPE_TYPES)}."
        )
    dims = compute_envelope_dims(complex_obj, shape, reach)
    envelope = cls(**dims)
    excluder = Excluder(complex_obj, probe=probe)
    return SurfaceSampler(envelope=envelope, excluder=excluder)

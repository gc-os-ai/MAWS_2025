"""
maws.space
==========

Sampling spaces for molecular conformational exploration.

This module defines geometric spaces from which MAWS samples configurations
(positions + orientations) during aptamer-ligand binding evaluation.

Classes
-------
Space : Base class for sampling spaces.
Box : Rectangular 3D sampling volume.
Cube : Cubic 3D sampling volume (special case of Box).
Sphere : Spherical 3D sampling volume.
SphericalShell : Hollow sphere sampling volume.
NAngles : N-dimensional torsion angle space.

Examples
--------
>>> import numpy as np
>>> from maws.space import NAngles

Create a 4-angle torsion space:

>>> angles = NAngles(4)
>>> sample = angles.generator()
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


class Envelope(Protocol):
    """Anything with `.generator() -> 7-element ndarray` is an envelope."""

    def generator(self) -> np.ndarray: ...


@dataclass(frozen=True)
class Cube:
    """Axis-aligned cube of side `width`, centered on `centre` (Å)."""

    width: float
    centre: np.ndarray

    def generator(self) -> np.ndarray:
        axis = np.random.uniform(-1, 1, 3)
        ax, ay, az = axis / np.linalg.norm(axis)
        rotation = np.random.uniform(0, np.pi)
        half = self.width / 2.0
        return np.array(
            [
                np.random.uniform(self.centre[0] - half, self.centre[0] + half),
                np.random.uniform(self.centre[1] - half, self.centre[1] + half),
                np.random.uniform(self.centre[2] - half, self.centre[2] + half),
                ax,
                ay,
                az,
                rotation,
            ]
        )


@dataclass(frozen=True)
class Sphere:
    """Sphere of given radius centered on `centre` (Å). Volume-correct sampling."""

    radius: float
    centre: np.ndarray

    def generator(self) -> np.ndarray:
        axis = np.random.uniform(-1, 1, 3)
        ax, ay, az = axis / np.linalg.norm(axis)
        # Volume-correct: r = R · u^(1/3) for u ~ U(0,1).
        u = np.random.uniform(0, 1)
        r = self.radius * u ** (1 / 3)
        phi = np.random.uniform(0, 2 * np.pi)
        # Uniform on sphere: cos(psi) ~ U(-1, 1).
        cos_psi = np.random.uniform(-1, 1)
        psi = np.arccos(cos_psi)
        rotation = np.random.uniform(0, 2 * np.pi)
        return np.array(
            [
                self.centre[0] + r * np.cos(phi) * np.sin(psi),
                self.centre[1] + r * np.sin(phi) * np.sin(psi),
                self.centre[2] + r * np.cos(psi),
                ax,
                ay,
                az,
                rotation,
            ]
        )


@dataclass(frozen=True)
class Shell:
    """Spherical shell between inner and outer radii, centered on `centre` (Å)."""

    inner: float
    outer: float
    centre: np.ndarray

    def generator(self) -> np.ndarray:
        axis = np.random.uniform(-1, 1, 3)
        ax, ay, az = axis / np.linalg.norm(axis)
        # Volume-correct shell: r = (u·(R_out^3 - R_in^3) + R_in^3)^(1/3).
        u = np.random.uniform(0, 1)
        r = (u * (self.outer**3 - self.inner**3) + self.inner**3) ** (1 / 3)
        phi = np.random.uniform(0, 2 * np.pi)
        cos_psi = np.random.uniform(-1, 1)
        psi = np.arccos(cos_psi)
        rotation = np.random.uniform(0, 2 * np.pi)
        return np.array(
            [
                self.centre[0] + r * np.cos(phi) * np.sin(psi),
                self.centre[1] + r * np.sin(phi) * np.sin(psi),
                self.centre[2] + r * np.cos(psi),
                ax,
                ay,
                az,
                rotation,
            ]
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
        for i, atom in enumerate(complex_obj.topology.atoms):
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


def compute_envelope_dims(complex_obj, shape: str, reach: float) -> dict:
    """
    Compute auto-sized envelope dimensions from the ligand atoms.

    Parameters
    ----------
    complex_obj
        Object with `.positions` (Quantity) and `.topology.atoms` (yielding
        atoms with `.element.mass`).
    shape : {"cube", "sphere", "shell"}
        Envelope shape.
    reach : float
        Aptamer reach beyond the surface (Å).

    Returns
    -------
    dict
        Kwargs suitable for the matching envelope dataclass:
          cube   → {"width", "centre"}
          sphere → {"radius", "centre"}
          shell  → {"inner", "outer", "centre"}
    """
    pos = np.asarray(nostrom(complex_obj.positions), dtype=float)

    def _mass(atom):
        m = atom.element.mass
        if hasattr(m, "value_in_unit"):
            return m.value_in_unit(unit.dalton)
        return float(m)

    masses = np.array([_mass(a) for a in complex_obj.topology.atoms], dtype=float)
    com = mass_weighted_center(pos, masses)
    dists = np.linalg.norm(pos - com, axis=1)
    R_max = float(dists.max())
    R_min = float(dists.min())

    if shape == "cube":
        return {"width": 2.0 * (R_max + reach), "centre": com}
    if shape == "sphere":
        return {"radius": R_max + reach, "centre": com}
    if shape == "shell":
        return {
            "inner": max(0.0, R_min - 5.0),
            "outer": R_max + reach,
            "centre": com,
        }
    raise ValueError(f"Unknown shape {shape!r}; expected one of cube / sphere / shell.")


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

    def generator(self) -> np.ndarray:
        for _ in range(self.max_rejections):
            sample = self.envelope.generator()
            if self.excluder.is_clear(sample[:3]):
                return sample
        raise SamplingError(
            f"Could not draw a clear point in {self.max_rejections} attempts. "
            f"Envelope may be fully buried — increase --reach, decrease --probe, "
            f"or check ligand size."
        )

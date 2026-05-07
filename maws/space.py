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

from dataclasses import dataclass
from typing import Protocol

import numpy as np

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

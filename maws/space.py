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
from openmm import unit

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


# ============================================================================
# Legacy Space + NAngles (still used by run.py / maws2023.py until migration).
# ============================================================================


class Space:
    """
    Base class representing a sampling space.

    A Space defines a region from which random samples can be drawn.
    Subclasses implement specific geometries (Box, Sphere, etc.).

    Parameters
    ----------
    membership_boolean_function : callable
        Function that returns True if a point is inside the space.
    units : openmm.unit, default=unit.angstroms
        Length units for the space.
    lower_bound : float, default=-1e6
        Lower bound for rejection sampling.
    upper_bound : float, default=1e6
        Upper bound for rejection sampling.

    Attributes
    ----------
    is_in : callable
        Membership test function.
    volume : float
        Volume of the space (for normalization).
    """

    def __init__(
        self,
        membership_boolean_function,
        units=unit.angstroms,
        lower_bound=-1e6,
        upper_bound=1e6,
    ):
        self.is_in = membership_boolean_function
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.units = units
        self.volume = 0

    def generator(self):
        """
        Generate a random sample from the space.

        Returns
        -------
        numpy.ndarray
            3D coordinates of a point inside the space.
        """
        result = None
        while result is None:
            candidate = np.random.uniform(self.lower_bound, self.upper_bound, 3)
            if self.is_in(candidate):
                result = candidate
        return result


class NAngles(Space):
    """
    N-dimensional torsion angle sampling space.

    Represents U(1)^N - N independent angles in [0, 2π). Used for sampling
    backbone torsion angles during conformational search.

    Parameters
    ----------
    number : int
        Number of angles to generate (N).

    Attributes
    ----------
    number : int
        Number of angles in the space.

    Examples
    --------
    >>> angles = NAngles(4)
    >>> sample = angles.generator()
    >>> len(sample)
    4
    >>> all(0 <= a < 2 * 3.14159 for a in sample)  # All in [0, 2π)
    True
    """

    def __init__(self, number):
        def boolean(position):
            return all(0 <= element and element <= 2 * np.pi for element in position)

        super().__init__(boolean)
        self.number = number

    def generator(self):
        """
        Generate N random angles in [0, 2π).

        Returns
        -------
        numpy.ndarray
            Array of `self.number` random angles in radians.
        """
        return np.array([np.random.uniform(0, 2 * np.pi) for i in range(self.number)])

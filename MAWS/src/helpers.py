"""
helpers.py – small unit-conversion & geometry helpers for MAWS.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import NDArray
from openmm import unit

__all__ = [
    # modern snake_case
    "to_angstrom",
    "from_angstrom",
    "to_kj",
    "from_kj",
    "angle",
    "directed_angle",
    # legacy camelCase aliases
    "angstrom",
    "nostrom",
    "kJ",
    "noJ",
]

# ---------------------------------------------------------------------------- #
# Unit helpers
# ---------------------------------------------------------------------------- #
def to_angstrom(vec: NDArray[np.floating] | Sequence[float]) -> unit.Quantity:
    """Return *vec* with OpenMM Å units attached."""
    return unit.Quantity(np.asarray(vec, dtype=float), unit.angstrom)


def from_angstrom(quantity: unit.Quantity) -> NDArray[np.floating]:
    """Strip OpenMM units, return **plain NumPy array** in Å."""
    return np.asarray(quantity.value_in_unit(unit.angstrom), dtype=float)


def to_kj(value: float | NDArray[np.floating]) -> unit.Quantity:
    """Attach kJ·mol⁻¹ units."""
    return unit.Quantity(np.asarray(value, dtype=float), unit.kilojoules_per_mole)


def from_kj(quantity: unit.Quantity) -> float | NDArray[np.floating]:
    """Return numerical value(s) in kJ·mol⁻¹ without units."""
    return quantity.value_in_unit(unit.kilojoules_per_mole)


# ---------------------------------------------------------------------------- #
# Geometry helpers
# ---------------------------------------------------------------------------- #
def _normalise(v: NDArray[np.floating]) -> NDArray[np.floating]:
    """Return the unit vector of *v*."""
    norm = np.linalg.norm(v)
    if norm == 0.0:
        raise ValueError("Zero-length vector cannot be normalised.")
    return v / norm


def angle(v1: Sequence[float], v2: Sequence[float]) -> float:
    """
    Smallest angle between two vectors in **radians** (0–π).

    Equivalent to ``arccos(dot(u,v))`` with numerical clipping.
    """
    u = _normalise(np.asarray(v1, dtype=float))
    v = _normalise(np.asarray(v2, dtype=float))
    cos_th = np.clip(np.dot(u, v), -1.0, 1.0)
    return float(np.arccos(cos_th))


def directed_angle(
    v1: Sequence[float],
    v2: Sequence[float],
    axis: Sequence[float],
) -> float:
    """
    Signed angle from *v1* to *v2* around *axis* (right-hand rule, radians).

    Range is (–π, π].
    """
    u1 = _normalise(np.asarray(v1, dtype=float))
    u2 = _normalise(np.asarray(v2, dtype=float))
    n = _normalise(np.asarray(axis, dtype=float))
    sin_th = np.dot(n, np.cross(u1, u2))
    cos_th = np.dot(u1, u2)
    return float(np.arctan2(sin_th, cos_th))


# ---------------------------------------------------------------------------- #
# Backward-compat CamelCase aliases
# ---------------------------------------------------------------------------- #
angstrom = to_angstrom
nostrom = from_angstrom
kJ = to_kj
noJ = from_kj

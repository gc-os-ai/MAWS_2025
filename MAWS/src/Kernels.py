"""
Kernels.py – high-performance numerical kernels compiled with Numba.

All public helpers return **plain NumPy arrays or floats** so that the hot
loops are numba-compilable.  Wrap the outputs with `openmm.unit.Quantity`
in user code if unit semantics are required.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np
from numba import njit
from openmm import unit  # only for callers who need to attach units

# --------------------------------------------------------------------------- #
# Constants / small helpers
# --------------------------------------------------------------------------- #
_EPS = 1.0e-50


@njit(cache=True, fastmath=True)
def _safe(x: float) -> float:
    """Prevent log(0) inside njit-compiled functions."""
    return x + _EPS


# --------------------------------------------------------------------------- #
# Kernels
# --------------------------------------------------------------------------- #
@njit(cache=True, fastmath=True)
def rotate_kernel(
    positions: np.ndarray,
    element: Tuple[int, int, int],  # (start, bond, end)
    axis: np.ndarray,
    angle: float,
) -> np.ndarray:
    """
    Rotate atoms *element[1] … element[2]-1* around *axis* by *angle* (rad).

    Parameters
    ----------
    positions
        (N, 3) Cartesian coordinates.
    element
        Tuple *(res_index, first_atom, last_atom_exclusive)*.
    axis
        3-vector, need not be normalised.
    angle
        Rotation angle in **radians**.

    Returns
    -------
    numpy.ndarray
        Rotated copy of *positions* (same shape).
    """
    # normalise axis ----------------------------------------------------- #
    axis_unit = axis / np.linalg.norm(axis)
    x, y, z = axis_unit
    phi_2 = angle * 0.5
    s = math.sin(phi_2)
    c = math.cos(phi_2)

    # quaternion → 3×3 rotation matrix ---------------------------------- #
    rot = np.empty((3, 3))
    rot[0, 0] = 2 * (x * x - 1) * s * s + 1
    rot[0, 1] = 2 * x * y * s * s - 2 * z * c * s
    rot[0, 2] = 2 * x * z * s * s + 2 * y * c * s
    rot[1, 0] = 2 * x * y * s * s + 2 * z * c * s
    rot[1, 1] = 2 * (y * y - 1) * s * s + 1
    rot[1, 2] = 2 * z * y * s * s - 2 * x * c * s
    rot[2, 0] = 2 * x * z * s * s - 2 * y * c * s
    rot[2, 1] = 2 * z * y * s * s + 2 * x * c * s
    rot[2, 2] = 2 * (z * z - 1) * s * s + 1

    # copy and centre ---------------------------------------------------- #
    pos = positions.copy()
    first, last = element[1], element[2]
    shift = -pos[first]

    # translate to origin, rotate, translate back
    for j in range(first, last):
        v = pos[j] + shift
        pos[j] = v @ rot.T - shift

    return pos


@njit(cache=True, fastmath=True)
def translate_kernel(
    positions: np.ndarray,
    element: Tuple[int, int, int],
    shift: np.ndarray,
) -> np.ndarray:
    """Translate atoms *element[0] … element[2]-1* by *shift* (3-vector)."""
    pos = positions.copy()
    start, end = element[0], element[2]
    for j in range(start, end):
        pos[j] += shift
    return pos


@njit(cache=True, fastmath=True)
def center_of_mass(positions: np.ndarray) -> np.ndarray:
    """Geometric centre (arithmetic mean) of positions."""
    n = positions.shape[0]
    total = np.zeros(3, dtype=positions.dtype)
    for i in range(n):
        total += positions[i]
    return total / n



@njit(cache=True, fastmath=True)
def radius(center: np.ndarray, positions: np.ndarray) -> float:
    """Maximum distance from *center* to any atom in *positions*."""
    diffs = positions - center
    r2 = (diffs**2).sum(axis=1)
    return math.sqrt(r2.max())


@njit(cache=True, fastmath=True)
def kl_divergence(sample: np.ndarray, reference: np.ndarray) -> float:
    """Kullback–Leibler divergence ∑ p log(p / q)."""
    return -(sample * np.log(sample / reference)).sum()


@njit(cache=True, fastmath=True)
def entropy(sample: np.ndarray) -> float:
    """Shannon entropy –∑ p log p for an un-normalised 1-D array."""
    p = sample / sample.sum()
    return -(p * np.log(_safe(p))).sum()


@njit(cache=True, fastmath=True)
def zps(energies: np.ndarray, beta: float = 1.0e-3):
    """
    Partition function Z, Boltzmann probabilities P, and entropy S
    for an array of energies (kJ mol⁻¹ or any consistent unit).
    """
    boltz = np.exp(-beta * energies)
    Z = boltz.sum()
    P = boltz / _safe(Z)
    S = entropy(P)
    return Z, P, S


@njit(cache=True, fastmath=True)
def S(energies: np.ndarray, beta: float = 1.0e-3) -> float:  # kept legacy name
    """Alias that returns the entropy only (compatibility shim)."""
    return zps(energies, beta)[2]


# --------------------------------------------------------------------------- #
# Back-compat camelCase aliases
# --------------------------------------------------------------------------- #
rotateKernel               = rotate_kernel
translateKernel            = translate_kernel
centerOfMass               = center_of_mass
kullbackLeiblerDivergenceKernel = kl_divergence
EntropyKernel              = entropy
ZPS                         = zps  # same as before

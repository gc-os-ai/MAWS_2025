# Kompilierte Versionen der rechenintensiven Operationen
# Compiled versions of computationally intensive operations.

from __future__ import annotations

import numpy as np
from numba import jit  # kept for compatibility with your code
from openmm import unit


def catch_zero(numeric):
    """
    Return numeric slightly nudged away from zero to avoid log/division issues.

    Parameters
    ----------
    numeric : array_like or float
        Value(s) that may contain zeros.

    Returns
    -------
    array_like or float
        `numeric` plus a tiny epsilon (1e-50) to avoid singularities.

    Notes
    -----
    This helper is used to stabilize `log` and divisions.
    """
    return numeric + 1e-50


@jit
def rotate_kernel(positions, element, axis, angle):
    """
    Rotate a contiguous slice of coordinates about an arbitrary axis through a pivot.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Cartesian coordinates in ångström **as plain floats** (no units).
    element : tuple of int
        Index triple controlling the affected region.
        By convention in this code:
          - `element[1]` is both the pivot index and the start index (inclusive),
          - `element[2]` is the end index (exclusive).
        All rows `j in [element[1], element[2])` are rotated about the pivot
        at `positions[element[1]]`.
    axis : array_like, shape (3,)
        Rotation axis vector (will be normalized internally).
    angle : float
        Rotation angle **in radians**.

    Returns
    -------
    openmm.unit.Quantity
        New coordinates with the same shape as `positions`, wrapped as
        `Quantity(..., unit.angstrom)`.

    Notes
    -----
    - The rotation matrix is built from the quaternion/axis–angle formula.
    - The implementation right-multiplies the row vector by the rotation matrix
      (`v @ R`). With the usual column-vector convention, this corresponds to
      applying `R^T`, i.e., a rotation by `-angle`. This preserves your original
      behavior.
    - Because this function returns an OpenMM `Quantity`, Numba will likely fall
      back to *object mode* (i.e., not nopython). For best performance, see the
      comment in the module doc or provide a unit-free kernel that returns a
      plain `ndarray` and wrap units outside.
    """
    x, y, z = np.asarray(axis) / (np.linalg.norm(np.asarray(axis)))
    phi_2 = angle / 2.0
    pos = np.array(positions[:])  # make a local copy
    shift_forward = -pos[element[1]]

    s = np.sin(phi_2)
    c = np.cos(phi_2)
    rot = np.array(
        [
            [
                2 * (np.power(x, 2) - 1) * np.power(s, 2) + 1,
                2 * x * y * np.power(s, 2) - 2 * z * c * s,
                2 * x * z * np.power(s, 2) + 2 * y * c * s,
            ],
            [
                2 * x * y * np.power(s, 2) + 2 * z * c * s,
                2 * (np.power(y, 2) - 1) * np.power(s, 2) + 1,
                2 * z * y * np.power(s, 2) - 2 * x * c * s,
            ],
            [
                2 * x * z * np.power(s, 2) - 2 * y * c * s,
                2 * z * y * np.power(s, 2) + 2 * x * c * s,
                2 * (np.power(z, 2) - 1) * np.power(s, 2) + 1,
            ],
        ]
    )

    # translate the affected block so the pivot sits at the origin
    for j in range(element[1], element[2]):
        pos[j] += shift_forward

    # rotate each selected point; note the right-multiplication v @ R
    for j in range(element[1], element[2]):
        roted = np.dot(pos[j], rot)
        pos[j] = roted - shift_forward  # translate back

    return unit.Quantity(pos, unit.angstrom)


@jit
def translate_kernel(positions, element, shift):
    """
    Translate a contiguous slice of coordinates by a fixed vector.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Cartesian coordinates (plain floats).
    element : tuple of int
        Index triple where `element[0]` is the start (inclusive) and
        `element[2]` is the end (exclusive). This mirrors your original code.
    shift : array_like, shape (3,)
        Translation vector in the same units as `positions`.

    Returns
    -------
    ndarray, shape (N, 3)
        Translated coordinates as a new array.
    """
    pos = np.array(positions)
    for j in range(element[0], element[2]):
        pos[j] += shift
    return pos


@jit
def center_of_mass(positions):
    """
    Compute the arithmetic center (centroid) of coordinates.

    Parameters
    ----------
    positions : array_like, shape (N, 3)
        Cartesian coordinates.

    Returns
    -------
    ndarray, shape (3,)
        Coordinate-wise mean.
    """
    return positions.sum(axis=0) / len(positions)


@jit
def radius(center, positions):
    """
    Maximum Euclidean distance from `center` to any point in `positions`.

    Parameters
    ----------
    center : array_like, shape (3,)
        Reference point.
    positions : array_like, shape (N, 3)
        Points to measure.

    Returns
    -------
    float
        `max_j ||positions[j] - center||_2`
    """
    return max(map(np.linalg.norm, np.asarray(positions) - np.asarray(center)))


@jit
def kl_divergence_kernel(sample, reference_sample):
    """
    Return **negative** Kullback–Leibler divergence of `sample` vs `reference_sample`.

    Parameters
    ----------
    sample : array_like, shape (M,)
        Nonnegative weights (often probabilities). Should sum to 1 for
        information-theoretic interpretation.
    reference_sample : array_like, shape (M,)
        Reference distribution; same shape as `sample`.

    Returns
    -------
    float
        `-(sum_i sample_i * log(sample_i / reference_i))`

    Notes
    -----
    - Standard KL is `sum p * log(p/q)` and is nonnegative.
      This function returns its **negative** (matching your original code).
    - Inputs are converted to `ndarray` and not normalized automatically.
    """
    return -(
        np.array(sample) * np.log(np.array(sample) / np.array(reference_sample))
    ).sum()


@jit
def entropy_kernel(sample):
    """
    Entropy relative to the uniform baseline (shifted by -log M).

    Parameters
    ----------
    sample : array_like, shape (M,)
        Nonnegative weights (often probabilities). Typically sum to 1.

    Returns
    -------
    float
        `-(sum_i p_i * log(p_i * M)) = H(p) - log M`, where `M=len(sample)`.

    Notes
    -----
    - This equals the usual Shannon entropy `H(p)` (nats) minus `log M`.
      It is 0 for the uniform distribution and negative otherwise.
    - `catchZero` is used to avoid `log(0)`.
    """
    return -(
        np.array(sample) * np.log(catch_zero(np.array(sample) * len(sample)))
    ).sum()


@jit
def ZPS(sample, beta=0.001):  # noqa: N802
    """
    Boltzmann partition Z, probabilities P, and (shifted) entropy S of P.

    Parameters
    ----------
    sample : array_like, shape (M,)
        Energy-like values `E_i`. Lower is better.
    beta : float, optional
        Inverse temperature. Default is 1e-3.

    Returns
    -------
    Z : float
        Partition function `Z = sum_i exp(-beta * E_i)`.
    P : ndarray, shape (M,)
        Boltzmann probabilities `P_i = exp(-beta * E_i) / Z`.
    S : float
        `EntropyKernel(P)`; note this is `H(P) - log M` (≤ 0).

    Notes
    -----
    This is a softmin over `sample`. As `beta` increases, P becomes peakier.
    """
    Z = np.exp(-beta * np.array(sample)).sum()
    P = np.exp(-beta * np.array(sample)) / catch_zero(Z)
    S = entropy_kernel(P)
    return Z, P, S


@jit
def S(sample, beta=0.001):  # noqa: N802
    """
    Convenience wrapper returning only the entropy from `ZPS`.

    Parameters
    ----------
    sample : array_like, shape (M,)
        Energy-like values.
    beta : float, optional
        Inverse temperature.

    Returns
    -------
    float
        `EntropyKernel(P)` computed on the Boltzmann probabilities of `sample`.
    """
    return ZPS(sample, beta)[2]

# helpers.py
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from openmm import unit
from openmm.unit import Quantity


def angstrom(array: ArrayLike) -> Quantity:
    """
    Attach Å units to a numeric array or vector.

    Parameters
    ----------
    array : array-like
        Numeric values interpreted as lengths in Å (unitless on input).

    Returns
    -------
    openmm.unit.Quantity
        The same values with units of Å (unit.angstrom).

    Notes
    -----
    This is the *only* place you should attach length units to raw arrays.
    Elsewhere, pass around unit-bearing vectors.
    """
    return np.asarray(array, dtype=float) * unit.angstrom


def nostrom(quantity: Quantity) -> np.ndarray:
    """
    Strip units from a length vector/array, returning pure Å as floats.

    Parameters
    ----------
    quantity : openmm.unit.Quantity
        Length(s) with units (must be convertible to Å).

    Returns
    -------
    numpy.ndarray
        The numeric values in Å, without units.

    Raises
    ------
    AttributeError
        If a unitless array is passed. Keep this strict to avoid silent mistakes.
    """
    return quantity.value_in_unit(unit.angstrom)


def kJ(array: ArrayLike) -> Quantity:  # noqa: N802
    """
    Attach kJ/mol to a numeric value or array.

    Parameters
    ----------
    array : array-like
        Unitless energy values.

    Returns
    -------
    openmm.unit.Quantity
        Values with units of kJ/mol.
    """
    return np.asarray(array, dtype=float) * unit.kilojoules_per_mole


def noJ(quantity: Quantity) -> float | np.ndarray:  # noqa: N802
    """
    Strip units from an energy, returning kJ/mol as plain numbers.

    Parameters
    ----------
    quantity : openmm.unit.Quantity
        Energy with units.

    Returns
    -------
    float or numpy.ndarray
        Numeric value(s) in kJ/mol.
    """
    return quantity.value_in_unit(unit.kilojoules_per_mole)


def angle(array1: ArrayLike, array2: ArrayLike) -> float:
    """
    Return the unsigned angle between two vectors (radians).

    Parameters
    ----------
    array1, array2 : array-like
        Unitless vectors. If you have unit-bearing vectors, strip first with
        :func:`nostrom`.

    Returns
    -------
    float
        Angle in radians (0..π).
    """
    a = np.asarray(array1, dtype=float)
    b = np.asarray(array2, dtype=float)
    return float(
        np.arccos(
            np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1.0, 1.0)
        )
    )


def directed_angle(array1: ArrayLike, array2: ArrayLike, axis: ArrayLike) -> float:
    """
    Return signed angle from array1 to array2 about 'axis' (radians).

    Parameters
    ----------
    array1, array2, axis : array-like
        Unitless vectors (use :func:`nostrom` beforehand if needed).

    Returns
    -------
    float
        Signed angle in radians (-π..π).
    """
    a1 = np.asarray(array1, dtype=float)
    a1 /= np.linalg.norm(a1)
    a2 = np.asarray(array2, dtype=float)
    a2 /= np.linalg.norm(a2)
    ax = np.asarray(axis, dtype=float)
    return float(np.arctan2(np.dot(ax, np.cross(a1, a2)), np.dot(a1, a2)))

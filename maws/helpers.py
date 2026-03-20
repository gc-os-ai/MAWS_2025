"""
maws.helpers
============

Unit conversion and vector utilities for MAWS.

This module provides functions to attach/strip OpenMM units from numeric arrays
and compute angles between vectors. All geometry in MAWS uses Ångströms (Å).

Functions
---------
angstrom : Attach Å units to raw numbers.
nostrom : Strip Å units to get raw floats.
kJ : Attach kJ/mol units to energies.
noJ : Strip kJ/mol units from energies.
angle : Compute unsigned angle between vectors.
directed_angle : Compute signed angle about an axis.

Examples
--------
>>> from maws.helpers import angstrom, nostrom
>>> import numpy as np
>>> vec = angstrom([1.0, 2.0, 3.0])
>>> raw = nostrom(vec)
>>> list(raw)
[1.0, 2.0, 3.0]
"""

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

    Examples
    --------
    >>> from maws.helpers import angstrom
    >>> vec = angstrom([1.0, 2.0, 3.0])
    >>> vec.unit == unit.angstrom
    True
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

    Examples
    --------
    >>> from maws.helpers import angstrom, nostrom
    >>> vec = angstrom([5.0, 10.0, 15.0])
    >>> raw = nostrom(vec)
    >>> list(raw)
    [5.0, 10.0, 15.0]
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

    Examples
    --------
    >>> from maws.helpers import kJ
    >>> energy = kJ(100.5)
    >>> energy.unit == unit.kilojoules_per_mole
    True
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

    Examples
    --------
    >>> from maws.helpers import kJ, noJ
    >>> energy = kJ(250.0)
    >>> float(noJ(energy))
    250.0
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

    Examples
    --------
    >>> from maws.helpers import angle
    >>> import numpy as np
    >>> angle([1, 0, 0], [0, 1, 0])  # 90 degrees
    1.5707963267948966
    >>> round(angle([1, 0, 0], [1, 0, 0]), 5)  # Same vector = 0
    0.0
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

    Examples
    --------
    >>> from maws.helpers import directed_angle
    >>> import numpy as np
    >>> # Rotate [1,0,0] to [0,1,0] about z-axis = +90 degrees
    >>> round(directed_angle([1, 0, 0], [0, 1, 0], [0, 0, 1]), 5)
    1.5708
    """
    a1 = np.asarray(array1, dtype=float)
    a1 /= np.linalg.norm(a1)
    a2 = np.asarray(array2, dtype=float)
    a2 /= np.linalg.norm(a2)
    ax = np.asarray(axis, dtype=float)
    return float(np.arctan2(np.dot(ax, np.cross(a1, a2)), np.dot(a1, a2)))


def center_of_mass(positions: ArrayLike) -> np.ndarray:
    """
    Compute the arithmetic center (centroid) of coordinates.

    Parameters
    ----------
    positions : array-like, shape (N, 3)
        Cartesian coordinates.

    Returns
    -------
    numpy.ndarray, shape (3,)
        Coordinate-wise mean (centroid).

    Examples
    --------
    >>> from maws.helpers import center_of_mass
    >>> import numpy as np
    >>> coords = np.array([[0, 0, 0], [2, 0, 0], [1, 2, 0]])
    >>> center_of_mass(coords)
    array([1.        , 0.66666667, 0.        ])
    """
    pos = np.asarray(positions, dtype=float)
    return pos.sum(axis=0) / len(pos)

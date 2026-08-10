"""
maws.geometry
=============

Small vector calculations shared by several parts of MAWS.

Everything here works on plain NumPy arrays of numbers in ångström and returns
plain numbers. There is no chemistry and no state.

Examples
--------
>>> import numpy as np
>>> round(angle_between(np.array([1.0, 0, 0]), np.array([0.0, 1, 0])), 6)
1.570796
"""

from __future__ import annotations

import numpy as np

from maws.errors import ConfigurationError

__all__ = ["angle_between", "centre_of_mass", "unit_vector"]


def unit_vector(vector: np.ndarray) -> np.ndarray:
    """Return a vector of length one pointing the same way as `vector`.

    Parameters
    ----------
    vector : numpy.ndarray
        Shape ``(3,)``. Any non-zero vector.

    Returns
    -------
    numpy.ndarray
        Shape ``(3,)``, with length 1.

    Raises
    ------
    maws.errors.ConfigurationError
        If `vector` has zero length, in which case it points nowhere.

    Examples
    --------
    >>> unit_vector(np.array([0.0, 3.0, 0.0]))
    array([0., 1., 0.])
    """
    length = float(np.linalg.norm(vector))
    if length == 0.0:
        raise ConfigurationError("a zero-length vector has no direction")
    return np.asarray(vector, dtype=np.float64) / length


def angle_between(first: np.ndarray, second: np.ndarray) -> float:
    """Return the angle between two directions, in radians.

    Parameters
    ----------
    first, second : numpy.ndarray
        Shape ``(3,)``. Their lengths are ignored; only direction matters.

    Returns
    -------
    float
        A value between 0 and π. Zero when the two point the same way, π when
        they point in opposite directions.

    Raises
    ------
    maws.errors.ConfigurationError
        If either vector has zero length.

    Notes
    -----
    The cosine is clipped to ``[-1, 1]`` before the arc-cosine is taken.
    Rounding in floating-point arithmetic can push it a hair outside that
    range for two nearly-parallel vectors, which would otherwise produce a
    not-a-number result.

    Examples
    --------
    >>> import numpy as np
    >>> round(angle_between(np.array([1.0, 0, 0]), np.array([-1.0, 0, 0])), 6)
    3.141593
    """
    cosine = float(np.dot(unit_vector(first), unit_vector(second)))
    return float(np.arccos(np.clip(cosine, -1.0, 1.0)))


def centre_of_mass(positions: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Return the balance point of a group of atoms.

    Heavier atoms pull the result towards themselves, so this is not the same
    as the plain average position.

    Parameters
    ----------
    positions : numpy.ndarray
        Shape ``(N, 3)``, in ångström.
    masses : numpy.ndarray
        Shape ``(N,)``, in daltons. One entry per atom.

    Returns
    -------
    numpy.ndarray
        Shape ``(3,)``, in ångström.

    Raises
    ------
    maws.errors.ConfigurationError
        If the two arrays disagree on how many atoms there are, or the masses
        add up to zero.

    Examples
    --------
    A light atom and a heavy one: the balance point sits near the heavy one.

    >>> import numpy as np
    >>> centre_of_mass(np.array([[0.0, 0, 0], [10.0, 0, 0]]), np.array([1.0, 9.0]))
    array([9., 0., 0.])
    """
    positions = np.asarray(positions, dtype=np.float64)
    masses = np.asarray(masses, dtype=np.float64)
    if positions.shape[0] != masses.shape[0]:
        raise ConfigurationError(
            f"got {positions.shape[0]} positions but {masses.shape[0]} masses"
        )
    total = float(masses.sum())
    if total == 0.0:
        raise ConfigurationError("total mass is zero, so there is no balance point")
    return (positions * masses[:, None]).sum(axis=0) / total

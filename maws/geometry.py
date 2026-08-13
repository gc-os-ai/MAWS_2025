"""
maws.geometry
=============

Small vector calculations shared by several parts of MAWS.

Three functions, each one line of arithmetic wrapped in a check that the input
makes sense: a direction of length one, the angle between two directions, and
the mass-weighted balance point of a group of atoms.

Everything here works on plain NumPy arrays and returns plain numbers.
Positions are in ångström and masses in daltons; angles come back in radians.
There is no chemistry and no state, so every function can be called with
numbers made up on the spot.

See Also
--------
maws.pose.rodrigues : Builds the rotation matrix these directions feed.

Examples
--------
>>> import numpy as np
>>> round(angle_between(np.array([1.0, 0, 0]), np.array([0.0, 1, 0])), 6)
1.570796
>>> unit_vector(np.array([0.0, 0.0, 5.0]))
array([0., 0., 1.])
"""

from __future__ import annotations

import numpy as np

from maws.errors import ConfigurationError

__all__ = ["angle_between", "centre_of_mass", "unit_vector"]


def unit_vector(vector: np.ndarray) -> np.ndarray:
    """Return a vector of length one pointing the same way as `vector`.

    Scaling a vector to length one throws away how long it was and keeps only
    where it points, which is what a rotation axis or a viewing direction
    needs.

    Parameters
    ----------
    vector : numpy.ndarray
        Shape ``(3,)``. Any non-zero vector, in whatever unit; the unit
        cancels.

    Returns
    -------
    numpy.ndarray
        Shape ``(3,)``, with length 1. A new array; `vector` is unchanged.

    Raises
    ------
    maws.errors.ConfigurationError
        If `vector` has zero length, in which case it points nowhere. In
        practice that means two atoms were found at the same position.

    See Also
    --------
    angle_between : Compares two directions, using this on each of them.

    Examples
    --------
    >>> import numpy as np
    >>> unit_vector(np.array([0.0, 3.0, 0.0]))
    array([0., 1., 0.])
    """
    length = float(np.linalg.norm(vector))
    if length == 0.0:
        raise ConfigurationError("a zero-length vector has no direction")
    return np.asarray(vector, dtype=np.float64) / length


def angle_between(first: np.ndarray, second: np.ndarray) -> float:
    r"""Return the angle between two directions, in radians.

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

    See Also
    --------
    unit_vector : What both arguments are reduced to first.

    Notes
    -----
    The angle is the arc-cosine of the dot product of the two directions:

    .. math::
        \theta = \arccos\!\left(
            \frac{\mathbf{a} \cdot \mathbf{b}}
                 {\lVert\mathbf{a}\rVert \, \lVert\mathbf{b}\rVert}
        \right)

    The cosine is clipped to ``[-1, 1]`` before the arc-cosine is taken.
    Rounding in floating-point arithmetic can push it a hair outside that
    range for two nearly-parallel vectors, which would otherwise produce a
    not-a-number result.

    The answer never exceeds π, so it says how far apart two directions are
    but not which way round they are. Use a signed torsion angle where the
    sense of the turn matters.

    Examples
    --------
    >>> import numpy as np
    >>> round(angle_between(np.array([1.0, 0, 0]), np.array([-1.0, 0, 0])), 6)
    3.141593

    Length makes no difference, only direction:

    >>> round(angle_between(np.array([2.0, 0, 0]), np.array([0.0, 0.5, 0])), 6)
    1.570796
    """
    cosine = float(np.dot(unit_vector(first), unit_vector(second)))
    return float(np.arccos(np.clip(cosine, -1.0, 1.0)))


def centre_of_mass(positions: np.ndarray, masses: np.ndarray) -> np.ndarray:
    r"""Return the balance point of a group of atoms.

    Heavier atoms pull the result towards themselves, so this is not the same
    as the plain average position.

    Parameters
    ----------
    positions : numpy.ndarray
        Shape ``(N, 3)`` positions in ångström, one row per atom.
    masses : numpy.ndarray
        Shape ``(N,)`` atomic masses in daltons, in the same order as
        `positions`. Passing every atom the same mass turns this into the plain
        average position.

    Returns
    -------
    numpy.ndarray
        Shape ``(3,)``, in ångström.

    Raises
    ------
    maws.errors.ConfigurationError
        If the two arrays disagree on how many atoms there are, or the masses
        add up to zero.

    See Also
    --------
    maws.pose.Pose.centroid : The unweighted average, which needs no masses.

    Notes
    -----
    .. math::
        \mathbf{c} = \frac{\sum_i m_i \, \mathbf{r}_i}{\sum_i m_i}

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

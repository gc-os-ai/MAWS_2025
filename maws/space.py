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

import numpy as np
from openmm import unit


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

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
>>> from maws.space import Cube, NAngles

Create a 20Å cube centered at the origin:

>>> cube = Cube(20.0, np.array([0.0, 0.0, 0.0]))
>>> sample = cube.generator()
>>> len(sample)  # [x, y, z, axis_x, axis_y, axis_z, angle]
7

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


class Box(Space):
    """
    A 3D rectangular box sampling space with rotation.

    Generates samples uniformly within the box, plus a random rotation
    axis and angle for molecular orientation.

    Parameters
    ----------
    x_width : float
        Width in the x-direction (Å).
    y_width : float
        Width in the y-direction (Å).
    z_width : float
        Width in the z-direction (Å).
    centre : array-like
        3D coordinates of the box center.
    units : openmm.unit, default=unit.angstroms
        Length units.

    Examples
    --------
    >>> import numpy as np
    >>> box = Box(10.0, 20.0, 15.0, np.array([0.0, 0.0, 0.0]))
    >>> sample = box.generator()
    >>> len(sample)  # [x, y, z, axis_x, axis_y, axis_z, rotation_angle]
    7
    """

    def __init__(self, x_width, y_width, z_width, centre, units=unit.angstroms):
        def boolean(position):
            x_in = (
                position[0] <= centre[0] + x_width / 2
                and position[0] >= centre[0] - x_width / 2
            )
            y_in = (
                position[1] <= centre[1] + y_width / 2
                and position[1] >= centre[1] - y_width / 2
            )
            z_in = (
                position[2] <= centre[2] + z_width / 2
                and position[2] >= centre[2] - z_width / 2
            )
            return x_in and y_in and z_in

        super().__init__(boolean, units=units)
        self.x_width = x_width
        self.y_width = y_width
        self.z_width = z_width
        self.centre = centre
        self.volume = x_width * y_width * z_width * (2 * np.pi) ** 3

    def generator(self):
        """
        Generate a random sample from the box.

        Returns
        -------
        numpy.ndarray
            7-element array: [x, y, z, axis_x, axis_y, axis_z, rotation_angle].
            The first 3 are position; next 3 are a random rotation axis;
            last is a rotation angle in [0, π].
        """
        axis = np.random.uniform(-1, 1, 3)
        x, y, z = axis / np.linalg.norm(axis)
        rotation_angle = np.random.uniform(0, np.pi)
        return np.array(
            [
                np.random.uniform(
                    self.centre[0] - self.x_width / 2,
                    self.centre[0] + self.x_width / 2,
                ),
                np.random.uniform(
                    self.centre[1] - self.y_width / 2,
                    self.centre[1] + self.y_width / 2,
                ),
                np.random.uniform(
                    self.centre[2] - self.z_width / 2,
                    self.centre[2] + self.z_width / 2,
                ),
                x,
                y,
                z,
                rotation_angle,
            ]
        )


class Cube(Box):
    """
    A 3D cubic sampling space (equal width in all dimensions).

    Parameters
    ----------
    width : float
        Side length in all directions (Å).
    centre : array-like
        3D coordinates of the cube center.
    units : openmm.unit, default=unit.angstroms
        Length units.

    Examples
    --------
    >>> import numpy as np
    >>> cube = Cube(20.0, np.array([0.0, 0.0, 0.0]))
    >>> sample = cube.generator()
    >>> len(sample)  # [x, y, z, axis_x, axis_y, axis_z, rotation_angle]
    7
    >>> -10.0 <= sample[0] <= 10.0  # Within ±width/2 of center
    True
    """

    def __init__(self, width, centre, units=unit.angstroms):
        super().__init__(width, width, width, centre, units=units)


class Sphere(Space):
    """
    A 3D spherical sampling space with rotation.

    Generates samples uniformly within the sphere, plus a random rotation
    axis and angle for molecular orientation.

    Parameters
    ----------
    radius : float
        Radius of the sphere (Å).
    centre : array-like
        3D coordinates of the sphere center.
    units : openmm.unit, default=unit.angstroms
        Length units.

    Examples
    --------
    >>> import numpy as np
    >>> sphere = Sphere(15.0, np.array([0.0, 0.0, 0.0]))
    >>> sample = sphere.generator()
    >>> len(sample)  # [x, y, z, axis_x, axis_y, axis_z, rotation_angle]
    7
    """

    def __init__(self, radius, centre, units=unit.angstroms):
        def boolean(position):
            return np.linalg.norm(position - centre) <= radius

        super().__init__(boolean, units=units)
        self.radius = radius
        self.centre = centre
        self.volume = 4.0 / 3.0 * np.pi * radius**3 * (2 * np.pi) ** 3

    def generator(self):
        """
        Generate a random sample from the sphere.

        Returns
        -------
        numpy.ndarray
            7-element array: [x, y, z, axis_x, axis_y, axis_z, rotation_angle].
        """
        axis = np.random.uniform(-1, 1, 3)
        x, y, z = axis / np.linalg.norm(axis)
        # r = np.random.uniform(0, self.radius)
        r = self.radius * np.random.uniform(0, 1) ** (
            1 / 3
        )  #  uniform-in-volume sampling, use u**(1/3)
        phi = np.random.uniform(0, 2 * np.pi)
        # psi = np.random.uniform(0, np.pi)
        # Uniform direction: cos(psi) uniform in [-1, 1]
        cos_psi = np.random.uniform(-1, 1)
        psi = np.arccos(cos_psi)
        rotation_angle = np.random.uniform(0, 2 * np.pi)
        result = np.array(
            [
                r * np.cos(phi) * np.sin(psi),
                r * np.sin(phi) * np.sin(psi),
                r * np.cos(psi),
                x,
                y,
                z,
                rotation_angle,
            ]
        )
        return result


class SphericalShell(Space):
    """
    A 3D spherical shell (hollow sphere) sampling space.

    Generates samples uniformly within the shell between inner and outer radii.

    Parameters
    ----------
    innerRadius : float
        Inner radius of the shell (Å).
    outerRadius : float
        Outer radius of the shell (Å).
    centre : array-like
        3D coordinates of the shell center.
    units : openmm.unit, default=unit.angstroms
        Length units.

    Examples
    --------
    >>> import numpy as np
    >>> shell = SphericalShell(5.0, 10.0, np.array([0.0, 0.0, 0.0]))
    >>> sample = shell.generator()
    >>> len(sample)
    7
    """

    def __init__(self, innerRadius, outerRadius, centre, units=unit.angstroms):
        def boolean(position):
            r = np.linalg.norm(position - centre)
            return innerRadius <= r <= outerRadius

        super().__init__(boolean, units=units)
        self.innerRadius = innerRadius
        self.outerRadius = outerRadius
        self.centre = centre
        self.volume = (
            4.0 / 3.0 * np.pi * outerRadius**3 * (2 * np.pi) ** 3
            - 4.0 / 3.0 * np.pi * innerRadius**3 * (2 * np.pi) ** 3
        )

    def generator(self):
        """
        Generate a random sample from the spherical shell.

        Returns
        -------
        numpy.ndarray
            7-element array: [x, y, z, axis_x, axis_y, axis_z, rotation_angle].
        """
        axis = np.random.uniform(-1, 1, 3)
        x, y, z = axis / np.linalg.norm(axis)
        r = np.random.uniform(self.innerRadius, self.outerRadius)
        phi = np.random.uniform(0, 2 * np.pi)
        psi = np.random.uniform(0, np.pi)
        rotation_angle = np.random.uniform(0, 2 * np.pi)
        result = np.array(
            [
                r * np.cos(phi) * np.sin(psi),
                r * np.sin(phi) * np.sin(psi),
                r * np.cos(psi),
                x,
                y,
                z,
                rotation_angle,
            ]
        )
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

"""spaces.py – modern random‐sampling utilities for MAWS.

All spaces expose a :py:meth:`generator` method that returns a NumPy array
containing **Cartesian coordinates** *(x, y, z)* followed by an axis
*(unit vector x, y, z)* and a rotation angle *θ* (rad).  The exact length of
the vector therefore depends on the shape:

* Box / Cube / Sphere / SphericalShell – 7 floats
* NAngles – N floats

The module is self-contained; OpenMM is imported only for its unit objects.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from typing import Callable, Sequence

import numpy as np
from numpy.typing import NDArray

from openmm import unit  # noqa: F401  (imported for side-effects / user code)

# Random-number generator (thread-safe, reproducible if needed)
_rng = np.random.default_rng()


# --------------------------------------------------------------------------- #
# Abstract base class
# --------------------------------------------------------------------------- #
class Space(ABC):
    """Abstract parent for all sampling spaces."""

    def __init__(
        self,
        membership_boolean_function: Callable[[NDArray[np.floating]], bool],
        *,
        lower_bound: float = -1.0e6,
        upper_bound: float = 1.0e6,
        units: unit = unit.angstrom,
    ) -> None:
        self._is_in = membership_boolean_function
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.units = units
        self.volume: float | None = None  # subclasses fill in

    # ------------------------------------------------------------------ #
    # Public helpers
    # ------------------------------------------------------------------ #
    def sample_point(self) -> NDArray[np.floating]:
        """Return one random point that lies inside the space."""
        while True:
            candidate: NDArray[np.floating] = _rng.uniform(
                self.lower_bound, self.upper_bound, 3
            )
            if self._is_in(candidate):
                return candidate

    @abstractmethod
    def generator(self) -> NDArray[np.floating]:
        """Return a random configuration (position + orientation)."""


# --------------------------------------------------------------------------- #
# Concrete shapes
# --------------------------------------------------------------------------- #
class Box(Space):
    """Axis-aligned rectangular box of size ``x×y×z`` centred at *centre*."""

    def __init__(
        self,
        x_width: float,
        y_width: float,
        z_width: float,
        centre: Sequence[float] = (0.0, 0.0, 0.0),
        *,
        units: unit = unit.angstrom,
    ) -> None:
        cx, cy, cz = centre

        def _inside(p: NDArray[np.floating]) -> bool:
            return (
                cx - x_width / 2 <= p[0] <= cx + x_width / 2
                and cy - y_width / 2 <= p[1] <= cy + y_width / 2
                and cz - z_width / 2 <= p[2] <= cz + z_width / 2
            )

        super().__init__(_inside, units=units)
        self.x_width, self.y_width, self.z_width = x_width, y_width, z_width
        self.centre = np.asarray(centre, dtype=float)
        self.volume = x_width * y_width * z_width  # no (2π)³ by default

    # -- orientation helper ------------------------------------------------ #
    @staticmethod
    def _random_axis_and_angle() -> tuple[NDArray[np.floating], float]:
        axis = _rng.uniform(-1.0, 1.0, 3)
        axis /= np.linalg.norm(axis)
        angle = _rng.uniform(0.0, math.pi)
        return axis, angle

    # ------------------------------------------------------------------ #
    # API
    # ------------------------------------------------------------------ #
    def generator(self) -> NDArray[np.floating]:
        x = _rng.uniform(
            self.centre[0] - self.x_width / 2, self.centre[0] + self.x_width / 2
        )
        y = _rng.uniform(
            self.centre[1] - self.y_width / 2, self.centre[1] + self.y_width / 2
        )
        z = _rng.uniform(
            self.centre[2] - self.z_width / 2, self.centre[2] + self.z_width / 2
        )
        axis, angle = self._random_axis_and_angle()
        return np.concatenate(([x, y, z], axis, [angle]))


class Cube(Box):
    """Cube of edge length *width* centred at *centre*."""

    def __init__(self, width: float, centre: Sequence[float] = (0.0, 0.0, 0.0)):
        super().__init__(width, width, width, centre)


class Sphere(Space):
    """Solid sphere of radius *r*."""

    def __init__(
        self, radius: float, centre: Sequence[float] = (0.0, 0.0, 0.0), *, units=unit.angstrom
    ):
        centre_arr = np.asarray(centre, dtype=float)

        def _inside(p: NDArray[np.floating]) -> bool:
            return np.linalg.norm(p - centre_arr) <= radius

        super().__init__(_inside, units=units)
        self.radius = radius
        self.centre = centre_arr
        self.volume = 4.0 / 3.0 * math.pi * radius**3

    def generator(self) -> NDArray[np.floating]:
        axis = _rng.uniform(-1.0, 1.0, 3)
        axis /= np.linalg.norm(axis)
        r = _rng.uniform(0.0, self.radius)
        phi = _rng.uniform(0.0, 2.0 * math.pi)
        theta = _rng.uniform(0.0, math.pi)
        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        angle = _rng.uniform(0.0, 2.0 * math.pi)
        return np.concatenate((self.centre + np.array([x, y, z]), axis, [angle]))


class SphericalShell(Space):
    """Hollow sphere with inner radius *r_in* and outer radius *r_out*."""

    def __init__(
        self,
        inner_radius: float,
        outer_radius: float,
        centre: Sequence[float] = (0.0, 0.0, 0.0),
        *,
        units=unit.angstrom,
    ):
        centre_arr = np.asarray(centre, dtype=float)

        def _inside(p: NDArray[np.floating]) -> bool:
            dist = np.linalg.norm(p - centre_arr)
            return inner_radius <= dist <= outer_radius

        super().__init__(_inside, units=units)
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.centre = centre_arr
        self.volume = 4.0 / 3.0 * math.pi * (outer_radius**3 - inner_radius**3)

    def generator(self) -> NDArray[np.floating]:
        axis = _rng.uniform(-1.0, 1.0, 3)
        axis /= np.linalg.norm(axis)
        r = _rng.uniform(self.inner_radius, self.outer_radius)
        phi = _rng.uniform(0.0, 2.0 * math.pi)
        theta = _rng.uniform(0.0, math.pi)
        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        angle = _rng.uniform(0.0, 2.0 * math.pi)
        return np.concatenate((self.centre + np.array([x, y, z]), axis, [angle]))


# --------------------------------------------------------------------------- #
# Angle space
# --------------------------------------------------------------------------- #
class NAngles(Space):
    """Direct product of *N* circles, i.e. ``U(1)^N``."""

    def __init__(self, number: int):
        def _inside(p: NDArray[np.floating]) -> bool:  # always true for angles
            return np.all((0.0 <= p) & (p <= 2.0 * math.pi))

        super().__init__(_inside)
        self.number = number
        self.volume = (2.0 * math.pi) ** number

    def generator(self) -> NDArray[np.floating]:
        return _rng.uniform(0.0, 2.0 * math.pi, self.number)


# # --------------------------------------------------------------------------- #
# # Quick smoke-test
# # --------------------------------------------------------------------------- #
# if __name__ == "__main__":
#     cube = Cube(width=10.0, centre=(0.0, 0.0, 0.0))
#     print("Cube sample:", cube.generator())

#     shell = SphericalShell(inner_radius=5.0, outer_radius=7.0)
#     print("Shell sample:", shell.generator())

#     angles = NAngles(4)
#     print("Angles:", angles.generator())

"""Where a sampled pose puts the aptamer.

``maws.space`` hands back an absolute point in the target's coordinate frame.
The aptamer's centre of mass has to end up *on* that point, whatever
orientation was drawn with it and wherever LEaP happened to build the strand.

Two ways that can go wrong, both pinned down here:

- treating the point as a displacement, so the strand lands at
  ``where_it_already_was + point``;
- turning the strand about one of its own atoms instead of its centre of
  mass, which drags the centre away from the point again.

See issue #48.
"""

from types import SimpleNamespace

import numpy as np
import openmm as mm
import pytest
from openmm import unit

from maws.complex import Complex

# Atoms 0-3 are the "aptamer", atoms 4-5 the rigid "target" that must not move.
ELEMENT = [0, 1, 4]
MASSES = [1.0, 12.0, 12.0, 16.0, 31.0, 31.0]
COORDS = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [10.0, 10.0, 10.0],
        [11.0, 10.0, 10.0],
    ]
)
# (1*[0,0,0] + 12*[1,0,0] + 12*[1,1,0] + 16*[0,1,0]) / 41. The heavy atoms
# pull it well off the plain average of the four coordinates, [0.5, 0.5, 0],
# so nothing here can pass by using the wrong kind of centre.
BUILT_CENTRE = np.array([24 / 41, 28 / 41, 0.0])

POINT = np.array([5.0, -3.0, 2.0])


def make_complex(coords=COORDS):
    """A Complex carrying positions and masses, and nothing else.

    Placement is coordinate arithmetic, so this needs neither LEaP nor a real
    OpenMM System - only something that answers ``.topology.atoms()`` with one
    mass per atom.
    """
    cpx = Complex()
    cpx.positions = [mm.Vec3(*row) for row in coords] * unit.angstrom
    cpx.topology = SimpleNamespace(
        atoms=lambda: [
            SimpleNamespace(element=SimpleNamespace(mass=m * unit.dalton))
            for m in MASSES
        ]
    )
    return cpx


class TestPlaceGlobal:
    """Tests for ``Complex.place_global``."""

    def test_centre_of_mass_lands_on_the_requested_point(self):
        """The element's centre of mass ends up exactly on the sampled point."""
        cpx = make_complex()

        cpx.place_global(
            ELEMENT,
            POINT * unit.angstrom,
            np.array([0.0, 0.0, 1.0]) * unit.angstrom,
            0.7,
        )

        after = np.asarray(cpx.positions.value_in_unit(unit.angstrom))
        masses = np.array(MASSES[0:4])
        centre = (after[0:4] * masses[:, None]).sum(axis=0) / masses.sum()
        assert np.allclose(centre, POINT, atol=1e-9)

    def test_landing_point_does_not_depend_on_where_the_element_started(self):
        """Moving the strand first, then placing it, gives the same result.

        This is the regression test for the point being used as a
        displacement: adding it to the current coordinates makes the outcome
        depend on where LEaP happened to build the strand.
        """
        axis = np.array([1.0, 1.0, 0.0]) * unit.angstrom

        here = make_complex()
        here.place_global(ELEMENT, POINT * unit.angstrom, axis, 0.9)

        moved = COORDS.copy()
        moved[0:4] += np.array([100.0, -50.0, 7.0])
        elsewhere = make_complex(moved)
        elsewhere.place_global(ELEMENT, POINT * unit.angstrom, axis, 0.9)

        assert np.allclose(
            here.positions.value_in_unit(unit.angstrom),
            elsewhere.positions.value_in_unit(unit.angstrom),
        )

    def test_orientation_does_not_change_where_the_centre_lands(self):
        """Two different turns place the centre of mass on the same point.

        This is the regression test for turning about the element's first
        atom. That moves the centre of mass by up to twice its distance from
        that atom, so the drawn angle would decide the final position.
        """
        axis = np.array([0.3, -0.6, 0.5]) * unit.angstrom
        masses = np.array(MASSES[0:4])

        for angle in (0.0, 2.4):
            cpx = make_complex()
            cpx.place_global(ELEMENT, POINT * unit.angstrom, axis, angle)

            after = np.asarray(cpx.positions.value_in_unit(unit.angstrom))
            centre = (after[0:4] * masses[:, None]).sum(axis=0) / masses.sum()
            assert np.allclose(centre, POINT, atol=1e-9), f"angle {angle}"

    def test_the_element_moves_as_a_rigid_body(self):
        """Placing changes no distance inside the element.

        A pose is a position and an orientation. It must not deform the
        molecule, so every internal distance survives untouched.
        """
        cpx = make_complex()

        cpx.place_global(
            ELEMENT,
            POINT * unit.angstrom,
            np.array([0.3, -0.6, 0.5]) * unit.angstrom,
            2.0,
        )

        after = np.asarray(cpx.positions.value_in_unit(unit.angstrom))
        before = np.linalg.norm(COORDS[0:4, None] - COORDS[None, 0:4], axis=-1)
        assert np.allclose(
            np.linalg.norm(after[0:4, None] - after[None, 0:4], axis=-1),
            before,
            atol=1e-9,
        )

    def test_atoms_outside_the_element_do_not_move(self):
        """The target the aptamer is being docked against stays put."""
        cpx = make_complex()

        cpx.place_global(
            ELEMENT,
            POINT * unit.angstrom,
            np.array([0.0, 0.0, 1.0]) * unit.angstrom,
            1.1,
        )

        after = np.asarray(cpx.positions.value_in_unit(unit.angstrom))
        assert np.array_equal(after[4:], COORDS[4:])

    def test_a_zero_angle_is_a_pure_translation(self):
        """With no turn, every atom of the element shifts by the same vector."""
        cpx = make_complex()

        cpx.place_global(
            ELEMENT,
            POINT * unit.angstrom,
            np.array([0.0, 0.0, 1.0]) * unit.angstrom,
            0.0,
        )

        after = np.asarray(cpx.positions.value_in_unit(unit.angstrom))
        assert np.allclose(after[0:4] - COORDS[0:4], POINT - BUILT_CENTRE, atol=1e-9)

    def test_placing_without_positions_raises(self):
        """A Complex that has not been built cannot be placed."""
        cpx = Complex()

        with pytest.raises(ValueError, match="positions"):
            cpx.place_global(
                ELEMENT,
                POINT * unit.angstrom,
                np.array([0.0, 0.0, 1.0]) * unit.angstrom,
                0.0,
            )


class TestElementCentre:
    """Tests for ``Complex.element_centre``."""

    def test_centre_is_mass_weighted(self):
        """Heavy atoms pull the centre away from the plain average."""
        assert np.allclose(make_complex().element_centre(ELEMENT), BUILT_CENTRE)

    def test_centre_covers_only_the_element(self):
        """Atoms outside the element are ignored."""
        # Atoms 0 and 1, masses 1 and 12: (12/13, 0, 0).
        assert np.allclose(
            make_complex().element_centre([0, 1, 2]), [12 / 13, 0.0, 0.0]
        )

    def test_centre_without_topology_raises(self):
        """Masses come from the topology, so there has to be one."""
        cpx = Complex()
        cpx.positions = [mm.Vec3(*row) for row in COORDS] * unit.angstrom

        with pytest.raises(ValueError, match="topology"):
            cpx.element_centre(ELEMENT)

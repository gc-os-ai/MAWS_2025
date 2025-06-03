"""
Unit tests for complex.py that do NOT spawn tleap / OpenMM.
"""

from unittest.mock import patch, MagicMock

import pytest

import MAWS.src.Structure as structure
from MAWS.src.Complex import Complex


@pytest.fixture
def dummy_structure() -> structure.Structure:
    """A minimal two-residue Structure object."""
    return structure.Structure(["ALA", "GLY"], residue_lengths=[3, 3])


def test_add_chain(dummy_structure):
    cplx = Complex()
    ch = cplx.add_chain("ALA GLY", dummy_structure)

    assert ch.length == 6
    assert len(cplx.chains) == 1
    assert ch.alias_sequence == "ALA GLY"          # no alias decorations


# patch order: Simulation → subprocess → Prmtop → Inpcrd
@patch("MAWS.src.Complex.app.Simulation")
@patch("MAWS.src.Complex.subprocess.run")
@patch("MAWS.src.Complex.app.AmberPrmtopFile")
@patch("MAWS.src.Complex.app.AmberInpcrdFile")
def test_build_mocked(
    mock_inp,        # last patched → first param
    mock_prm,
    mock_sub,
    mock_sim,        # first patched → last param
    dummy_structure,
    tmp_path,
):
    # Fake OpenMM objects
    mock_prm.return_value.topology = MagicMock()
    mock_inp.return_value.positions = [0] * 6

    cplx = Complex()
    cplx.add_chain("ALA GLY", dummy_structure)
    cplx.build("x", workdir=tmp_path)

    # Assertions: our mocks were used exactly once
    mock_sub.assert_called_once()      # tleap invoked
    mock_prm.assert_called_once()
    mock_inp.assert_called_once()
    mock_sim.assert_called_once()      # Simulation ctor patched


def test_rotate_global_pivot(dummy_structure):
    import numpy as np
    from openmm import unit

    cplx = Complex()
    cplx.positions = [
        unit.Quantity(np.array([0.0, 0.0, 0.0]), unit.angstrom),
        unit.Quantity(np.array([1.0, 0.0, 0.0]), unit.angstrom),
        unit.Quantity(np.array([2.0, 0.0, 0.0]), unit.angstrom),
    ]

    # rotate atoms 1..2 around z-axis passing through atom 1
    cplx.rotate_global([0, 1, 3], unit.Quantity(np.array([0.0, 0.0, 1.0]), unit.angstrom), np.pi / 2)

    assert np.allclose(cplx.positions[0].value_in_unit(unit.angstrom), [0.0, 0.0, 0.0])
    assert np.allclose(cplx.positions[1].value_in_unit(unit.angstrom), [1.0, 0.0, 0.0])
    assert not np.allclose(cplx.positions[2].value_in_unit(unit.angstrom), [2.0, 0.0, 0.0])

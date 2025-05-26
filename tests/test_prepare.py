"""
Unit tests for prepare.py that *mock* subprocess calls so we don't need AmberTools
during CI runs.
"""

from pathlib import Path
from unittest.mock import patch, call

import numpy as np
import pytest
from openmm.app import PDBFile

import MAWS.src.Prepare as prepare


@pytest.fixture
def dummy_pdb(tmp_path: Path) -> Path:
    """Write a 3-atom PDB to a temp dir and return the path."""
    text = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.008   1.422   0.000  1.00  0.00           C
TER
END
"""
    path = tmp_path / "dummy.pdb"
    path.write_text(text)
    return path


@patch("MAWS.src.Prepare._run")
def test_make_lib_builds_commands(mock_run, dummy_pdb):
    """make_lib should call antechamber, parmchk, tleap when parameterized=False."""
    atoms = prepare.make_lib(
        dummy_pdb,
        residue_name="DMY",
        charges="bcc",
        atom_type="gaff",
        parameterized=False,
    )

    # Three atoms in the dummy file
    assert atoms == 3

    # antechamber, parmchk, tleap called exactly once each
    commands = [c.args[0] for c in mock_run.call_args_list]
    assert any("antechamber" in str(cmd) for cmd in commands)
    assert any("parmchk" in str(cmd) for cmd in commands)
    assert any("tleap" in str(cmd) for cmd in commands)


@patch("MAWS.src.Prepare._run")
def test_toggle_hydrogens_calls_reduce(mock_run, dummy_pdb):
    prepare.toggle_hydrogens(dummy_pdb, keep=True)
    mock_run.assert_has_calls(
        [
            call(f"reduce -Trim {dummy_pdb} > {dummy_pdb}"),
            call(f"reduce -Build {dummy_pdb} > {dummy_pdb}"),
        ]
    )

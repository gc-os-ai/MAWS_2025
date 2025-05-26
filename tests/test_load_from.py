"""
Tests for load_from.py that do *not* depend on AmberTools or any external files.
"""

from pathlib import Path
from textwrap import dedent
from unittest.mock import patch

import numpy as np
import pytest
from openmm.app import PDBFile

import MAWS.src.LoadFrom as load_from


@pytest.fixture
def sample_xml(tmp_path: Path) -> Path:
    text = dedent(
        """\
        <root>
          <residuePath>./params</residuePath>

          <alias>
            <element name="ALA">
              <alone>ALA</alone>
              <start>ALA-</start>
              <middle>-x-</middle>
              <end>-ALA</end>
            </element>
          </alias>

          <residue name="ALA" length="10">
            <rotation>
              <start>0</start>
              <bond>1</bond>
              <end>2</end>
            </rotation>
            <backbone>
              <start>0</start>
              <middle_pre>1</middle_pre>
              <bond>2</bond>
              <middle_post>3</middle_post>
              <end>4</end>
            </backbone>
            <append>
              <newFirstAtom>0</newFirstAtom>
              <oldLastAtom>-1</oldLastAtom>
              <bondLength>1.6</bondLength>
            </append>
            <prepend>
              <newLastAtom>-1</newLastAtom>
              <oldFirstAtom>0</oldFirstAtom>
              <bondLength>1.6</bondLength>
            </prepend>
          </residue>
        </root>
        """
    )
    p = tmp_path / "test.xml"
    p.write_text(text)
    return p


def test_parse_xml_structure(sample_xml):
    names, lengths, rotations, backbones, connects, rpath, aliases = load_from._parse_xml_structure(
        sample_xml
    )
    assert names == ["ALA"]
    assert lengths == [10]
    assert rotations[0][:3] == ("ALA", 0, 1)
    assert isinstance(rpath, str) and rpath.endswith("params")
    assert aliases[0][0] == "ALA"


@pytest.fixture
def sample_pdb(tmp_path: Path) -> Path:
    pdb = dedent(
        """\
        ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N
        ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C
        ATOM      3  C   GLY A   1       2.008   1.422   0.000  1.00  0.00           C
        TER
        END
        """
    )
    p = tmp_path / "dummy.pdb"
    p.write_text(pdb)
    return p


def test_parse_pdb_structure(sample_pdb):
    names, lengths = load_from._parse_pdb_structure(sample_pdb)
    assert names == ["GLY"]
    assert lengths == [3]


@patch("MAWS.src.LoadFrom.structure.Structure")
def test_xml_structure_builds(mock_structure, sample_xml):
    load_from.xml_structure(sample_xml)
    mock_structure.assert_called_once()


@patch("MAWS.src.LoadFrom.structure.Structure")
def test_pdb_structure_builds(mock_structure, sample_pdb):
    load_from.pdb_structure(sample_pdb)
    mock_structure.assert_called_once()

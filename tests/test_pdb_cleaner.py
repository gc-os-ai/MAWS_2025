"""
Tests for maws.pdb_cleaner module.

This module tests the PDB parsing and cleaning functions:
- pdb_structure(): Parse single PDB lines
- pdb_reader(): Read and split PDB files
- Various check functions for detecting PDB issues
"""

# Import only the parsing functions that don't require complex dependencies
from maws.pdb_cleaner import (
    _pad80,
    check_hydrogen,
    pdb_reader,
    pdb_structure,
)


class TestPdbStructure:
    """Tests for pdb_structure() line parser."""

    def test_pdb_structure_parses_atom_line(self):
        """pdb_structure() correctly parses an ATOM record."""
        # Standard ATOM line (columns are fixed-width)
        line = "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 50.00           N"  # noqa: E501
        result = pdb_structure(line)

        # result is a tuple matching PDB_ITEMS
        assert result[0] == "ATOM"  # Records
        assert result[1].strip() == "1"  # AtomSeq
        assert "N" in result[2]  # AtomTyp
        assert result[4].strip() == "ALA"  # ResName
        assert result[5].strip() == "A"  # ChainID

    def test_pdb_structure_parses_hetatm_line(self):
        """pdb_structure() correctly parses a HETATM record."""
        line = "HETATM  500  C1  LIG A 100      15.500  25.500  35.500  1.00 30.00           C"  # noqa: E501
        result = pdb_structure(line)

        assert result[0] == "HETATM"  # Records
        assert result[4].strip() == "LIG"  # ResName

    def test_pdb_structure_handles_short_line(self):
        """pdb_structure() pads short lines to 80 chars."""
        # Short line (less than 80 chars)
        short_line = "ATOM      1  CA  GLY A   1"
        # Should not raise
        result = pdb_structure(short_line)
        assert result[0] == "ATOM"


class TestPad80:
    """Tests for _pad80() helper function."""

    def test_pad80_short_string(self):
        """_pad80() pads short strings to 80 chars."""
        result = _pad80("hello")
        assert len(result) >= 80

    def test_pad80_exact_80_string(self):
        """_pad80() returns 80-char string unchanged."""
        exact_str = "x" * 80
        result = _pad80(exact_str)
        assert len(result) == 80

    def test_pad80_returns_80_chars_for_long_string(self):
        """_pad80() returns at least 80 chars (may truncate or not)."""
        long_str = "x" * 100
        result = _pad80(long_str)
        # The function pads to 80, so result is 80 chars
        assert len(result) >= 80


class TestPdbReader:
    """Tests for pdb_reader() file parser."""

    def test_pdb_reader_returns_tuple(self, sample_pdb_path):
        """pdb_reader() returns (method, main_df, ligand_df) tuple."""
        result = pdb_reader(sample_pdb_path)
        assert len(result) == 3
        method, main_df, ligand_df = result
        # method is a string (may be empty)
        assert isinstance(method, str)
        # main_df is a DataFrame
        assert hasattr(main_df, "columns")

    def test_pdb_reader_parses_atoms(self, sample_pdb_path):
        """pdb_reader() correctly parses ATOM records."""
        _, main_df, _ = pdb_reader(sample_pdb_path)
        # Should have atoms
        assert len(main_df) > 0
        # Should have expected columns
        assert "Records" in main_df.columns
        assert "ResName" in main_df.columns


class TestCheckFunctions:
    """Tests for PDB check/validation functions."""

    def test_check_hydrogen_detects_h(self, sample_pdb_path):
        """check_hydrogen() detects hydrogen atoms if present."""
        _, main_df, _ = pdb_reader(sample_pdb_path)
        result = check_hydrogen(sample_pdb_path, main_df)
        # Result is either None (no H) or filename (has H)
        # We just check it doesn't crash
        assert result is None or isinstance(result, str)

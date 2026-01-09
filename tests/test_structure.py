"""
Tests for maws.structure module.

This module tests the Structure class which manages residue templates
and topology rules for nucleic acids. Tests cover:
- Structure initialization with residue names/lengths
- Alias translation (sequence alias → canonical residue names)
- Index resolution (positive/negative index normalization)
- Torsion retrieval
- Bond connectivity for append/prepend operations
"""

import pytest

from maws.structure import Structure


class TestStructureInit:
    """Tests for Structure instantiation."""

    def test_structure_with_residue_names(self):
        """Structure can be created with just residue names."""
        struct = Structure(["A", "B", "C"])
        # Structure stores names in residue_names list
        assert "A" in struct.residue_names
        assert "B" in struct.residue_names
        assert "C" in struct.residue_names

    def test_structure_with_lengths(self):
        """Structure stores residue lengths correctly."""
        struct = Structure(["A", "B"], residue_length=[10, 20])
        # Structure stores lengths in residue_length dict
        assert struct.residue_length["A"] == 10
        assert struct.residue_length["B"] == 20

    def test_structure_with_alias(self):
        """Structure stores alias mappings correctly."""
        # Alias format: [name, alone, start, middle, end]
        aliases = [
            ["G", "GN", "G5", "G", "G3"],
            ["A", "AN", "A5", "A", "A3"],
        ]
        struct = Structure(
            ["G", "A", "GN", "AN", "G5", "A5", "G3", "A3"],
            alias=aliases,
        )
        # Alias stored as [alone, start, middle, end]
        assert struct.alias["G"] == ["GN", "G5", "G", "G3"]


class TestTranslate:
    """Tests for Structure.translate() method."""

    def test_translate_single_residue(self):
        """Single residue translates to 'alone' form."""
        struct = Structure(
            ["GN", "G"],
            residue_length=[33, 34],
            alias=[["G", "GN", "G5", "G", "G3"]],
        )
        # Single 'G' should become 'GN' (alone) - returns string
        result = struct.translate("G")
        assert result == "GN"

    def test_translate_multi_residue(self):
        """Multi-residue sequence uses start/middle/end forms."""
        struct = Structure(
            ["GN", "AN", "G", "A", "G5", "A5", "G3", "A3"],
            residue_length=[33, 32, 34, 33, 32, 31, 35, 34],
            alias=[
                ["G", "GN", "G5", "G", "G3"],
                ["A", "AN", "A5", "A", "A3"],
            ],
        )
        # 'G A G' → 'G5 A G3' (start, middle, end) - returns space-separated string
        result = struct.translate("G A G")
        assert result == "G5 A G3"

    def test_translate_two_residues(self):
        """Two-residue sequence uses start/end forms (no middle)."""
        struct = Structure(
            ["GN", "AN", "G", "A", "G5", "A5", "G3", "A3"],
            residue_length=[33, 32, 34, 33, 32, 31, 35, 34],
            alias=[
                ["G", "GN", "G5", "G", "G3"],
                ["A", "AN", "A5", "A", "A3"],
            ],
        )
        # 'G A' → 'G5 A3' (start, end) - returns space-separated string
        result = struct.translate("G A")
        assert result == "G5 A3"


class TestResolveIndex:
    """Tests for Structure.resolve_index() method."""

    @pytest.fixture
    def simple_structure(self):
        """Structure with known residue lengths."""
        return Structure(["RES"], residue_length=[10])

    def test_resolve_positive_index(self, simple_structure):
        """Positive index returns unchanged."""
        result = simple_structure.resolve_index("RES", 3)
        assert result == 3

    def test_resolve_negative_index(self, simple_structure):
        """Negative index resolves from end."""
        # -1 with length 10 → 10 - 1 = 9
        result = simple_structure.resolve_index("RES", -1)
        assert result == 9

    def test_resolve_index_unknown_residue(self, simple_structure):
        """Unknown residue raises ValueError."""
        with pytest.raises(ValueError):
            simple_structure.resolve_index("UNKNOWN", 0)


class TestTorsions:
    """Tests for Structure.torsions() method."""

    def test_torsions_returns_rotation_triples(self):
        """torsions() returns normalized rotation triples."""
        struct = Structure(
            ["RES"],
            residue_length=[20],
            rotating_elements=[
                ("RES", 0, 1, None),  # (start, bond, end)
                ("RES", 5, 10, -5),  # end=-5 → 20-5=15
            ],
        )
        torsions = struct.torsions("RES")
        # Should have 2 torsions
        assert len(torsions) == 2
        # First: (0, 1, None)
        assert torsions[0] == (0, 1, None)
        # Second: end=-5 normalized to 15
        assert torsions[1] == (5, 10, 15)


class TestConnectivity:
    """Tests for append_bond() and prepend_bond() methods."""

    @pytest.fixture
    def structure_with_connect(self):
        """Structure with connectivity info."""
        # connect format: [[append_first, append_last], [prepend_last, prepend_first],
        # append_len, prepend_len]
        return Structure(
            ["RES"],
            residue_length=[20],
            connect=[[[0, -1], [-2, 0], 1.6, 1.6]],
        )

    def test_append_bond_returns_connection_info(self, structure_with_connect):
        """append_bond() returns first atom, last atom, bond length."""
        # With prev_residue_length provided, old_last is resolved
        result = structure_with_connect.append_bond("RES", prev_residue_length=20)
        # Returns (first_atom_idx, last_atom_idx, bond_length)
        assert len(result) == 3
        # first atom = 0, last = -1 → 19 (with prev_length=20), bond_length = 1.6
        assert result[0] == 0
        assert result[1] == 19  # -1 + 20 = 19
        assert result[2] == pytest.approx(1.6)

    def test_append_bond_without_prev_length(self, structure_with_connect):
        """append_bond() returns raw index when prev_residue_length not provided."""
        result = structure_with_connect.append_bond("RES")
        # Without prev_residue_length, old_last stays as -1
        assert result[0] == 0
        assert result[1] == -1  # Raw index
        assert result[2] == pytest.approx(1.6)

    def test_prepend_bond_returns_connection_info(self, structure_with_connect):
        """prepend_bond() returns connection info for prepending."""
        result = structure_with_connect.prepend_bond("RES")
        assert len(result) == 3
        # new_last = -2 → 18 (resolved), old_first = 0, bond_length = 1.6
        assert result[0] == 18
        assert result[1] == 0
        assert result[2] == pytest.approx(1.6)

"""
Tests for maws.rna_structure and maws.dna_structure modules.

These modules provide pre-defined Structure objects for RNA and DNA
nucleic acids. These are UNIT tests - no external dependencies required
(load_rna_structure and load_dna_structure are pure Python).

Tests verify that the loaded structures have the expected residues,
lengths, and mappings.
"""

from maws.dna_structure import (
    RESIDUE_LENGTH as DNA_RESIDUE_LENGTH,
)
from maws.dna_structure import (
    RESIDUE_NAMES as DNA_RESIDUE_NAMES,
)
from maws.dna_structure import (
    load_dna_structure,
)
from maws.rna_structure import (
    RESIDUE_LENGTH as RNA_RESIDUE_LENGTH,
)
from maws.rna_structure import (
    RESIDUE_NAMES as RNA_RESIDUE_NAMES,
)
from maws.rna_structure import (
    load_rna_structure,
)


class TestRNAStructureData:
    """Tests for RNA structure data tables."""

    def test_rna_has_16_residues(self):
        """RNA structure defines 16 residue types."""
        assert len(RNA_RESIDUE_NAMES) == 16

    def test_rna_residue_lengths_match_names(self):
        """RNA residue lengths list matches residue names count."""
        assert len(RNA_RESIDUE_LENGTH) == len(RNA_RESIDUE_NAMES)

    def test_rna_contains_expected_bases(self):
        """RNA structure contains G, A, U, C variants."""
        # Should have standard bases
        assert "G" in RNA_RESIDUE_NAMES
        assert "A" in RNA_RESIDUE_NAMES
        assert "U" in RNA_RESIDUE_NAMES
        assert "C" in RNA_RESIDUE_NAMES
        # Should have 5' and 3' variants
        assert "G5" in RNA_RESIDUE_NAMES
        assert "A3" in RNA_RESIDUE_NAMES


class TestDNAStructureData:
    """Tests for DNA structure data tables."""

    def test_dna_has_16_residues(self):
        """DNA structure defines 16 residue types."""
        assert len(DNA_RESIDUE_NAMES) == 16

    def test_dna_residue_lengths_match_names(self):
        """DNA residue lengths list matches residue names count."""
        assert len(DNA_RESIDUE_LENGTH) == len(DNA_RESIDUE_NAMES)

    def test_dna_contains_expected_bases(self):
        """DNA structure contains DG, DA, DT, DC variants."""
        # Should have standard DNA bases (with D prefix)
        assert "DG" in DNA_RESIDUE_NAMES
        assert "DA" in DNA_RESIDUE_NAMES
        assert "DT" in DNA_RESIDUE_NAMES  # Thymine instead of Uracil
        assert "DC" in DNA_RESIDUE_NAMES


class TestLoadRNAStructure:
    """Tests for load_rna_structure() factory function."""

    def test_load_rna_structure_returns_structure(self):
        """load_rna_structure() returns a Structure object."""
        from maws.structure import Structure

        rna = load_rna_structure()
        assert isinstance(rna, Structure)

    def test_load_rna_structure_has_residue_names(self):
        """Loaded RNA structure has expected residue names."""
        rna = load_rna_structure()
        # Check some residue names exist in the list
        assert "G" in rna.residue_names
        assert "A" in rna.residue_names
        assert "GN" in rna.residue_names

    def test_load_rna_structure_has_lengths(self):
        """Loaded RNA structure has residue lengths defined."""
        rna = load_rna_structure()
        # Standard G has 34 atoms (stored in residue_length dict)
        assert rna.residue_length["G"] == 34

    def test_load_rna_structure_can_translate(self):
        """Loaded RNA structure can translate alias sequences."""
        rna = load_rna_structure()
        # Single G → GN (alone form) - translate returns string
        result = rna.translate("G")
        assert result == "GN"


class TestLoadDNAStructure:
    """Tests for load_dna_structure() factory function."""

    def test_load_dna_structure_returns_structure(self):
        """load_dna_structure() returns a Structure object."""
        from maws.structure import Structure

        dna = load_dna_structure()
        assert isinstance(dna, Structure)

    def test_load_dna_structure_has_residue_names(self):
        """Loaded DNA structure has expected residue names."""
        dna = load_dna_structure()
        # Check residue names in the list
        assert "DG" in dna.residue_names
        assert "DA" in dna.residue_names
        assert "DT" in dna.residue_names

    def test_load_dna_structure_has_lengths(self):
        """Loaded DNA structure has residue lengths defined."""
        dna = load_dna_structure()
        # DG has 33 atoms (stored in residue_length dict)
        assert dna.residue_length["DG"] == 33

    def test_load_dna_structure_can_translate(self):
        """Loaded DNA structure can translate alias sequences."""
        dna = load_dna_structure()
        # Single G → DGN (alone form) - translate returns string
        result = dna.translate("G")
        assert result == "DGN"

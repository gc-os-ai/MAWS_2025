"""
Tests for the maws.run Python API.
"""

import pytest

from maws.run import MAWSConfig, MAWSResult


class TestMAWSConfig:
    """Tests for the MAWSConfig dataclass."""

    def test_config_defaults(self):
        """MAWSConfig has sensible defaults."""
        config = MAWSConfig(pdb_path="test.pdb", num_nucleotides=15)
        assert config.pdb_path == "test.pdb"
        assert config.num_nucleotides == 15
        assert config.name == "MAWS_aptamer"
        assert config.aptamer_type == "RNA"
        assert config.molecule_type == "protein"
        assert config.beta == 0.01
        assert config.first_chunk_size == 5000
        assert config.second_chunk_size == 5000
        assert config.clean_pdb is False
        assert config.verbose is True

    def test_config_custom_values(self):
        """MAWSConfig accepts custom values."""
        config = MAWSConfig(
            pdb_path="data/ligand.pdb",
            num_nucleotides=20,
            name="my_aptamer",
            aptamer_type="DNA",
            molecule_type="organic",
            beta=0.02,
            first_chunk_size=100,
            second_chunk_size=50,
        )
        assert config.num_nucleotides == 20
        assert config.aptamer_type == "DNA"
        assert config.molecule_type == "organic"
        assert config.beta == 0.02

    def test_config_aptamer_types(self):
        """MAWSConfig accepts valid aptamer types."""
        rna_config = MAWSConfig(
            pdb_path="test.pdb", num_nucleotides=10, aptamer_type="RNA"
        )
        dna_config = MAWSConfig(
            pdb_path="test.pdb", num_nucleotides=10, aptamer_type="DNA"
        )
        assert rna_config.aptamer_type == "RNA"
        assert dna_config.aptamer_type == "DNA"

    def test_config_molecule_types(self):
        """MAWSConfig accepts valid molecule types."""
        for mol_type in ["protein", "organic", "lipid"]:
            config = MAWSConfig(
                pdb_path="test.pdb", num_nucleotides=10, molecule_type=mol_type
            )
            assert config.molecule_type == mol_type


class TestMAWSResult:
    """Tests for the MAWSResult dataclass."""

    def test_result_attributes(self):
        """MAWSResult has expected attributes."""
        config = MAWSConfig(pdb_path="test.pdb", num_nucleotides=10)
        # Create a minimal result (without a real Complex)
        result = MAWSResult(
            sequence="G A U C",
            energy=-100.0,
            entropy=0.5,
            complex=None,  # type: ignore
            config=config,
            pdb_path="test_RESULT.pdb",
        )
        assert result.sequence == "G A U C"
        assert result.energy == -100.0
        assert result.entropy == 0.5
        assert result.config is config


class TestRunMAWS:
    """Integration tests for run_maws (require AmberTools/OpenMM)."""

    @pytest.fixture
    def test_pdb_path(self):
        """Path to a test PDB file."""
        import os

        pdb = "data/1BRQ.pdb"
        if os.path.exists(pdb):
            return pdb
        pytest.skip("Test PDB not available")

    @pytest.mark.slow
    def test_run_maws_minimal(self, test_pdb_path):
        """run_maws executes with minimal configuration.

        Note: This test requires GPU/OpenMM and may be slow.
        Run with: pytest -m slow
        """
        pytest.skip("Skipping - requires GPU and is slow")
        from maws.run import run_maws

        config = MAWSConfig(
            pdb_path=test_pdb_path,
            name="test_api",
            num_nucleotides=1,  # Very short for testing
            first_chunk_size=5,
            second_chunk_size=5,
            molecule_type="protein",
            clean_pdb=True,
            remove_h=True,
            drop_hetatm=True,
            verbose=False,
        )

        result = run_maws(config)

        assert result.sequence is not None
        assert len(result.sequence.split()) == 2  # 1 + 1 nucleotides
        assert result.complex is not None
        assert result.pdb_path.endswith("_RESULT.pdb")

"""
Unit tests for maws.chain module.

These tests verify the Chain class's pure Python functionality:
- sequence management (create, append, prepend)
- index and offset tracking
- no external dependencies (AmberTools/OpenMM not required)

Note: Rotation tests require a real Complex with positions, so they are
in test_chain_complex.py as integration tests.
"""

import pytest

from maws.chain import Chain
from maws.structure import Structure


class MockComplex:
    """Minimal mock of Complex for Chain unit tests.

    Chain only uses Complex for:
    - storing reference (self.complex = Complex)
    - iterating chains (Complex.chains) in update_chains()
    - calling rotate_element() which we don't test here
    """

    def __init__(self):
        self.chains: list[Chain] = []

    def rotate_element(self, element, angle, reverse=False):
        """Mock rotation - does nothing in unit tests."""
        pass


@pytest.fixture
def simple_structure():
    """Create a simple Structure for testing Chain."""
    # Define residues with known lengths
    residues = ["A", "B", "C"]
    lengths = [10, 15, 20]
    alias = [
        ["A", "A", "A", "A", "A"],  # [name, alone, start, middle, end]
        ["B", "B", "B", "B", "B"],
        ["C", "C", "C", "C", "C"],
    ]
    return Structure(residues, residue_length=lengths, alias=alias)


@pytest.fixture
def mock_complex():
    """Create a mock Complex for Chain tests."""
    return MockComplex()


class TestChainInit:
    """Tests for Chain initialization."""

    def test_chain_empty_init(self, mock_complex, simple_structure):
        """Chain can be created with empty sequence."""
        chain = Chain(mock_complex, simple_structure, sequence=None, start=0, ID=0)
        assert chain.id == 0
        assert chain.start == 0
        assert chain.length == 0
        assert chain.sequence == ""

    def test_chain_with_sequence(self, mock_complex, simple_structure):
        """Chain computes length and offsets from sequence."""
        chain = Chain(mock_complex, simple_structure, sequence="A B", start=0, ID=0)
        # Length should be A(10) + B(15) = 25
        assert chain.length == 25
        assert chain.alias_sequence == "A B"
        # residues_start[0] = 0, residues_start[1] = 10
        assert chain.residues_start == [0, 10]

    def test_chain_sequence_array(self, mock_complex, simple_structure):
        """Chain splits sequence into arrays."""
        chain = Chain(mock_complex, simple_structure, sequence="A B C", start=0, ID=0)
        assert chain.sequence_array == ["A", "B", "C"]
        assert chain.alias_sequence_array == ["A", "B", "C"]


class TestChainSequenceOperations:
    """Tests for Chain sequence modification methods."""

    def test_create_sequence(self, mock_complex, simple_structure):
        """create_sequence() overwrites entire sequence."""
        mock_complex.chains = []
        chain = Chain(mock_complex, simple_structure, sequence=None, start=0, ID=0)
        mock_complex.chains.append(chain)

        chain.create_sequence("A B")

        assert chain.alias_sequence == "A B"
        assert chain.length == 25  # 10 + 15

    def test_create_sequence_unknown_residue(self, mock_complex, simple_structure):
        """create_sequence() raises for unknown residues."""
        mock_complex.chains = []
        chain = Chain(mock_complex, simple_structure, sequence=None, start=0, ID=0)
        mock_complex.chains.append(chain)

        # Unknown residue raises KeyError (from translate) or ValueError
        with pytest.raises((KeyError, ValueError)):
            chain.create_sequence("UNKNOWN")

    def test_append_sequence(self, mock_complex, simple_structure):
        """append_sequence() adds to right (3') end."""
        mock_complex.chains = []
        chain = Chain(mock_complex, simple_structure, sequence="A", start=0, ID=0)
        mock_complex.chains.append(chain)

        chain.append_sequence("B")

        assert chain.alias_sequence_array == ["A", "B"]
        assert chain.length == 25  # 10 + 15
        assert len(chain.append_history) > 0  # Records appended residues

    def test_prepend_sequence(self, mock_complex, simple_structure):
        """prepend_sequence() adds to left (5') end."""
        mock_complex.chains = []
        chain = Chain(mock_complex, simple_structure, sequence="B", start=0, ID=0)
        mock_complex.chains.append(chain)

        chain.prepend_sequence("A")

        assert chain.alias_sequence_array == ["A", "B"]
        assert chain.length == 25
        assert len(chain.prepend_history) > 0  # Records prepended residues


class TestChainUpdateChains:
    """Tests for Chain.update_chains() method."""

    def test_update_chains_single(self, mock_complex, simple_structure):
        """update_chains() recomputes offsets for single chain."""
        mock_complex.chains = []
        chain = Chain(mock_complex, simple_structure, sequence="A", start=0, ID=0)
        mock_complex.chains.append(chain)

        # Manually change sequence array
        chain.sequence_array = ["A", "B"]
        chain.update_chains()

        assert chain.length == 25
        assert chain.residues_start == [0, 10]

    def test_update_chains_updates_downstream(self, mock_complex, simple_structure):
        """update_chains() shifts downstream chain starts."""
        mock_complex.chains = []
        chain1 = Chain(mock_complex, simple_structure, sequence="A", start=0, ID=0)
        chain2 = Chain(mock_complex, simple_structure, sequence="B", start=10, ID=1)
        mock_complex.chains = [chain1, chain2]

        # Expand chain1
        chain1.sequence_array = ["A", "C"]  # Now length = 10 + 20 = 30
        chain1.update_chains()

        # chain2.start should shift by +20 (chain1 grew by 20)
        assert chain2.start == 30


class TestChainHistory:
    """Tests for Chain history tracking."""

    def test_append_history(self, mock_complex, simple_structure):
        """Appending records history correctly."""
        mock_complex.chains = []
        chain = Chain(mock_complex, simple_structure, sequence="A", start=0, ID=0)
        mock_complex.chains.append(chain)

        old_length = chain.length
        chain.append_sequence("B")

        assert chain.length_history == old_length
        assert "B" in chain.sequence_array

    def test_prepend_history(self, mock_complex, simple_structure):
        """Prepending records history correctly."""
        mock_complex.chains = []
        chain = Chain(mock_complex, simple_structure, sequence="B", start=0, ID=0)
        mock_complex.chains.append(chain)

        old_length = chain.length
        chain.prepend_sequence("A")

        assert chain.length_history == old_length
        assert chain.start_history != chain.start  # Start shifted

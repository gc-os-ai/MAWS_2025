"""
Pytest configuration for MAWS test suite.

This file provides shared fixtures and markers for all tests:
- `integration` marker: Tests requiring OpenMM/AmberTools
- `slow` marker: Tests that take more than a few seconds
- Common fixtures for test data paths
"""

import os

import pytest


def pytest_configure(config):
    """Register custom markers to avoid warnings."""
    config.addinivalue_line(
        "markers", "integration: marks tests requiring OpenMM/AmberTools"
    )
    config.addinivalue_line("markers", "slow: marks tests as slow-running")


# ---------------------------------------------------------------------------
# Shared Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def data_dir():
    """Path to the test data directory."""
    return os.path.join(os.path.dirname(__file__), "..", "data")


@pytest.fixture
def sample_pdb_path(data_dir):
    """Path to the 1BRQ.pdb test file."""
    pdb_path = os.path.join(data_dir, "1BRQ.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip("Test PDB file not found: data/1BRQ.pdb")
    return pdb_path


@pytest.fixture
def cleaned_pdb_path(data_dir):
    """Path to the pre-cleaned 1BRQ PDB file."""
    pdb_path = os.path.join(data_dir, "1BRQ_cleaned.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip("Cleaned PDB file not found: data/1BRQ_cleaned.pdb")
    return pdb_path


# ---------------------------------------------------------------------------
# Synthetic Complex-like fixtures (no OpenMM build) for surface-sampler tests
# ---------------------------------------------------------------------------

import numpy as np  # noqa: E402
from openmm import unit  # noqa: E402


class _SyntheticAtom:
    def __init__(self, symbol):
        self.element = type("E", (), {"symbol": symbol, "mass": 12.0 * unit.dalton})()


class _SyntheticTopology:
    def __init__(self, symbols):
        self._atoms = [_SyntheticAtom(s) for s in symbols]

    def atoms(self):
        # Match OpenMM Topology.atoms(), which is a generator method.
        return iter(self._atoms)


class SyntheticComplex:
    """Duck-typed Complex for surface-sampling tests (no OpenMM build)."""

    def __init__(self, positions_angstrom, symbols):
        self.positions = np.asarray(positions_angstrom, dtype=float) * unit.angstrom
        self.topology = _SyntheticTopology(symbols)


@pytest.fixture
def synthetic_two_carbon_complex():
    """Two C atoms 10 Å apart along the x-axis."""
    return SyntheticComplex(
        positions_angstrom=[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
        symbols=["C", "C"],
    )


@pytest.fixture
def synthetic_octahedron_complex():
    """Six C atoms in a unit-axis octahedron at distance 5 Å from origin."""
    return SyntheticComplex(
        positions_angstrom=[
            [5.0, 0.0, 0.0],
            [-5.0, 0.0, 0.0],
            [0.0, 5.0, 0.0],
            [0.0, -5.0, 0.0],
            [0.0, 0.0, 5.0],
            [0.0, 0.0, -5.0],
        ],
        symbols=["C"] * 6,
    )

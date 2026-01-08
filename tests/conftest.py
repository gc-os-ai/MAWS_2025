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

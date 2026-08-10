"""
Shared fixtures for the MAWS test suite.

The suite is in three tiers, one directory each:

``tests/unit``
    Pure functions and value types. Nothing is built, nothing is scored, and
    no random numbers are drawn. Every assertion can be checked by hand.
``tests/algorithm``
    The search, the samplers, the reporters and the command-line program, run
    end to end against the stand-in builder and the stand-in scorer. No
    AmberTools, no OpenMM, milliseconds per test.
``tests/integration``
    Anything needing AmberTools or OpenMM installed. Marked ``integration``
    and skipped by default; run them with ``pytest -m integration``.

The fixtures below serve the middle tier. They build a small design out of
:class:`~maws.build.FakeBuilder` grids, which have the right chains, the right
atom counts and real turnable bonds, but no chemistry.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from maws.build import FakeBuilder
from maws.forcefield import ForceField
from maws.libraries import rna
from maws.topology import Assembly, BuiltSystem

SMALL_TARGET_ATOMS = 20
"""Atom count of the stand-in target used throughout the suite.

Large enough that its grid is genuinely three-dimensional, small enough that a
distance between any two of its atoms can be worked out by hand.
"""


@pytest.fixture
def rng() -> np.random.Generator:
    """A random generator with a fixed seed, so every test repeats exactly."""
    return np.random.default_rng(20250810)


@pytest.fixture
def data_dir() -> Path:
    """The repository's ``data`` directory, which holds example structures."""
    return Path(__file__).resolve().parent.parent / "data"


@pytest.fixture
def sample_pdb_path(data_dir: Path) -> Path:
    """A real downloaded protein structure, skipping the test if it is absent."""
    path = data_dir / "1BRQ.pdb"
    if not path.exists():
        pytest.skip(f"example structure not found: {path}")
    return path


@pytest.fixture
def forcefield() -> ForceField:
    """The usual RNA-against-protein setup."""
    return ForceField.for_target("RNA", "protein")


@pytest.fixture
def builder() -> FakeBuilder:
    """A builder that puts atoms on a grid instead of running AmberTools."""
    return FakeBuilder()


@pytest.fixture
def empty_system(builder: FakeBuilder, forcefield: ForceField) -> BuiltSystem:
    """A target with no aptamer yet, which is what a search starts from."""
    return builder.build(
        Assembly().with_aptamer(rna()).with_ligand_stub(SMALL_TARGET_ATOMS),
        forcefield,
    )


@pytest.fixture
def one_residue_system(builder: FakeBuilder, forcefield: ForceField) -> BuiltSystem:
    """A target with a single-guanine strand beside it."""
    return builder.build(
        Assembly().with_aptamer(rna(), "G").with_ligand_stub(SMALL_TARGET_ATOMS),
        forcefield,
    )


@pytest.fixture
def two_residue_system(builder: FakeBuilder, forcefield: ForceField) -> BuiltSystem:
    """A target with a two-nucleotide strand beside it."""
    return builder.build(
        Assembly().with_aptamer(rna(), "G A").with_ligand_stub(SMALL_TARGET_ATOMS),
        forcefield,
    )

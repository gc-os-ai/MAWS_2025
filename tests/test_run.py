"""
Integration test for maws.run.

These tests require AmberTools/OpenMM and external binaries (e.g., tleap).
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import pytest

from maws.run import MawsResult, MawsRunner

ShapeT = Literal["cube", "sphere", "shell"]


@pytest.mark.integration
def test_maws_runner_smoke(tmp_path: Path) -> None:
    pdb = Path("data/1BRQ.pdb")
    if not pdb.exists():
        pytest.skip("Test PDB not available (data/1BRQ.pdb).")

    runner = MawsRunner(
        num_nucleotides=1,
        aptamer_type="RNA",
        molecule_type="protein",
        first_chunck_size=2,
        second_chunck_size=2,
        clean_pdb=True,
        remove_h=True,
        drop_hetatm=False,
        verbose=False,
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    result = runner.run(pdb=pdb, name="test_api", output_pdb=out_dir)

    assert isinstance(result, MawsResult)
    assert isinstance(result.sequence, str)
    assert result.sequence.strip() != ""

    # not NaN
    assert result.energy == result.energy
    assert result.entropy == result.entropy

    assert result.pdb_path is not None
    assert result.pdb_path.endswith("test_api_RESULT.pdb")
    assert Path(result.pdb_path).exists()
    assert Path(result.pdb_path).stat().st_size > 0


@pytest.mark.integration
@pytest.mark.parametrize("shape", ["cube", "sphere", "shell"])
def test_maws_runner_each_shape(tmp_path: Path, shape: str) -> None:
    """End-to-end: each surface-aware envelope shape produces a result."""
    pdb = Path("data/1BRQ.pdb")
    if not pdb.exists():
        pytest.skip("Test PDB not available (data/1BRQ.pdb).")

    runner = MawsRunner(
        num_nucleotides=1,
        aptamer_type="RNA",
        molecule_type="protein",
        first_chunck_size=2,
        second_chunck_size=2,
        clean_pdb=True,
        remove_h=True,
        drop_hetatm=False,
        verbose=False,
        shape=cast(ShapeT, shape),
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    result = runner.run(pdb=pdb, name=f"test_{shape}", output_pdb=out_dir)
    assert isinstance(result, MawsResult)
    assert isinstance(result.sequence, str) and result.sequence.strip() != ""


def test_runner_rejects_negative_reach() -> None:
    """MawsRunner raises ValueError on negative reach (no integration setup needed)."""
    with pytest.raises(ValueError, match="reach must be >= 0"):
        MawsRunner(
            num_nucleotides=1,
            aptamer_type="RNA",
            molecule_type="protein",
            reach=-1.0,
        )


def test_runner_rejects_negative_probe() -> None:
    """MawsRunner raises ValueError on negative probe (no integration setup needed)."""
    with pytest.raises(ValueError, match="probe must be >= 0"):
        MawsRunner(
            num_nucleotides=1,
            aptamer_type="RNA",
            molecule_type="protein",
            probe=-1.0,
        )

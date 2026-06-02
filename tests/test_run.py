"""
Integration test for maws.run.

These tests require AmberTools/OpenMM and external binaries (e.g., tleap).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from maws.run import MawsResult, MawsRunner


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


def test_runner_rejects_negative_salt_conc() -> None:
    """MawsRunner raises ValueError on negative salt_conc (no integration setup)."""
    with pytest.raises(ValueError, match="salt_conc must be >= 0"):
        MawsRunner(
            num_nucleotides=1,
            aptamer_type="RNA",
            molecule_type="protein",
            salt_conc=-0.1,
        )


def test_runner_default_salt_conc() -> None:
    """MawsRunner defaults salt_conc to physiological 0.15 mol/L."""
    runner = MawsRunner(
        num_nucleotides=1,
        aptamer_type="RNA",
        molecule_type="protein",
    )
    assert runner.salt_conc == 0.15


def test_runner_accepts_zero_salt_conc() -> None:
    """salt_conc=0.0 (documented unscreened mode) is accepted, not rejected."""
    runner = MawsRunner(
        num_nucleotides=1,
        aptamer_type="RNA",
        molecule_type="protein",
        salt_conc=0.0,
    )
    assert runner.salt_conc == 0.0


def test_run_threads_salt_conc_into_every_complex(monkeypatch) -> None:
    """run() passes salt_conc into every Complex construction (wiring coverage).

    Locks the end-to-end seam without LEaP/OpenMM: a recording double captures
    constructor kwargs and a sentinel halts run() right after the builds.
    """
    import maws.run as run_mod

    constructed: list[dict] = []

    class _RecordingComplex:
        def __init__(self, **kwargs):
            constructed.append(kwargs)

        def add_chain(self, *args, **kwargs):
            pass

        def add_chain_from_pdb(self, *args, **kwargs):
            pass

        def build(self, *args, **kwargs):
            pass

    class _StopHereError(Exception):
        pass

    def _stop(*args, **kwargs):
        raise _StopHereError

    monkeypatch.setattr(run_mod, "Complex", _RecordingComplex)
    monkeypatch.setattr(
        run_mod, "resolve_pdb_path", lambda *a, **k: ("lig.pdb", "lig.pdb")
    )
    monkeypatch.setattr(run_mod.space, "make_sampler", _stop)

    runner = MawsRunner(
        num_nucleotides=1,
        aptamer_type="RNA",
        molecule_type="protein",
        salt_conc=0.3,
    )

    with pytest.raises(_StopHereError):
        runner.run(pdb="lig.pdb")

    # Both the template complex and the ligand-only complex must be built.
    assert len(constructed) == 2
    assert all(kwargs.get("salt_conc") == 0.3 for kwargs in constructed)

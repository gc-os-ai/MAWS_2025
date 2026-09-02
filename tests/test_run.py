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
        first_chunk_size=2,
        second_chunk_size=2,
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


class TestSelectBeam:
    """Tests for select_beam, which decides what survives a search step.

    The beam is what makes a wrong early choice recoverable. Keeping one
    candidate is the greedy search EFBA specifies; keeping more is the
    departure described in `MawsRunner`.
    """

    @staticmethod
    def _candidate(entropy, sequence, total=None):
        from maws.run import Candidate

        return Candidate(
            entropy=entropy,
            total=entropy if total is None else total,
            energy=0.0,
            sequence=sequence,
            positions=[],
        )

    def test_candidates_are_ranked_by_the_running_total(self) -> None:
        """A lucky step does not let a weak lineage displace a strong one.

        EFBA's entropy is extensive, so an aptamer's score is the sum over
        every nucleotide in it. Beam members carry different histories, so
        ranking their children on the current step alone would weigh a
        strong lineage against a weak one as though they were equal.
        """
        from maws.run import select_beam

        lucky_step = self._candidate(-0.9, "AA", total=-1.0)
        strong_line = self._candidate(-0.2, "GG", total=-2.0)
        beam = select_beam([lucky_step, strong_line], width=1)
        assert [c.sequence for c in beam] == ["GG"]

    def test_the_lowest_scoring_candidate_comes_first(self) -> None:
        """Candidates are ordered by score, lowest first."""
        from maws.run import select_beam

        beam = select_beam(
            [
                self._candidate(-0.2, "A"),
                self._candidate(-0.9, "G"),
                self._candidate(-0.5, "C"),
            ],
            width=3,
        )
        assert [c.sequence for c in beam] == ["G", "C", "A"]

    def test_a_width_of_one_keeps_only_the_winner(self) -> None:
        """Width 1 is the greedy search EFBA specifies."""
        from maws.run import select_beam

        beam = select_beam(
            [self._candidate(-0.2, "A"), self._candidate(-0.9, "G")], width=1
        )
        assert [c.sequence for c in beam] == ["G"]

    def test_a_width_wider_than_the_field_keeps_everything(self) -> None:
        """Asking for more candidates than exist returns the ones that do."""
        from maws.run import select_beam

        beam = select_beam([self._candidate(-0.2, "A")], width=5)
        assert len(beam) == 1

    def test_an_empty_field_gives_an_empty_beam(self) -> None:
        """A step where every candidate was unviable returns nothing.

        `draw_clear_torsions` raises when a nucleotide cannot be bent clear
        of the target, which happens when the strand's growing end is buried.
        The search skips that candidate rather than dying, so a step can end
        with fewer candidates than it started with, or none at all.
        """
        from maws.run import select_beam

        assert select_beam([], width=3) == []

    def test_tied_scores_do_not_compare_the_rest_of_the_candidate(self) -> None:
        """A tie is broken without touching `positions`.

        Positions are arrays. Ordering two candidates by comparing whole
        records would reach them on a tie, and comparing arrays raises
        rather than returning an order.
        """
        import numpy as np

        from maws.run import Candidate, select_beam

        tied = [
            Candidate(-0.5, -0.5, 0.0, "G", np.zeros((3, 3))),
            Candidate(-0.5, -0.5, 0.0, "A", np.ones((3, 3))),
        ]
        assert len(select_beam(tied, width=2)) == 2


def test_runner_defaults_to_the_greedy_search_efba_specifies() -> None:
    """`beam` defaults to 1, so a default run follows the published method."""
    runner = MawsRunner(num_nucleotides=1, aptamer_type="RNA", molecule_type="protein")
    assert runner.beam == 1


@pytest.mark.parametrize("beam", [0, -1])
def test_runner_rejects_a_beam_narrower_than_one(beam: int) -> None:
    """A search has to carry at least one candidate between steps."""
    with pytest.raises(ValueError, match="beam must be >= 1"):
        MawsRunner(
            num_nucleotides=1,
            aptamer_type="RNA",
            molecule_type="protein",
            beam=beam,
        )


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


@pytest.mark.parametrize("num_nucleotides", [0, -1])
def test_runner_rejects_non_positive_num_nucleotides(num_nucleotides: int) -> None:
    """MawsRunner rejects num_nucleotides <= 0 with a message matching the check."""
    with pytest.raises(
        ValueError,
        match=f"num_nucleotides must be greater than 0, got {num_nucleotides}",
    ):
        MawsRunner(
            num_nucleotides=num_nucleotides,
            aptamer_type="RNA",
            molecule_type="protein",
        )


@pytest.mark.parametrize(
    ("first", "second"), [(0, 5000), (5000, 0), (-1, 5000), (5000, -1)]
)
def test_runner_rejects_non_positive_chunk_size(first: int, second: int) -> None:
    """MawsRunner raises ValueError on a chunk size that is not positive."""
    with pytest.raises(ValueError, match="Chunk size must be greater than 0"):
        MawsRunner(
            num_nucleotides=1,
            aptamer_type="RNA",
            molecule_type="protein",
            first_chunk_size=first,
            second_chunk_size=second,
        )


def test_runner_chunk_size_kwargs_reach_attributes() -> None:
    """Chunk-size kwargs are spelled 'chunk' and land on the matching attributes."""
    runner = MawsRunner(
        num_nucleotides=1,
        aptamer_type="RNA",
        molecule_type="protein",
        first_chunk_size=7,
        second_chunk_size=11,
    )
    assert runner.first_chunk_size == 7
    assert runner.second_chunk_size == 11


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

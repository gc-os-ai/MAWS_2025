"""
Tests for the three ways into MAWS, and for what they report.

:func:`maws.api.design`, :class:`maws.AptamerDesigner` and the ``maws`` command
are three faces of one run. These tests check they agree with each other and
with :func:`maws.search.grow_aptamer`, and that the reporters turn events into
the right messages and files.

Nothing here builds a real structure. The design paths that would need
AmberTools are exercised through :func:`maws.api.collect`, which takes an
already-running search, so the same result-shaping code is covered without the
external programs.
"""

from __future__ import annotations

import io
import json
import logging
from pathlib import Path

import numpy as np
import pytest

from maws import __version__
from maws._designer import AptamerDesigner, _BaseEstimator
from maws.api import MawsResult, collect
from maws.cli import build_parser, load_config, main, resolve_settings
from maws.errors import ConfigurationError
from maws.reporting import JsonlReporter, LogReporter, ProgressReporter
from maws.search import (
    CandidateScored,
    SearchFinished,
    StepCompleted,
    StepStarted,
    grow_aptamer,
    stub_energy,
)
from maws.values import NucleotideSequence


@pytest.fixture
def events(empty_system, builder):
    """A short but complete run of the search, as a list of events."""
    return list(
        grow_aptamer(
            empty_system,
            energy=stub_energy(),
            builder=builder,
            n_nucleotides=2,
            first_samples=4,
            samples=4,
            rng=np.random.default_rng(0),
        )
    )


class TestMawsResult:
    """The record a finished run hands back."""

    def test_length_counts_the_nucleotides_in_the_strand(self):
        """A strand is written with spaces between its nucleotides."""
        assert MawsResult("G A U", -1421.9, -0.072, steps=3).length == 3

    def test_printing_a_result_shows_the_strand_and_both_numbers(self):
        """One line carries everything needed to compare two runs."""
        text = str(MawsResult("G A U", -1421.9, -0.072, steps=3))
        assert "G A U" in text
        assert "-1421.90" in text
        assert "-0.072000" in text

    def test_a_finished_run_is_marked_successful(self):
        """Success is the default, so only failures have to say so."""
        assert MawsResult("G", 0.0, 0.0, steps=1).success

    def test_a_result_cannot_be_changed_after_it_is_made(self):
        """Results are records of what happened, not working state."""
        result = MawsResult("G", 0.0, 0.0, steps=1)
        with pytest.raises(AttributeError):
            result.sequence = "A"


class TestCollect:
    """Reducing a run's events to a single answer."""

    def test_the_answer_comes_from_the_final_event(self, events):
        """The strand reported is the one the search settled on."""
        result = collect(iter(events))
        assert result.sequence == str(events[-1].winner.sequence)

    def test_the_score_and_energy_come_from_the_same_candidate(self, events):
        """Both numbers describe the winning strand, not different ones."""
        result = collect(iter(events))
        assert result.energy == events[-1].winner.energy
        assert result.entropy == events[-1].winner.entropy

    def test_the_final_structure_is_kept(self, events):
        """The result carries enough to write the strand out afterwards."""
        result = collect(iter(events))
        assert result.system is not None
        assert result.pose.n_atoms == result.system.n_atoms

    def test_every_event_is_offered_to_the_reporter(self, events):
        """A reporter sees the whole run, not only its conclusion."""
        seen = []
        collect(iter(events), on_event=seen.append)
        assert len(seen) == len(events)

    def test_a_run_that_stops_early_is_marked_unsuccessful(self, events):
        """A partial run must be distinguishable from a finished one."""
        truncated = [e for e in events if not isinstance(e, SearchFinished)]
        result = collect(iter(truncated))
        assert not result.success
        assert "stopped before it finished" in result.message


class TestProgressReporter:
    """The reporter that prints one line per step."""

    def test_the_first_step_says_the_strand_is_empty(self):
        """There is nothing to grow from until a nucleotide is chosen."""
        stream = io.StringIO()
        ProgressReporter(stream)(
            StepStarted(step=1, total=3, sequence=NucleotideSequence(()))
        )
        assert "empty strand" in stream.getvalue()

    def test_a_later_step_says_what_it_is_growing_from(self):
        """Progress is only readable if it names the strand so far."""
        stream = io.StringIO()
        ProgressReporter(stream)(
            StepStarted(step=2, total=3, sequence=NucleotideSequence(("G", "A")))
        )
        assert "Step 2/3: growing from G A" in stream.getvalue()

    def test_scored_candidates_are_not_printed(self, events):
        """This reporter is for watching progress, not for the full record."""
        stream = io.StringIO()
        reporter = ProgressReporter(stream)
        for event in events:
            reporter(event)
        printed = stream.getvalue().splitlines()
        assert len(printed) == 3  # two steps started, one search finished

    def test_the_end_of_the_run_reports_the_strand_and_its_numbers(self, events):
        """The last line is the answer, so it carries everything."""
        stream = io.StringIO()
        ProgressReporter(stream)(events[-1])
        assert "Done after 2 steps" in stream.getvalue()


class TestLogReporter:
    """The reporter that writes a record of the run to disk."""

    def test_it_writes_both_of_its_files(self, tmp_path, events):
        """A run leaves behind a score record and a structure record."""
        with LogReporter("demo", tmp_path) as reporter:
            for event in events:
                reporter(event)
        names = sorted(path.name for path in tmp_path.iterdir())
        assert names == ["demo_scores.log", "demo_steps.pdb"]

    def test_the_score_file_holds_one_line_per_candidate(self, tmp_path, events):
        """The record explains how the search chose what it chose."""
        with LogReporter("demo", tmp_path) as reporter:
            for event in events:
                reporter(event)
        lines = (tmp_path / "demo_scores.log").read_text().splitlines()
        scored = sum(isinstance(e, CandidateScored) for e in events)
        assert len(lines) == scored

    def test_each_score_line_names_the_strand_and_both_numbers(self, tmp_path, events):
        """A line of the record is readable on its own."""
        with LogReporter("demo", tmp_path) as reporter:
            for event in events:
                reporter(event)
        first = (tmp_path / "demo_scores.log").read_text().splitlines()[0]
        assert "SEQUENCE:" in first
        assert "SCORE:" in first
        assert "ENERGY:" in first

    def test_it_creates_the_directory_it_writes_into(self, tmp_path):
        """Naming an output directory is enough; it need not exist."""
        target = tmp_path / "runs" / "today"
        with LogReporter("demo", target):
            pass
        assert target.is_dir()

    def test_closing_twice_is_harmless(self, tmp_path):
        """Leaving the context after an explicit close must not fail."""
        reporter = LogReporter("demo", tmp_path)
        reporter.close()
        reporter.close()

    def test_progress_goes_to_the_logger_it_was_given(self, tmp_path, caplog):
        """Where messages go is the caller's choice, not this class's."""
        logger = logging.getLogger("maws.test.reporter")
        with caplog.at_level(logging.INFO, logger="maws.test.reporter"):
            with LogReporter("demo", tmp_path, logger=logger) as reporter:
                reporter(StepStarted(step=1, total=1, sequence=NucleotideSequence(())))
        assert "Step 1/1" in caplog.text


class TestJsonlReporter:
    """The reporter that writes one JSON object per event."""

    def test_every_event_becomes_one_line(self, events):
        """One line per event is what makes the output streamable."""
        stream = io.StringIO()
        reporter = JsonlReporter(stream)
        for event in events:
            reporter(event)
        assert len(stream.getvalue().splitlines()) == len(events)

    def test_every_line_is_valid_json_naming_its_event(self, events):
        """Another program has to be able to tell the events apart."""
        stream = io.StringIO()
        reporter = JsonlReporter(stream)
        for event in events:
            reporter(event)
        kinds = [json.loads(line)["event"] for line in stream.getvalue().splitlines()]
        assert kinds[0] == "StepStarted"
        assert kinds[-1] == "SearchFinished"

    def test_a_scored_candidate_carries_its_numbers(self, events):
        """The machine-readable record has to hold the scores too."""
        stream = io.StringIO()
        reporter = JsonlReporter(stream)
        for event in events:
            if isinstance(event, CandidateScored):
                reporter(event)
                break
        record = json.loads(stream.getvalue())
        assert {"sequence", "token", "direction", "score", "energy"} <= set(record)

    def test_structures_are_left_out(self, events):
        """Atom positions do not belong in JSON; the PDB file holds those."""
        stream = io.StringIO()
        completed = next(e for e in events if isinstance(e, StepCompleted))
        JsonlReporter(stream)(completed)
        assert "pose" not in json.loads(stream.getvalue())


class TestEstimatorConventions:
    """The settings interface :class:`AptamerDesigner` follows."""

    def test_building_one_does_no_work(self):
        """Nothing is read, run or validated until fit is called."""
        AptamerDesigner(n_nucleotides=-5, molecule="nonsense")

    def test_settings_are_stored_exactly_as_given(self):
        """Each argument is kept under its own name and left alone."""
        designer = AptamerDesigner(n_nucleotides=7, n_samples=42)
        assert designer.n_nucleotides == 7
        assert designer.n_samples == 42

    def test_get_params_lists_every_setting(self):
        """The settings can be read back out, which is what makes copying work."""
        params = AptamerDesigner().get_params()
        assert params["n_nucleotides"] == 15
        assert "random_state" in params

    def test_an_object_can_be_rebuilt_from_its_own_settings(self):
        """Taking one apart and putting it back gives an equal configuration."""
        original = AptamerDesigner(n_nucleotides=8, aptamer="DNA")
        copy = AptamerDesigner(**original.get_params())
        assert copy.get_params() == original.get_params()

    def test_set_params_changes_one_setting_and_returns_the_object(self):
        """Chaining a change onto a copy is the usual way to vary a run."""
        designer = AptamerDesigner()
        assert designer.set_params(n_nucleotides=3) is designer
        assert designer.n_nucleotides == 3

    def test_an_unknown_setting_is_rejected_by_name(self):
        """A typo in a setting name is caught, and the real names are listed."""
        with pytest.raises(ConfigurationError, match="has no setting 'ntides'"):
            AptamerDesigner().set_params(ntides=5)

    def test_the_repr_shows_the_settings(self):
        """Printing a designer shows how it was configured."""
        assert "n_nucleotides=3" in repr(AptamerDesigner(n_nucleotides=3))

    def test_results_are_only_available_after_fitting(self):
        """Asking for a score before running says so plainly."""
        with pytest.raises(ConfigurationError, match="not been fitted"):
            AptamerDesigner().score()

    def test_settings_are_discovered_from_the_constructor(self):
        """The list of settings is derived, so it cannot fall out of step."""
        assert "n_samples" in _BaseEstimator._param_names.__func__(AptamerDesigner)

    @pytest.mark.parametrize(
        ("setting", "value", "message"),
        [
            ("n_nucleotides", 0, "at least 1"),
            ("n_samples", 0, "at least 1"),
            ("beta", -1.0, "must not be negative"),
        ],
    )
    def test_settings_are_checked_when_the_run_starts(self, setting, value, message):
        """Checking happens at fit, which is where the settings are first used."""
        designer = AptamerDesigner(**{setting: value})
        with pytest.raises(ConfigurationError, match=message):
            designer.fit(["target.pdb"])

    def test_fitting_without_a_target_is_rejected(self):
        """There is nothing to design against."""
        with pytest.raises(ConfigurationError, match="at least one target"):
            AptamerDesigner().fit([])


class TestCommandLineParsing:
    """Turning command-line arguments into settings."""

    def test_the_design_subcommand_accepts_its_main_flags(self):
        """The flags a user reaches for first are all there."""
        args = build_parser().parse_args(
            ["design", "--target", "t.pdb", "--length", "4", "--aptamer", "DNA"]
        )
        assert (args.target, args.length, args.aptamer) == ("t.pdb", 4, "DNA")

    @pytest.mark.parametrize("command", ["design", "prepare", "inspect", "clean"])
    def test_every_subcommand_is_reachable(self, command):
        """All four are wired up and give a --help without failing."""
        with pytest.raises(SystemExit) as exit_info:
            build_parser().parse_args([command, "--help"])
        assert exit_info.value.code == 0

    def test_a_subcommand_is_required(self):
        """Running the bare command is an error, not a silent no-op."""
        with pytest.raises(SystemExit):
            build_parser().parse_args([])

    def test_verbose_and_quiet_cannot_both_be_given(self):
        """The two contradict each other, so argparse refuses both."""
        with pytest.raises(SystemExit):
            build_parser().parse_args(["design", "--target", "t.pdb", "-v", "-q"])

    def test_the_version_flag_reports_the_package_version(self, capsys):
        """``maws --version`` prints the same version the package reports."""
        with pytest.raises(SystemExit):
            main(["--version"])
        assert __version__ in capsys.readouterr().out


class TestConfigFile:
    """Reading settings from a file instead of typing them all."""

    def _write(self, tmp_path: Path, text: str) -> Path:
        """Write a config file and return its path.

        Parameters
        ----------
        tmp_path : pathlib.Path
            Directory to write into.
        text : str
            The file's contents.

        Returns
        -------
        pathlib.Path
            Where it was written.
        """
        path = tmp_path / "maws.toml"
        path.write_text(text)
        return path

    def test_settings_are_read_from_every_section(self, tmp_path):
        """Sections group settings for readability and are then merged."""
        path = self._write(
            tmp_path,
            "[design]\nlength = 20\n[sampling]\nseed = 42\n[physics]\n"
            "salt_conc = 0.3\n",
        )
        settings = load_config(path)
        assert settings == {"length": 20, "seed": 42, "salt_conc": 0.3}

    def test_hyphens_in_names_become_underscores(self, tmp_path):
        """A file may spell a setting the way the flag does."""
        path = self._write(tmp_path, "[design]\nclean-pdb = true\n")
        assert load_config(path) == {"clean_pdb": True}

    def test_unknown_sections_are_ignored(self, tmp_path):
        """A file may carry notes or settings meant for something else."""
        path = self._write(tmp_path, "[notes]\nauthor = 'me'\n[design]\nlength = 5\n")
        assert load_config(path) == {"length": 5}

    def test_a_missing_file_is_reported_by_name(self, tmp_path):
        """A mistyped path should not look like an empty configuration."""
        with pytest.raises(ConfigurationError, match="cannot read config file"):
            load_config(tmp_path / "absent.toml")

    def test_a_malformed_file_is_reported_as_such(self, tmp_path):
        """Broken TOML is a different problem from a missing file."""
        path = self._write(tmp_path, "this is not toml [[[")
        with pytest.raises(ConfigurationError, match="not valid TOML"):
            load_config(path)

    def test_a_flag_typed_on_the_command_line_wins(self, tmp_path):
        """The file supplies defaults; what you type overrides them."""
        path = self._write(tmp_path, "[design]\nlength = 20\n")
        args = build_parser().parse_args(
            ["design", "--target", "t.pdb", "--length", "25", "--config", str(path)]
        )
        settled = resolve_settings(args, args.design_defaults)
        assert settled["length"] == 25

    def test_the_file_fills_in_what_was_not_typed(self, tmp_path):
        """That is the point of a config file: fewer flags to type."""
        path = self._write(tmp_path, "[design]\nlength = 20\n")
        args = build_parser().parse_args(
            ["design", "--target", "t.pdb", "--config", str(path)]
        )
        settled = resolve_settings(args, args.design_defaults)
        assert settled["length"] == 20

    def test_settings_fall_back_to_the_built_in_defaults(self):
        """With no file and no flag, the documented default applies."""
        args = build_parser().parse_args(["design", "--target", "t.pdb"])
        settled = resolve_settings(args, args.design_defaults)
        assert settled["length"] == 15

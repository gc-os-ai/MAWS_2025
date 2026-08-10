"""
maws.reporting
==============

Turning a search's progress into messages and files.

:func:`maws.search.grow_aptamer` reports what it is doing as a stream of
events. It does not print anything and does not write anything: what to do with
each event is left to whoever is reading them. This module holds the readers.

:class:`LogReporter`
    Prints progress and writes two files: a record of every candidate and its
    score, and the best structure from each step.
:class:`ProgressReporter`
    Prints one line per step and nothing else.
:class:`JsonlReporter`
    Writes one JSON object per event, for feeding into another program.

Each is a callable taking one event, so any of them can be passed straight to
:func:`maws.api.design` as ``on_event``, and writing another takes a few lines.

Examples
--------
>>> from maws.search import StepStarted
>>> from maws.values import NucleotideSequence
>>> reporter = ProgressReporter()
>>> reporter(StepStarted(step=1, total=3, sequence=NucleotideSequence(())))
Step 1/3: starting from an empty strand
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from types import TracebackType
from typing import IO, Any

from maws.search import (
    CandidateScored,
    SearchFinished,
    StepCompleted,
    StepEvent,
    StepStarted,
)

__all__ = ["JsonlReporter", "LogReporter", "ProgressReporter"]


class ProgressReporter:
    """Prints one line per step, and one at the end.

    Parameters
    ----------
    stream : io.TextIOBase, optional
        Where to print. Defaults to standard output.

    See Also
    --------
    LogReporter : Also writes files and reports every candidate.

    Examples
    --------
    >>> from maws.search import StepStarted
    >>> from maws.values import NucleotideSequence
    >>> ProgressReporter()(
    ...     StepStarted(step=2, total=5, sequence=NucleotideSequence(("G", "A")))
    ... )
    Step 2/5: growing from G A
    """

    __slots__ = ("_stream",)

    def __init__(self, stream: IO[str] | None = None) -> None:
        self._stream = stream if stream is not None else sys.stdout

    def __call__(self, event: StepEvent) -> None:
        """Print a line if this event calls for one.

        Parameters
        ----------
        event : maws.search.StepEvent
            The event to report. Events other than the start of a step and the
            end of the search produce no output.
        """
        match event:
            case StepStarted(step=step, total=total, sequence=sequence):
                start = (
                    f"growing from {sequence}"
                    if sequence
                    else ("starting from an empty strand")
                )
                print(f"Step {step}/{total}: {start}", file=self._stream)
            case SearchFinished(winner=winner, steps=steps):
                print(
                    f"Done after {steps} steps: {winner.sequence} "
                    f"(E={winner.energy:.2f} kJ/mol, S={winner.entropy:.6f})",
                    file=self._stream,
                )


class LogReporter:
    """Prints progress and writes a record of the whole run to disk.

    Two files are written into `directory`:

    ``<job>_scores.log``
        One line per candidate: its sequence, its score, and its best energy.
        This is the record of how the search chose what it chose.
    ``<job>_steps.pdb``
        The best structure from each step, one model per step, so the strand
        can be watched growing in a molecular viewer.

    Parameters
    ----------
    job : str
        Stem for the filenames.
    directory : str or pathlib.Path, default="."
        Where to write them. Created if it does not exist.
    logger : logging.Logger, optional
        Where progress messages go. Defaults to this module's logger.

    See Also
    --------
    ProgressReporter : Prints progress without writing anything.

    Notes
    -----
    Holds two files open, so use it as a context manager, or call
    :meth:`close` when the search has finished.

    Examples
    --------
    >>> import tempfile
    >>> from maws.search import StepStarted
    >>> from maws.values import NucleotideSequence
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     with LogReporter("demo", tmp) as reporter:
    ...         reporter(StepStarted(step=1, total=1, sequence=NucleotideSequence(())))
    ...     sorted(p.name for p in Path(tmp).iterdir())
    ['demo_scores.log', 'demo_steps.pdb']
    """

    __slots__ = ("_logger", "_scores", "_steps")

    def __init__(
        self,
        job: str,
        directory: str | Path = ".",
        *,
        logger: logging.Logger | None = None,
    ) -> None:
        base = Path(directory)
        base.mkdir(parents=True, exist_ok=True)
        self._logger = logger if logger is not None else logging.getLogger(__name__)
        self._scores = (base / f"{job}_scores.log").open("w")
        self._steps = (base / f"{job}_steps.pdb").open("w")

    def __enter__(self) -> LogReporter:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def close(self) -> None:
        """Close both output files.

        Safe to call more than once.
        """
        for handle in (self._scores, self._steps):
            if not handle.closed:
                handle.close()

    def __call__(self, event: StepEvent) -> None:
        """Report one event.

        Parameters
        ----------
        event : maws.search.StepEvent
            The event to report.
        """
        match event:
            case StepStarted(step=step, total=total, sequence=sequence):
                self._logger.info(
                    "Step %d/%d: %s",
                    step,
                    total,
                    f"growing from {sequence}" if sequence else "empty strand",
                )
            case CandidateScored(candidate=candidate):
                self._logger.info(
                    "  %s %s -> score=%.6f energy=%.2f",
                    candidate.token,
                    candidate.direction,
                    candidate.entropy,
                    candidate.energy,
                )
                self._scores.write(
                    f"SEQUENCE: {candidate.sequence}  "
                    f"SCORE: {candidate.entropy}  "
                    f"ENERGY: {candidate.energy}\n"
                )
            case StepCompleted(step=step, winner=winner):
                self._logger.info("Step %d complete. Best: %s", step, winner.sequence)
                write_model(self._steps, winner.pose, winner.system, model=step)
            case SearchFinished(winner=winner, steps=steps):
                self._logger.info(
                    "Done after %d steps: %s (E=%.2f kJ/mol)",
                    steps,
                    winner.sequence,
                    winner.energy,
                )


class JsonlReporter:
    """Writes one JSON object per event, one per line.

    Suitable for piping into another program, or for keeping a machine-readable
    record of a run alongside the human-readable one.

    Parameters
    ----------
    stream : io.TextIOBase, optional
        Where to write. Defaults to standard output.

    Examples
    --------
    >>> import io
    >>> from maws.search import StepStarted
    >>> from maws.values import NucleotideSequence
    >>> buffer = io.StringIO()
    >>> JsonlReporter(buffer)(
    ...     StepStarted(step=1, total=2, sequence=NucleotideSequence(("G",)))
    ... )
    >>> buffer.getvalue().strip()
    '{"event": "StepStarted", "step": 1, "total": 2, "sequence": "G"}'
    """

    __slots__ = ("_stream",)

    def __init__(self, stream: IO[str] | None = None) -> None:
        self._stream = stream if stream is not None else sys.stdout

    def __call__(self, event: StepEvent) -> None:
        """Write one event as a line of JSON.

        Parameters
        ----------
        event : maws.search.StepEvent
            The event to write.
        """
        self._stream.write(json.dumps(_as_dict(event)) + "\n")


def _as_dict(event: StepEvent) -> dict[str, Any]:
    """Flatten one event into plain JSON-serialisable values.

    Parameters
    ----------
    event : maws.search.StepEvent
        The event to flatten.

    Returns
    -------
    dict
        Always carries an ``"event"`` key naming the kind. Structures and atom
        positions are left out: they do not survive JSON, and the ``.pdb`` file
        :class:`LogReporter` writes is the right place for them.
    """
    match event:
        case StepStarted(step=step, total=total, sequence=sequence):
            return {
                "event": "StepStarted",
                "step": step,
                "total": total,
                "sequence": str(sequence),
            }
        case CandidateScored(step=step, candidate=candidate):
            return {
                "event": "CandidateScored",
                "step": step,
                "sequence": str(candidate.sequence),
                "token": candidate.token,
                "direction": candidate.direction,
                "score": candidate.entropy,
                "energy": candidate.energy,
            }
        case StepCompleted(step=step, winner=winner):
            return {
                "event": "StepCompleted",
                "step": step,
                "sequence": str(winner.sequence),
                "score": winner.entropy,
                "energy": winner.energy,
            }
        case SearchFinished(winner=winner, steps=steps):
            return {
                "event": "SearchFinished",
                "steps": steps,
                "sequence": str(winner.sequence),
                "score": winner.entropy,
                "energy": winner.energy,
            }
    return {"event": type(event).__name__}


def write_model(stream: IO[str], pose: Any, system: Any, *, model: int = 1) -> None:
    """Write one structure to an open PDB file as a numbered model.

    Parameters
    ----------
    stream : io.TextIOBase
        An open file to append to.
    pose : maws.pose.Pose
        The atom positions to write.
    system : maws.topology.BuiltSystem
        The structure they belong to, which supplies the atom names and
        residues.
    model : int, default=1
        The model number to label this structure with. A viewer shows numbered
        models as frames, so writing one per step animates the strand growing.

    Notes
    -----
    Does nothing when `system` was not built by AmberTools, since there is then
    no list of atoms and residues to write. That is the case for the stand-in
    builder used in tests, so a test may pass a reporter without arranging for
    a real structure.
    """
    amber = getattr(system, "amber", None)
    if amber is None:
        return

    from openmm import app

    app.PDBFile.writeModel(
        amber.topology, pose.to_openmm(), file=stream, modelIndex=model
    )

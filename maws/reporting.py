"""
maws.reporting
==============

Turning a search's progress into messages and files.

:func:`maws.search.grow_aptamer` grows an *aptamer* — a short strand of RNA or
DNA designed to stick to a target molecule — one nucleotide at a time, and
reports what it is doing as a stream of events. It does not print anything and
does not write anything: what to do with each event is left to whoever is
reading them. This module holds the readers.

:class:`ProgressReporter`
    Prints one line per step and nothing else.
:class:`LogReporter`
    Prints progress and writes two files: a record of every candidate and its
    score, and the best structure from each step.
:class:`JsonlReporter`
    Writes one JSON object per event, for feeding into another program.

Each is a callable taking one event, so any of them can be passed straight to
:func:`maws.api.design` as ``on_event``, and writing another takes a few lines.
None of them returns anything or changes the event it is given, so several can
be run over the same search without interfering.

Structures are written as PDB files: the plain-text format molecular viewers
read, one line per atom giving its name, its residue and its position in
ångström. Several structures can share one file as numbered *models*, which a
viewer shows as the frames of an animation.

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
    """An event reader that prints one line per step, and one at the end.

    The quietest of the three readers: it ignores every candidate that is
    tried and reports only which step has begun and what the search settled
    on. Nothing is written to disk.

    Parameters
    ----------
    stream : io.TextIOBase, optional
        Where to print. Defaults to standard output. Pass an
        :class:`io.StringIO` to collect the lines instead of showing them.

    See Also
    --------
    LogReporter : Also writes files and reports every candidate.
    JsonlReporter : Reports every event, as JSON rather than prose.

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
        r"""Print a line if this event calls for one.

        Parameters
        ----------
        event : maws.search.StepEvent
            The event to report. Events other than the start of a step and the
            end of the search produce no output.

        Examples
        --------
        The start of a step is announced; an individual candidate's score is
        not:

        >>> import io
        >>> from maws.search import CandidateScored, StepStarted
        >>> from maws.values import NucleotideSequence
        >>> buffer = io.StringIO()
        >>> reporter = ProgressReporter(buffer)
        >>> reporter(StepStarted(step=3, total=8, sequence=NucleotideSequence(("G",))))
        >>> reporter(CandidateScored(step=3, candidate=None))
        >>> buffer.getvalue()
        'Step 3/8: growing from G\n'
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
                    f"(E={winner.energy:.2f} kJ/mol, score={winner.score:.2f})",
                    file=self._stream,
                )


class LogReporter:
    """An event reader that keeps a record of the whole run on disk.

    Two files are written into `directory`:

    ``<job>_scores.log``
        One line per candidate: its sequence, its score, and its best energy in
        kJ/mol. This is the record of how the search chose what it chose.
    ``<job>_steps.pdb``
        The best structure from each step, one model per step, so the strand
        can be watched growing in a molecular viewer.

    Progress messages go through the :mod:`logging` module rather than being
    printed, so a program embedding MAWS decides where they end up and at what
    level of detail.

    Parameters
    ----------
    job : str
        Stem the two filenames are built from. A run called ``"demo"`` writes
        ``demo_scores.log`` and ``demo_steps.pdb``.
    directory : str or pathlib.Path, default="."
        Where to write them. Created if it does not exist. Files already there
        under these names are overwritten.
    logger : logging.Logger, optional
        Where progress messages go. Defaults to this module's own logger, which
        a program can configure by name.

    See Also
    --------
    ProgressReporter : Prints progress without writing anything.
    write_model : Appends one structure to the ``.pdb`` file.

    Notes
    -----
    .. warning::
        Both files are opened as soon as the reporter is made and stay open.
        Use it as a context manager, or call :meth:`close` when the search has
        finished, or the last few lines may never reach the disk.

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

        Safe to call more than once, so a program that closes explicitly and
        also uses the reporter as a context manager comes to no harm.

        Examples
        --------
        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as tmp:
        ...     reporter = LogReporter("demo", tmp)
        ...     reporter.close()
        ...     reporter.close()  # closing twice is not an error
        ...     (Path(tmp) / "demo_scores.log").exists()
        True
        """
        for handle in (self._scores, self._steps):
            if not handle.closed:
                handle.close()

    def __call__(self, event: StepEvent) -> None:
        """Report one event, to the log and where applicable to a file.

        Every kind of event produces a log message. A scored candidate is also
        appended to the ``.log`` file, and a completed step appends its winning
        structure to the ``.pdb`` file.

        Parameters
        ----------
        event : maws.search.StepEvent
            The event to report.

        Examples
        --------
        >>> import tempfile
        >>> from pathlib import Path
        >>> from maws.search import StepStarted
        >>> from maws.values import NucleotideSequence
        >>> with tempfile.TemporaryDirectory() as tmp:
        ...     with LogReporter("demo", tmp) as reporter:
        ...         reporter(
        ...             StepStarted(step=1, total=2, sequence=NucleotideSequence(()))
        ...         )
        ...     (Path(tmp) / "demo_scores.log").read_text()
        ''
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
                    candidate.score,
                    candidate.energy,
                )
                self._scores.write(
                    f"SEQUENCE: {candidate.sequence}  "
                    f"SCORE: {candidate.score}  "
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
    """An event reader that writes one JSON object per event, one per line.

    Suitable for piping into another program, or for keeping a
    machine-readable record of a run alongside the human-readable one. Every
    event is written, and each line stands on its own, so a reader can act on
    the search while it is still running.

    Parameters
    ----------
    stream : io.TextIOBase, optional
        Where to write. Defaults to standard output. Each line is written as
        soon as its event arrives, but the stream is not flushed, so a
        long-running pipe may lag behind the search.

    See Also
    --------
    ProgressReporter : The same events as prose, for a person to read.

    Notes
    -----
    Structures and atom positions are left out of the JSON: they do not fit
    into it usefully. Use :class:`LogReporter` when the coordinates are wanted
    as well.

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
        r"""Write one event as a line of JSON.

        Parameters
        ----------
        event : maws.search.StepEvent
            The event to write. Every kind produces a line; an unrecognised one
            produces a line carrying only the event name.

        Examples
        --------
        >>> import io
        >>> from maws.search import StepEvent
        >>> buffer = io.StringIO()
        >>> JsonlReporter(buffer)(StepEvent())
        >>> buffer.getvalue()
        '{"event": "StepEvent"}\n'
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
                "score": candidate.score,
                "energy": candidate.energy,
            }
        case StepCompleted(step=step, winner=winner):
            return {
                "event": "StepCompleted",
                "step": step,
                "sequence": str(winner.sequence),
                "score": winner.score,
                "energy": winner.energy,
            }
        case SearchFinished(winner=winner, steps=steps):
            return {
                "event": "SearchFinished",
                "steps": steps,
                "sequence": str(winner.sequence),
                "score": winner.score,
                "energy": winner.energy,
            }
    return {"event": type(event).__name__}


def write_model(stream: IO[str], pose: Any, system: Any, *, model: int = 1) -> None:
    """write_model(stream, pose, system, *, model=1) -> None

    Write one structure to an open PDB file as a numbered model.

    Parameters
    ----------
    stream : io.TextIOBase
        An open text file to append to. Left open afterwards, so several
        structures can be added to the same file.
    pose : maws.pose.Pose
        The atom positions to write.
    system : maws.topology.BuiltSystem
        The structure they belong to, which supplies the atom names and
        residues. The positions are not checked against it.
    model : int, default=1
        The number to label this structure with. A viewer shows numbered models
        as frames, so writing one per step animates the strand growing. Reusing
        a number does not overwrite anything; it just confuses the viewer.

    See Also
    --------
    LogReporter : Calls this once per completed step.

    Notes
    -----
    .. note::
        Does nothing at all when `system` has no AmberTools output attached,
        since there is then no list of atoms and residues to write against. The
        grid-based builder used in tests produces such a system, so a test may
        pass a reporter without arranging for a real structure.

    Examples
    --------
    >>> from pathlib import Path  # doctest: +SKIP
    >>> with Path("frames.pdb").open("w") as handle:  # doctest: +SKIP
    ...     write_model(handle, system.pose, system, model=1)
    ...     write_model(handle, grown.pose, grown, model=2)
    >>> sum(
    ...     1 for line in Path("frames.pdb").open() if line.startswith("MODEL")
    ... )  # doctest: +SKIP
    2
    """
    amber = getattr(system, "amber", None)
    if amber is None:
        return

    from openmm import app

    app.PDBFile.writeModel(
        amber.topology, pose.to_openmm(), file=stream, modelIndex=model
    )

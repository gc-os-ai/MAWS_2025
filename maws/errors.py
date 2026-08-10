"""
maws.errors
===========

Every exception MAWS raises on purpose, and what each one asks you to do.

MAWS designs an *aptamer*: a short strand of RNA or DNA that folds up against a
target molecule and sticks to it. Doing that means reading files a person wrote
and running programs that are installed separately from this package, so a
failure comes from one of two places, and the type of the exception says which:

1. *The run was described wrongly* — :class:`ConfigurationError`. An argument,
   a sequence, or a file path has to change. Running the same thing again, on
   any machine, fails the same way.
2. *An outside program broke* — :class:`ToolchainError`. The description is
   fine; a program MAWS hands work to is missing or exited with an error. The
   very same run may succeed on a machine where that program is installed.

The outside programs are AmberTools. ``tleap`` — the program itself is called
LEaP — turns a sequence plus a set of parameter files into a built structure.
``antechamber`` and ``parmchk2`` work out parameters for a molecule the force
field does not already cover; a *force field* is the collection of numbers
describing how the atoms of a molecule push and pull on each other.

Every exception here inherits from :class:`MawsError`, so ``except MawsError``
catches everything this package raises deliberately and nothing else.

Examples
--------
>>> from maws.errors import ConfigurationError, MawsError
>>> issubclass(ConfigurationError, MawsError)
True
>>> try:
...     raise ConfigurationError("sequence 'X' is not in this library")
... except MawsError as exc:
...     print(exc)
sequence 'X' is not in this library
"""

from __future__ import annotations

from pathlib import Path

__all__ = [
    "BuildError",
    "ConfigurationError",
    "MawsError",
    "MissingLibrary",
    "SamplingError",
    "ToolchainError",
]


class MawsError(Exception):
    """Base class for every error MAWS raises deliberately.

    Catching this catches all of them. Anything else that escapes from MAWS,
    a bare ``TypeError`, is a bug in the package rather than a failure
    it knows how to explain.
    """


class ConfigurationError(MawsError):
    """The run was described wrongly, and the description can be corrected.

    Raised for a nucleotide the residue library does not know, residue tables
    whose rows do not line up, two chains given the same name, and an index
    that names no residue, atom or bond. The message states what was asked for
    and lists what was there instead, so the fix is usually one edit to the
    arguments.

    Examples
    --------
    >>> from maws.libraries import rna
    >>> try:
    ...     rna().alias("X")
    ... except ConfigurationError as exc:
    ...     print(str(exc).split(".")[0])
    token 'X' is not in this library
    """


class ToolchainError(MawsError):
    """An outside program is missing, or ran and failed.

    Raised when ``tleap``, ``antechamber`` or ``parmchk2`` cannot be found on
    ``PATH``, or is found and exits with a non-zero status. Nothing about the
    run description is at fault, so changing the sequence or the arguments will
    not help; installing AmberTools, or putting it on ``PATH``, will.
    """


class BuildError(ToolchainError):
    """LEaP ran to the end but left behind nothing that can be simulated.

    LEaP reports many problems by printing a complaint and carrying on, so an
    exit status of zero does not mean it produced a structure. This is raised
    when the ``.prmtop`` file it was told to write — the file holding the
    parameters and the list of which atoms are bonded to which — is missing or
    empty afterwards.
    """


class SamplingError(MawsError):
    """A sampler gave up before finding a position for the aptamer.

    Samplers propose a position near the target and throw away any proposal
    that would push the aptamer into the target's atoms. After a fixed number
    of consecutive rejections they stop rather than loop forever, which is the
    point at which this is raised.

    A sampler that never succeeds usually means the region it draws from lies
    entirely inside the target. Widening that region (a larger ``reach``),
    shrinking the clearance demanded around each target atom (a smaller
    ``probe``), or checking that the target really is the size you expect are
    the three things worth trying.
    """


class MissingLibrary(ConfigurationError, FileNotFoundError):  # noqa: N818
    """One or more residue library files are not on disk.

    LEaP only knows what atoms a residue has if it can read a ``.lib`` file
    describing that residue. This is raised before LEaP is started, naming
    every residue whose file is absent, rather than letting LEaP fail one
    residue at a time.

    Inherits from :class:`FileNotFoundError` as well as
    :class:`ConfigurationError`, so ``except FileNotFoundError`` catches it
    too.

    Parameters
    ----------
    directory : pathlib.Path
        The directory that was searched. Quoted in the message so the fix does
        not require guessing where MAWS looked.
    residues : sequence of str
        Names of the residues whose ``.lib`` files were not found. Kept on the
        exception as `residues`, so code that catches it can report or
        regenerate exactly those.

    Examples
    --------
    >>> from pathlib import Path
    >>> err = MissingLibrary(Path("/params"), ["LIG", "GN"])
    >>> err.residues
    ['LIG', 'GN']
    >>> "missing from /params: LIG, GN" in str(err)
    True
    """

    def __init__(self, directory: Path, residues: list[str]) -> None:
        self.directory = directory
        self.residues = list(residues)
        super().__init__(
            f"{len(self.residues)} LEaP library file(s) missing from "
            f"{directory}: {', '.join(self.residues)}. "
            f"Run `maws prepare` to generate them."
        )

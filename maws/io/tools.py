"""
maws.io.tools
=============

Running the AmberTools programs, and reporting their failures usefully.

AmberTools is a suite of programs installed separately from this package.
MAWS calls three of them: ``antechamber`` and ``parmchk2`` work out the *force
field* — the numbers describing how a molecule's atoms pull on and push away
each other — for a molecule nothing already describes, and ``tleap`` turns a
sequence plus those numbers into a built structure. The program ``tleap`` is
the command-line form of a modelling program called LEaP, and the two names are
used interchangeably.

They are ordinary command-line programs, so the two things that go wrong are
"it is not installed" and "it ran and failed". :func:`find_tool` reports the
first, :func:`run_tool` reports both, and each raises a
:class:`~maws.errors.ToolchainError` naming the program and what to do about
it. Nothing else in MAWS starts another program.

See Also
--------
maws.errors.ToolchainError : What both functions raise.
maws.build.LeapBuilder : The main user of ``tleap``.

Examples
--------
>>> find_tool("python").endswith("python")
True
>>> print(run_tool(["python", "-c", "pass"]))  # nothing to return on success
None
"""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from maws.errors import ToolchainError

__all__ = ["find_tool", "run_tool"]


def find_tool(name: str) -> str:
    """Return the full path to an installed program.

    Parameters
    ----------
    name : str
        Program to look for, e.g. ``"tleap"``. Looked up on ``PATH`` exactly as
        a shell would, so which one is found depends on the active environment.

    Returns
    -------
    str
        The full path to it, which is what should be run: resolving the name
        once here means a later failure cannot be blamed on the wrong copy of
        the program having been picked up.

    Raises
    ------
    maws.errors.ToolchainError
        If it is not on ``PATH``. The message names the program and says where
        it comes from, since the usual cause is that AmberTools is not
        installed or its environment is not active.

    See Also
    --------
    run_tool : Finds a program and runs it.

    Examples
    --------
    >>> find_tool("python").endswith("python")
    True

    A program that does not exist is reported rather than guessed at:

    >>> from maws.errors import ToolchainError
    >>> try:
    ...     find_tool("definitely-not-a-real-program")
    ... except ToolchainError as exc:
    ...     print(str(exc).split(".")[0])
    definitely-not-a-real-program is not on PATH
    """
    found = shutil.which(name)
    if found is None:
        raise ToolchainError(
            f"{name} is not on PATH. It is part of AmberTools, which is "
            f"installed separately from this package; check that AmberTools "
            f"is installed and its environment is active."
        )
    return found


def run_tool(command: Sequence[str], *, cwd: str | Path | None = None) -> None:
    """run_tool(command, *, cwd=None) -> None

    Run one program and wait for it, raising if it fails.

    Parameters
    ----------
    command : sequence of str
        The program and its arguments. The first entry is looked up on
        ``PATH``; the rest are passed through untouched. No shell is involved,
        so quoting, globbing and redirection have no effect here.
    cwd : str or pathlib.Path, optional
        Directory to run in. Defaults to the current working directory. The
        AmberTools programs scatter intermediate files beside wherever they
        run, so this is normally a temporary directory that is thrown away
        afterwards.

    Returns
    -------
    None
        Nothing. A program that fails raises instead of returning a status.

    Raises
    ------
    maws.errors.ToolchainError
        If the program is not on ``PATH``, or exits with a non-zero status.

    See Also
    --------
    find_tool : Does the ``PATH`` lookup, and can be used on its own to check
        that a program is available before starting any work.

    Notes
    -----
    .. note::
        Output is not captured, so whatever the program prints goes straight to
        the terminal. That is deliberate: ``tleap`` explains what it could not
        find in its own output, and hiding it would leave the error message
        raised here with nothing useful to say.

    Examples
    --------
    >>> print(run_tool(["python", "-c", "pass"]))
    None

    A non-zero exit status becomes an exception:

    >>> from maws.errors import ToolchainError
    >>> try:
    ...     run_tool(["python", "-c", "raise SystemExit(3)"])
    ... except ToolchainError as exc:
    ...     print(str(exc).split(".")[0])
    python exited with status 3
    """
    executable = find_tool(command[0])
    try:
        subprocess.run([executable, *command[1:]], cwd=cwd, check=True)
    except subprocess.CalledProcessError as exc:
        raise ToolchainError(
            f"{command[0]} exited with status {exc.returncode}. Its own "
            f"output above normally says what it could not do."
        ) from exc

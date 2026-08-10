"""
maws.io.tools
=============

Running the AmberTools programs, and reporting their failures usefully.

AmberTools is a suite of programs installed separately from this package.
MAWS calls three of them: ``antechamber`` and ``parmchk2`` work out the
parameters describing how a molecule's atoms pull on each other, and ``tleap``
turns a sequence plus those parameters into a built structure.

They are ordinary command-line programs, so the two things that go wrong are
"it is not installed" and "it ran and failed". :func:`run_tool` turns both into
a :class:`~maws.errors.ToolchainError` that says which program and what to do.

Examples
--------
>>> run_tool(["python", "-c", "pass"])
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
        Program to look for, e.g. ``"tleap"``. Searched for on ``PATH``.

    Returns
    -------
    str
        The full path to it.

    Raises
    ------
    maws.errors.ToolchainError
        If it is not on ``PATH``. The message names the program and says where
        it comes from, since the usual cause is that AmberTools is not
        installed or its environment is not active.

    Examples
    --------
    >>> find_tool("python").endswith("python")
    True
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
    """Run one program and wait for it, raising if it fails.

    Parameters
    ----------
    command : sequence of str
        The program and its arguments. The first entry is looked up on
        ``PATH``.
    cwd : str or pathlib.Path, optional
        Directory to run in. Defaults to the current working directory. The
        AmberTools programs write intermediate files beside wherever they run,
        so this is normally a temporary directory.

    Raises
    ------
    maws.errors.ToolchainError
        If the program is not on ``PATH``, or exits with a non-zero status.

    Notes
    -----
    Output is not captured, so whatever the program prints goes straight to
    the terminal. That is deliberate: ``tleap`` explains what it could not find
    in its own output, and hiding it would leave the error message here with
    nothing useful to say.

    Examples
    --------
    >>> run_tool(["python", "-c", "pass"])
    """
    executable = find_tool(command[0])
    try:
        subprocess.run([executable, *command[1:]], cwd=cwd, check=True)
    except subprocess.CalledProcessError as exc:
        raise ToolchainError(
            f"{command[0]} exited with status {exc.returncode}. Its own "
            f"output above normally says what it could not do."
        ) from exc

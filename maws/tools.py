"""
maws.tools
==========

Utility functions for executing external AmberTools programs.

This module provides thin wrappers for finding and running external executables
required by MAWS (tleap, antechamber, parmchk2, etc.).

Functions
---------
find_exe : Locate an executable on PATH.
run : Execute a command with subprocess.

Examples
--------
>>> from maws.tools import find_exe
>>> # Only works if tleap is installed
>>> # tleap_path = find_exe('tleap')
"""

import shutil
import subprocess
from pathlib import Path


class ExecError(RuntimeError):
    """
    Exception raised when a required executable is not found.

    Raised by :func:`find_exe` when the requested program is not on PATH.
    """

    pass


def find_exe(name: str) -> str:
    """
    Locate an executable on the system PATH.

    Parameters
    ----------
    name : str
        Name of the executable to find (e.g., 'tleap', 'antechamber').

    Returns
    -------
    str
        Full path to the executable.

    Raises
    ------
    ExecError
        If the executable is not found on PATH.

    """
    exe = shutil.which(name)
    if not exe:
        raise ExecError(
            f"{name} not found on PATH. Install AmberTools and ensure it's on PATH."
        )
    return exe


def run(cmd: list[str], cwd: str | Path | None = None) -> None:
    """
    Execute a command using subprocess.

    Parameters
    ----------
    cmd : list[str]
        Command and arguments to execute.
    cwd : str or Path, optional
        Working directory for the command.

    Raises
    ------
    subprocess.CalledProcessError
        If the command exits with non-zero status.

    Examples
    --------
    >>> from maws.tools import run
    >>> run(["echo", "hello"])  # Runs silently
    """
    subprocess.run(cmd, cwd=cwd, check=True)

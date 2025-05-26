"""
prepare.py – utilities to parametrize small molecules / peptides with AmberTools
and to toggle hydrogens using `reduce`.

The heavy lifting (`antechamber`, `parmchk`, `tleap`, `reduce`) happens in
external binaries; this module only builds command lines and runs them with
sub-processes.  All paths are handled by ``pathlib.Path``.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Optional

from openmm import app


# --------------------------------------------------------------------------- #
# _private helpers
# --------------------------------------------------------------------------- #
def _run(cmd: str | list[str], *, cwd: Optional[Path] = None) -> None:
    """Wrapper around subprocess.run with sane defaults."""
    subprocess.run(cmd, cwd=cwd, check=True, shell=isinstance(cmd, str))


def _write_file(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def make_lib(
    file_path: str | Path,
    residue_name: str,
    *,
    connect0: Optional[str] = None,
    connect1: Optional[str] = None,
    charges: str = "bcc",
    atom_type: str = "gaff",
    force_field: str = "leaprc.protein.ff14SB",
    parameterized: bool = False,
) -> int:
    """
    Parametrise *file_path* and write ``<residue_name>.lib`` for tleap.

    Parameters
    ----------
    file_path
        Input PDB or MOL2 file.
    residue_name
        Name under which the residue is saved in the library.
    connect0, connect1
        Atom names to mark as *head* / *tail* if you want peptide-like
        connectivity.
    charges, atom_type, force_field
        Controls passed through to *antechamber* / *tleap*.
    parameterized
        If *True*, assume ``file_path`` is already parametrised (GAFF atoms &
        charges present) and skip antechamber/parmchk.

    Returns
    -------
    int
        Atom count of the original structure (via OpenMM).
    """
    in_path = Path(file_path).expanduser().resolve()
    name = in_path.stem
    ext = in_path.suffix.lstrip(".").lower()

    if ext not in {"pdb", "mol2"}:
        raise ValueError("file_path must be .pdb or .mol2")

    lib_base = in_path.with_name(residue_name)  # <dir>/<residue_name>
    leap_in = in_path.with_suffix(".in")        # temporary input to tleap

    # ------------------------------------------------------------------ #
    # 1.  Build tleap script
    # ------------------------------------------------------------------ #
    if parameterized:
        leap_script = f"""
source leaprc.gaff
source {force_field}
{residue_name} = loadpdb {in_path.name}
saveoff {residue_name} {lib_base}.lib
savepdb {residue_name} {lib_base}_tmp.pdb
quit
"""
    else:
        leap_script = f"""
source {force_field}
source leaprc.gaff
{residue_name} = loadmol2 {name}.mol2
loadamberparams {lib_base}.frcmod
"""

        if connect0 and connect1:
            leap_script += f"""
set {residue_name} head {residue_name}.1.{connect0}
set {residue_name} tail {residue_name}.1.{connect1}
set {residue_name}.1 connect0 {residue_name}.head
set {residue_name}.1 connect1 {residue_name}.tail
"""

        leap_script += f"""
check {residue_name}
saveoff {residue_name} {lib_base}.lib
savepdb {residue_name} {lib_base}_tmp.pdb
quit
"""

    _write_file(leap_in, leap_script.strip())

    # ------------------------------------------------------------------ #
    # 2.  Run external AmberTools commands
    # ------------------------------------------------------------------ #
    if not parameterized:
        _run(
            (
                "antechamber",
                "-i",
                str(in_path),
                "-fi",
                ext,
                "-nc",
                "2.02",
                "-o",
                f"{name}.mol2",
                "-fo",
                "mol2",
                "-c",
                charges,
                "-rn",
                residue_name,
                "-at",
                atom_type,
            )
        )
        _run(
            (
                "parmchk",
                "-i",
                f"{name}.mol2",
                "-f",
                "mol2",
                "-o",
                f"{lib_base}.frcmod",
            )
        )

    # _run(f"tleap -f {leap_in}")
    workdir = in_path.parent          # <-- directory with the PDB
    _run(f"tleap -f {leap_in.name}", cwd=workdir)

    # ------------------------------------------------------------------ #
    # 3.  Atom count for caller
    # ------------------------------------------------------------------ #
    pdb = app.PDBFile(str(in_path))
    atom_count = sum(1 for _ in pdb.topology.atoms())

    # ------------------------------------------------------------------ #
    # 4.  House-keeping
    # ------------------------------------------------------------------ #
    leap_in.unlink(missing_ok=True)
    Path("leap.log").unlink(missing_ok=True)
    Path(f"{lib_base}_tmp.pdb").unlink(missing_ok=True)

    if not parameterized:
        Path(f"{name}.mol2").unlink(missing_ok=True)
        Path(f"{lib_base}.frcmod").unlink(missing_ok=True)

    return atom_count


def toggle_hydrogens(path: str | Path, *, keep: bool = True) -> None:
    """
    Use the **reduce** tool to strip (and optionally rebuild) hydrogens
    in-place.

    Parameters
    ----------
    path
        PDB file to operate on.
    keep
        If *True* (default) hydrogens are re-added after trimming.  If *False*
        the PDB remains without hydrogens.
    """
    p = Path(path).expanduser().resolve()
    _run(f"reduce -Trim {p} > {p}")
    if keep:
        _run(f"reduce -Build {p} > {p}")


__all__ = ["make_lib", "toggle_hydrogens"]

"""
Thin, testable wrappers around AmberTools executables for ligand/fragment prep,
plus a pure-Python hydrogen toggle using OpenMM.

Functions
---------
makeLib
    Parameterize a molecule (optionally pre-parameterized) and produce
    Amber OFF library (.lib) and frcmod (if needed) using an isolated temp dir.

toggleHydrogens
    Add or remove hydrogens in-place on a PDB using OpenMM's Modeller.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from openmm import app
from openmm.app import ForceField, Modeller, PDBFile

from maws.tools import find_exe, run


def make_lib(
    file_path: str | Path,
    residue_name: str,
    connect0: str | None = None,
    connect1: str | None = None,
    charges: str = "bcc",
    atom_type: str = "gaff",
    force_field_aptamer: str = "leaprc.RNA.OL3",
    force_field_ligand: str = "leaprc.protein.ff19SB",
    parameterized: bool = False,
) -> int:
    """
    Generate Amber OFF library (.lib) and (when needed) a .frcmod for a residue.

    This function wraps `antechamber`, `parmchk2`, and `tleap` behind a clean
    Python API. It works entirely in a temporary directory and moves only the
    final artifacts next to the input file.

    Parameters
    ----------
    file_path
        Path to the input structure (e.g., .pdb, .mol2, .sdf). If
        `parameterized=True`, this must be a PDB already carrying the coordinates
        you want used by LEaP.
    residue_name
        Name of the residue created inside LEaP (used to name .lib/.frcmod).
    connect0, connect1
        Optional atom names (within the new residue) used to define polymer
        head/tail connectivity for LEaP (`set head/tail`, `connect0/connect1`).
        Both must be provided to activate connectivity directives.
    charges
        Charge model for antechamber (e.g., 'bcc' for AM1-BCC). Ignored if
        `parameterized=True`.
    atom_type
        Antechamber atom type set (e.g., 'gaff' or 'gaff2'). Ignored if
        `parameterized=True`.
    force_field_aptamer
        LEaP “source” line for the aptamer/nucleic acid FF (e.g., RNA.OL3/DNA.OL21).
    force_field_ligand
        LEaP “source” line for the ligand/protein FF (e.g., ff19SB or gaff2).
    parameterized
        If True, skip antechamber/parmchk2 and `loadpdb` instead of `loadmol2`.
        Useful when the input is already parameterized or you only need a
        polymerizable building block with coordinates.

    Returns
    -------
    int
        Atom count of the temporary PDB produced by LEaP (matches legacy behavior).

    Side Effects
    ------------
    Writes the following next to the input file:
      - `<residue_name>.lib`
      - `<residue_name>.frcmod` (only when `parameterized` is False)

    Raises
    ------
    ExecError
        If required executables (`antechamber`, `parmchk2`, `tleap`) are not on PATH.
    CalledProcessError
        If any external command exits with a non-zero status.

    Notes
    -----
    - No `conda run` is used; tools must be discoverable via `PATH`.
    - Work is done in a temp dir; only final artifacts are moved.
    """
    src = Path(file_path).resolve()
    name, ext = src.stem, src.suffix[1:].lower()
    out_base = src.parent / residue_name

    with tempfile.TemporaryDirectory() as td:
        w = Path(td)

        if not parameterized:
            # antechamber: input may be pdb/mol2/sdf...; output is temp/{name}.mol2
            run(
                [
                    find_exe("antechamber"),
                    "-i",
                    str(src),
                    "-fi",
                    ext,
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
                ],
                cwd=w,
            )

            # parmchk2 → temp/{residue_name}.frcmod
            run(
                [
                    find_exe("parmchk2"),
                    "-i",
                    f"{name}.mol2",
                    "-f",
                    "mol2",
                    "-o",
                    f"{residue_name}.frcmod",
                ],
                cwd=w,
            )
        else:
            # Ensure LEaP can find the PDB in the temp dir
            shutil.copy2(src, w / f"{name}.pdb")

        # Build LEaP input
        lines = [
            f"source {force_field_aptamer}",
            f"source {force_field_ligand}",
        ]
        if parameterized:
            lines.append(f"{residue_name} = loadpdb {name}.pdb")
        else:
            lines += [
                f"{residue_name} = loadmol2 {name}.mol2",
                f"loadamberparams {residue_name}.frcmod",
            ]

        if connect0 and connect1:
            lines += [
                f"set {residue_name} head {residue_name}.1.{connect0}",
                f"set {residue_name} tail {residue_name}.1.{connect1}",
                f"set {residue_name}.1 connect0 {residue_name}.head",
                f"set {residue_name}.1 connect1 {residue_name}.tail",
            ]

        lines += [
            f"check {residue_name}",
            f"saveoff {residue_name} {residue_name}.lib",
            f"savepdb {residue_name} {residue_name}_tmp.pdb",
            "quit",
        ]

        (w / "leap.in").write_text("\n".join(lines))
        run([find_exe("tleap"), "-f", "leap.in"], cwd=w)

        # Move outputs next to the input
        shutil.move(w / f"{residue_name}.lib", out_base.with_suffix(".lib"))
        if not parameterized:
            shutil.move(w / f"{residue_name}.frcmod", out_base.with_suffix(".frcmod"))

        # Report atom count as before
        pdb = app.PDBFile(str(w / f"{residue_name}_tmp.pdb"))
        length = sum(1 for _ in pdb.topology.atoms())

    return length


def toggle_hydrogens(path: str, add: bool = True, ph: float = 7.0) -> None:
    """
    Add or strip hydrogens in-place on a PDB using OpenMM Modeller.

    Parameters
    ----------
    path
        Path to a PDB file to modify in-place.
    add
        If True, add hydrogens according to the provided force field and pH.
        If False, remove all hydrogens.
    ph
        Target pH used by `Modeller.addHydrogens()`.

    Notes
    -----
    - Uses `amber19-all.xml` and `amber19/tip3pfb.xml` when adding hydrogens.
      Ensure OpenMM’s forcefield files are installed/available.
    """
    pdb = PDBFile(path)
    modeller = Modeller(pdb.topology, pdb.positions)
    if add:
        ff = ForceField("amber19-all.xml", "amber19/tip3pfb.xml")
        modeller.addHydrogens(ff, pH=ph)
    else:
        hydrogens = [
            a for a in modeller.getTopology().atoms() if a.element.symbol == "H"
        ]
        modeller.delete(hydrogens)
    with open(path, "w") as f:
        PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(), f)

"""
maws.io.prepare
===============

Getting a target molecule ready for the program that builds structures.

MAWS designs an *aptamer*: a short strand of RNA or DNA that folds up against a
target molecule and sticks to it. Before anything can be built or scored, the
target has to be described in terms the modelling program understands. That
description has two parts, and this module produces both.

A *force field* is the collection of numbers saying how strongly each pair of
atoms in a molecule pulls on or pushes away the other. Standard force fields
already cover proteins and nucleic acids, but not an arbitrary small organic
molecule, so for those the numbers have to be worked out first. Three
AmberTools programs, installed separately from this package, do that work:

``antechamber``
    Reads a structure file, works out what chemical type each atom is, and
    assigns it a partial charge. Writes a ``mol2`` file — a text format
    carrying atoms, bonds, types and charges together.
``parmchk2``
    Reads that ``mol2`` file and writes a ``.frcmod`` file supplying whatever
    parameters the force field is missing for those atom types.
``tleap``
    The command-line form of LEaP, the modelling program that assembles
    molecules. Given those two pieces it writes a ``.lib`` file: a residue
    definition recording the molecule's atoms, their names, types, charges and
    bonds, so that later runs can refer to the whole molecule by a single
    residue name.

:func:`make_lib` runs whichever of the three are needed and leaves the ``.lib``
and ``.frcmod`` files behind. :func:`toggle_hydrogens` is separate and needs
none of them: it adds or strips the hydrogen atoms of a PDB file, the
plain-text format that lists one atom per line with its position in ångström.

.. note::
    Every function here needs AmberTools or OpenMM installed and a real
    structure file to work on, so none of the examples run on their own.

See Also
--------
maws.io.tools.run_tool : Starts each of the three programs.
maws.build.LeapBuilder : Calls :func:`make_lib` and then builds with the result.

Examples
--------
Put the hydrogens back into a protein structure, then describe it as one
residue called ``LIG``:

>>> toggle_hydrogens("target.pdb")  # doctest: +SKIP
>>> make_lib("target.pdb", "LIG", parameterized=True)  # doctest: +SKIP
1462
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from openmm import app
from openmm.app import ForceField, Modeller, PDBFile

from maws.io.tools import run_tool


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
    output_dir: str | Path | None = None,
) -> int:
    """make_lib(file_path, residue_name, connect0=None, connect1=None,
                charges="bcc", atom_type="gaff",
                force_field_aptamer="leaprc.RNA.OL3",
                force_field_ligand="leaprc.protein.ff19SB",
                parameterized=False, output_dir=None) -> int

    Describe a molecule as a single named residue, and write that description.

    A ``.lib`` file is always written, holding the residue definition LEaP will
    later load. A ``.frcmod`` file, supplying force field parameters LEaP does
    not otherwise have, is written as well whenever `parameterized` is False.
    Both are named after `residue_name`, so a second call with the same name
    overwrites them.

    All the intermediate files go into a temporary directory that is thrown
    away, so nothing is left behind except those one or two results.

    Parameters
    ----------
    file_path : str or pathlib.Path
        The structure to read. Its suffix tells ``antechamber`` what format to
        expect, so ``.pdb``, ``.mol2`` and ``.sdf`` all work. When
        `parameterized` is True this must be a PDB file, and the positions in
        it are the ones LEaP will use.
    residue_name : str
        What to call the molecule from now on. Everything afterwards refers to
        it by this name, including the sequence handed to ``tleap``, and it
        names the two output files. Keep it to the three or four upper-case
        characters a PDB file has room for, e.g. ``"LIG"``.
    connect0, connect1 : str, optional
        Atom names within the new residue to use as its two joining points, so
        that copies of it can be chained together. Both must be given, or
        neither: with only one, no joining points are set at all and the
        residue can only stand alone.
    charges : str, default="bcc"
        Which method ``antechamber`` uses to work out each atom's partial
        charge. ``"bcc"`` is the usual choice: fast, and accurate enough for
        comparing shapes of the same molecule. Not used when `parameterized`
        is True, because then no charges are being worked out.
    atom_type : str, default="gaff"
        Which family of chemical atom types ``antechamber`` assigns. This has
        to agree with `force_field_ligand`, so a `force_field_ligand`
        mentioning ``gaff2`` overrides whatever is passed here with ``"gaff2"``
        rather than producing a structure that cannot be built. Not used when
        `parameterized` is True.
    force_field_aptamer : str, default="leaprc.RNA.OL3"
        Name LEaP is given to load the nucleic acid parameters. It matters even
        though no aptamer is involved yet, because the residue definition must
        be written against the same parameters the design run will use.
    force_field_ligand : str, default="leaprc.protein.ff19SB"
        Name LEaP is given to load the target's own parameters. Must suit what
        the target actually is: ``"leaprc.protein.ff19SB"`` for a protein,
        ``"leaprc.gaff2"`` for a small organic molecule.
    parameterized : bool, default=False
        Whether `force_field_ligand` already covers this molecule. True skips
        ``antechamber`` and ``parmchk2`` entirely and hands the PDB file
        straight to LEaP, which is much faster and is the right choice for a
        protein. False works the parameters out first, which is the only
        option for a molecule no standard force field describes.
    output_dir : str or pathlib.Path, optional
        Where to put the results. Defaults to the directory `file_path` is in.

    Returns
    -------
    int
        How many atoms the residue has, counted from the structure LEaP built.
        This is the number MAWS needs to know which stretch of a design's atom
        array belongs to this molecule, and it may differ from the atom count
        of `file_path` when hydrogens were added or removed along the way.

    Raises
    ------
    maws.errors.ToolchainError
        If AmberTools is not installed, or one of its three programs exits with
        an error.
    FileNotFoundError
        If ``tleap`` exits cleanly but writes no ``.lib`` file. It reports many
        problems by printing a complaint and carrying on, so a clean exit is
        not on its own evidence that anything was produced.

    See Also
    --------
    toggle_hydrogens : Adds or removes hydrogens before this is run.
    maws.io.tools.run_tool : Starts each program and reports its failures.

    Notes
    -----
    .. note::
        When a PDB file's parameters are being worked out, its ``CONECT``
        records — the lines stating which atoms are bonded — are stripped
        before ``antechamber`` sees it, and the bonding is worked out from the
        atom positions instead. Programs that write PDB files often encode a
        double or aromatic bond as a repeated ``CONECT`` entry, which
        ``antechamber`` reads as two separate bonds between the same pair of
        atoms, and LEaP then refuses to build the result.

    Examples
    --------
    A protein target, whose parameters the force field already covers:

    >>> n_atoms = make_lib(  # doctest: +SKIP
    ...     "target.pdb",
    ...     "LIG",
    ...     force_field_ligand="leaprc.protein.ff19SB",
    ...     parameterized=True,
    ... )
    >>> n_atoms  # doctest: +SKIP
    1462

    A small organic molecule, whose parameters have to be worked out:

    >>> make_lib(  # doctest: +SKIP
    ...     "ligand.mol2",
    ...     "LIG",
    ...     force_field_ligand="leaprc.gaff2",
    ...     atom_type="gaff2",
    ...     parameterized=False,
    ... )
    42
    """
    src = Path(file_path).resolve()
    name, ext = src.stem, src.suffix[1:].lower()
    base = Path(output_dir) if output_dir is not None else src.parent
    base.mkdir(parents=True, exist_ok=True)
    out_base = base / residue_name

    with tempfile.TemporaryDirectory() as td:
        w = Path(td)

        if not parameterized:
            # Match antechamber/parmchk2 to the ligand force field. maws sources
            # leaprc.gaff2 for organics, so atoms must be typed with gaff2; the
            # default "gaff" produces an empty frcmod and missing parameters,
            # making saveamberparm fail.
            ante_at = "gaff2" if "gaff2" in force_field_ligand else atom_type
            parmchk_set = "gaff2" if ante_at == "gaff2" else "gaff"

            # antechamber input. Strip PDB CONECT records first: tools like RDKit
            # encode double/aromatic bonds as repeated CONECT entries, which
            # antechamber turns into duplicate bonds that LEaP rejects ("cannot
            # add bond"). Let antechamber perceive connectivity from geometry.
            ante_in = str(src)
            if ext == "pdb":
                stripped = w / f"{name}_noconect.pdb"
                stripped.write_text(
                    "".join(
                        line
                        for line in src.read_text().splitlines(keepends=True)
                        if not line.startswith("CONECT")
                    )
                )
                ante_in = str(stripped)

            # -fi is taken from the input's own suffix rather than assumed, so
            # a .mol2 or .sdf target works without a separate code path.
            run_tool(
                [
                    "antechamber",
                    "-i",
                    ante_in,
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
                    ante_at,
                ],
                cwd=w,
            )

            # -s must name the same atom-type family antechamber just used, or
            # parmchk2 fills the gaps from the wrong table.
            run_tool(
                [
                    "parmchk2",
                    "-i",
                    f"{name}.mol2",
                    "-f",
                    "mol2",
                    "-o",
                    f"{residue_name}.frcmod",
                    "-s",
                    parmchk_set,
                ],
                cwd=w,
            )
        else:
            # tleap runs with the temporary directory as its working directory,
            # so every file it is told to load has to be sitting in there.
            shutil.copy2(src, w / f"{name}.pdb")

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

        # Both joining points or neither: a residue with only one end defined
        # cannot be chained, so setting one alone would only mislead.
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
        run_tool(["tleap", "-f", "leap.in"], cwd=w)

        # Only these one or two files are worth keeping; the temporary
        # directory and everything else in it disappears on the way out.
        shutil.move(w / f"{residue_name}.lib", out_base.with_suffix(".lib"))
        if not parameterized:
            shutil.move(w / f"{residue_name}.frcmod", out_base.with_suffix(".frcmod"))

        # Count atoms from what LEaP actually built rather than from the input
        # file: LEaP may add or drop hydrogens while assembling the residue.
        pdb = app.PDBFile(str(w / f"{residue_name}_tmp.pdb"))
        length = sum(1 for _ in pdb.topology.atoms())

    return length


def toggle_hydrogens(path: str, add: bool = True, ph: float = 7.0) -> None:
    """toggle_hydrogens(path, add=True, ph=7.0) -> None

    Add every missing hydrogen atom to a PDB file, or remove them all.

    Structures determined experimentally often have no hydrogens in them at
    all, because the method used cannot see atoms that light. A force field
    needs them, so they have to be put back before the molecule can be scored.
    Which ones belong depends on how acidic the surroundings are, which is what
    `ph` says.

    Parameters
    ----------
    path : str
        The PDB file to rewrite. It is overwritten in place, so keep a copy if
        the original is wanted afterwards.
    add : bool, default=True
        True puts the missing hydrogens in; False takes every hydrogen out.
        Removing them is useful when a file's hydrogens are unreliable and it
        is easier to strip them and start again.
    ph : float, default=7.0
        How acidic the surroundings are taken to be. This decides which of the
        chemical groups that can gain or lose a hydrogen actually carry one, so
        it changes the molecule's charge as well as its atom count. 7.0 is
        neutral, close to the conditions inside a cell. Not used when `add` is
        False.

    Returns
    -------
    None
        Nothing. The file named by `path` is rewritten.

    Raises
    ------
    OSError
        If `path` cannot be read or cannot be written back.

    See Also
    --------
    make_lib : What is normally run next, once the hydrogens are settled.

    Notes
    -----
    .. warning::
        Which hydrogens can be added is decided from the ``amber19-all.xml``
        parameters that ship with OpenMM, and those describe proteins and
        nucleic acids only. A residue they do not cover keeps whatever
        hydrogens it already had.

    Examples
    --------
    >>> from openmm.app import PDBFile  # doctest: +SKIP
    >>> toggle_hydrogens("target.pdb")  # doctest: +SKIP
    >>> sum(  # doctest: +SKIP
    ...     1
    ...     for atom in PDBFile("target.pdb").topology.atoms()
    ...     if atom.element.symbol == "H"
    ... )
    1204
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

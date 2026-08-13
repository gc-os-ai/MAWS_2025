"""
maws.io.pdb_cleaner
===================

Tidying a downloaded structure file so the modelling program can read it.

A PDB file is the standard way a measured molecular structure is published: one
line per atom, each carrying the atom's name, which residue it belongs to,
which chain, and where it sits in ångström. Files as published record what the
experiment saw, which is not always a single, complete, unambiguous molecule.
LEaP, the program that builds structures here, wants exactly that, and refuses
files that are not.

Four things in a published file get in its way:

*several models*
    A structure measured by NMR is published as twenty or so complete copies of
    the molecule, one per plausible solution. Read all at once they become
    twenty molecules sitting on top of each other.
*alternate positions*
    Where the experiment could not decide between two positions for a side
    chain, both are recorded, marked ``A`` and ``B`` and each given the
    fraction of the time it is thought to be there. Read at once they become
    one residue with two side chains.
*waters and other extras*
    Crystallographic waters, buffer molecules and ions the experiment happened
    to catch. Usually unwanted, and there are often hundreds of them.
*several chains*
    A structure may hold a whole complex when only one protein of it is the
    target.

:func:`clean_pdb_file` deals with those four and nothing else, and
:func:`resolve_pdb_path` is the thin wrapper the rest of MAWS calls.

What it will not do
-------------------
Silence matters more here than anywhere else in MAWS, because a damaged
structure still builds, still scores, and still produces a strand. So:

*nothing is reordered*
    LEaP works out which atoms are bonded from the order the lines come in.
    Moving a line rewrites the molecule's chemistry. Every step here keeps the
    file in the order it arrived.
*no residue is deleted from the middle of a chain*
    Deleting one joins its neighbours to each other, inventing a bond that is
    not there. Residues a standard force field does not know by name are
    renamed to the nearest one it does, never dropped.
*no chain is silently joined to another*
    The ``TER`` line that separates two chains is always written.
*nothing fails quietly*
    :func:`clean_pdb_file` returns a count of everything it removed, and
    refuses outright to write a file that has lost more than it should have.

See Also
--------
maws.io.prepare.make_lib : What the cleaned file is handed to next.

Examples
--------
>>> from pathlib import Path
>>> import tempfile
>>> pdb = (
...     "ATOM      1  N   ALA A   1      11.10  22.20  33.30  1.00  0.00  N\\n"
...     "HETATM    2  O   HOH A   2      44.40  55.50  66.60  1.00  0.00  O\\n"
...     "END\\n"
... )
>>> with tempfile.TemporaryDirectory() as tmp:
...     source = Path(tmp) / "target.pdb"
...     _ = source.write_text(pdb)
...     cleaned, report = clean_pdb_file(source, Path(tmp) / "clean.pdb")
...     report.waters_removed, report.atoms_kept
(1, 1)
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from maws.errors import ConfigurationError

__all__ = [
    "CleaningReport",
    "clean_pdb_file",
    "resolve_pdb_path",
]

WATER_NAMES = frozenset({"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3"})
"""Residue names meaning "a water molecule the experiment happened to see".

Several names for the same thing are in circulation, because different
programs write different ones.
"""

RENAMED_RESIDUES = {
    "MSE": "MET",
    "SEC": "CYS",
    "PYL": "LYS",
}
"""Residues to rename rather than remove, and what to call them instead.

These are the standard amino acids with one atom swapped for a heavier one of
the same group, used to make a crystal structure easier to solve. A protein
force field has no entry for them, so LEaP would refuse the file. Renaming
them to the residue they stand in for keeps the chain whole and changes one
atom of it; deleting them would tear a hole in the middle of the chain, and
LEaP would then bond the two residues either side of the hole together.
"""

RENAMED_ATOMS = {
    ("MSE", "SE"): ("SD", "S"),
    ("SEC", "SE"): ("SG", "S"),
}
"""The one atom each renamed residue carries that its stand-in does not.

Keyed by the residue's *original* name and the atom's name; the value is the
atom name and element symbol to write instead. Selenium takes the place of the
sulfur in methionine and in cysteine, so it becomes that sulfur.
"""

MINIMUM_KEPT = 0.5
"""How much of a structure must survive cleaning, as a fraction of its atoms.

Not counting waters, hydrogens and chains that were deliberately dropped —
only what was meant to be kept. Anything below this means something went wrong
that was not asked for, and the file is refused rather than handed on.
"""


@dataclass(frozen=True, slots=True)
class CleaningReport:
    """What cleaning a file removed, counted.

    Parameters
    ----------
    atoms_read : int
        How many atom lines the original file held, across all its models.
    atoms_kept : int
        How many were written out.
    models_dropped : int
        How many whole copies of the molecule were left behind. Zero for an
        ordinary crystal structure; nineteen or so for an NMR one.
    waters_removed : int
        How many water molecules' atoms were dropped.
    heteroatoms_removed : int
        How many atoms were dropped for being marked ``HETATM``, not counting
        waters. Zero unless that was asked for.
    hydrogens_removed : int
        How many hydrogen atoms were dropped. Zero unless that was asked for.
    chains_dropped : tuple of str
        Which chains were left out, by their letter.
    conformers_resolved : int
        How many residues had more than one recorded position, and so had one
        of them chosen.
    residues_renamed : int
        How many residues were renamed to a standard one, as described by
        :data:`RENAMED_RESIDUES`.

    See Also
    --------
    clean_pdb_file : Produces one of these.

    Notes
    -----
    The counts do not have to add up to ``atoms_read - atoms_kept``: an atom
    can be dropped for more than one reason, and only the first that applies is
    counted.
    """

    atoms_read: int
    atoms_kept: int
    models_dropped: int
    waters_removed: int
    heteroatoms_removed: int
    hydrogens_removed: int
    chains_dropped: tuple[str, ...]
    conformers_resolved: int
    residues_renamed: int

    def summary(self) -> str:
        """Return a one-line description of what was removed.

        Returns
        -------
        str
            Suitable for a log line. Mentions only what actually happened, so
            a file that needed nothing doing to it says so.

        Examples
        --------
        >>> CleaningReport(10, 10, 0, 0, 0, 0, (), 0, 0).summary()
        'kept all 10 atoms'
        """
        parts = []
        if self.models_dropped:
            parts.append(f"{self.models_dropped} further models")
        if self.waters_removed:
            parts.append(f"{self.waters_removed} water atoms")
        if self.heteroatoms_removed:
            parts.append(f"{self.heteroatoms_removed} heteroatoms")
        if self.hydrogens_removed:
            parts.append(f"{self.hydrogens_removed} hydrogens")
        if self.chains_dropped:
            parts.append(f"chains {', '.join(self.chains_dropped)}")
        if self.conformers_resolved:
            parts.append(f"{self.conformers_resolved} alternate positions")
        if self.residues_renamed:
            parts.append(f"{self.residues_renamed} residues renamed")
        if not parts:
            return f"kept all {self.atoms_kept} atoms"
        return (
            f"kept {self.atoms_kept} of {self.atoms_read} atoms; "
            f"removed {', '.join(parts)}"
        )


@dataclass(frozen=True, slots=True)
class _Atom:
    """One atom line of a PDB file, with the fields this module needs.

    Parameters
    ----------
    record : {"ATOM", "HETATM"}
        Which kind of line it came from. ``HETATM`` marks anything that is not
        part of a standard protein or nucleic acid chain.
    name : str
        The atom's name within its residue, e.g. ``"CA"``.
    alt_loc : str
        Which recorded position this is, when the experiment found more than
        one. A single space means there was only one.
    residue : str
        The residue's name, e.g. ``"ALA"``.
    chain : str
        Which chain the residue belongs to. Case matters: ``"a"`` and ``"A"``
        are different chains.
    seq_num : str
        The residue's number within its chain, as written. Kept as text
        because published files are not consistent about padding it.
    insertion : str
        A letter distinguishing residues inserted between two numbered ones. A
        single space when there is none.
    occupancy : float
        What fraction of the time the atom is thought to be here, between 0
        and 1. Only meaningful alongside `alt_loc`.
    element : str
        The chemical symbol, right-aligned in two characters as the format
        requires. May be blank in older files.
    line : str
        The whole original line, padded to 80 characters. Written back out
        with only the fields this module changes replaced, so nothing else
        about the atom is disturbed.
    """

    record: str
    name: str
    alt_loc: str
    residue: str
    chain: str
    seq_num: str
    insertion: str
    occupancy: float
    element: str
    line: str

    @property
    def residue_key(self) -> tuple[str, str, str]:
        """tuple of str : What identifies the residue this atom belongs to."""
        return (self.chain, self.seq_num, self.insertion)

    @property
    def is_hydrogen(self) -> bool:
        """bool : Whether this is a hydrogen or a deuterium.

        Read from the element column when the file fills it in. Older files
        leave it blank, and then it is read from the atom's name instead —
        which starts with H or D for a hydrogen, except that the name may be
        shifted one character right when it belongs to a heavy atom with a
        two-letter symbol.
        """
        if self.element.strip():
            return self.element.strip().upper() in ("H", "D")
        stripped = self.name.strip()
        if not stripped:
            return False
        first = stripped[0]
        if first.isdigit():
            stripped = stripped[1:]
        return bool(stripped) and stripped[0] in ("H", "D")

    def rewritten(
        self, *, alt_loc: str, seq_num: int, serial: int, name: str, element: str
    ) -> str:
        """Return this atom's line with the numbered fields replaced.

        Parameters
        ----------
        alt_loc : str
            What to put in the alternate-position column, one character.
        seq_num : int
            The residue number to write.
        serial : int
            The atom number to write.
        name : str
            The atom name to write, already positioned within its four
            characters.
        element : str
            The chemical symbol to write, right-aligned in two characters.

        Returns
        -------
        str
            An 80-character line, with a newline.

        Notes
        -----
        Every other column is carried across from the original untouched,
        including the positions, so cleaning a file never moves an atom.

        The insertion-code column is always cleared, because the residues are
        renumbered so that nothing needs it. Leaving it filled in would give
        two different residues the same number and letter.
        """
        line = self.line
        return (
            f"{line[:6]}{serial:>5}{line[11]}{name:<4}{alt_loc}{line[17:21]}"
            f"{line[21]}{seq_num:>4} {line[27:76]}{element:>2}{line[78:80]}\n"
        )


def _parse_atom(line: str) -> _Atom:
    """Read one ATOM or HETATM line into its fields.

    Parameters
    ----------
    line : str
        The line, with or without its newline.

    Returns
    -------
    _Atom
        Its fields, and the line itself padded to 80 characters.

    Notes
    -----
    The PDB format is fixed-column, so every field is read by its position
    rather than by splitting on spaces. Splitting would go wrong the moment two
    columns ran together, which happens for four-character atom names and for
    coordinates past 100 ångström.

    An occupancy that is missing or unreadable is taken as 1.0, which is what a
    file means when it does not say.
    """
    padded = line.rstrip("\n").ljust(80)
    try:
        occupancy = float(padded[54:60])
    except ValueError:
        occupancy = 1.0
    return _Atom(
        record=padded[:6].strip(),
        name=padded[12:16],
        alt_loc=padded[16],
        residue=padded[17:20].strip(),
        chain=padded[21],
        seq_num=padded[22:26].strip(),
        insertion=padded[26],
        occupancy=occupancy,
        element=padded[76:78],
        line=padded,
    )


def _first_model(lines: list[str]) -> tuple[list[str], int]:
    """Return the lines of the first model, and how many were left behind.

    Parameters
    ----------
    lines : list of str
        Every line of the file.

    Returns
    -------
    kept : list of str
        The lines up to the end of the first model, or all of them when the
        file has no models marked.
    dropped : int
        How many further models there were.

    Notes
    -----
    A structure measured by NMR is published as a series of complete copies of
    the molecule, each wrapped in ``MODEL`` and ``ENDMDL`` lines and each a
    plausible answer. They are not different parts of one molecule, so reading
    them all gives twenty molecules sitting on top of each other, with the same
    residue numbers repeated twenty times over.
    """
    if not any(line.startswith("MODEL ") for line in lines):
        return lines, 0

    kept: list[str] = []
    depth = 0
    finished = False
    dropped = 0
    for line in lines:
        if line.startswith("MODEL "):
            depth += 1
            if depth > 1:
                dropped += 1
            continue
        if line.startswith("ENDMDL"):
            if depth == 1:
                finished = True
            continue
        if not finished:
            kept.append(line)
    return kept, dropped


def _chosen_conformers(atoms: list[_Atom]) -> dict[tuple[str, str, str], str]:
    """Pick one recorded position per residue that has more than one.

    Parameters
    ----------
    atoms : list of _Atom
        Every atom of the structure.

    Returns
    -------
    dict
        Residue key to the alternate-position letter to keep. Residues with
        only one recorded position are left out.

    Notes
    -----
    The choice is made once for a whole residue, not atom by atom. Choosing per
    atom would take the backbone from one recorded position and part of the
    side chain from the other, giving a shape that was never observed and
    typically one where the side chain runs into itself. It can also leave one
    residue number carrying two different residue names, when the experiment
    could not decide which amino acid was there.

    The position kept is the one whose atoms carry the most occupancy in total
    — that is, the one the experiment thought was there most of the time.
    """
    weights: dict[tuple[str, str, str], Counter] = {}
    for atom in atoms:
        if atom.alt_loc == " ":
            continue
        weights.setdefault(atom.residue_key, Counter())[atom.alt_loc] += atom.occupancy
    return {
        key: max(sorted(counter), key=lambda letter: counter[letter])
        for key, counter in weights.items()
    }


def _selected_chains(atoms: list[_Atom], keep: str) -> set[str]:
    """Work out which chains to keep.

    Parameters
    ----------
    atoms : list of _Atom
        Every atom of the structure.
    keep : str
        ``"all"`` for every chain; ``"one"`` for the longest; otherwise a
        string of chain letters to keep, e.g. ``"AB"``.

    Returns
    -------
    set of str
        The chain letters to keep.

    Raises
    ------
    maws.errors.ConfigurationError
        If `keep` names no chain the file actually holds. Left as an error
        rather than a warning because it would otherwise produce an empty
        file.

    Notes
    -----
    Chain letters are compared exactly. Published structures of large
    assemblies do use both ``A`` and ``a`` for different chains, so treating
    them as the same would quietly keep two chains where one was asked for.

    "The longest" is measured in residues rather than atoms. Counting atoms
    lets a short chain win by having more of them, which happens whenever one
    chain is a nucleic acid and another a protein.
    """
    present = {atom.chain for atom in atoms}
    if keep == "all":
        return present
    if keep == "one":
        sizes: dict[str, set] = {}
        for atom in atoms:
            sizes.setdefault(atom.chain, set()).add(atom.residue_key)
        longest = max(sorted(sizes), key=lambda chain: len(sizes[chain]))
        return {longest}

    wanted = set(keep)
    usable = wanted & present
    if not usable:
        raise ConfigurationError(
            f"keep_chains={keep!r} names no chain in this file. It holds "
            f"{', '.join(sorted(present))}. Chain letters are case-sensitive."
        )
    return usable


def clean_pdb_file(
    source: str | Path,
    destination: str | Path,
    *,
    keep_chains: str = "all",
    remove_hydrogens: bool = False,
    drop_heteroatoms: bool = False,
) -> tuple[Path, CleaningReport]:
    """Write a tidied copy of a PDB file, and say what was removed.

    Parameters
    ----------
    source : str or pathlib.Path
        The file to read.
    destination : str or pathlib.Path
        Where to write the tidied copy. Overwritten if it exists.
    keep_chains : str, default="all"
        Which chains to keep. ``"all"`` keeps every one; ``"one"`` keeps the
        chain with the most residues and drops the rest; anything else is read
        as the letters to keep, so ``"AB"`` keeps chains A and B. Letters are
        matched exactly, upper and lower case being different chains.
    remove_hydrogens : bool, default=False
        Whether to drop hydrogen and deuterium atoms. Worth doing when the
        hydrogens are to be added back at the right places for the force
        field being used, which is not always where the experiment put them.
    drop_heteroatoms : bool, default=False
        Whether to drop everything marked ``HETATM``: bound ligands,
        cofactors, metal ions and sugars. Waters are dropped either way. Take
        care with this one — a metal ion held inside a protein is part of what
        holds its shape, and removing it changes both the charge and the fold.

    Returns
    -------
    path : pathlib.Path
        Where the tidied copy was written, the same as `destination`.
    report : CleaningReport
        What was removed, counted.

    Raises
    ------
    maws.errors.ConfigurationError
        If the file holds no atoms, if `keep_chains` names no chain it holds,
        or if less than :data:`MINIMUM_KEPT` of what was meant to be kept
        survived. The last of these means something went wrong that was not
        asked for, and passing the file on would hide it.

    See Also
    --------
    resolve_pdb_path : The wrapper the rest of MAWS calls.
    CleaningReport : What comes back alongside the file.

    Notes
    -----
    The lines are written in the order they were read. Nothing is sorted at any
    point, because LEaP decides which atoms are bonded from that order.

    Residues are renumbered from 1 within each chain and the insertion-code
    column is cleared. That is what makes it safe to keep residues carrying an
    insertion code — the letters that told two residues with the same number
    apart are no longer needed once they have different numbers. The
    alternative, dropping those residues, would leave a hole in the middle of a
    chain and LEaP would join the residues either side of it.

    A ``TER`` line is written at the end of every chain, whether or not the
    original had one. Without it LEaP joins two chains into a single molecule.

    Examples
    --------
    >>> from pathlib import Path
    >>> import tempfile
    >>> pdb = (
    ...     "ATOM      1  N   ALA A   1      11.10  22.20  33.30  1.00  0.00  N\\n"
    ...     "ATOM      2  N   GLY B   1      44.40  55.50  66.60  1.00  0.00  N\\n"
    ...     "END\\n"
    ... )
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     source = Path(tmp) / "target.pdb"
    ...     _ = source.write_text(pdb)
    ...     _, report = clean_pdb_file(source, Path(tmp) / "clean.pdb", keep_chains="A")
    ...     report.chains_dropped, report.atoms_kept
    (('B',), 1)
    """
    source = Path(source)
    destination = Path(destination)

    lines, models_dropped = _first_model(source.read_text().splitlines())
    atoms = [
        _parse_atom(line) for line in lines if line.startswith(("ATOM  ", "HETATM"))
    ]
    if not atoms:
        raise ConfigurationError(
            f"{source} holds no atom records. It may be a header-only entry, "
            f"or not a PDB file at all."
        )

    chains = _selected_chains(atoms, keep_chains)
    conformers = _chosen_conformers(atoms)

    kept: list[_Atom] = []
    removed = Counter()
    for atom in atoms:
        if atom.residue in WATER_NAMES:
            removed["water"] += 1
        elif atom.chain not in chains:
            removed["chain"] += 1
        elif drop_heteroatoms and atom.record == "HETATM":
            removed["hetatm"] += 1
        elif remove_hydrogens and atom.is_hydrogen:
            removed["hydrogen"] += 1
        elif atom.alt_loc not in (" ", conformers.get(atom.residue_key, " ")):
            removed["conformer"] += 1
        else:
            kept.append(atom)

    meant_to_keep = len(atoms) - sum(
        removed[reason] for reason in ("water", "chain", "hetatm", "hydrogen")
    )
    if len(kept) < MINIMUM_KEPT * meant_to_keep:
        raise ConfigurationError(
            f"cleaning {source} would keep only {len(kept)} of the "
            f"{meant_to_keep} atoms that were meant to survive it. Something "
            f"is wrong with the file rather than with the settings."
        )

    renamed = _write_cleaned(destination, kept)
    return destination, CleaningReport(
        atoms_read=len(atoms),
        atoms_kept=len(kept),
        models_dropped=models_dropped,
        waters_removed=removed["water"],
        heteroatoms_removed=removed["hetatm"],
        hydrogens_removed=removed["hydrogen"],
        chains_dropped=tuple(sorted({atom.chain for atom in atoms} - chains)),
        conformers_resolved=len(conformers),
        residues_renamed=renamed,
    )


def _write_cleaned(destination: Path, atoms: list[_Atom]) -> int:
    """Write the atoms out as a PDB file, renumbered and separated by chain.

    Parameters
    ----------
    destination : pathlib.Path
        Where to write. Overwritten if it exists.
    atoms : list of _Atom
        The atoms to write, already in the order they should appear.

    Returns
    -------
    int
        How many residues were renamed to a standard one on the way out.

    Notes
    -----
    Atom serials are renumbered from 1 and residues from 1 within each chain,
    so that nothing is repeated and nothing needs an insertion code. A ``TER``
    line closes each chain and takes the next serial after the atom before it,
    as the format requires.
    """
    out: list[str] = []
    renamed_residues: set[tuple[str, str, str]] = set()
    serial = 0
    residue_number = 0
    previous_residue: tuple[str, str, str] | None = None
    previous_chain: str | None = None
    previous_atom: _Atom | None = None

    for atom in atoms:
        if previous_chain is not None and atom.chain != previous_chain:
            serial += 1
            out.append(_ter_line(previous_atom, serial, residue_number))
            residue_number = 0
            previous_residue = None
        if atom.residue_key != previous_residue:
            residue_number += 1
            previous_residue = atom.residue_key

        name, element = RENAMED_ATOMS.get(
            (atom.residue, atom.name.strip()), (atom.name, atom.element)
        )
        if name is not atom.name:
            name = f" {name:<3}"
        if atom.residue in RENAMED_RESIDUES:
            renamed_residues.add(atom.residue_key)

        serial += 1
        line = atom.rewritten(
            alt_loc=" ",
            seq_num=residue_number,
            serial=serial,
            name=name,
            element=element,
        )
        replacement = RENAMED_RESIDUES.get(atom.residue)
        if replacement is not None:
            line = f"{line[:17]}{replacement:>3}{line[20:]}"
        out.append(line)
        previous_chain = atom.chain
        previous_atom = atom

    if previous_atom is not None:
        serial += 1
        out.append(_ter_line(previous_atom, serial, residue_number))
    out.append("END\n")

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(out))
    return len(renamed_residues)


def _ter_line(atom: _Atom, serial: int, residue_number: int) -> str:
    """Return the ``TER`` line closing the chain that `atom` ends.

    Parameters
    ----------
    atom : _Atom
        The last atom of the chain.
    serial : int
        The atom number to write, which is one past that atom's.
    residue_number : int
        The renumbered residue the chain ends on.

    Returns
    -------
    str
        An 80-character line, with a newline.
    """
    residue = RENAMED_RESIDUES.get(atom.residue, atom.residue)
    return (
        f"TER   {serial:>5}      {residue:>3} {atom.chain}{residue_number:>4}"
    ).ljust(80) + "\n"


def resolve_pdb_path(
    pdb_path: str | Path,
    molecule_type: str,
    *,
    clean_pdb: bool,
    keep_chains: str = "all",
    remove_h: bool = False,
    drop_hetatm: bool = False,
    logger: logging.Logger | None = None,
) -> tuple[str, str]:
    """Return the target file to use, cleaning it first if asked.

    Parameters
    ----------
    pdb_path : str or pathlib.Path
        The target file as the caller gave it.
    molecule_type : str
        What the target is: ``"protein"``, ``"organic"`` or ``"lipid"``.
    clean_pdb : bool
        Whether to tidy the file at all. When False the file is used exactly
        as given.
    keep_chains : str, default="all"
        Which chains to keep. See :func:`clean_pdb_file`.
    remove_h : bool, default=False
        Whether to drop hydrogen atoms. See :func:`clean_pdb_file`.
    drop_hetatm : bool, default=False
        Whether to drop bound ligands, cofactors and ions. See
        :func:`clean_pdb_file`.
    logger : logging.Logger, optional
        Where to report what was removed. Nothing is reported when this is
        left out.

    Returns
    -------
    path : str
        The file to hand on, which is the tidied copy when one was made.
    original : str
        The file as given, so a caller can report both.

    Raises
    ------
    maws.errors.ConfigurationError
        If cleaning was asked for and could not be done. It is not attempted
        and then abandoned: a run that believes its target was tidied and is
        quietly working from the untidied file is worse than one that stops.

    See Also
    --------
    clean_pdb_file : Does the work.

    Notes
    -----
    Cleaning is skipped for anything that is not a protein. A small organic
    molecule is a single residue with a name of its own, so the tidying here
    has nothing to do to it, and dropping its heteroatoms would drop the
    molecule.
    """
    original = str(pdb_path)
    if not clean_pdb:
        return original, original
    if molecule_type != "protein":
        if logger is not None:
            logger.info(
                "not cleaning %s: cleaning applies to proteins, and this target is %s",
                original,
                molecule_type,
            )
        return original, original

    source = Path(pdb_path)
    destination = source.with_name(f"{source.stem}_cleaned.pdb")
    cleaned, report = clean_pdb_file(
        source,
        destination,
        keep_chains=keep_chains,
        remove_hydrogens=remove_h,
        drop_heteroatoms=drop_hetatm,
    )
    if logger is not None:
        logger.info("cleaned %s -> %s: %s", original, cleaned, report.summary())
    return str(cleaned), original

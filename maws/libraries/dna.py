"""
maws.libraries.dna
==================

DNA residue tables, and the ``dna()`` factory that validates them.

A DNA strand is a chain of *residues* — single nucleotides, written ``G``,
``A``, ``T`` and ``C``. Its two ends are chemically different and are called 5'
and 3'; a sequence is written 5'-to-3', left to right. A residue in the middle
of a strand, at one of its ends, or on its own is a slightly different molecule
in each case, so each written nucleotide corresponds to four residue names.
Those names all begin with ``D``, so the token ``G`` maps to ``DGN``, ``DG5``,
``DG`` or ``DG3``.

The tables here describe each of those sixteen residues: how many atoms it has,
which of its bonds can be turned, and which atoms join it to its neighbours.
The numbers are atom positions within one residue, counted from its first atom,
and they follow the order in which LEaP lays the atoms out. LEaP is the
molecular modelling program that builds the actual structure, run as the
command ``tleap``.

This module holds data and nothing else: no file access, no OpenMM, no LEaP.
Everything here can be read and checked without anything installed.

Reading the tables
------------------
``RESIDUE_NAMES`` and ``RESIDUE_LENGTH`` are positionally aligned: entry *i* of
one describes the same residue as entry *i* of the other. ``ROTATIONS`` and
``ALIASES`` are keyed by name instead. That mix is exactly why
:meth:`maws.values.ResidueLibrary.from_tables` checks the alignment once, on
load, rather than letting a misaligned row surface later as a LEaP complaint
about a missing parameter.

Rotation rows are ``(residue, start, bond, end_or_None)``, describing a
*torsion*: a bond that can be turned, swinging everything on one side of it
around while the other side stays put. ``start`` and ``bond`` are the two atoms
of that bond. An ``end_or_None`` of ``None`` means "rotate the rest of the
chain" and becomes a :class:`~maws.values.BackboneTorsion`; a number bounds the
moving atoms inside the residue and becomes a
:class:`~maws.values.LocalTorsion`. Negative indices count back from the
residue's last atom and are resolved on load.

Alias rows are ``[token, alone, start, middle, end]``: the four residue names
one written token takes, depending on where it sits in the chain.

See Also
--------
maws.libraries.rna : The same tables for RNA.
maws.values.ResidueLibrary : What the tables are turned into.

Examples
--------
>>> dna()["DG"].n_atoms
33
>>> dna().alias("T").end
'DT3'
"""

from functools import cache

from maws.values import ResidueLibrary

__all__ = [
    "ALIASES",
    "CONNECT",
    "RESIDUE_LENGTH",
    "RESIDUE_NAMES",
    "ROTATIONS",
    "dna",
]

# The sixteen residue names LEaP knows for DNA: four nucleotides in each of
# four positions — alone (N suffix), middle (no suffix), 5' end, 3' end.
# Positionally aligned with RESIDUE_LENGTH.
RESIDUE_NAMES: list[str] = [
    "DGN",
    "DAN",
    "DTN",
    "DCN",
    "DG",
    "DA",
    "DT",
    "DC",
    "DG5",
    "DA5",
    "DT5",
    "DC5",
    "DG3",
    "DA3",
    "DT3",
    "DC3",
]

# How many atoms each residue has, entry i describing RESIDUE_NAMES[i]. The
# counts differ between positions because an end residue carries an extra group
# that a middle one does not.
RESIDUE_LENGTH: list[int] = [
    32,
    31,
    31,
    29,
    33,
    32,
    32,
    30,
    31,
    30,
    30,
    28,
    34,
    33,
    33,
    31,
]

# Turnable bonds, four per residue: (residue, start, bond, end_or_None), where
# start and bond are the two atoms of the bond and end_or_None bounds the atoms
# that swing — None meaning the rest of the strand. Indices count from the
# residue's own first atom, negative ones back from its last.
ROTATIONS: list[tuple[str, int, int, int | None]] = [
    ("DGN", 0, 1, None),
    ("DGN", 1, 2, None),
    ("DGN", 8, 10, -7),
    ("DGN", -4, -2, None),
    ("DAN", 0, 1, None),
    ("DAN", 1, 2, None),
    ("DAN", 8, 10, -7),
    ("DAN", -4, -2, None),
    ("DTN", 0, 1, None),
    ("DTN", 1, 2, None),
    ("DTN", 8, 10, -7),
    ("DTN", -4, -2, None),
    ("DCN", 0, 1, None),
    ("DCN", 1, 2, None),
    ("DCN", 8, 10, -7),
    ("DCN", -4, -2, None),
    ("DG", 0, 3, None),
    ("DG", 3, 4, None),
    ("DG", 10, 12, -6),
    ("DG", -6, -1, None),
    ("DA", 0, 3, None),
    ("DA", 3, 4, None),
    ("DA", 10, 12, -6),
    ("DA", -6, -1, None),
    ("DT", 0, 3, None),
    ("DT", 3, 4, None),
    ("DT", 10, 12, -6),
    ("DT", -6, -1, None),
    ("DC", 0, 3, None),
    ("DC", 3, 4, None),
    ("DC", 10, 12, -6),
    ("DC", -6, -1, None),
    ("DG5", 0, 1, None),
    ("DG5", 1, 2, None),
    ("DG5", 8, 10, -6),
    ("DG5", -6, -1, None),
    ("DA5", 0, 1, None),
    ("DA5", 1, 2, None),
    ("DA5", 8, 10, -6),
    ("DA5", -6, -1, None),
    ("DT5", 0, 1, None),
    ("DT5", 1, 2, None),
    ("DT5", 8, 10, -6),
    ("DT5", -6, -1, None),
    ("DC5", 0, 1, None),
    ("DC5", 1, 2, None),
    ("DC5", 8, 10, -6),
    ("DC5", -6, -1, None),
    ("DG3", 0, 3, None),
    ("DG3", 3, 4, None),
    ("DG3", 10, 12, -7),
    ("DG3", -4, -2, None),
    ("DA3", 0, 3, None),
    ("DA3", 3, 4, None),
    ("DA3", 10, 12, -7),
    ("DA3", -4, -2, None),
    ("DT3", 0, 3, None),
    ("DT3", 3, 4, None),
    ("DT3", 10, 12, -7),
    ("DT3", -4, -2, None),
    ("DC3", 0, 3, None),
    ("DC3", 3, 4, None),
    ("DC3", 10, 12, -7),
    ("DC3", -4, -2, None),
]

# Which atoms bond this residue to its neighbours, one row per residue, aligned
# with RESIDUE_NAMES: [[own_head, other_head], [own_tail, other_tail],
# head_len, tail_len]. "head" is the join towards the 5' end and "tail" the one
# towards the 3' end; the two lengths are in ångström. Every DNA residue joins
# the same way, so one row is repeated rather than written out sixteen times.
CONNECT: list[list] = [[[0, -1], [-2, 0], 1.6, 1.6] for _ in RESIDUE_NAMES]

# What each written token becomes, by where it sits in the strand:
# [token, alone, start, middle, end]. A residue name is a valid token too, so
# a sequence may be written as G A T C or as DG DA DT DC.
ALIASES: list[list[str]] = [
    ["DCN", "DCN", "DC5", "DC", "DC3"],
    ["A", "DAN", "DA5", "DA", "DA3"],
    ["C", "DCN", "DC5", "DC", "DC3"],
    ["DTN", "DTN", "DT5", "DT", "DT3"],
    ["G", "DGN", "DG5", "DG", "DG3"],
    ["DG3", "DG3", "DG", "DG", "DG3"],
    ["DG", "DG", "DG", "DG", "DG"],
    ["DAN", "DAN", "DA5", "DA", "DA3"],
    ["DA3", "DA3", "DA", "DA", "DA3"],
    ["DGN", "DGN", "DG5", "DG", "DG3"],
    ["DC", "DC", "DC", "DC", "DC"],
    ["DA", "DA", "DA", "DA", "DA"],
    ["DA5", "DA5", "DA5", "DA", "DA"],
    ["T", "DTN", "DT5", "DT", "DT3"],
    ["DG5", "DG5", "DG5", "DG", "DG"],
    ["DT3", "DT3", "DT", "DT", "DT3"],
    ["DT", "DT", "DT", "DT", "DT"],
    ["DC5", "DC5", "DC5", "DC", "DC"],
    ["DC3", "DC3", "DC", "DC", "DC3"],
    ["DT5", "DT5", "DT5", "DT", "DT"],
]


@cache
def dna() -> ResidueLibrary:
    """Return the DNA residue library, built from the tables in this module.

    Everything MAWS needs to design a DNA aptamer, other than the physics: the
    sixteen residues, their atom counts, their turnable bonds and the four
    names each written nucleotide can take.

    Returns
    -------
    maws.values.ResidueLibrary
        DNA templates keyed by residue name, together with the alias table that
        maps written tokens (``G A T C``) onto them. The same object every
        time, so do not try to change it.

    Raises
    ------
    maws.errors.ConfigurationError
        If the tables in this module do not line up, or a turnable bond names
        an atom its residue does not have. That would be an error in this
        file, so it surfaces on the first call rather than during a design run.

    See Also
    --------
    maws.libraries.rna : The same thing for RNA.
    maws.values.ResidueLibrary.from_tables : Does the checking and building.

    Notes
    -----
    Cached, because the tables are constant and the library never changes, so
    one instance per process is safe to share and the checking is paid for
    once.

    Examples
    --------
    >>> lib = dna()
    >>> len(lib)
    16
    >>> lib["DGN"].n_atoms
    32

    A written ``T`` becomes a different residue at each position:

    >>> lib.alias("T")
    AliasSet(alone='DTN', start='DT5', middle='DT', end='DT3')

    The same object comes back every time:

    >>> dna() is dna()
    True
    """
    return ResidueLibrary.from_tables(
        names=RESIDUE_NAMES,
        lengths=RESIDUE_LENGTH,
        rotations=ROTATIONS,
        connect=CONNECT,
        aliases=ALIASES,
    )

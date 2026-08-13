"""
maws.libraries.rna
==================

RNA residue tables, and the ``rna()`` factory that validates them.

An RNA strand is a chain of *residues* — single nucleotides, written ``G``,
``A``, ``U`` and ``C``. Its two ends are chemically different and are called 5'
and 3'; a sequence is written 5'-to-3', left to right. A residue in the middle
of a strand, at one of its ends, or on its own is a slightly different molecule
in each case, so each written nucleotide corresponds to four residue names.

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
maws.libraries.dna : The same tables for DNA.
maws.values.ResidueLibrary : What the tables are turned into.

Examples
--------
>>> rna()["G"].n_atoms
34
>>> rna().alias("U").end
'U3'
"""

from functools import cache

from maws.values import ResidueLibrary

__all__ = [
    "ALIASES",
    "CONNECT",
    "RESIDUE_LENGTH",
    "RESIDUE_NAMES",
    "ROTATIONS",
    "rna",
]

# The sixteen residue names LEaP knows for RNA: four nucleotides in each of
# four positions — alone (N suffix), middle (no suffix), 5' end, 3' end.
# Positionally aligned with RESIDUE_LENGTH.
RESIDUE_NAMES: list[str] = [
    "GN",
    "AN",
    "UN",
    "CN",
    "G",
    "A",
    "U",
    "C",
    "G5",
    "A5",
    "U5",
    "C5",
    "G3",
    "A3",
    "U3",
    "C3",
]

# How many atoms each residue has, entry i describing RESIDUE_NAMES[i]. The
# counts differ between positions because an end residue carries an extra group
# that a middle one does not.
RESIDUE_LENGTH: list[int] = [
    33,
    32,
    29,
    30,
    34,
    33,
    30,
    31,
    32,
    31,
    28,
    29,
    35,
    34,
    31,
    32,
]

# Turnable bonds, four per residue: (residue, start, bond, end_or_None), where
# start and bond are the two atoms of the bond and end_or_None bounds the atoms
# that swing — None meaning the rest of the strand. Indices count from the
# residue's own first atom, negative ones back from its last.
ROTATIONS: list[tuple[str, int, int, int | None]] = [
    ("GN", 0, 1, None),
    ("GN", 1, 2, None),
    ("GN", 8, 10, -8),
    ("GN", -8, -2, None),
    ("AN", 0, 1, None),
    ("AN", 1, 2, None),
    ("AN", 8, 10, -8),
    ("AN", -8, -2, None),
    ("UN", 0, 1, None),
    ("UN", 1, 2, None),
    ("UN", 8, 10, -8),
    ("UN", -8, -2, None),
    ("CN", 0, 1, None),
    ("CN", 1, 2, None),
    ("CN", 8, 10, -8),
    ("CN", -8, -2, None),
    ("G", 0, 3, None),
    ("G", 3, 4, None),
    ("G", 10, 12, -7),
    ("G", -7, -1, None),
    ("A", 0, 3, None),
    ("A", 3, 4, None),
    ("A", 10, 12, -7),
    ("A", -7, -1, None),
    ("U", 0, 3, None),
    ("U", 3, 4, None),
    ("U", 10, 12, -7),
    ("U", -7, -1, None),
    ("C", 0, 3, None),
    ("C", 3, 4, None),
    ("C", 10, 12, -7),
    ("C", -7, -1, None),
    ("G5", 0, 1, None),
    ("G5", 1, 2, None),
    ("G5", 8, 10, -7),
    ("G5", -7, -1, None),
    ("A5", 0, 1, None),
    ("A5", 1, 2, None),
    ("A5", 8, 10, -7),
    ("A5", -7, -1, None),
    ("U5", 0, 1, None),
    ("U5", 1, 2, None),
    ("U5", 8, 10, -7),
    ("U5", -7, -1, None),
    ("C5", 0, 1, None),
    ("C5", 1, 2, None),
    ("C5", 8, 10, -7),
    ("C5", -7, -1, None),
    ("G3", 0, 3, None),
    ("G3", 3, 4, None),
    ("G3", 10, 12, -8),
    ("G3", -8, -2, None),
    ("A3", 0, 3, None),
    ("A3", 3, 4, None),
    ("A3", 10, 12, -8),
    ("A3", -8, -2, None),
    ("U3", 0, 3, None),
    ("U3", 3, 4, None),
    ("U3", 10, 12, -8),
    ("U3", -8, -2, None),
    ("C3", 0, 3, None),
    ("C3", 3, 4, None),
    ("C3", 10, 12, -8),
    ("C3", -8, -2, None),
]

# Which atoms bond this residue to its neighbours, one row per residue, aligned
# with RESIDUE_NAMES: [[own_head, other_head], [own_tail, other_tail],
# head_len, tail_len]. "head" is the join towards the 5' end and "tail" the one
# towards the 3' end; the two lengths are in ångström. Every RNA residue joins
# the same way, so one row is repeated rather than written out sixteen times.
CONNECT: list[list] = [[[0, -1], [-2, 0], 1.6, 1.6] for _ in RESIDUE_NAMES]

# What each written token becomes, by where it sits in the strand:
# [token, alone, start, middle, end]. A residue name is a valid token too, so
# a sequence may be written either way.
ALIASES: list[list[str]] = [
    ["CN", "CN", "C5", "C", "C3"],
    ["A", "AN", "A5", "A", "A3"],
    ["C", "CN", "C5", "C", "C3"],
    ["UN", "UN", "U5", "U", "U3"],
    ["G", "GN", "G5", "G", "G3"],
    ["G3", "G3", "G", "G", "G3"],
    ["AN", "AN", "A5", "A", "A3"],
    ["A3", "A3", "A", "A", "A3"],
    ["GN", "GN", "G5", "G", "G3"],
    ["A5", "A5", "A5", "A", "A"],
    ["U", "UN", "U5", "U", "U3"],
    ["G5", "G5", "G5", "G", "G"],
    ["U3", "U3", "U", "U", "U3"],
    ["C5", "C5", "C5", "C", "C"],
    ["C3", "C3", "C", "C", "C3"],
    ["U5", "U5", "U5", "U", "U"],
]


@cache
def rna() -> ResidueLibrary:
    """Return the RNA residue library, built from the tables in this module.

    Everything MAWS needs to design an RNA aptamer, other than the physics: the
    sixteen residues, their atom counts, their turnable bonds and the four
    names each written nucleotide can take.

    Returns
    -------
    maws.values.ResidueLibrary
        RNA templates keyed by residue name, together with the alias table that
        maps written tokens (``G A U C``) onto them. The same object every
        time, so do not try to change it.

    Raises
    ------
    maws.errors.ConfigurationError
        If the tables in this module do not line up, or a turnable bond names
        an atom its residue does not have. That would be an error in this
        file, so it surfaces on the first call rather than during a design run.

    See Also
    --------
    maws.libraries.dna : The same thing for DNA.
    maws.values.ResidueLibrary.from_tables : Does the checking and building.

    Notes
    -----
    Cached, because the tables are constant and the library never changes, so
    one instance per process is safe to share and the checking is paid for
    once.

    Examples
    --------
    >>> lib = rna()
    >>> len(lib)
    16
    >>> lib["GN"].n_atoms
    33

    A lone G is a different residue from a G at the 5' end:

    >>> lib.alias("G").alone, lib.alias("G").start
    ('GN', 'G5')

    The same object comes back every time:

    >>> rna() is rna()
    True
    """
    return ResidueLibrary.from_tables(
        names=RESIDUE_NAMES,
        lengths=RESIDUE_LENGTH,
        rotations=ROTATIONS,
        connect=CONNECT,
        aliases=ALIASES,
    )

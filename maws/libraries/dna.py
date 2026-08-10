"""
maws.libraries.dna
==================

DNA residue tables, and the ``dna()`` factory that validates them.

This module holds data and nothing else: no file access, no OpenMM, no LEaP.
The tables are transcribed from the original ``DNA.xml`` and carry the same
numbers the legacy ``maws.dna_structure`` module did.

Reading the tables
------------------
``RESIDUE_NAMES`` and ``RESIDUE_LENGTH`` are positionally aligned: entry *i* of
one describes the same residue as entry *i* of the other. ``ROTATIONS`` and
``ALIASES`` are keyed by name instead. That mix is exactly why
:meth:`maws.values.ResidueLibrary.from_tables` checks the alignment once, on
load, rather than letting a misaligned row surface later as a LEaP complaint
about a missing parameter.

Rotation rows are ``(residue, start, bond, end_or_None)``. An ``end_or_None`` of
``None`` means "rotate the rest of the chain" and becomes a
:class:`~maws.values.BackboneTorsion`; a number bounds the moving atoms inside
the residue and becomes a :class:`~maws.values.LocalTorsion`. Negative indices
count back from the residue's last atom and are resolved on load.

Alias rows are ``[token, alone, start, middle, end]``: the four canonical LEaP
names one written token takes, depending on where it sits in the chain.
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

# Canonical residue names. Positionally aligned with RESIDUE_LENGTH.
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

# Atom count per residue, aligned with RESIDUE_NAMES.
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

# Rotatable bonds: (residue, start, bond, end_or_None).
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

# Polymer connectivity, one row per residue, aligned with RESIDUE_NAMES:
# [[own_head, other_head], [own_tail, other_tail], head_len, tail_len].
# Every residue connects the same way, so the rows are built by repetition.
CONNECT: list[list] = [[[0, -1], [-2, 0], 1.6, 1.6] for _ in RESIDUE_NAMES]

# Alias rows: [token, alone, start, middle, end].
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
    """Return the validated DNA residue library.

    Returns
    -------
    maws.values.ResidueLibrary
        DNA templates keyed by canonical residue name, together with the
        alias table that maps written tokens (``G A T C``) onto them.

    Notes
    -----
    Cached, because the tables are constant and the library is immutable, so
    one instance per process is safe to share.

    This replaces the legacy ``load_dna_structure()``, which also read
    ``.lib`` files from disk inside its constructor and so could not be called
    in a test without AmberTools output already on disk.

    Examples
    --------
    >>> lib = dna()
    >>> len(lib)
    16
    >>> lib["DGN"].n_atoms
    32
    """
    return ResidueLibrary.from_tables(
        names=RESIDUE_NAMES,
        lengths=RESIDUE_LENGTH,
        rotations=ROTATIONS,
        connect=CONNECT,
        aliases=ALIASES,
    )

"""
maws.libraries.rna
==================

RNA residue tables, and the ``rna()`` factory that validates them.

This module holds data and nothing else: no file access, no OpenMM, no LEaP.
The tables are transcribed from the original ``RNA.xml`` and carry the same
numbers the legacy ``maws.rna_structure`` module did.

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
    "rna",
]

# Canonical residue names. Positionally aligned with RESIDUE_LENGTH.
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

# Atom count per residue, aligned with RESIDUE_NAMES.
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

# Rotatable bonds: (residue, start, bond, end_or_None).
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

# Polymer connectivity, one row per residue, aligned with RESIDUE_NAMES:
# [[own_head, other_head], [own_tail, other_tail], head_len, tail_len].
# Every residue connects the same way, so the rows are built by repetition.
CONNECT: list[list] = [[[0, -1], [-2, 0], 1.6, 1.6] for _ in RESIDUE_NAMES]

# Alias rows: [token, alone, start, middle, end].
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
    """Return the validated RNA residue library.

    Returns
    -------
    maws.values.ResidueLibrary
        RNA templates keyed by canonical residue name, together with the
        alias table that maps written tokens (``G A U C``) onto them.

    Notes
    -----
    Cached, because the tables are constant and the library is immutable, so
    one instance per process is safe to share.

    This replaces the legacy ``load_rna_structure()``, which also read
    ``.lib`` files from disk inside its constructor and so could not be called
    in a test without AmberTools output already on disk.

    Examples
    --------
    >>> lib = rna()
    >>> len(lib)
    16
    >>> lib["GN"].n_atoms
    33
    """
    return ResidueLibrary.from_tables(
        names=RESIDUE_NAMES,
        lengths=RESIDUE_LENGTH,
        rotations=ROTATIONS,
        connect=CONNECT,
        aliases=ALIASES,
    )

"""
MAWS residue template and topology rules.

This module defines the **static chemistry** used by MAWS:
- Residue metadata (atom counts, alias mapping for LEaP names)
- Polymer connectivity rules for appending/prepending residues
- Rotating elements (torsion definitions) used at runtime
- Backbone anchor indices (normalized at construction)

Design
------
`Structure` is a lightweight container that holds *per-residue rules*.
Runtime code (e.g., :class:`maws.chain.Chain` and :class:`maws.complex.Complex`)
consults these rules to:
  - translate human-friendly alias sequences into canonical LEaP names,
  - calculate residue lengths and index offsets,
  - decide which atom triples define rotatable elements,
  - find which atoms should be joined when growing a polymer chain.

Indexing conventions
--------------------
- Atom indices are **0-based** inside each residue.
- Negative indices mean “from the end” (Python style), e.g. -1 → last atom.
- Backbone indices are normalized to **non-negative** at construction since
  downstream code uses them directly without further normalization.
- Torsion (rotation) indices are **not** normalized at construction: negative
  indices are kept and normalized later at runtime in :meth:`Structure.torsions`.

LEaP resources
--------------
If `residue_path` is given (or `""` for current directory), `Structure` creates
a deterministic `init_string` that loads `<name>.lib` and `<name>.frcmod` for
each residue, suitable to paste into a LEaP input file.

Examples
--------
Build a minimal single-residue template:

>>> s = Structure(
...     residue_names=["X"],
...     residue_length=[3],
...     rotating_elements=[("X", 0, 1, None)],  # torsion about bond (0,1)
...     backbone_elements=[("X", 0, 1, 2, 1, 2)],  # connection anchors
...     connect=[[[0, -1], [-2, 0], 1.6, 1.6]],  # polymer rules
...     alias=[["X", "X", "X", "X", "X"]],  # alone/start/middle/end
... )
>>> s.resolve_index("X", -1)
2
>>> s.torsions("X")
[(0, 1, None)]
>>> s.translate("X X X")
'X X X'  # first=start, middle=middle, last=end

See Also
--------
maws.chain.Chain
    Consumes `Structure` to manage per-chain sequences and geometry.
maws.complex.Complex
    Owns coordinates and uses `Chain` to apply rotations/translations.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import TypeAlias

# Type aliases to make shapes explicit
ResidueName: TypeAlias = str
AtomIndex: TypeAlias = (
    int  # 0-based; negative values mean "from the end" (Python-style)
)

# Rotation specification at construction:
# (residue, start_atom_idx, bond_atom_idx, end_atom_idx_or_None)
RotationSpec: TypeAlias = tuple[ResidueName, AtomIndex, AtomIndex, AtomIndex | None]

# Backbone specification at construction (raw indices, may be negative):
# (residue, start, middle_pre, bond, middle_post, end)
BackboneSpec: TypeAlias = tuple[
    ResidueName, AtomIndex, AtomIndex, AtomIndex, AtomIndex, AtomIndex
]

# Stored backbone entry after normalization to non-negative indices:
# [[start, middle_pre, bond], [middle_post, end]]
BackboneEntry: TypeAlias = list[list[int]]

# A single rotation triple stored internally as a list: [start, bond, end_or_None]
RotTripleStored: TypeAlias = list[int | None]

# The per-residue rotation list:
# either [[start, bond, end], ...] OR the sentinel [None] meaning “not defined yet”.
RotListStored: TypeAlias = list[RotTripleStored]  # OR the sentinel: [None]

# Alias entry stored per residue: [alone, start, middle, end]
AliasEntry: TypeAlias = list[str]

# Connectivity entry stored per residue:
# [[append_first, append_last],
# [prepend_last, prepend_first],
# append_bond_len, prepend_bond_len]
# Atom indices may be negative (Python-style). Bond lengths are in Å.
ConnectEntry: TypeAlias = list[list[int] | float]


class Structure:
    """
    Container for residue templates and per-residue topology rules.

    Parameters
    ----------
    residue_names : Sequence[str]
        Residue template names in **fixed order** (e.g. `["A", "C", "G", "U"]` or
        `["DG", ...]`).
        This order is used to align other per-residue arrays (lengths,
        connectivity, etc.).
    residue_length : Sequence[int], optional
        Number of atoms per residue, aligned 1:1 with `residue_names`.
        Used for:
          * normalizing negative indices in backbone specs,
          * resolving indices in helpers like :meth:`resolve_index`,
          * computing chain atom counts.
    rotating_elements : Sequence[RotationSpec], optional
        Torsion definitions: triples `(residue, start, bond, end_or_None)`.
        Negative indices are allowed and will be normalized *later* by :meth:`torsions`.
    backbone_elements : Sequence[BackboneSpec], optional
        `(residue, start, middle_pre, bond, middle_post, end)`; indices may be negative.
        These are normalized to **non-negative** indices at construction and stored as
        `[[start, middle_pre, bond], [middle_post, end]]`.
    connect : Sequence[ConnectEntry], optional
        Per-residue polymer connectivity rules:
        `[[append_first, old_last], [new_last, old_first], append_len, prepend_len]`.
        Indices may be negative. Bond lengths are floats in Å.
        If omitted, a generic default `[[0, -1], [-2, 0], 1.6, 1.6]` is used.
    residue_path : str | None, optional
        Directory containing `<name>.lib` and `<name>.frcmod` for each residue.
        If `""`, use `"."`. If `None`, no LEaP `init_string` is produced.
    alias : Sequence[Sequence[str]] | None, optional
        Alias mapping rows `[residue_name, alone, start, middle, end]`.
        Internally stored per residue as `[alone, start, middle, end]`.
        If omitted, defaults to identity mapping for all residues.

    Attributes
    ----------
    init_string : str
        Deterministic LEaP bootstrap (one `loadoff` + `loadamberparams` pair per
        residue),or empty if `residue_path is None`.
    residue_length : dict[str, int]
        Map residue → atom count.
    connect : dict[str, ConnectEntry]
        Map residue → connectivity entry used by chain growth logic.
    alias : dict[str, list[str]]
        Map residue → `[alone, start, middle, end]` LEaP names.
    rotating_elements : dict[str, RotListStored]
        Map residue → rotation triples or sentinel `[None]`.
    backbone_elements : dict[str, BackboneEntry]
        Map residue → normalized backbone anchor indices.

    Notes
    -----
    - Indices are 0-based. Negative values mean “from the end”.
    - Backbone indices are normalized at construction; torsion indices are not.
    - `init_string` is generated only when `residue_path` is not `None`.
    """

    def __init__(
        self,
        residue_names: Sequence[ResidueName],
        residue_length: Sequence[int] | None = None,
        rotating_elements: Sequence[RotationSpec] | None = None,
        backbone_elements: Sequence[BackboneSpec] | None = None,
        connect: Sequence[ConnectEntry] | None = None,
        residue_path: str | None = None,
        alias: Sequence[Sequence[str]] | None = None,
    ):
        # Ordered list of residue names used to align all per-residue arrays
        self.residue_names: list[ResidueName] = list(residue_names)

        # Path for LEaP resources; if None, do not emit LEaP commands
        self.residue_path: str | None = residue_path

        # LEaP bootstrap string (pairs of loadoff/loadamberparams)
        self.init_string: str = ""
        if self.residue_path is not None:
            base = self.residue_path if self.residue_path != "" else "."
            for name in self.residue_names:
                self.init_string += (
                    f"loadoff {base}/{name}.lib\nloadamberparams {base}/{name}.frcmod\n"
                )

        # Map residue -> atom count (used to normalize negatives and report lengths)
        self.residue_length: dict[ResidueName, int] = {}
        if residue_length:
            for idx, res in enumerate(self.residue_names):
                self.residue_length[res] = int(residue_length[idx])

        # Map residue -> connectivity entry
        self.connect: dict[ResidueName, ConnectEntry] = {}
        if connect:
            for idx, res in enumerate(self.residue_names):
                self.connect[res] = list(connect[idx])
        else:
            default_conn: ConnectEntry = [[0, -1], [-2, 0], 1.6, 1.6]
            for res in self.residue_names:
                self.connect[res] = list(default_conn)

        # Map residue -> alias entry [alone, start, middle, end]
        # Initialize to identity mapping for known residues, then overlay user entries
        self.alias: dict[ResidueName, AliasEntry] = {
            res: [res, res, res, res] for res in self.residue_names
        }
        if alias:
            for row in alias:
                # incoming rows are [name, alone, start, middle, end]
                if len(row) != 5:
                    continue  # ignore malformed rows
                name = row[0]
                self.alias[name] = [row[1], row[2], row[3], row[4]]

        # Map residue -> rotation list.
        # Start each residue with a sentinel [None] meaning “no rotations defined yet”.
        self.rotating_elements: dict[ResidueName, RotListStored] = {}
        for name in self.residue_names:
            self.rotating_elements[name] = [None]  # sentinel

        if rotating_elements:
            for residue, start, bond, end in rotating_elements:
                if self.rotating_elements[residue] == [None]:
                    self.rotating_elements[residue] = [[start, bond, end]]
                elif self.rotating_elements.get(residue) is None:
                    raise ValueError(
                        "Residue does not exist! CANNOT assign rotability!"
                    )
                else:
                    self.rotating_elements[residue].append([start, bond, end])

        # Map residue -> normalized backbone entry:
        # [[start, middle_pre, bond], [middle_post, end]] with all indices >= 0
        self.backbone_elements: dict[ResidueName, BackboneEntry] = {}
        if backbone_elements:
            for residue, start, middle_pre, bond, middle_post, end in backbone_elements:
                if (
                    residue not in self.residue_length
                    or self.residue_length[residue] <= 0
                ):
                    raise ValueError(
                        f"Backbone specified for {residue!r} but its length is not set."
                    )
                L = self.residue_length[residue]

                def norm(i: int, *, L: int = L) -> int:
                    # Convert negative indices to absolute positions
                    return i + L if i < 0 else i

                self.backbone_elements[residue] = [
                    [norm(start), norm(middle_pre), norm(bond)],
                    [norm(middle_post), norm(end)],
                ]

    # -------------------------------------------------------------------------
    # Public helpers used by Chain/Complex
    # -------------------------------------------------------------------------

    def add_rotation(
        self,
        residue_name: ResidueName,
        rotations: tuple[int, int, int | None] | Iterable[tuple[int, int, int | None]],
    ) -> dict[ResidueName, RotListStored]:
        """
        Add new rotating-element definitions for a residue.

        Parameters
        ----------
        residue_name : str
            Residue to augment.
        rotations : (int, int, int | None) or Iterable[tuple[int, int, int | None]]
            Either a single triple `[start, bond, end_or_None]`, or an iterable
            of such triples. Negative indices are allowed and normalized later.
        basestring : type, optional
            Ignored; kept for API compatibility with older code.

        Returns
        -------
        dict[str, RotListStored]
            Updated :attr:`rotating_elements`.

        Raises
        ------
        ValueError
            If any rotation is not a 3-tuple/list.
        """
        # Ensure the entry exists and is a list (not the sentinel)
        if self.rotating_elements.get(residue_name) == [None]:
            self.rotating_elements[residue_name] = []

        def is_triplet(x: object) -> bool:
            return isinstance(x, list | tuple) and len(x) == 3

        if is_triplet(rotations):
            s, b, e = rotations  # type: ignore[index]
            self.rotating_elements[residue_name].append(
                [int(s), int(b), None if e is None else int(e)]
            )
        elif isinstance(rotations, Iterable):
            for rot in rotations:  # type: ignore[assignment]
                if not is_triplet(rot):
                    raise ValueError(
                        "Each rotation must be a [start, bond, end] triple."
                    )
                s, b, e = rot  # type: ignore[misc]
                self.rotating_elements[residue_name].append(
                    [int(s), int(b), None if e is None else int(e)]
                )
        else:
            raise ValueError(
                "rotations must be a triple or an iterable of triples "
                "[start, bond, end]."
            )

        return self.rotating_elements

    def translate(self, sequence: str) -> str:
        """
        Translate an **alias** sequence (space-separated) into a **canonical**
        LEaP residue-name sequence using per-residue alias mapping.

        Rules
        -----
        - Single residue → ``alone``
        - First residue → ``start``
        - Middle residues → ``middle``
        - Last residue → ``end``

        Parameters
        ----------
        sequence : str
            Space-separated alias tokens, e.g. ``"G A U"``.

        Returns
        -------
        str
            Space-separated canonical residue names ready for LEaP.

        Examples
        --------
        If ``alias['A'] = [A, A5, A, A3]`` then
        - ``'A'`` → ``'A'``
        - ``'A B C'`` → ``'A5 B C3'``
        """
        sequence_array = sequence.split()  # robust to multiple whitespaces
        if len(sequence_array) == 1:
            return self.alias[sequence_array[0]][0]  # alone
        first = self.alias[sequence_array[0]][1]
        middles = [self.alias[name][2] for name in sequence_array][1:-1]
        last = self.alias[sequence_array[-1]][3]
        return " ".join([first] + middles + [last])

    def resolve_index(self, residue: ResidueName, i: AtomIndex) -> int:
        """
        Normalize an atom index for a given residue.

        Parameters
        ----------
        residue : str
            Residue name to resolve against.
        i : int
            0-based index; negative values mean “from the end”.

        Returns
        -------
        int
            Absolute non-negative index within the residue.

        Raises
        ------
        ValueError
            If residue is unknown or its length is not set.
        """
        if residue not in self.residue_length or self.residue_length[residue] <= 0:
            raise ValueError(f"Unknown residue or length not set: {residue!r}")
        L = self.residue_length[residue]
        return i + L if i < 0 else i

    def append_bond(
        self, residue: ResidueName, prev_residue_length: int | None = None
    ) -> tuple[int, int, float]:
        """
        Return connection info when **appending** this residue to the right (3') end.

        Parameters
        ----------
        residue : str
            Residue being appended (new rightmost residue).
        prev_residue_length : int | None, optional
            Atom count of the *previous* (left) residue. If provided, the
            returned ``old_atom_idx`` is resolved to an absolute index; if not,
            it may be negative (relative to the previous residue).

        Returns
        -------
        tuple[int, int, float]
            ``(new_atom_idx, old_atom_idx, bond_length_Å)``

        Raises
        ------
        ValueError
            If connectivity for `residue` is missing.
        """
        try:
            append_pair, _, append_len, _ = self.connect[
                residue
            ]  # [[new_first, old_last], ...]
        except KeyError as e:
            raise ValueError(f"No connectivity entry for residue {residue!r}") from e

        new_first, old_last = append_pair
        new_idx = self.resolve_index(residue, int(new_first))

        if prev_residue_length is not None:
            old_idx = (
                old_last + int(prev_residue_length) if old_last < 0 else int(old_last)
            )
        else:
            old_idx = int(old_last)

        return new_idx, old_idx, float(append_len)

    def prepend_bond(
        self, residue: ResidueName, next_residue_length: int | None = None
    ) -> tuple[int, int, float]:
        """
        Return connection info when **prepending** this residue to the left (5') end.

        Parameters
        ----------
        residue : str
            Residue being prepended (new leftmost residue).
        next_residue_length : int | None, optional
            Atom count of the *next* (right) residue. If provided, the returned
            ``old_atom_idx`` is resolved to an absolute index; if not, it may be
            negative (relative to the next residue).

        Returns
        -------
        tuple[int, int, float]
            ``(new_atom_idx, old_atom_idx, bond_length_Å)``

        Raises
        ------
        ValueError
            If connectivity for `residue` is missing.
        """
        try:
            _, prepend_pair, _, prepend_len = self.connect[
                residue
            ]  # [[...], [new_last, old_first], ...]
        except KeyError as e:
            raise ValueError(f"No connectivity entry for residue {residue!r}") from e

        new_last, old_first = prepend_pair
        new_idx = self.resolve_index(residue, int(new_last))

        if next_residue_length is not None:
            old_idx = (
                old_first + int(next_residue_length)
                if old_first < 0
                else int(old_first)
            )
        else:
            old_idx = int(old_first)

        return new_idx, old_idx, float(prepend_len)

    def torsions(self, residue: ResidueName) -> list[tuple[int, int, int | None]]:
        """
        Return all rotation triples for a residue with indices normalized.

        Parameters
        ----------
        residue : str
            Residue name to query.

        Returns
        -------
        list[tuple[int, int, int | None]]
            Each tuple is ``(start, bond, end_or_None)`` with indices normalized
            to absolute 0-based; if `end` was `None`, it remains `None`.

        Raises
        ------
        ValueError
            If the residue is unknown or its length is not set.
        """
        if residue not in self.rotating_elements:
            raise ValueError(f"Unknown residue {residue!r}")
        triples = self.rotating_elements[residue]
        if triples == [None]:  # sentinel for "no rotations defined"
            return []

        if residue not in self.residue_length or self.residue_length[residue] <= 0:
            raise ValueError(f"Length not set for residue {residue!r}")
        L = self.residue_length[residue]

        def norm(x: int | None) -> int | None:
            if x is None:
                return None
            return x + L if x < 0 else x

        out: list[tuple[int, int, int | None]] = []
        for t in triples:
            s = int(t[0])
            b = int(t[1])
            e = t[2] if (t[2] is None) else int(t[2])
            out.append((norm(s), norm(b), norm(e)))  # type: ignore[arg-type]
        return out

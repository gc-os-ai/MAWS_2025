"""
structure.py

Defines residue-based chemistry used by MAWS:
- residue metadata (length, alias mapping)
- polymer connectivity rules
- rotating elements (handled dynamically at runtime)
- backbone anchor indices (normalized here)

Notes
-----
- LEaP `init_string` is generated deterministically from `residue_path` and `residue_names`.
- Rotating elements may use negative indices (relative to residue length). We keep
  them as-is because normalization is performed later during rotation.
- Backbone elements are normalized here to positive indices because later code
  uses them directly without further normalization.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import TypeAlias

# ── Type aliases to make shapes explicit ─────────────────────────────────────
ResidueName: TypeAlias = str
AtomIndex: TypeAlias = (
    int  # 0-based; negative values mean "from the end" (Python-style)
)

# Rotation specification provided at construction time
# (residue, start_atom_idx, bond_atom_idx, end_atom_idx_or_None)
RotationSpec: TypeAlias = tuple[ResidueName, AtomIndex, AtomIndex, AtomIndex | None]

# Backbone specification provided at construction time
# (residue, start, middle_pre, bond, middle_post, end) — raw indices (may be negative)
BackboneSpec: TypeAlias = tuple[
    ResidueName, AtomIndex, AtomIndex, AtomIndex, AtomIndex, AtomIndex
]

# Stored backbone entry after normalization to positive indices:
# [[start, middle_pre, bond], [middle_post, end]]
BackboneEntry: TypeAlias = list[list[int]]

# A single rotation triple stored internally as a list: [start, bond, end_or_None]
RotTripleStored: TypeAlias = list[int | None]

# The per-residue rotation list.
# IMPORTANT: we use a sentinel value [None] to mean "not defined yet".
# Otherwise it is a list of [start, bond, end] triples.
RotListStored: TypeAlias = list[RotTripleStored]  # OR the sentinel: [None]

# Alias entry stored per residue: [alone, start, middle, end]
AliasEntry: TypeAlias = list[str]

# Connectivity entry stored per residue:
# [[append_first, append_last], [prepend_last, prepend_first], append_bond_len, prepend_bond_len]
# Atom indices may be negative (Python-style). Bond lengths are in Å.
ConnectEntry: TypeAlias = list[list[int] | float]


class Structure:
    """
    Container for residue templates and per-residue topology rules.

    Parameters
    ----------
    residue_names : Sequence[str]
        Names of residue templates in a fixed order (e.g., ["A", "C", "G", "U"] or ["DG", ...]).
        This order is used to align other per-residue arrays (lengths, connectivity, etc.).
    residue_length : Sequence[int], optional
        Number of atoms per residue, aligned 1:1 with `residue_names`.
        (Used to normalize negative indices in backbone specs.)
    rotating_elements : Sequence[RotationSpec], optional
        Triples ``(residue, start, bond, end_or_None)``. Negative indices are allowed
        and are interpreted later at runtime (not normalized here).
    backbone_elements : Sequence[BackboneSpec], optional
        Tuples ``(residue, start, middle_pre, bond, middle_post, end)`` possibly
        containing negative indices, which are normalized here to non-negative indices
        using the residue's length (from `residue_length`).
    connect : Sequence[ConnectEntry], optional
        Per-residue connectivity:
        ``[[append_first, append_last], [prepend_last, prepend_first], append_len, prepend_len]``.
        Atom indices can be negative; bond lengths are floats in Å.
    residue_path : str | None
        Directory containing ``<name>.lib`` and ``<name>.frcmod`` for each residue.
        If "", use the current directory ``"."``. If ``None``, no LEaP `init_string` is produced.
    alias : Sequence[Sequence[str]] | None
        Alias mapping entries of the form ``[residue_name, alone, start, middle, end]``.
        Internally stored per residue as ``[alone, start, middle, end]``.

    Notes
    -----
    - Indices are 0-based. Negative values mean "from the end" (Python-style).
    - Backbone indices are normalized at construction time; rotation indices are not.
    - `init_string` is generated deterministically from `residue_path` and `residue_names`.

    Examples
    --------
    Build from the bundled RNA templates:

    >>> from RNA import rna_structure
    >>> s = rna_structure()
    >>> s.residue_length["A"]
    33
    >>> s.torsions("A")[:2]
    [(0, 3, None), (3, 4, None)]
    >>> s.append_bond("A", prev_residue_length=33)
    (0, 31, 1.6)
    >>> s.translate("G A U")
    'G5 A U3'

    Minimal manual construction:

    >>> s = Structure(
    ...     residue_names=["X"],
    ...     residue_length=[3],
    ...     rotating_elements=[("X", 0, 1, None)],
    ...     backbone_elements=[("X", 0, 1, 2, 1, 2)],
    ...     connect=[[[0, -1], [-2, 0], 1.6, 1.6]],
    ...     alias=[["X", "X", "X", "X", "X"]],
    ... )
    >>> s.resolve_index("X", -1)
    2
    >>> s.torsions("X")
    [(0, 1, None)]
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
        # Ordered list of residue template names used throughout this object.
        self.residue_names: list[ResidueName] = list(residue_names)

        # Path for LEaP resources; if None, we do not emit LEaP commands.
        # If "", treat as current directory ".".
        self.residue_path: str | None = residue_path

        # LEaP bootstrap string (lines like `loadoff X.lib` and `loadamberparams X.frcmod`)
        self.init_string: str = ""
        if self.residue_path is not None:
            base = self.residue_path if self.residue_path != "" else "."
            for name in self.residue_names:
                self.init_string += (
                    f"loadoff {base}/{name}.lib\nloadamberparams {base}/{name}.frcmod\n"
                )

        # Map residue -> atom count (used to normalize negative indices for backbone)
        self.residue_length: dict[ResidueName, int] = {}
        if residue_length:
            for idx, res in enumerate(self.residue_names):
                self.residue_length[res] = int(residue_length[idx])

        # Map residue -> connectivity entry.
        # Default if not provided: [[0, -1], [-2, 0], 1.6, 1.6]
        self.connect: dict[ResidueName, ConnectEntry] = {}
        if connect:
            for idx, res in enumerate(self.residue_names):
                self.connect[res] = list(connect[idx])
        else:
            default_conn: ConnectEntry = [[0, -1], [-2, 0], 1.6, 1.6]
            for res in self.residue_names:
                self.connect[res] = list(default_conn)

        # Map residue -> alias entry [alone, start, middle, end]
        # Initialize to identity mapping for known residues.
        self.alias: dict[ResidueName, AliasEntry] = {
            res: [res, res, res, res] for res in self.residue_names
        }
        if alias:
            # incoming entries are [name, alone, start, middle, end]
            for elem in alias:
                name = elem[0]
                if len(elem) != 5:
                    continue  # silently ignore malformed rows
                self.alias[name] = [elem[1], elem[2], elem[3], elem[4]]

        # Map residue -> rotation list.
        # Start each residue with a sentinel [None] meaning "no rotations defined yet".
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
        # [[start, middle_pre, bond], [middle_post, end]] with all indices >= 0.
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

                def norm(i: int) -> int:
                    # Convert negative indices to positive positions relative to residue length.
                    # Example: -1 -> L-1, -2 -> L-2
                    return i + L if i < 0 else i

                revised_start = norm(start)
                revised_middle_pre = norm(middle_pre)
                revised_bond = norm(bond)
                revised_middle_post = norm(middle_post)
                revised_end = norm(end)
                self.backbone_elements[residue] = [
                    [revised_start, revised_middle_pre, revised_bond],
                    [revised_middle_post, revised_end],
                ]

    # -------------------------------------------------------------------------

    def add_rotation(
        self,
        residue_name: ResidueName,
        rotations: tuple[int, int, int | None] | Iterable[tuple[int, int, int | None]],
        basestring: type = str,  # kept for backward-compat with old API that passed `str`
    ) -> dict[ResidueName, RotListStored]:
        """
        Add new rotating-element definitions for a residue.

        Parameters
        ----------
        residue_name : str
            Residue to augment.
        rotations :
            Either a single triple [start, bond, end_or_None] or an iterable of such triples.
            Negative indices are allowed and interpreted relative to residue length at runtime.
        basestring :
            Ignored; kept for backward compatibility with older APIs that passed `str`.

        Returns
        -------
        dict
            Updated `self.rotating_elements`.
        """
        # Ensure the entry exists and is a list
        if self.rotating_elements.get(residue_name) == [None]:
            self.rotating_elements[residue_name] = []

        # Helper: validates a 3-tuple/list
        def is_triplet(x: object) -> bool:
            return isinstance(x, (list, tuple)) and len(x) == 3

        if is_triplet(rotations):
            start, bond, end = rotations  # type: ignore[index]
            self.rotating_elements[residue_name].append(
                [int(start), int(bond), None if end is None else int(end)]
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
                "rotations must be a triple or an iterable of triples [start, bond, end]."
            )

        return self.rotating_elements

    def translate(self, sequence: str) -> str:
        """
        Translate an alias sequence (space-separated) into a residue-name sequence
        using the per-residue alias mapping.

        Examples
        --------
        If alias['A'] = [A, A5, A, A3], then:
            'A'           -> 'A'
            'A B C'       -> 'A5 B C3' (first uses start, last uses end, middle use middle)

        Returns
        -------
        str
            Space-separated residue-name sequence for LEaP.
        """
        sequence_array = sequence.split()  # robust to multiple whitespaces
        if len(sequence_array) == 1:
            # Single residue: use "alone" form.
            return self.alias[sequence_array[0]][0]
        else:
            # First -> "start", last -> "end", middle -> "middle".
            first = self.alias[sequence_array[0]][1]
            middles = [self.alias[name][2] for name in sequence_array][1:-1]
            last = self.alias[sequence_array[-1]][3]
            return " ".join([first] + middles + [last])

    # ---------- Convenience helpers (keep existing API untouched) ----------

    def resolve_index(self, residue: ResidueName, i: AtomIndex) -> int:
        """
        Normalize an atom index for `residue`: negative values become absolute 0-based.
        Raises if residue is unknown or its length is not set.
        """
        if residue not in self.residue_length or self.residue_length[residue] <= 0:
            raise ValueError(f"Unknown residue or length not set: {residue!r}")
        L = self.residue_length[residue]
        return i + L if i < 0 else i

    def append_bond(
        self,
        residue: ResidueName,
        prev_residue_length: int | None = None,
    ) -> tuple[int, int, float]:
        """
        Return the (new_atom_idx, old_atom_idx, bond_length) used when APPENDING this residue
        to the right (3') end of an existing chain.

        Notes
        -----
        - new_atom_idx is resolved (absolute index within the NEW residue).
        - old_atom_idx is:
            * resolved if `prev_residue_length` is provided,
            * otherwise returned as stored (may be negative, i.e., relative to the PREVIOUS residue).
        """
        try:
            append_pair, _, append_len, _ = self.connect[
                residue
            ]  # [[new_first, old_last], [..], append_len, ..]
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
        self,
        residue: ResidueName,
        next_residue_length: int | None = None,
    ) -> tuple[int, int, float]:
        """
        Symmetric helper for PREPENDING this residue to the left (5') end of an existing chain.

        Returns (new_atom_idx, old_atom_idx, bond_length) where:
        - new_atom_idx is resolved within THIS residue.
        - old_atom_idx is for the FIRST residue of the existing chain:
            resolved if `next_residue_length` is provided, else returned as stored (may be negative).
        """
        try:
            _, prepend_pair, _, prepend_len = self.connect[
                residue
            ]  # [[..], [old_last, old_first], .., prepend_len]
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
        Return all rotation triples for `residue` with indices normalized to absolute 0-based.
        Each item is (start, bond, end_or_None). If end is None, it stays None.
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
            # t is [start, bond, end_or_None]
            s = int(t[0])
            b = int(t[1])
            e = t[2] if (t[2] is None) else int(t[2])
            out.append((norm(s), norm(b), norm(e)))  # type: ignore[arg-type]
        return out

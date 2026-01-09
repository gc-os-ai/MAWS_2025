"""
maws.structure
==============

Residue templates and topology rules for MAWS.

Defines per-residue metadata: atom counts, alias mapping, connectivity rules,
and torsion definitions used by Chain and Complex.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import TypeAlias

# Type aliases
ResidueName: TypeAlias = str
AtomIndex: TypeAlias = int  # 0-based; negative = from end

RotationSpec: TypeAlias = tuple[ResidueName, AtomIndex, AtomIndex, AtomIndex | None]
BackboneSpec: TypeAlias = tuple[
    ResidueName, AtomIndex, AtomIndex, AtomIndex, AtomIndex, AtomIndex
]
BackboneEntry: TypeAlias = list[list[int]]
RotTripleStored: TypeAlias = list[int | None]
RotListStored: TypeAlias = list[RotTripleStored]
AliasEntry: TypeAlias = list[str]
ConnectEntry: TypeAlias = list[list[int] | float]


class Structure:
    """
    Container for residue templates and per-residue topology rules.

    Parameters
    ----------
    residue_names : Sequence[str]
        Residue template names in fixed order.
    residue_length : Sequence[int], optional
        Atom count per residue.
    rotating_elements : Sequence[RotationSpec], optional
        Torsion definitions: (residue, start, bond, end_or_None).
    backbone_elements : Sequence[BackboneSpec], optional
        Connection anchors for polymer growth.
    connect : Sequence[ConnectEntry], optional
        Polymer connectivity rules.
    residue_path : str | None, optional
        Directory containing .lib/.frcmod files for LEaP.
    alias : Sequence[Sequence[str]] | None, optional
        Alias mapping rows [residue, alone, start, middle, end].
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

"""
maws.chain
==========

A Chain is a lightweight view into a Complex's coordinate array.

Each Chain tracks its slice of atoms and provides methods for sequence
manipulation and per-residue rotations. Geometry is delegated to Complex.
"""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # Only for static type checkers; avoids runtime circular import
    from maws.complex import Complex


class Chain:
    """
    A single polymer chain inside a Complex.

    Manages sequence state and provides methods for per-residue rotations.
    Geometry operations are delegated to the owning Complex.

    Attributes
    ----------
    id : int
        Chain index within Complex.chains.
    start : int
        Global atom index where this chain begins.
    length : int
        Total atoms in this chain.
    structure : Structure
        Residue templates and topology rules.
    sequence : str
        Space-separated canonical residue names for LEaP.
    alias_sequence : str
        Space-separated alias tokens (human-readable).
    """

    def __init__(
        self,
        Complex: Complex,
        Structure,
        sequence: str | None = None,
        start: int = 0,
        ID: int = 0,
    ):
        """
        Parameters
        ----------
        Complex : Complex
            Owning complex instance. Geometry is delegated to this object.
        Structure : Structure
            Residue/template bank used to translate aliases, determine residue
            lengths, rotating elements, and connectivity.
        sequence : str or None, optional
            Initial **alias** sequence (space-separated). If provided, it is
            immediately translated to canonical residue names via
            ``Structure.translate`` and all per-residue offsets and lengths are
            computed. If ``None`` or empty, the chain starts empty and can be
            configured later with :meth:`create_sequence`, :meth:`append_sequence`,
            or :meth:`prepend_sequence`.
        start : int, default=0
            Global atom index where this chain begins in the complex.
        ID : int, default=0
            Chain identifier.

        Notes
        -----
        This constructor initializes sequence bookkeeping only. To materialize a
        topology/coordinate set you must call :meth:`Complex.build` (or
        :meth:`Complex.rebuild` after edits).
        """
        self.id = ID
        self.start = start
        self.start_history = start
        self.complex = Complex

        self.residues_start: list[int] = []
        self.length = 0
        self.length_history = self.length
        self.element = [self.start, self.start + 1, self.start + self.length]

        self.structure = Structure
        self.alias_sequence = ""
        self.sequence = ""
        self.sequence_array: list[str] = []
        self.alias_sequence_array: list[str] = []
        self.append_history: list[str] = []
        self.prepend_history: list[str] = []

        # Optional eager initialization from an alias sequence
        if sequence:
            self.alias_sequence = sequence
            self.sequence = self.structure.translate(self.alias_sequence)
            self.sequence_array = self.sequence.split(" ")
            self.alias_sequence_array = self.alias_sequence.split(" ")
            self.length = sum(
                map(self.structure.residue_length.__getitem__, self.sequence_array)
            )
            self.length_history = self.length

            tally = 0
            for residue in self.sequence_array:
                self.residues_start.append(tally)
                tally += self.structure.residue_length[residue]

            self.element = [self.start, self.start + 1, self.start + self.length]

    def update_chains(self):
        """
        Recompute this chain's length and per-residue offsets, then update
        **downstream** chains' global starts to keep indices consistent.

        Notes
        -----
        When this chain's sequence changes, its atom count (``length``) can grow
        or shrink. Since all atoms of all chains live in one flattened list
        (``Complex.positions``), any chain that starts **at or after** this
        chain must have its ``start`` shifted by the delta
        ``(new_length - old_length)``. This method performs:

        1. Recompute ``length`` from ``sequence_array`` and rebuild
           ``residues_start``.
        2. For every chain in ``Complex.chains`` with ``chain.start >= self.start``,
           shift ``chain.start`` by the delta and refresh its ``element`` triple.
        3. Counter-adjust this chain's own ``start`` so that only *later* chains
           move (net effect: this chain's ``start`` remains stable, but later
           chains shift).
        """
        old_length = self.length
        self.length = sum(
            map(self.structure.residue_length.__getitem__, self.sequence_array)
        )

        self.residues_start = []
        tally = 0
        for residue in self.sequence_array:
            self.residues_start.append(tally)
            tally += self.structure.residue_length[residue]

        self.element = [self.start, self.start + 1, self.start + self.length]
        start = copy.deepcopy(self.start)  # int; deepcopy is harmless and explicit

        for chain in self.complex.chains:
            chain.start_history = chain.start
            if chain.start >= start:
                chain.start += self.length - old_length
                chain.element = [
                    chain.start,
                    chain.start + 1,
                    chain.start + chain.length,
                ]

        # Counter-adjust "self" so only later chains are effectively shifted
        self.start -= self.length - old_length
        self.element = [self.start, self.start + 1, self.start + self.length]

    def create_sequence(self, sequence: str):
        """
        Overwrite the entire chain sequence (aliases), validate, and update
        lengths/offsets. **Must call** :meth:`Complex.build` (or
        :meth:`Complex.rebuild`) afterwards to realize the change.

        Parameters
        ----------
        sequence : str
            New alias sequence (space-separated), e.g. ``"G A U"``.

        Raises
        ------
        ValueError
            If any translated residue name is unknown to the provided
            :class:`Structure`.

        Notes
        -----
        - Translates aliases via :meth:`Structure.translate`.
        - Recomputes ``sequence_array``, ``residues_start``, and ``length``.
        - Snapshots ``start_history`` and ``length_history`` for later rebuild
          coordinate mapping.
        """
        alias_sequence_array = sequence.split(" ")
        sequence_array = self.structure.translate(sequence).split(" ")
        for name in sequence_array:
            if name not in self.structure.residue_names:
                raise ValueError("Residue not defined! CANNOT create sequence!")

        self.alias_sequence = sequence
        self.sequence = self.structure.translate(self.alias_sequence)
        self.alias_sequence_array = alias_sequence_array
        self.sequence_array = sequence_array

        self.update_chains()
        self.start_history = self.start
        self.length_history = self.length
        self.sequence_array_history = self.sequence_array  # kept for compatibility

    def append_sequence(self, sequence: str):
        """
        Append one alias residue to the **right** (3′) end and record history
        for :meth:`Complex.rebuild`.

        Parameters
        ----------
        sequence : str
            Alias token to append (e.g., ``"A"`` or a ligand alias).

        Notes
        -----
        Internally calls :meth:`create_sequence` with the new alias added,
        then:
        - sets ``append_history`` to the canonical residues newly appended;
        - clears ``prepend_history``;
        - snapshots old ``length`` and ``start`` for rebuild mapping.
        """
        start_history = self.start
        old_length = self.length
        prev_len = len(self.sequence_array)

        self.create_sequence(" ".join(self.alias_sequence_array[:] + [sequence]))

        self.length_history = old_length
        self.start_history = start_history
        self.prepend_history = []
        self.append_history = self.sequence_array[prev_len:]

    def prepend_sequence(self, sequence: str):
        """
        Prepend one alias residue to the **left** (5′) end and record history
        for :meth:`Complex.rebuild`.

        Parameters
        ----------
        sequence : str
            Alias token to prepend.

        Notes
        -----
        Internally calls :meth:`create_sequence` with the alias inserted at the
        front, then:
        - sets ``prepend_history`` to the canonical residues newly prepended;
        - snapshots old ``length``/``start`` for rebuild mapping.
        """
        old_length = self.length
        prev_len = len(self.sequence_array)

        self.create_sequence(" ".join([sequence] + self.alias_sequence_array[:]))

        self.length_history = old_length
        self.start_history = self.start + self.length - old_length
        self.prepend_history = self.sequence_array[
            : len(self.sequence_array) - prev_len
        ]

    def rotate_element(self, element, angle: float, reverse: bool = False):
        """Turn one bond of this chain by `angle` radians.

        Parameters
        ----------
        element : sequence of int
            ``[first, second]`` atom indices naming the bond to turn, counted
            from the start of this chain. Longer sequences are accepted and
            the extra entries ignored.
        angle : float
            How far to turn, in radians.
        reverse : bool, default=False
            Turn the part of the molecule joined to `first` rather than the
            part joined to `second`.

        Raises
        ------
        ValueError
            If fewer than two indices are given.

        See Also
        --------
        rotate_in_residue : Names a bond by residue and torsion instead.
        maws.complex.Complex.rotate_element : Does the work.

        Notes
        -----
        Indices are counted from the start of this chain; they are shifted to
        whole-Complex indices here.
        """
        if len(element) < 2:
            raise ValueError(
                f"A bond needs two atom indices; got {len(element)}: {list(element)}"
            )
        self.complex.rotate_element(
            [element[0] + self.start, element[1] + self.start],
            angle,
            reverse=reverse,
        )

    def rotate_in_residue(
        self,
        residue_index: int,
        residue_element_index: int,
        angle: float,
        reverse: bool = False,
    ):
        """Turn one of a residue's named bonds by `angle` radians.

        Parameters
        ----------
        residue_index : int
            Which residue of this chain, counting from zero. Negative values
            count back from the last residue.
        residue_element_index : int
            Which of that residue's turnable bonds, counting from zero, in the
            order its :class:`~maws.structure.Structure` lists them.
        angle : float
            How far to turn, in radians.
        reverse : bool, default=False
            Turn the part of the molecule joined to the bond's first atom
            rather than the part joined to its second.

        Raises
        ------
        IndexError
            If either index is out of range for this chain or residue.

        See Also
        --------
        maws.structure.Structure.torsions : Lists a residue's turnable bonds.

        Notes
        -----
        A residue's bonds are stored as atom indices counted from the start of
        that residue, with negative values counting back from its end. Both
        are resolved here and shifted to indices counted across the whole
        chain.
        """
        if residue_index < 0:
            residue_index += len(self.sequence_array)

        residue_name = self.sequence_array[residue_index]
        residue_length = self.structure.residue_length[residue_name]
        offset = self.residues_start[residue_index]
        bond = self.structure.rotating_elements[residue_name][residue_element_index]

        # Build a new list rather than resolving in place: rotating_elements
        # holds one entry per residue type, shared by every chain using it.
        self.rotate_element(
            [(atom + residue_length if atom < 0 else atom) + offset for atom in bond],
            angle,
            reverse=reverse,
        )

    # ---- Compatibility helpers---------------------------

    def rotate_historic_element(self, historic_element, angle: float):
        """
        Rotate using an element triple that was recorded against an **older**
        ``start`` (e.g., before a prepend/append).

        Parameters
        ----------
        historic_element : list[int | None]
            Element triple based on ``start_history`` rather than the current
            ``start``.
        angle : float
            Rotation angle in radians.
        """
        if historic_element[2]:
            self.rotate_element(
                [
                    historic_element[0] + self.start_history - self.start,
                    historic_element[1] + self.start_history - self.start,
                    historic_element[2] + self.start_history - self.start,
                ],
                angle,
            )
        else:
            self.rotate_element(
                [
                    historic_element[0] + self.start_history - self.start,
                    historic_element[0] + self.start_history - self.start,
                    None,
                ],
                angle,
            )

    def rotate_in_historic_residue(
        self, historic_index: int, element_index: int, angle: float
    ):
        """
        Rotate a torsion using a residue index captured **before** a prepend.

        Parameters
        ----------
        historic_index : int
            Residue index based on the old sequence (prior to prepends).
        element_index : int
            Torsion index within that residue type.
        angle : float
            Rotation angle in radians.

        Notes
        -----
        If residues were prepended, all old indices shift right by
        ``len(prepend_history)``. This method adjusts the index and forwards to
        :meth:`rotate_in_residue`.
        """
        offset = len(self.prepend_history)
        self.rotate_in_residue(historic_index + offset, element_index, angle)

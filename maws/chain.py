# maws/chain.py
from __future__ import annotations

import copy
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # Only for static type checkers; avoids runtime circular import
    from .complex import Complex


class Chain:
    """
    Representation of a single polymer chain inside a :class:`~maws.complex.Complex`.

    This class keeps sequence-related state (alias and canonical forms), per-residue
    atom-offsets, and convenience helpers for applying rotations/translations to
    either the whole chain or to sub-elements (torsions) defined in the
    associated :class:`~maws.structure.Structure`.  **It does not own coordinates**;
    geometric operations are delegated to the owning :class:`~maws.complex.Complex`,
    which holds ``positions`` and the OpenMM machinery.

    Attributes
    ----------
    id : int
        Stable identifier for this chain within ``Complex.chains`` (0-based).
    start : int
        Global atom index (0-based) where this chain begins inside
        ``Complex.positions`` (the complex stores all atoms from all chains in
        one flat list).
    start_history : int
        Snapshot of the previous ``start`` used by :meth:`Complex.rebuild`
        to splice old and new coordinate blocks correctly.
    complex : Complex
        Back-reference to the owning complex. All geometry is ultimately executed
        on the complex since it owns ``positions``, ``topology``, etc.
    residues_start : list[int]
        Per-residue atom offsets **relative to** ``start``.  For residue lengths
        ``[33, 31, 32]`` this becomes ``[0, 33, 64]``.
    length : int
        Total number of atoms in this chain (sum over per-residue atom counts
        from ``Structure.residue_length``).
    length_history : int
        Snapshot of the previous ``length`` used during :meth:`Complex.rebuild`.
    element : list[int]
        Convenience triple ``[start_atom, bond_atom, end_atom_exclusive]`` used
        for whole-chain transforms (``end`` is exclusive, like slicing).
    structure : Structure
        Template bank providing residue metadata and topology rules:
        ``residue_length``, ``alias``, ``rotating_elements``, and ``connect``.
    alias_sequence : str
        Space-separated alias tokens (human-facing), e.g. ``"G A U"``.
    sequence : str
        Space-separated canonical residue names used by LEaP, e.g. ``"G5 A U3"``.
    sequence_array : list[str]
        Tokenized canonical sequence (``sequence.split()``).
    alias_sequence_array : list[str]
        Tokenized alias sequence (``alias_sequence.split()``).
    append_history : list[str]
        Canonical residues appended on the right (3′) in the most recent edit.
    prepend_history : list[str]
        Canonical residues prepended on the left (5′) in the most recent edit.
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
        """
        Rotate a chain-local element by ``angle`` radians.

        Parameters
        ----------
        element : list[int | None]
            Triple ``[start, bond, end_or_None]`` **relative to this chain**.
            - ``start`` : first atom index (inclusive)
            - ``bond``  : atom forming the rotation axis with ``start``
            - ``end``   : one-past-last atom index (exclusive). If ``None``,
              the element extends to the end of the chain.
        angle : float
            Rotation angle in radians.
        reverse : bool, default=False
            If ``True`` and ``end`` is ``None``, rotate the **complement** of
            the specified range (the part *not* selected).

        Notes
        -----
        This method converts chain-local indices to **global** atom indices by
        adding ``self.start``, then delegates to :meth:`Complex.rotate_element`.
        """
        revised_element = element[:]
        rev = reverse
        if rev:
            # Special handling: encode "complement" rotation into indices
            if revised_element[2] is None:
                revised_element[2] = 0
            else:
                revised_element[2] = revised_element[1]
        rev = False

        if len(revised_element) == 3 and revised_element[2] is not None:
            revised_element = [idx + self.start for idx in revised_element]
            self.complex.rotate_element(revised_element, angle, reverse=rev)
        elif len(revised_element) == 3 and revised_element[2] is None:
            revised_element = [
                revised_element[0] + self.start,
                revised_element[1] + self.start,
                self.length + self.start,
            ]
            self.complex.rotate_element(revised_element, angle, reverse=rev)
        else:
            raise ValueError("Rotable element contains too many or too few components!")

    def rotate_in_residue(
        self,
        residue_index: int,
        residue_element_index: int,
        angle: float,
        reverse: bool = False,
    ):
        """
        Rotate one of the template-defined torsions **inside a specific residue**.

        Parameters
        ----------
        residue_index : int
            Index of the residue within this chain. Negative values count from
            the end (Python-style).
        residue_element_index : int
            Which torsion to rotate for that residue **type**; an index into
            ``Structure.rotating_elements[resname]``.
        angle : float
            Rotation angle in radians.
        reverse : bool, default=False
            If ``True`` and the torsion has ``end=None``, rotate the complement.

        Notes
        -----
        The torsion triple is stored **residue-local** in the :class:`Structure`.
        This method normalizes any negative atom indices using the residue's
        length, translates them to **chain-local** using ``residues_start``, and
        then calls :meth:`rotate_element` to perform the actual rotation.
        """
        rev = reverse
        revised_residue_index = residue_index
        if residue_index < 0:
            revised_residue_index += len(self.sequence_array)

        element = self.structure.rotating_elements[
            self.sequence_array[revised_residue_index]
        ][residue_element_index]

        # Normalize possibly negative indices within the residue
        for i in range(len(element)):
            if element[i] and element[i] < 0:
                element[i] += self.structure.residue_length[
                    self.sequence_array[revised_residue_index]
                ]

            if element[2] is None:
                revised_element = [
                    element[0] + self.residues_start[revised_residue_index],
                    element[1] + self.residues_start[revised_residue_index],
                    None,
                ]
            elif element[2] == 0:
                revised_element = [
                    element[0] + self.residues_start[revised_residue_index],
                    element[1] + self.residues_start[revised_residue_index],
                    element[2],
                ]
            else:
                revised_element = [
                    element[0] + self.residues_start[revised_residue_index],
                    element[1] + self.residues_start[revised_residue_index],
                    element[2] + self.residues_start[revised_residue_index],
                ]
                rev = False

            self.rotate_element(revised_element, angle, reverse=rev)

    # ---- Compatibility helpers (historic indices) ---------------------------

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

    # ---- Whole-chain transforms (delegated) --------------------------------

    def rotate_global(self, axis, angle: float):
        """
        Rotate the **entire chain** around a given axis by ``angle`` radians.

        Parameters
        ----------
        axis : array-like or openmm.Vec3
            Rotation axis (direction only matters; magnitude is ignored).
        angle : float
            Rotation angle in radians.

        Notes
        -----
        Delegates to :meth:`Complex.rotate_global` using this chain's
        ``element`` triple.
        """
        self.complex.rotate_global(self.element, axis, angle)

    def translate_global(self, shift):
        """
        Translate the **entire chain** by a 3-vector.

        Parameters
        ----------
        shift : array-like or openmm.unit.Quantity
            Displacement vector. If given as a quantity, its unit should be
            compatible with Ångström.

        Notes
        -----
        Delegates to :meth:`Complex.translate_global` using this chain's
        ``element`` triple.
        """
        self.complex.translate_global(self.element, shift)

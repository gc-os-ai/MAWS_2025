# maws/chain.py
from __future__ import annotations

import copy
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # Only for type checkers; avoids runtime circular import
    from .complex import Complex


class Chain:
    """
    A polymer chain backed by a Structure template. Maintains sequence state,
    residue start indices, and provides rotation/translation helpers that
    delegate to the parent Complex.
    """

    def __init__(
        self,
        Complex: Complex,
        Structure,
        sequence: str | None = None,
        start: int = 0,
        ID: int = 0,
    ):
        self.id = ID
        self.start = start  # global atom index (coming from Complex.positions?)
        self.start_history = start  # Snapshot of the old start used when remapping coordinates during Complex.rebuild(). It’s a single integer, not an array. The mapping logic uses this number to splice old vs new atom blocks correctly.
        self.complex = Complex  # No. It’s a reference to the owning Complex instance. Chain calls back into it (e.g., self.complex.rotate_element(...)) because the Complex owns positions, topology, etc.
        self.residues_start: list[
            int
        ] = []  # Per-residue atom offsets inside this chain: the first atom index of each residue, relative to self.start. Example: if residue lengths are [33, 31, 32], then residues_start becomes [0, 33, 64].
        self.length = 0  # Total atom count in this chain (sum of per-residue atom counts from Structure.residue_length).
        self.length_history = self.length  # Snapshot of the old atom count, used to splice old coordinates into the rebuilt system. (still unclear to me).
        self.element = [self.start, self.start + 1, self.start + self.length]
        # A rotation/translation “window”: [start_atom, bond_atom, end_atom_exclusive].
        # start_atom = first atom index of the element
        # bond_atom = the pivot partner when constructing the rotation axis
        # end_atom_exclusive = one past the last atom (like Python slicing)
        # For whole-chain ops, this targets all atoms of the chain.

        self.structure = Structure  # The Structure object, is basically all the information is coming from.t stores residue names, atom counts, alias rules, rotating elements, and connectivity. Chain uses it to translate aliases, compute lengths, and look up rotation elements.
        self.alias_sequence = ""  # Yes, a string like "G A U" using alias tokens (start/middle/end aware). Structure.translate(...) converts this into the canonical residue-name sequence for LEaP.
        self.sequence = ""  # The translated sequence string using canonical residue names (e.g., "G5 A U3"). This is what LEaP expects.
        self.sequence_array: list[str] = []
        self.alias_sequence_array: list[
            str
        ] = []  # The alias tokens split. We keep both alias and canonical forms (and their tokenized lists) because: you append/prepend with aliases (human-friendly),but builds/lengths/rotations rely on canonical names and per-residue metadata.
        self.append_history: list[str] = []
        self.prepend_history: list[
            str
        ] = []  # Which residues were added on the right/left in the most recent change (canonical names).

        if sequence:  # Because Chain can be created empty ('') and filled later, or with an initial sequence.
            self.alias_sequence = sequence
            self.sequence = self.structure.translate(
                self.alias_sequence
            )  # Convert alias tokens to canonical residue names (start/middle/end forms applied).
            self.sequence_array = self.sequence.split(" ")
            self.alias_sequence_array = self.alias_sequence.split(" ")
            self.length = sum(
                map(self.structure.residue_length.__getitem__, self.sequence_array)
            )  # Compute total atoms by summing each residue’s atom count from the template.
            self.length_history = self.length
            tally = 0
            for residue in self.sequence_array:
                self.residues_start.append(tally)
                tally += self.structure.residue_length[residue]
            self.element = [self.start, self.start + 1, self.start + self.length]

    def update_chains(self):
        """Recompute this chain's length/starts and update downstream chains' starts."""  ## I have no clue what does this mean?
        length = self.length
        self.length = sum(
            map(self.structure.residue_length.__getitem__, self.sequence_array)
        )  # I have no clue what line does this code do??
        self.residues_start = []
        tally = 0
        for residue in self.sequence_array:
            self.residues_start.append(tally)
            tally += self.structure.residue_length[residue]
        self.element = [self.start, self.start + 1, self.start + self.length]
        start = copy.deepcopy(self.start)  # Why are we even creating deep copy?

        for chain in self.complex.chains:
            chain.start_history = chain.start
            if chain.start >= start:
                chain.start += self.length - length
                chain.start_history += 0
                chain.element = [
                    chain.start,
                    chain.start + 1,
                    chain.start + chain.length,
                ]

        self.start -= self.length - length
        self.start_history -= 0
        self.element = [self.start, self.start + 1, self.start + self.length]

    def create_sequence(
        self, sequence: str
    ):  # why overwrite this , I don't understnd the usecse?
        """Overwrite the chain's sequence (rebuild Complex afterward)."""
        alias_sequence_array = sequence.split(" ")
        sequence_array = self.structure.translate(sequence).split(" ")
        for letter in sequence_array:
            if letter not in self.structure.residue_names:
                raise ValueError("Residue not defined! CANNOT create sequence!")
        self.alias_sequence = sequence
        self.sequence = self.structure.translate(self.alias_sequence)
        self.alias_sequence_array = alias_sequence_array
        self.sequence_array = sequence_array
        self.update_chains()
        self.start_history = self.start
        self.length_history = self.length
        self.sequence_array_history = self.sequence_array

    def append_sequence(self, sequence: str):  # What does this mean?
        """Append a residue (alias) to the right side and track history for rebuild."""
        start_history = self.start
        length = self.length
        seq_ar_length = len(self.sequence_array)
        self.create_sequence(" ".join(self.alias_sequence_array[:] + [sequence]))
        self.length_history = length
        self.start_history = start_history
        self.prepend_history = []
        self.append_history = self.sequence_array[seq_ar_length:]

    def prepend_sequence(self, sequence: str):
        """Prepend a residue (alias) to the left side and track history for rebuild."""
        length = self.length
        seq_ar_length = len(self.sequence_array)
        self.create_sequence(" ".join([sequence] + self.alias_sequence_array[:]))
        self.length_history = length
        self.start_history = self.start + self.length - length
        self.prepend_history = self.sequence_array[
            : len(self.sequence_array) - seq_ar_length
        ]

    def rotate_element(
        self, element, angle: float, reverse: bool = False
    ):  # What are we rotating arrays, but how are we doing all this ?
        """
        Rotate the specified element [start, bond, end] by angle (radians).
        If reverse=True and end is None, rotate the complement instead.
        """
        revised_element = element[:]
        rev = reverse
        if rev:
            if revised_element[2] is None:
                revised_element[2] = 0
            else:
                revised_element[2] = revised_element[1]
        rev = False

        if len(revised_element) == 3 and revised_element[2] is not None:
            revised_element = [index + self.start for index in revised_element]
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
        """Rotate one of the template-defined rotating elements inside a residue."""
        rev = reverse
        revised_residue_index = residue_index
        if residue_index < 0:
            revised_residue_index += len(self.sequence_array)
        element = self.structure.rotating_elements[
            self.sequence_array[revised_residue_index]
        ][residue_element_index]
        # normalize negative indices
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

    # deprecated compatibility helpers
    def rotate_historic_element(self, historic_element, angle: float):
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
        offset = len(self.prepend_history)
        self.rotate_in_residue(historic_index + offset, element_index, angle)

    def rotate_global(self, axis, angle: float):
        """Rotate the entire chain around the given axis by angle (radians)."""
        self.complex.rotate_global(self.element, axis, angle)

    def translate_global(self, shift):
        """Translate the entire chain by the given 3-vector (OpenMM Quantity)."""
        self.complex.translate_global(self.element, shift)

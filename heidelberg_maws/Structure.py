"""
structure.py

Defines residue-based chemistry used by MAWS:
- residue metadata (length, alias mapping)
- polymer connectivity rules
- rotating elements (handled dynamically at runtime)
- backbone anchor indices (normalized here)

Notes
-----
- LEaP init_string is generated deterministically from `residue_path` and names.
- Rotating elements may use negative indices (relative to residue length). We keep
  them as-is because normalization is performed later during rotation.
- Backbone elements are normalized here to positive indices because later code
  uses them directly without further normalization.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple


class Structure(object):
    """
    Container for residue templates and per-residue topology rules.

    Parameters
    ----------
    residue_names : list[str]
        Names of residue templates (e.g., ["A", "C", "G", "U"] or ["LIG"]).
    residue_length : list[int], optional
        Number of atoms per residue, aligned to `residue_names`.
    rotating_elements : list[tuple], optional
        Triples of (residue_name, start, bond, end_or_None). Negative indices are
        allowed and are interpreted relative to residue length at runtime.
    backbone_elements : list[tuple], optional
        Tuples (residue_name, start, middle_pre, bond, middle_post, end). Negative
        indices are normalized here to positive, since later code uses them directly.
    connect : list[list], optional
        Per-residue connectivity: [[append_first, append_last], [prepend_last, prepend_first],
        append_bond_length, prepend_bond_length].
    residue_path : str | None
        Directory containing `<name>.lib` and `<name>.frcmod` for each residue. If "",
        use current directory. If None, no LEaP init_string is produced.
    alias : list[list[str]] | None
        Alias mapping entries of the form [residue_name, alone, start, middle, end].
        If omitted, defaults to [name, name, name, name] for each residue.
    """

    def __init__(
        self,
        residue_names: Sequence[str],
        residue_length: Optional[Sequence[int]] = None,
        rotating_elements: Optional[Sequence[Tuple[str, int, int, Optional[int]]]] = None,
        backbone_elements: Optional[Sequence[Tuple[str, int, int, int, int, int]]] = None,
        connect: Optional[Sequence[Sequence]] = None,
        residue_path: Optional[str] = None,
        alias: Optional[Sequence[Sequence[str]]] = None,
    ):
        self.residue_names: List[str] = list(residue_names)
        self.residue_path: Optional[str] = residue_path

        # ---- LEaP init string (loadoff/loadamberparams) ---------------------
        self.init_string: str = ""
        if self.residue_path is not None:
            base = self.residue_path if self.residue_path != "" else "."
            for name in self.residue_names:
                self.init_string += (
                    f"loadoff {base}/{name}.lib\n"
                    f"loadamberparams {base}/{name}.frcmod\n"
                )

        # ---- residue length mapping -----------------------------------------
        self.residue_length: Dict[str, int] = defaultdict(lambda: 0)
        if residue_length:
            for idx, res in enumerate(self.residue_names):
                self.residue_length[res] = int(residue_length[idx])

        # ---- polymer connectivity rules -------------------------------------
        # default: head/tail indices and default bond lengths
        self.connect: Dict[str, list] = defaultdict(lambda: [[0, -1], [-2, 0], 1.6, 1.6])
        if connect:
            for idx, res in enumerate(self.residue_names):
                self.connect[res] = list(connect[idx])

        # ---- alias mapping (alone/start/middle/end) --------------------------
        self.alias: Dict[str, List[str]] = defaultdict(lambda: None)
        if self.residue_names:
            for res in self.residue_names:
                self.alias[res] = [res] * 4
        if alias:
            # alias entries like: [name, alone, start, middle, end]
            for elem in alias:
                self.alias[elem[0]] = elem[1:]

        # ---- rotating elements (kept with possible negative indices) ---------
        # default placeholder per residue
        self.rotating_elements: Dict[str, list] = defaultdict(lambda: None)
        for name in self.residue_names:
            self.rotating_elements[name] = [None]

        if rotating_elements:
            for residue, start, bond, end in rotating_elements:
                if self.rotating_elements[residue] == [None]:
                    self.rotating_elements[residue] = [[start, bond, end]]
                elif self.rotating_elements[residue] is None:
                    raise ValueError('Residue does not exist! CANNOT assign rotability!')
                else:
                    self.rotating_elements[residue].append([start, bond, end])

        # ---- backbone elements (normalized to positive indices) --------------
        # Stored as: {residue: [[start, middle_pre, bond], [middle_post, end]]}
        self.backbone_elements: Dict[str, list] = defaultdict(lambda: None)
        if backbone_elements:
            for residue, start, middle_pre, bond, middle_post, end in backbone_elements:
                # normalize negatives relative to residue length
                L = self.residue_length[residue]
                def norm(i: int) -> int:
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

    def add_rotation(self, residue_name: str, rotations, basestring=str):
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
        if self.rotating_elements[residue_name] == [None]:
            self.rotating_elements[residue_name] = []

        # Accept a single triple or a list of triples
        def is_triplet(x):
            return isinstance(x, (list, tuple)) and len(x) == 3

        if is_triplet(rotations):
            self.rotating_elements[residue_name].append(list(rotations))
        elif isinstance(rotations, (list, tuple)) and rotations and is_triplet(rotations[0]):
            for rot in rotations:
                if not is_triplet(rot):
                    raise ValueError("Each rotation must be a [start, bond, end] triple.")
                self.rotating_elements[residue_name].append(list(rot))
        else:
            raise ValueError("rotations must be a triple or a list of triples [start, bond, end].")

        return self.rotating_elements

    def translate(self, sequence: str) -> str:
        """
        Translate an alias sequence (with spaces) into a residue-name sequence
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
        sequence_array = sequence.split(' ')
        if len(sequence_array) == 1:
            return self.alias[sequence_array[0]][0]
        else:
            first = self.alias[sequence_array[0]][1]
            middles = [self.alias[name][2] for name in sequence_array][1:-1]
            last = self.alias[sequence_array[-1]][3]
            return " ".join([first] + middles + [last])

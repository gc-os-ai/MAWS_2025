"""structure.py – modern-Python rewrite of the original MAWS Structure helper."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

# --------------------------------------------------------------------------- #
# Type aliases
# --------------------------------------------------------------------------- #
RotTriplet   = Tuple[int, int, Optional[int]]
BackboneSpec = Tuple[int, int, int, int, int]
ConnectRule  = Tuple[Tuple[int, int], Tuple[int, int], float, float]
Alias        = Tuple[str, str, str, str]


# --------------------------------------------------------------------------- #
# Helper
# --------------------------------------------------------------------------- #
def _normalise(idx: int, length: int) -> int:
    """Convert negative indices to non-negative Pythonic positions."""
    return idx + length if idx < 0 else idx


# --------------------------------------------------------------------------- #
# Dataclass
# --------------------------------------------------------------------------- #
@dataclass
class Structure:
    """
    Lightweight container for residue metadata used by MAWS.
    """

    # Required ------------------------------------------------------------ #
    residue_names: Sequence[str]

    # Optional / advanced -------------------------------------------------- #
    residue_lengths: Mapping[str, int] | Sequence[int] = field(default_factory=dict)
    rotating_elements: Dict[str, List[RotTriplet]]     = field(default_factory=dict)
    backbone_elements: Dict[str, List[BackboneSpec]]   = field(default_factory=dict)
    connect: Dict[str, ConnectRule]                    = field(default_factory=dict)
    alias: Dict[str, Alias]                            = field(default_factory=dict)
    residue_path: str | Path | None                    = None

    # Derived -------------------------------------------------------------- #
    init_string: str                                   = field(init=False, repr=False)

    # ------------------------------------------------------------------ #
    # Life-cycle
    # ------------------------------------------------------------------ #
    def __post_init__(self) -> None:
        # --- Coerce residue_path ---------------------------------------- #
        if isinstance(self.residue_path, str):
            self.residue_path = Path(self.residue_path)

        # --- Coerce residue_lengths ------------------------------------- #
        if isinstance(self.residue_lengths, Sequence):
            # Turn [10, 12] into {"A": 10, "B": 12}
            self.residue_lengths = dict(zip(self.residue_names, self.residue_lengths))

        # --- Coerce rotating_elements (list ➝ dict) --------------------- #
        if isinstance(self.rotating_elements, list):
            rot_dict: Dict[str, List[RotTriplet]] = {n: [] for n in self.residue_names}
            for res, s, b, e in self.rotating_elements:
                rot_dict.setdefault(res, []).append((s, b, e))
            self.rotating_elements = rot_dict
        # --- Coerce backbone_elements (list ➝ dict) --------------------- #
        if isinstance(self.backbone_elements, list):
            bb_dict: Dict[str, List[BackboneSpec]] = {n: [] for n in self.residue_names}
            for res, *vals in self.backbone_elements:
                bb_dict.setdefault(res, []).append(tuple(vals))  # type: ignore[arg-type]
            self.backbone_elements = bb_dict

        # --- Coerce connect (list aligned with residue_names) ----------- #
        if isinstance(self.connect, list):
            if len(self.connect) != len(self.residue_names):
                raise ValueError("connect list length must match residue_names")
            self.connect = dict(zip(self.residue_names, self.connect))

        # --- Coerce alias (list ➝ dict) --------------------------------- #
        if isinstance(self.alias, list):
            self.alias = {elem[0]: tuple(elem[1:5]) for elem in self.alias}

        # 1.  Build LEaP initialisation fragment ------------------------- #
        self.init_string = self._build_leap_script()

        # 2.  Ensure every residue has dictionary slots ------------------ #
        for name in self.residue_names:
            self.rotating_elements.setdefault(name, [])
            self.backbone_elements.setdefault(name, [])
            self.connect.setdefault(name, ((0, -1), (-2, 0), 1.6, 1.6))
            # self.alias.setdefault(name, (name, f"{name}-", "-x-", f"-{name}"))
            
            # legacy MAWS default: no decoration unless user overrides
            self.alias.setdefault(name, (name, name, name, name))

        # 3.  Normalise negative indices in rotating_elements ------------ #
        for res, torsions in self.rotating_elements.items():
            length = self.residue_lengths.get(res, 0)
            self.rotating_elements[res] = [
                (
                    _normalise(start, length),
                    _normalise(bond,  length),
                    _normalise(end,   length) if end is not None else None,
                )
                for start, bond, end in torsions
            ]

        # 4.  Normalise negative indices in backbone_elements ------------ #
        for res, elems in self.backbone_elements.items():
            length = self.residue_lengths.get(res, 0)
            self.backbone_elements[res] = [
                (
                    _normalise(s,     length),
                    _normalise(mp,    length),
                    _normalise(b,     length),
                    _normalise(mpost, length),
                    _normalise(e,     length),
                )
                for s, mp, b, mpost, e in elems
            ]

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #
    def add_rotation(self, residue: str, rotation: str | Sequence[str]) -> None:
        """Append one or more torsion labels to *residue*."""
        if isinstance(rotation, str):
            rotation = [rotation]
        elif not all(isinstance(r, str) for r in rotation):
            raise TypeError("All rotation labels must be strings.")
        self.rotating_elements[residue].extend(rotation)

    def translate(self, sequence: str) -> str:
        """Convert a space-separated residue string to the aliased printable form."""
        tokens = sequence.split()
        if not tokens:
            raise ValueError("Empty sequence supplied.")

        if len(tokens) == 1:
            return self.alias[tokens[0]][0]

        prefix  = self.alias[tokens[0]][1]
        infixes = (self.alias[t][2] for t in tokens[1:-1])
        suffix  = self.alias[tokens[-1]][3]
        return " ".join([prefix, *infixes, suffix])

    # ------------------------------------------------------------------ #
    # Private helpers
    # ------------------------------------------------------------------ #
    def _build_leap_script(self) -> str:
        """Return the text block that LEaP should execute."""
        if self.residue_path is None:
            return ""

        lines: List[str] = []
        for name in self.residue_names:
            base = self.residue_path / name
            lines.append(f"loadoff {base.with_suffix('.lib')}")
            # lines.append(f"loadamberparams {base.with_suffix('.frcmod')}")
            frc = base.with_suffix('.frcmod')
            if frc.exists():                       # ← only if the file exists
                lines.append(f"loadamberparams {frc}")
        return "\n".join(lines)


# # --------------------------------------------------------------------------- #
# # Demo / quick smoke-test
# # --------------------------------------------------------------------------- #
# if __name__ == "__main__":
#     s = Structure(
#         residue_names   = ["A", "B"],
#         residue_lengths = [10, 12],               # list is fine, now coerced
#         rotating_elements = {"A": [(0, 1, 2)]},
#         residue_path    = "./params",             # str is fine, now coerced
#     )
#     print(s.translate("A B A"))
#     print(s.init_string)

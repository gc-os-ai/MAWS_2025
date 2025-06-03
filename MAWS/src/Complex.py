"""
complex.py – container classes *Chain* and *Complex* used by the MAWS pipeline.

* 100 % Python 3.10+
* no more `simtk` namespace
* uses pathlib, snake_case, type hints
* keeps backward-compat camel-case aliases where old code relies on them
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np
import openmm as mm                  # low-level API
from openmm import app, unit

from MAWS.src.helpers   import from_angstrom
from MAWS.src.Prepare   import make_lib
import MAWS.src.Structure as _st


# ──────────────────────────────────────────────────────────────────────────────
# Chain
# ──────────────────────────────────────────────────────────────────────────────
class Chain:
    """Single polymer chain inside a :class:`Complex`."""

    def __init__(
        self,
        parent: "Complex",
        st: _st.Structure,
        sequence: str | None = None,
        *,
        start: int = 0,
        cid:   int = 0,
    ) -> None:
        self.parent:   Complex         = parent
        self.structure: _st.Structure  = st

        self.id     = cid
        self.start  = start            # index of first atom
        self.length = 0                # atom count (filled later)

        self.residue_offsets: List[int] = []      # atom index of each residue
        self.alias_sequence = self.sequence = ""
        self.alias_tokens:   List[str]   = []
        self.sequence_tokens: List[str]   = []

        # history (needed by old MAWS code – harmless if unused)
        self._start_prev   = start
        self._length_prev  = 0
        self._prepend_hist: List[str] = []
        self._append_hist: List[str] = []

        # whole-chain element used for global rigid moves
        self.element = [start, start + 1, start + self.length]

        if sequence:
            self._set_sequence(sequence)

    # ────────────────────────────────────────────────────────────────────
    # internal helpers
    # ────────────────────────────────────────────────────────────────────
    def _set_sequence(self, seq: str) -> None:
        """Populate all sequence-related bookkeeping fields."""
        self.alias_sequence  = seq
        self.alias_tokens    = seq.split()
        self.sequence        = self.structure.translate(seq)
        self.sequence_tokens = self.sequence.split()

        # length & residue offsets
        self.length = sum(self.structure.residue_lengths[r] for r in self.sequence_tokens)
        self.residue_offsets.clear()
        offset = 0
        for res in self.sequence_tokens:
            self.residue_offsets.append(offset)
            offset += self.structure.residue_lengths[res]

        # element & history
        self.element        = [self.start, self.start + 1, self.start + self.length]
        self._length_prev   = self.length
        self._start_prev    = self.start

    def _refresh_downstream_starts(self, delta: int) -> None:
        """Shift starts of *later* chains when our own length changed."""
        for ch in self.parent.chains:
            if ch is self:
                continue
            if ch.start >= self.start:
                ch.start   += delta
                ch.element = [x + delta for x in ch.element]

    # ────────────────────────────────────────────────────────────────────
    # MAWS-0.x compatibility helpers (rotate_element / rotate_in_residue)
    # ────────────────────────────────────────────────────────────────────
    def rotate_element(
        self,
        element: list[int | None],
        angle:   float,
        *,
        reverse: bool = False,
    ) -> None:
        """Rotate arbitrary element `[start, bond, end]` inside this chain."""
        elem = element[:]
        if elem[2] is None:                # open-ended range
            elem[2] = self.length
        if elem[2] <= elem[0]:             # legacy “reverse” trick
            elem[1], elem[2] = elem[2], elem[1]

        abs_elem = [self.start + x for x in elem]
        axis_vec = self.parent.positions[abs_elem[1]] - self.parent.positions[abs_elem[0]]

        self.parent.rotate_global(abs_elem, axis_vec, angle, reverse=reverse)

    def rotate_in_residue(
        self,
        res_idx:   int,
        rotor_idx: int,
        angle:     float,
        *,
        reverse:   bool = False,
    ) -> None:
        """Rotate *rotor_idx* of residue *res_idx* (legacy API)."""
        if res_idx < 0:
            res_idx += len(self.sequence_tokens)

        rotor_tbl = getattr(self.structure, "rotating_elements", None)
        if rotor_tbl is None:
            raise AttributeError("Structure lacks 'rotating_elements' table")

        try:
            rel_elem = rotor_tbl[self.sequence_tokens[res_idx]][rotor_idx]
        except (KeyError, IndexError) as exc:
            raise ValueError(f"Residue {self.sequence_tokens[res_idx]!r} has no rotor #{rotor_idx}") from exc

        res_start = self.start + self.residue_offsets[res_idx]
        abs_elem: list[int | None] = []
        for x in rel_elem:
            if x is None:
                abs_elem.append(None)
            elif x >= 0:
                abs_elem.append(res_start + x)
            else:                          # negative index → wrap from end
                abs_elem.append(res_start + self.structure.residue_lengths[self.sequence_tokens[res_idx]] + x)

        self.rotate_element(abs_elem, angle, reverse=reverse)

    # camel-case aliases for un-refactored scripts
    rotateElement   = rotate_element
    rotateInResidue = rotate_in_residue

    # ────────────────────────────────────────────────────────────────────
    # public sequence API
    # ────────────────────────────────────────────────────────────────────
    def create_sequence(self, seq: str) -> None:
        old_start, old_len      = self.start, self.length
        self._length_prev, self._start_prev = old_len, old_start
        self._set_sequence(seq)
        self._refresh_downstream_starts(self.length - old_len)

    def append_sequence(self, seq: str) -> None:
        self._append_hist = seq.split()
        self.create_sequence(self.alias_sequence + " " + seq)

    def prepend_sequence(self, seq: str) -> None:
        self._prepend_hist = seq.split()
        old_start = self.start
        self.create_sequence(seq + " " + self.alias_sequence)
        self.start   = old_start - sum(self.structure.residue_lengths[r] for r in self._prepend_hist)
        self.element = [self.start, self.start + 1, self.start + self.length]

    # rigid-body helpers
    def rotate_global(self, axis: unit.Quantity, angle: float) -> None:
        self.parent.rotate_global(self.element, axis, angle)

    def translate_global(self, shift: unit.Quantity) -> None:
        self.parent.translate_global(self.element, shift)

    def __repr__(self) -> str:            # for debugging
        return f"<Chain id={self.id} start={self.start} len={self.length} '{self.alias_sequence}'>"


# ──────────────────────────────────────────────────────────────────────────────
# Complex
# ──────────────────────────────────────────────────────────────────────────────
class Complex:
    """A set of chains plus an (optional) OpenMM Simulation used for scoring."""

    _TEMPLATE_PREAMBLE = "source {force}\nsource leaprc.gaff\n"

    def __init__(self, force_field: str = "leaprc.protein.ff14SB"):
        self._force_field = force_field
        self.chains: List[Chain] = []

        # OpenMM objects (set after build)
        self.prmtop:   app.AmberPrmtopFile | None = None
        self.inpcrd:   app.AmberInpcrdFile | None = None
        self.positions: List[unit.Quantity] | None = None
        self.topology: app.Topology | None       = None
        self.system:   mm.System | None          = None
        self.integrator: mm.Integrator | None    = None
        self.simulation: app.Simulation | None   = None

    # ────────────────────────────────────────────────────────────────────
    # chain management
    # ────────────────────────────────────────────────────────────────────
    def add_chain(self, seq: str, st: _st.Structure) -> Chain:
        start = sum(ch.length for ch in self.chains)
        ch    = Chain(self, st, sequence=seq, start=start, cid=len(self.chains))
        self.chains.append(ch)
        return ch

    def add_chain_from_pdb(
        self,
        pdb_path:      str | Path,
        protein_name:  str,
        *,
        lib_name:      str = "PDB",
        parameterized: bool = False,
    ) -> Chain:
        length = make_lib(pdb_path, lib_name,
                          parameterized=parameterized,
                          force_field=self._force_field)
        st = _st.Structure([lib_name],
                           residue_lengths=[length],
                           residue_path=Path(pdb_path).parent)
        return self.add_chain(lib_name, st)

    # ────────────────────────────────────────────────────────────────────
    # (re)building with tleap / OpenMM
    # ────────────────────────────────────────────────────────────────────
    def build(self, name: str, *, workdir: str | Path = ".") -> None:
        if not self.chains:
            raise ValueError("Cannot build an empty Complex")

        workdir = Path(workdir)
        leap_in = workdir / f"{name}.in"

        script: list[str] = [self._TEMPLATE_PREAMBLE.format(force=self._force_field)]
        script.extend(ch.structure.init_string for ch in self.chains)

        for idx, ch in enumerate(self.chains):
            if ch.sequence:
                script.append(f"CHAIN{idx} = sequence {{{ch.sequence}}}")

        union = " ".join(f"CHAIN{idx}" for idx, ch in enumerate(self.chains) if ch.sequence)
        script += [
            f"UNION = combine {{{union}}}",
            f"saveamberparm UNION {name}.prmtop {name}.inpcrd",
            "quit",
        ]
        leap_in.write_text("\n".join(script))

        subprocess.run(f"tleap -f {leap_in}", shell=True, check=True, cwd=workdir)

        # ---- OpenMM ---------------------------------------------------
        self.prmtop   = app.AmberPrmtopFile(str(workdir / f"{name}.prmtop"))
        self.inpcrd   = app.AmberInpcrdFile(str(workdir / f"{name}.inpcrd"))
        self.topology = self.prmtop.topology
        self.positions = list(self.inpcrd.positions)

        self.integrator = mm.LangevinIntegrator(
            300 * unit.kelvin, 1.0 / unit.picosecond, 0.00005 * unit.picoseconds
        )
        self.system = self.prmtop.createSystem(
            nonbondedCutoff=5 * unit.angstrom,
            nonbondedMethod=app.NoCutoff,
            constraints=None,
            implicitSolvent=app.OBC1,
        )
        self.simulation = app.Simulation(
            self.topology, self.system, self.integrator, mm.Platform.getPlatformByName("CPU")
        )

    # ❱❱❱ NEW – skinny back-compat wrapper
    def rebuild(self, name: str, *, workdir: str | Path = ".", exclusion: list[Chain] | None = None) -> None:
        """
        Legacy entry-point retained for MAWS-0.x scripts.

        Current MAWS2 workflow immediately overwrites positions after each
        rebuild, so a plain rebuild == build is sufficient.
        """
        self.build(name, workdir=workdir)

    # ────────────────────────────────────────────────────────────────────
    # geometry transforms
    # ────────────────────────────────────────────────────────────────────
    def rotate_global(
        self,
        element: Sequence[int],
        axis:    unit.Quantity,
        angle:   float,
        *,
        reverse: bool = False,
    ) -> None:
        if self.positions is None:
            raise RuntimeError("Call .build() first")

        pos   = self.positions

        # rotation always happens around the "bond" atom (element[1])
        pivot = pos[element[1]]

        if reverse:
            start_idx, end_idx = element[0], element[1]
        else:
            start_idx, end_idx = element[1], element[2]

        axis_np   = from_angstrom(axis)
        axis_np  /= np.linalg.norm(axis_np)
        x, y, z   = axis_np
        s, c      = np.sin(angle / 2.0), np.cos(angle / 2.0)

        rot = np.array([
            [2*(x*x-1)*s*s+1,  2*x*y*s*s-2*z*c*s, 2*x*z*s*s+2*y*c*s],
            [2*x*y*s*s+2*z*c*s,2*(y*y-1)*s*s+1,   2*z*y*s*s-2*x*c*s],
            [2*x*z*s*s-2*y*c*s,2*z*y*s*s+2*x*c*s, 2*(z*z-1)*s*s+1 ],
        ])

        for i in range(start_idx, end_idx):
            vec_np = from_angstrom(pos[i] - pivot)
            pos[i] = unit.Quantity(np.dot(vec_np, rot), unit.angstrom) + pivot

    def translate_global(self, element: Sequence[int], shift: unit.Quantity) -> None:
        if self.positions is None:
            raise RuntimeError("Call .build() first")
        for i in range(element[0], element[2]):
            self.positions[i] += shift

    # ────────────────────────────────────────────────────────────────────
    # energetics
    # ────────────────────────────────────────────────────────────────────
    def get_energy(self) -> Tuple[float, List[unit.Quantity]]:
        if self.positions is None or self.simulation is None:
            raise RuntimeError("Call .build() before scoring")
        self.simulation.context.setPositions(self.positions)
        state = self.simulation.context.getState(getEnergy=True)
        e_val = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        return e_val, self.positions

    def minimize(self, max_iter: int = 200) -> float:
        if self.simulation is None:
            raise RuntimeError("Call .build() before minimisation")
        self.simulation.context.setPositions(self.positions)
        self.simulation.minimizeEnergy(maxIterations=max_iter)
        self.positions = self.simulation.context.getPositions(asNumpy=True)
        return self.get_energy()[0]

    # ────────────────────────────────────────────────────────────────────
    def __repr__(self) -> str:           # pragma: no cover
        return f"<Complex chains={len(self.chains)}>"

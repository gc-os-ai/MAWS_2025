"""
complex.py

Core data structures for MAWS:
- Chain: a polymer chain built from residue templates with rotation utilities.
- Complex: a collection of chains; responsible for building AMBER prmtop/inpcrd
           via LEaP, caching builds, and running energies/MD via OpenMM.

Notes
-----
- External executables are invoked via tools.run/find_exe (no 'conda run').
- LEaP builds are cached in .maws_cache/ keyed by the exact input (FF sources,
  residue init strings, and sequences). Subsequent identical builds reuse cached
  prmtop/inpcrd files (no LEaP invocation).
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import List, Optional

import numpy as np
from openmm import unit
from openmm import app
import openmm as mm

from Prepare import makeLib
from helpers import angle as ang
from helpers import directed_angle as d_ang
from helpers import (angstrom, nostrom, kJ, noJ)
from Structure import *
from tools import find_exe, run


## Represents a molecule chain, comprising multiple residues
class Chain(object):
    """
    A polymer chain backed by a Structure template. Maintains sequence state,
    residue start indices, and provides rotation/translation helpers that
    delegate to the parent Complex.
    """

    def __init__(self, Complex, Structure, sequence: Optional[str] = None, start: int = 0, ID: int = 0):
        self.id = ID
        self.start = start
        self.start_history = start
        self.complex = Complex
        self.residues_start: List[int] = []
        self.length = 0
        self.length_history = self.length
        self.element = [self.start, self.start + 1, self.start + self.length]
        self.structure = Structure
        self.alias_sequence = ''
        self.sequence = ''
        self.sequence_array: List[str] = []
        self.alias_sequence_array: List[str] = []
        self.append_history: List[str] = []
        self.prepend_history: List[str] = []

        if sequence:
            self.alias_sequence = sequence
            self.sequence = self.structure.translate(self.alias_sequence)
            self.sequence_array = self.sequence.split(' ')
            self.alias_sequence_array = self.alias_sequence.split(' ')
            self.length = sum(map(self.structure.residue_length.__getitem__, self.sequence_array))
            self.length_history = self.length
            tally = 0
            for residue in self.sequence_array:
                self.residues_start.append(tally)
                tally += self.structure.residue_length[residue]
            self.element = [self.start, self.start + 1, self.start + self.length]

    def update_chains(self):
        """
        Recompute this chain's length/starts and update downstream chains' starts.
        Users should not call this directly; it's invoked by sequence mutators.
        """
        length = self.length
        self.length = sum(map(self.structure.residue_length.__getitem__, self.sequence_array))
        self.residues_start = []
        tally = 0
        for residue in self.sequence_array:
            self.residues_start.append(tally)
            tally += self.structure.residue_length[residue]
        self.element = [self.start, self.start + 1, self.start + self.length]
        start = copy.deepcopy(self.start)

        for chain in self.complex.chains:
            chain.start_history = chain.start
            if chain.start >= start:
                chain.start += self.length - length
                chain.start_history += 0
                chain.element = [chain.start, chain.start + 1, chain.start + chain.length]

        self.start -= self.length - length
        self.start_history -= 0
        self.element = [self.start, self.start + 1, self.start + self.length]

    def create_sequence(self, sequence: str):
        """
        Overwrite the chain's sequence.

        After calling this, you must call Complex.build()/rebuild() to realize
        the new topology/coordinates.
        """
        alias_sequence_array = sequence.split(' ')
        sequence_array = self.structure.translate(sequence).split(' ')
        for letter in sequence_array:
            if letter not in self.structure.residue_names:
                raise ValueError('Residue not defined! CANNOT create sequence!')
        self.alias_sequence = sequence
        self.sequence = self.structure.translate(self.alias_sequence)
        self.alias_sequence_array = alias_sequence_array
        self.sequence_array = sequence_array
        self.update_chains()
        self.start_history = self.start
        self.length_history = self.length
        self.sequence_array_history = self.sequence_array

    def append_sequence(self, sequence: str):
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
        self.prepend_history = self.sequence_array[:len(self.sequence_array) - seq_ar_length]

    def rotate_element(self, element, angle: float, reverse: bool = False):
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
            revised_element = [revised_element[0] + self.start,
                               revised_element[1] + self.start,
                               self.length + self.start]
            self.complex.rotate_element(revised_element, angle, reverse=rev)
        else:
            raise ValueError('Rotable element contains too many or too few components!')

    def rotate_in_residue(self, residue_index: int, residue_element_index: int, angle: float, reverse: bool = False):
        """
        Rotate one of the pre-defined rotating elements of a residue in this chain.

        Parameters
        ----------
        residue_index : int
            Index into the chain's residue list; negative indices count from the end.
        residue_element_index : int
            Index into Structure.rotating_elements[residue].
        angle : float
            Rotation angle in radians.
        reverse : bool
            If True and the element has end=None, rotate the complement.
        """
        rev = reverse
        revised_residue_index = residue_index
        if residue_index < 0:
            revised_residue_index += len(self.sequence_array)
        element = self.structure.rotating_elements[self.sequence_array[revised_residue_index]][residue_element_index]
        # normalize possibly negative element indices
        for i in range(len(element)):
            if element[i] and element[i] < 0:
                element[i] += self.structure.residue_length[self.sequence_array[revised_residue_index]]

            if element[2] is None:
                revised_element = [element[0] + self.residues_start[revised_residue_index],
                                   element[1] + self.residues_start[revised_residue_index], None]
            elif element[2] == 0:
                revised_element = [element[0] + self.residues_start[revised_residue_index],
                                   element[1] + self.residues_start[revised_residue_index],
                                   element[2]]
            else:
                revised_element = [element[0] + self.residues_start[revised_residue_index],
                                   element[1] + self.residues_start[revised_residue_index],
                                   element[2] + self.residues_start[revised_residue_index]]
                rev = False
            self.rotate_element(revised_element, angle, reverse=rev)

    # deprecated helpers kept for compatibility
    def rotate_historic_element(self, historic_element, angle: float):
        if historic_element[2]:
            self.rotate_element([historic_element[0] + self.start_history - self.start,
                                 historic_element[1] + self.start_history - self.start,
                                 historic_element[2] + self.start_history - self.start], angle)
        else:
            self.rotate_element([historic_element[0] + self.start_history - self.start,
                                 historic_element[0] + self.start_history - self.start,
                                 None], angle)

    def rotate_in_historic_residue(self, historic_index: int, element_index: int, angle: float):
        offset = len(self.prepend_history)
        self.rotate_in_residue(historic_index + offset, element_index, angle)

    def rotate_global(self, axis, angle: float):
        """Rotate the entire chain around the given axis by angle (radians)."""
        self.complex.rotate_global(self.element, axis, angle)

    def translate_global(self, shift):
        """Translate the entire chain by the given 3-vector (OpenMM Quantity)."""
        self.complex.translate_global(self.element, shift)


## Represents a complex containing multiple molecule chains.
class Complex(object):
    """
    Collection of Chain objects that can be built via LEaP into AMBER prmtop/inpcrd
    and simulated/queried via OpenMM.

    Parameters
    ----------
    force_field_aptamer : str
        LEaP 'source' line for the aptamer/nucleic acid force field.
    force_field_ligand : str
        LEaP 'source' line for the ligand/protein/small-molecule force field.
    conda_env : str
        Ignored (kept for backward compatibility).
    """

    def __init__(self, force_field_aptamer: str = "leaprc.RNA.OL3",
                 force_field_ligand: str = "leaprc.protein.ff19SB",
                 conda_env: str = "maws_p3"):
        self.build_string = f"""
                            source {force_field_aptamer}
                            source {force_field_ligand}
                            """
        self.prmtop = None
        self.inpcrd = None
        self.positions = None
        self.topology = None
        self.chains: List[Chain] = []
        self.system = None
        self.integrator = None
        self.simulation = None
        self.cenv = conda_env  # legacy, unused

    def add_chain(self, sequence: str, structure):
        """
        Append a Chain with the given sequence and Structure template.
        """
        if self.chains:
            start = sum([chain.length for chain in self.chains])
            chainID = len(self.chains)
        else:
            start = 0
            chainID = 0
        self.chains.append(Chain(self, structure, sequence=sequence, start=start, ID=chainID))

    def add_chain_from_PDB(self, pdb_path: str, force_field_aptamer: str, force_field_ligand: str,
                           structure=None, pdb_name: str = 'LIG', parameterized: bool = False):
        """
        Create a one-residue Structure from a PDB (or pre-parameterized building block),
        generate its .lib/.frcmod via makeLib, and add as a Chain.
        """
        # makeLib no longer accepts conda_env; tools must be on PATH
        length = makeLib(pdb_path, pdb_name,
                         force_field_aptamer=force_field_aptamer,
                         force_field_ligand=force_field_ligand,
                         parameterized=parameterized)
        path = str(Path(pdb_path).resolve().parent)
        structure = Structure([pdb_name], residue_length=[length], residue_path=path)
        self.add_chain(pdb_name, structure)

    # ---- LEaP build + cache -------------------------------------------------

    def _build_cache_key(self) -> str:
        """
        Compute a deterministic cache key from the build inputs:
        - current build preamble (FF sources),
        - each chain's residue init_string (loadoff/loadamberparams),
        - sequences of all chains.
        """
        payload = {
            "build": " ".join(self.build_string.split()),  # normalize whitespace
            "inits": [("".join(ch.structure.init_string.split())) for ch in self.chains],
            "seqs": [ch.sequence for ch in self.chains],
        }
        return hashlib.sha1(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def build(self, target_path: str = "", file_name: str = "out"):
        """
        Materialize the current set of chains into AMBER prmtop/inpcrd using LEaP.

        - If an identical build was done before, reuses cached files from .maws_cache/.
        - Otherwise writes a temporary LEaP input and runs `tleap -f ...`.

        Parameters
        ----------
        target_path : str
            Path prefix where LEaP would write outputs (kept for compatibility).
            Cache files are stored under .maws_cache regardless.
        file_name : str
            Base name for output files (prmtop/inpcrd). Ignored by cache when reusing.
        """
        if not self.chains:
            raise ValueError('Empty Complex! CANNOT build!')

        # Assemble LEaP input
        build_string_base = self.build_string
        leap_str = [self.build_string]
        for chain in self.chains:
            leap_str.append(chain.structure.init_string)
        for index, chain in enumerate(self.chains):
            if chain.sequence:
                leap_str.append(f"CHAIN{index} = sequence {{{chain.sequence}}}")
        chain_names = [f"CHAIN{idx}" for idx, ch in enumerate(self.chains) if ch.sequence]
        chain_string = " ".join(chain_names)
        leap_str.append(f"UNION = combine {{{chain_string}}}")
        out_prefix = f"{target_path}{file_name}"
        leap_str.append(f"saveamberparm UNION {out_prefix}.prmtop {out_prefix}.inpcrd")
        leap_str.append("quit")
        leap_input = "\n".join(leap_str)

        # Caching
        cache_dir = Path(".maws_cache")
        cache_dir.mkdir(exist_ok=True)
        key = self._build_cache_key()
        cache_prm = cache_dir / f"{key}.prmtop"
        cache_crd = cache_dir / f"{key}.inpcrd"

        # Try cache
        if not (cache_prm.exists() and cache_crd.exists()):
            # Write LEaP input near the intended output for easier debugging
            in_file = Path(f"{target_path}{file_name}.in")
            in_file.write_text(leap_input)
            # Invoke tleap directly (no conda run)
            run([find_exe("tleap"), "-f", str(in_file)])
            # Move results into cache
            produced_prm = Path(f"{out_prefix}.prmtop")
            produced_crd = Path(f"{out_prefix}.inpcrd")
            if not (produced_prm.exists() and produced_crd.exists()):
                raise RuntimeError("LEaP did not produce expected .prmtop/.inpcrd outputs.")
            produced_prm.replace(cache_prm)
            produced_crd.replace(cache_crd)

        # Load from cache
        self.build_string = build_string_base  # restore original
        self.prmtop = app.AmberPrmtopFile(str(cache_prm))
        self.inpcrd = app.AmberInpcrdFile(str(cache_crd))
        self.topology = self.prmtop.topology
        self.positions = self.inpcrd.positions
        self.integrator = mm.LangevinIntegrator(300. * unit.kelvin, 1. / unit.picosecond, 0.002 * unit.picoseconds)
        self.system = self.prmtop.createSystem(nonbondedCutoff=5 * unit.angstrom,
                                               nonbondedMethod=app.NoCutoff,
                                               constraints=None, implicitSolvent=app.OBC1)
        self.simulation = app.Simulation(self.topology, self.system, self.integrator)

    # ---- Geometry utilities -------------------------------------------------

    def rebuild(self, target_path: str = "", file_name: str = "out", exclusion: list = []):
        """
        Rebuild prmtop/inpcrd/topology/positions after sequence changes, attempting to
        preserve coordinates of atoms outside the modified regions.
        """
        old_positions = self.positions[:]
        self.build(target_path=target_path, file_name=file_name)

        for index, chain in enumerate(self.chains):
            if chain in exclusion:
                continue

            pre_positions = self.positions[chain.start:chain.start_history]
            chain_positions = old_positions[chain.start:chain.start + chain.length_history]
            post_positions = self.positions[chain.start_history + chain.length_history:chain.start + chain.length]

            if len(pre_positions) != 0 and chain.prepend_history:
                # ---- fix prepended atoms ----
                pre_positions = self.positions[chain.start:chain.start_history + 1]
                pre_vector = self.positions[chain.start_history + chain.structure.connect[chain.prepend_history[-1]][1][0]] - self.positions[chain.start_history + 1]
                old_pre_vector = old_positions[chain.start] - old_positions[chain.start + 1]
                angle = -ang(nostrom(pre_vector), nostrom(old_pre_vector))
                axis = np.cross(np.asarray(nostrom(pre_vector)), np.asarray(nostrom(old_pre_vector)))
                if all(axis == np.zeros(3)):
                    axis = np.array([1., 0., 0.])
                    angle = 0
                else:
                    axis /= np.linalg.norm(axis)
                x, y, z = axis
                phi_2 = angle / 2.
                pos = pre_positions[:]
                shift_forward = mm.Vec3(0, 0, 0) * unit.angstroms - pos[-1 + chain.structure.connect[chain.prepend_history[-1]][1][0]]
                s = np.math.sin(phi_2)
                c = np.math.cos(phi_2)
                rot = np.array([[2 * (np.power(x, 2) - 1) * np.power(s, 2) + 1, 2 * x * y * np.power(s, 2) - 2 * z * c * s,
                                 2 * x * z * np.power(s, 2) + 2 * y * c * s],
                                [2 * x * y * np.power(s, 2) + 2 * z * c * s, 2 * (np.power(y, 2) - 1) * np.power(s, 2) + 1,
                                 2 * z * y * np.power(s, 2) - 2 * x * c * s],
                                [2 * x * z * np.power(s, 2) - 2 * y * c * s, 2 * z * y * np.power(s, 2) + 2 * x * c * s,
                                 2 * (np.power(z, 2) - 1) * np.power(s, 2) + 1]])

                for j in range(0, len(pos)):
                    pos[j] += shift_forward

                # bond-length correction for the new connection
                shift_back = chain_positions[chain.structure.connect[chain.sequence_array[len(chain.prepend_history)]][1][1]]
                pre_bond_shift = (chain.structure.connect[chain.prepend_history[-1]][2]) * old_pre_vector / np.linalg.norm(
                    np.asarray(nostrom(old_pre_vector))) - old_pre_vector

                for j in range(0, len(pos)):
                    roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)), rot)
                    pos[j] = mm.Vec3(roted[0], roted[1], roted[2]) * unit.angstrom
                    pos[j] += shift_back + pre_bond_shift

                pre_positions = pos[:]
                chain_positions[0] += pre_bond_shift

                self.positions = self.positions[:chain.start] + pre_positions[:] + chain_positions[1:] + self.positions[chain.start + chain.length:]

            if len(post_positions) != 0 and chain.append_history:
                # ---- fix appended atoms ----
                post_positions = self.positions[chain.start_history + chain.length_history - 1:chain.start_history + chain.length]
                post_vector = self.positions[chain.start_history + chain.length_history - 1] - self.positions[chain.start_history + chain.length_history - 2]
                old_post_vector = old_positions[chain.start_history + chain.length_history - 1] - old_positions[chain.start_history + chain.length_history - 2]
                angle = -ang(nostrom(post_vector), nostrom(old_post_vector))
                axis = np.cross(np.asarray(nostrom(post_vector)), np.asarray(nostrom(old_post_vector)))
                if all(axis == np.zeros(3)):
                    axis = np.array([1., 0., 0.])
                    angle = 0.
                else:
                    axis /= np.linalg.norm(axis)
                x, y, z = axis
                phi_2 = angle / 2.
                pos = post_positions[:]
                shift_forward = mm.Vec3(0, 0, 0) * unit.angstroms - pos[chain.structure.connect[chain.append_history[0]][0][0]]
                s = np.math.sin(phi_2)
                c = np.math.cos(phi_2)
                rot = np.array([[2 * (np.power(x, 2) - 1) * np.power(s, 2) + 1, 2 * x * y * np.power(s, 2) - 2 * z * c * s,
                                 2 * x * z * np.power(s, 2) + 2 * y * c * s],
                                [2 * x * y * np.power(s, 2) + 2 * z * c * s, 2 * (np.power(y, 2) - 1) * np.power(s, 2) + 1,
                                 2 * z * y * np.power(s, 2) - 2 * x * c * s],
                                [2 * x * z * np.power(s, 2) - 2 * y * c * s, 2 * z * y * np.power(s, 2) + 2 * x * c * s,
                                 2 * (np.power(z, 2) - 1) * np.power(s, 2) + 1]])

                for j in range(0, len(pos)):
                    pos[j] += shift_forward

                post_bond_shift = (chain.structure.connect[chain.append_history[0]][2]) * old_post_vector / np.linalg.norm(
                    np.asarray(nostrom(old_post_vector))) - old_post_vector
                shift_back = chain_positions[chain.structure.connect[chain.sequence_array[-len(chain.append_history)]][0][1]]

                for pos_idx, pos_elem in enumerate(pos):
                    roted = np.dot(np.array(pos_elem.value_in_unit(unit.angstrom)), rot)
                    pos[pos_idx] = mm.Vec3(roted[0], roted[1], roted[2]) * unit.angstrom
                    pos[pos_idx] += shift_back + post_bond_shift

                post_positions = pos[:]
                chain_positions[-1] += post_bond_shift
                self.positions = self.positions[:chain.start] + chain_positions[:-1] + post_positions[:] + self.positions[chain.start + chain.length:]

            if not (chain.append_history or chain.prepend_history):
                self.positions = (self.positions[:chain.start]
                                  + old_positions[chain.start_history:chain.start_history + chain.length_history]
                                  + self.positions[chain.start + chain.length:])

    def rotate_element(self, element, angle: float, reverse: bool = False):
        """Rotate a contiguous element [start, bond, end] by angle (radians)."""
        revised_element = element[:]
        if not self.positions:
            raise ValueError('This Complex contains no positions! You CANNOT rotate!')
        pos = self.positions[:]
        vec_a = (pos[revised_element[1]] - pos[revised_element[0]])
        if revised_element[2] <= revised_element[0]:
            revised_element_1 = revised_element[1]
            revised_element[1] = revised_element[2]
            revised_element[2] = revised_element_1
        self.rotate_global(revised_element, vec_a, angle, reverse=reverse, glob=False)

    def rotate_global(self, element, axis, angle: float, reverse: bool = False, glob: bool = True):
        """Rotate either whole chain (glob=True) or sub-element (glob=False) around axis by angle."""
        if not self.positions:
            raise ValueError('This Complex contains no positions! You CANNOT rotate!')
        x, y, z = np.asarray(nostrom(axis)) / (np.linalg.norm(np.asarray(nostrom(axis))))
        phi_2 = angle / 2.
        pos = self.positions[:]
        starting_index = 0 if glob else 1
        shift_forward = mm.Vec3(0, 0, 0) * unit.angstroms - pos[element[2] if reverse else element[starting_index]]
        s = np.math.sin(phi_2)
        c = np.math.cos(phi_2)
        rot = np.array([[2 * (np.power(x, 2) - 1) * np.power(s, 2) + 1, 2 * x * y * np.power(s, 2) - 2 * z * c * s, 2 * x * z * np.power(s, 2) + 2 * y * c * s],
                        [2 * x * y * np.power(s, 2) + 2 * z * c * s, 2 * (np.power(y, 2) - 1) * np.power(s, 2) + 1, 2 * z * y * np.power(s, 2) - 2 * x * c * s],
                        [2 * x * z * np.power(s, 2) - 2 * y * c * s, 2 * z * y * np.power(s, 2) + 2 * x * c * s, 2 * (np.power(z, 2) - 1) * np.power(s, 2) + 1]])

        for j in range(element[starting_index], element[2]):
            pos[j] += shift_forward
        for j in range(element[starting_index], element[2]):
            roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)), rot)
            pos[j] = mm.Vec3(roted[0], roted[1], roted[2]) * unit.angstrom
            pos[j] -= shift_forward

        self.positions = pos[:]

    def translate_global(self, element, shift):
        """Translate either whole chain (glob) or a sub-element by the given shift (OpenMM Quantity)."""
        if not self.positions:
            raise ValueError('This Complex contains no positions! You CANNOT translate!')
        vec_shift = shift
        pos = self.positions[:]
        for j in range(element[0], element[2]):
            pos[j] += vec_shift
        self.positions = pos[:]

    # ---- Energy and dynamics -----------------------------------------------

    def get_energy(self):
        """Return (potential_energy_kJ_per_mol, positions)."""
        self.simulation.context.setPositions(self.positions)
        state = self.simulation.context.getState(getPositions=True, getEnergy=True, groups=1)
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        return free_E, self.positions

    def minimize(self, max_iterations: int = 100):
        """Local energy minimization (OpenMM). Returns final potential energy (kJ/mol)."""
        self.simulation.context.setPositions(self.positions)
        self.simulation.minimizeEnergy(maxIterations=max_iterations)
        state = self.simulation.context.getState(getPositions=True, getEnergy=True, groups=1)
        self.positions = state.getPositions()
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        return free_E

    def step(self, number_of_steps: int):
        """Run MD steps and return (potential_energy_kJ_per_mol, positions)."""
        self.simulation.step(number_of_steps)
        self.positions = self.simulation.context.getPositions()
        return self.get_energy()

    def rigid_minimize(self, max_iterations: int = 100, max_step_iterations: int = 100):
        """
        Experimental: random rotations on residues followed by local minimization
        to explore lower-energy conformations.
        """
        energy = None
        for i in range(max_iterations):
            for chain in self.chains:
                for idx, residue in enumerate(chain.sequence_array):
                    for j in range(max_step_iterations):
                        positions = self.positions[:]
                        chain.rotate_in_residue(idx,
                                                np.random.choice([elem for elem in range(len(chain.structure.rotating_elements[residue]))]),
                                                np.random.uniform(-np.math.pi, np.math.pi))
                        free_E = self.get_energy()[0]
                        if free_E < energy or energy is None:
                            energy = free_E
                            self.positions = positions[:]

    def pert_min(self, size: float = 1e-1, iterations: int = 50):
        """Chain-wriggling heuristic: small random coordinate kicks followed by minimization."""
        for repeat in range(iterations):
            for i in range(len(self.positions)):
                self.positions[i] += np.random.uniform(-size, size, 3) * unit.angstrom
            self.minimize()

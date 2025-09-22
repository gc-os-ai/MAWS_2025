# maws/complex.py
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import openmm as mm
from helpers import angle as ang
from helpers import nostrom
from openmm import app, unit
from tools import find_exe, run

from maws.chain import Chain  # runtime import is one-way: Complex -> Chain
from maws.prepare import makeLib
from maws.structure import *  # Structure, etc.


class Complex:
    """
    Collection of Chain objects that can be built via LEaP into AMBER prmtop/inpcrd
    and simulated/queried via OpenMM.
    """

    def __init__(
        self,
        force_field_aptamer: str = "leaprc.RNA.OL3",
        force_field_ligand: str = "leaprc.protein.ff19SB",
        conda_env: str = "maws_p3",
    ):
        self.build_string = f"""
                            source {force_field_aptamer}
                            source {force_field_ligand}
                            """
        self.prmtop = None
        self.inpcrd = None
        self.positions = None
        self.topology = None
        self.chains: list[Chain] = []
        self.system = None
        self.integrator = None
        self.simulation = None
        self.cenv = conda_env  # legacy, unused

    def add_chain(self, sequence: str, structure):
        """Append a Chain with the given sequence and Structure template."""
        if self.chains:
            start = sum(chain.length for chain in self.chains)
            chainID = len(self.chains)
        else:
            start = 0
            chainID = 0
        self.chains.append(
            Chain(self, structure, sequence=sequence, start=start, ID=chainID)
        )

    def add_chain_from_PDB(
        self,
        pdb_path: str,
        force_field_aptamer: str,
        force_field_ligand: str,
        structure=None,
        pdb_name: str = "LIG",
        parameterized: bool = False,
    ):
        """
        Create a one-residue Structure from a PDB (or pre-parameterized building block),
        generate its .lib/.frcmod via makeLib, and add as a Chain.
        """
        length = makeLib(
            pdb_path,
            pdb_name,
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            parameterized=parameterized,
        )
        path = str(Path(pdb_path).resolve().parent)
        structure = Structure([pdb_name], residue_length=[length], residue_path=path)
        self.add_chain(pdb_name, structure)

    # ---- LEaP build + cache -------------------------------------------------

    def _build_cache_key(self) -> str:
        payload = {
            "build": " ".join(self.build_string.split()),
            "inits": [
                ("".join(ch.structure.init_string.split())) for ch in self.chains
            ],
            "seqs": [ch.sequence for ch in self.chains],
        }
        return hashlib.sha1(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def build(self, target_path: str = "", file_name: str = "out"):
        if not self.chains:
            raise ValueError("Empty Complex! CANNOT build!")

        build_string_base = self.build_string
        leap_str = [self.build_string]
        for chain in self.chains:
            leap_str.append(chain.structure.init_string)
        for index, chain in enumerate(self.chains):
            if chain.sequence:
                leap_str.append(f"CHAIN{index} = sequence {{{chain.sequence}}}")
        chain_names = [
            f"CHAIN{idx}" for idx, ch in enumerate(self.chains) if ch.sequence
        ]
        chain_string = " ".join(chain_names)
        leap_str.append(f"UNION = combine {{{chain_string}}}")
        out_prefix = f"{target_path}{file_name}"
        leap_str.append(f"saveamberparm UNION {out_prefix}.prmtop {out_prefix}.inpcrd")
        leap_str.append("quit")
        leap_input = "\n".join(leap_str)

        cache_dir = Path(".maws_cache")
        cache_dir.mkdir(exist_ok=True)
        key = self._build_cache_key()
        cache_prm = cache_dir / f"{key}.prmtop"
        cache_crd = cache_dir / f"{key}.inpcrd"

        if not (cache_prm.exists() and cache_crd.exists()):
            in_file = Path(f"{target_path}{file_name}.in")
            in_file.write_text(leap_input)
            run([find_exe("tleap"), "-f", str(in_file)])
            produced_prm = Path(f"{out_prefix}.prmtop")
            produced_crd = Path(f"{out_prefix}.inpcrd")
            if not (produced_prm.exists() and produced_crd.exists()):
                raise RuntimeError(
                    "LEaP did not produce expected .prmtop/.inpcrd outputs."
                )
            produced_prm.replace(cache_prm)
            produced_crd.replace(cache_crd)

        self.build_string = build_string_base
        self.prmtop = app.AmberPrmtopFile(str(cache_prm))
        self.inpcrd = app.AmberInpcrdFile(str(cache_crd))
        self.topology = self.prmtop.topology
        self.positions = self.inpcrd.positions
        self.integrator = mm.LangevinIntegrator(
            300.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
        )
        self.system = self.prmtop.createSystem(
            nonbondedCutoff=5 * unit.angstrom,
            nonbondedMethod=app.NoCutoff,
            constraints=None,
            implicitSolvent=app.OBC1,
        )
        self.simulation = app.Simulation(self.topology, self.system, self.integrator)

    # ---- Geometry utilities -------------------------------------------------

    def rebuild(
        self, target_path: str = "", file_name: str = "out", exclusion: list = []
    ):
        old_positions = self.positions[:]
        self.build(target_path=target_path, file_name=file_name)

        for index, chain in enumerate(self.chains):
            if chain in exclusion:
                continue

            pre_positions = self.positions[chain.start : chain.start_history]
            chain_positions = old_positions[
                chain.start : chain.start + chain.length_history
            ]
            post_positions = self.positions[
                chain.start_history + chain.length_history : chain.start + chain.length
            ]

            if len(pre_positions) != 0 and chain.prepend_history:
                pre_positions = self.positions[chain.start : chain.start_history + 1]
                pre_vector = (
                    self.positions[
                        chain.start_history
                        + chain.structure.connect[chain.prepend_history[-1]][1][0]
                    ]
                    - self.positions[chain.start_history + 1]
                )
                old_pre_vector = (
                    old_positions[chain.start] - old_positions[chain.start + 1]
                )
                angle = -ang(nostrom(pre_vector), nostrom(old_pre_vector))
                axis = np.cross(
                    np.asarray(nostrom(pre_vector)), np.asarray(nostrom(old_pre_vector))
                )
                if all(axis == np.zeros(3)):
                    axis = np.array([1.0, 0.0, 0.0])
                    angle = 0
                else:
                    axis /= np.linalg.norm(axis)
                x, y, z = axis
                phi_2 = angle / 2.0
                pos = pre_positions[:]
                shift_forward = (
                    mm.Vec3(0, 0, 0) * unit.angstroms
                    - pos[-1 + chain.structure.connect[chain.prepend_history[-1]][1][0]]
                )
                s = np.math.sin(phi_2)
                c = np.math.cos(phi_2)
                rot = np.array(
                    [
                        [
                            2 * (np.power(x, 2) - 1) * np.power(s, 2) + 1,
                            2 * x * y * np.power(s, 2) - 2 * z * c * s,
                            2 * x * z * np.power(s, 2) + 2 * y * c * s,
                        ],
                        [
                            2 * x * y * np.power(s, 2) + 2 * z * c * s,
                            2 * (np.power(y, 2) - 1) * np.power(s, 2) + 1,
                            2 * z * y * np.power(s, 2) - 2 * x * c * s,
                        ],
                        [
                            2 * x * z * np.power(s, 2) - 2 * y * c * s,
                            2 * z * y * np.power(s, 2) + 2 * x * c * s,
                            2 * (np.power(z, 2) - 1) * np.power(s, 2) + 1,
                        ],
                    ]
                )
                for j in range(0, len(pos)):
                    pos[j] += shift_forward
                shift_back = chain_positions[
                    chain.structure.connect[
                        chain.sequence_array[len(chain.prepend_history)]
                    ][1][1]
                ]
                pre_bond_shift = (
                    (chain.structure.connect[chain.prepend_history[-1]][2])
                    * old_pre_vector
                    / np.linalg.norm(np.asarray(nostrom(old_pre_vector)))
                    - old_pre_vector
                )
                for j in range(0, len(pos)):
                    roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)), rot)
                    pos[j] = mm.Vec3(roted[0], roted[1], roted[2]) * unit.angstrom
                    pos[j] += shift_back + pre_bond_shift
                pre_positions = pos[:]
                chain_positions[0] += pre_bond_shift
                self.positions = (
                    self.positions[: chain.start]
                    + pre_positions[:]
                    + chain_positions[1:]
                    + self.positions[chain.start + chain.length :]
                )

            if len(post_positions) != 0 and chain.append_history:
                post_positions = self.positions[
                    chain.start_history + chain.length_history - 1 : chain.start_history
                    + chain.length
                ]
                post_vector = (
                    self.positions[chain.start_history + chain.length_history - 1]
                    - self.positions[chain.start_history + chain.length_history - 2]
                )
                old_post_vector = (
                    old_positions[chain.start_history + chain.length_history - 1]
                    - old_positions[chain.start_history + chain.length_history - 2]
                )
                angle = -ang(nostrom(post_vector), nostrom(old_post_vector))
                axis = np.cross(
                    np.asarray(nostrom(post_vector)),
                    np.asarray(nostrom(old_post_vector)),
                )
                if all(axis == np.zeros(3)):
                    axis = np.array([1.0, 0.0, 0.0])
                    angle = 0.0
                else:
                    axis /= np.linalg.norm(axis)
                x, y, z = axis
                phi_2 = angle / 2.0
                pos = post_positions[:]
                shift_forward = (
                    mm.Vec3(0, 0, 0) * unit.angstroms
                    - pos[chain.structure.connect[chain.append_history[0]][0][0]]
                )
                s = np.math.sin(phi_2)
                c = np.math.cos(phi_2)
                rot = np.array(
                    [
                        [
                            2 * (np.power(x, 2) - 1) * np.power(s, 2) + 1,
                            2 * x * y * np.power(s, 2) - 2 * z * c * s,
                            2 * x * z * np.power(s, 2) + 2 * y * c * s,
                        ],
                        [
                            2 * x * y * np.power(s, 2) + 2 * z * c * s,
                            2 * (np.power(y, 2) - 1) * np.power(s, 2) + 1,
                            2 * z * y * np.power(s, 2) - 2 * x * c * s,
                        ],
                        [
                            2 * x * z * np.power(s, 2) - 2 * y * c * s,
                            2 * z * y * np.power(s, 2) + 2 * x * c * s,
                            2 * (np.power(z, 2) - 1) * np.power(s, 2) + 1,
                        ],
                    ]
                )
                for j in range(0, len(pos)):
                    pos[j] += shift_forward
                post_bond_shift = (
                    (chain.structure.connect[chain.append_history[0]][2])
                    * old_post_vector
                    / np.linalg.norm(np.asarray(nostrom(old_post_vector)))
                    - old_post_vector
                )
                shift_back = chain_positions[
                    chain.structure.connect[
                        chain.sequence_array[-len(chain.append_history)]
                    ][0][1]
                ]
                for pos_idx, pos_elem in enumerate(pos):
                    roted = np.dot(np.array(pos_elem.value_in_unit(unit.angstrom)), rot)
                    pos[pos_idx] = mm.Vec3(roted[0], roted[1], roted[2]) * unit.angstrom
                    pos[pos_idx] += shift_back + post_bond_shift
                post_positions = pos[:]
                chain_positions[-1] += post_bond_shift
                self.positions = (
                    self.positions[: chain.start]
                    + chain_positions[:-1]
                    + post_positions[:]
                    + self.positions[chain.start + chain.length :]
                )

            if not (chain.append_history or chain.prepend_history):
                self.positions = (
                    self.positions[: chain.start]
                    + old_positions[
                        chain.start_history : chain.start_history + chain.length_history
                    ]
                    + self.positions[chain.start + chain.length :]
                )

    def rotate_element(self, element, angle: float, reverse: bool = False):
        revised_element = element[:]
        if not self.positions:
            raise ValueError("This Complex contains no positions! You CANNOT rotate!")
        pos = self.positions[:]
        vec_a = pos[revised_element[1]] - pos[revised_element[0]]
        if revised_element[2] <= revised_element[0]:
            revised_element_1 = revised_element[1]
            revised_element[1] = revised_element[2]
            revised_element[2] = revised_element_1
        self.rotate_global(revised_element, vec_a, angle, reverse=reverse, glob=False)

    def rotate_global(
        self, element, axis, angle: float, reverse: bool = False, glob: bool = True
    ):
        if not self.positions:
            raise ValueError("This Complex contains no positions! You CANNOT rotate!")
        x, y, z = np.asarray(nostrom(axis)) / (
            np.linalg.norm(np.asarray(nostrom(axis)))
        )
        phi_2 = angle / 2.0
        pos = self.positions[:]
        starting_index = 0 if glob else 1
        shift_forward = (
            mm.Vec3(0, 0, 0) * unit.angstroms
            - pos[element[2] if reverse else element[starting_index]]
        )
        s = np.math.sin(phi_2)
        c = np.math.cos(phi_2)
        rot = np.array(
            [
                [
                    2 * (np.power(x, 2) - 1) * np.power(s, 2) + 1,
                    2 * x * y * np.power(s, 2) - 2 * z * c * s,
                    2 * x * z * np.power(s, 2) + 2 * y * c * s,
                ],
                [
                    2 * x * y * np.power(s, 2) + 2 * z * c * s,
                    2 * (np.power(y, 2) - 1) * np.power(s, 2) + 1,
                    2 * z * y * np.power(s, 2) - 2 * x * c * s,
                ],
                [
                    2 * x * z * np.power(s, 2) - 2 * y * c * s,
                    2 * z * y * np.power(s, 2) + 2 * x * c * s,
                    2 * (np.power(z, 2) - 1) * np.power(s, 2) + 1,
                ],
            ]
        )
        for j in range(element[starting_index], element[2]):
            pos[j] += shift_forward
        for j in range(element[starting_index], element[2]):
            roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)), rot)
            pos[j] = mm.Vec3(roted[0], roted[1], roted[2]) * unit.angstrom
            pos[j] -= shift_forward
        self.positions = pos[:]

    def translate_global(self, element, shift):
        if not self.positions:
            raise ValueError(
                "This Complex contains no positions! You CANNOT translate!"
            )
        vec_shift = shift
        pos = self.positions[:]
        for j in range(element[0], element[2]):
            pos[j] += vec_shift
        self.positions = pos[:]

    # ---- Energy and dynamics -----------------------------------------------

    def get_energy(self):
        self.simulation.context.setPositions(self.positions)
        state = self.simulation.context.getState(
            getPositions=True, getEnergy=True, groups=1
        )
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        return free_E, self.positions

    def minimize(self, max_iterations: int = 100):
        self.simulation.context.setPositions(self.positions)
        self.simulation.minimizeEnergy(maxIterations=max_iterations)
        state = self.simulation.context.getState(
            getPositions=True, getEnergy=True, groups=1
        )
        self.positions = state.getPositions()
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        return free_E

    def step(self, number_of_steps: int):
        self.simulation.step(number_of_steps)
        self.positions = self.simulation.context.getPositions()
        return self.get_energy()

    def rigid_minimize(self, max_iterations: int = 100, max_step_iterations: int = 100):
        energy = None
        for i in range(max_iterations):
            for chain in self.chains:
                for idx, residue in enumerate(chain.sequence_array):
                    for j in range(max_step_iterations):
                        positions = self.positions[:]
                        chain.rotate_in_residue(
                            idx,
                            np.random.choice(
                                [
                                    elem
                                    for elem in range(
                                        len(chain.structure.rotating_elements[residue])
                                    )
                                ]
                            ),
                            np.random.uniform(-np.math.pi, np.math.pi),
                        )
                        free_E = self.get_energy()[0]
                        if free_E < energy or energy is None:
                            energy = free_E
                            self.positions = positions[:]

    def pert_min(self, size: float = 1e-1, iterations: int = 50):
        for repeat in range(iterations):
            for i in range(len(self.positions)):
                self.positions[i] += np.random.uniform(-size, size, 3) * unit.angstrom
            self.minimize()

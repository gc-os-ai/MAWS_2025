# maws/complex.py
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit

from maws.chain import Chain
from maws.helpers import angle as ang
from maws.helpers import nostrom
from maws.prepare import make_lib
from maws.structure import Structure
from maws.tools import find_exe, run


class Complex:
    """
    Container for one or more :class:`~maws.chain.Chain` objects that can be
    materialized via AmberTools LEaP (AMBER ``.prmtop``/``.inpcrd``) and
    simulated or queried via OpenMM.

    The :class:`Complex` owns:
      - the list of chains (topology composition),
      - the *coordinates* and OpenMM objects (``positions``, ``topology``,
        ``system``, ``integrator``, and ``simulation``),
      - the LEaP build preamble (``build_string``) and a deterministic
        caching mechanism that avoids repeated LEaP runs for identical inputs.

    Notes
    -----
    - All atoms from all chains are stored as a **single flat list** in
      :attr:`positions`. Each chain keeps a ``start`` offset and a ``length``
      so it can address its slice of atoms.
    - LEaP builds are cached in ``.maws_cache/``. The cache key is a SHA1
      hash over:
        * normalized ``build_string`` (force-field ``source`` lines),
        * each chain's LEaP ``init_string`` (``loadoff``/``loadamberparams``),
        * each chain's *canonical* sequence (after alias translation).
    - External tools must be discoverable on ``PATH``.

    Attributes
    ----------
    build_string : str
        LEaP preamble, typically a pair of ``source`` lines for aptamer and ligand FFs.
    prmtop : app.AmberPrmtopFile | None
        Loaded AMBER topology after :meth:`build`.
    inpcrd : app.AmberInpcrdFile | None
        Loaded AMBER coordinates after :meth:`build`.
    topology : app.Topology | None
        OpenMM topology object loaded from the AMBER files.
    positions : list[mm.Vec3] | None
        Coordinates for the whole complex (Å). Updated by geometry ops and MD.
    chains : list[Chain]
        Chain collection in build order (their ``start`` indices are monotonic).
    system : mm.System | None
        OpenMM System created by :meth:`build`.
    integrator : mm.Integrator | None
        OpenMM integrator (defaults to Langevin).
    simulation : app.Simulation | None
        OpenMM Simulation wrapper with context.
    """

    def __init__(
        self,
        force_field_aptamer: str = "leaprc.RNA.OL3",
        force_field_ligand: str = "leaprc.protein.ff19SB",
    ):
        """
        Parameters
        ----------
        force_field_aptamer : str, default="leaprc.RNA.OL3"
            LEaP ``source`` line for the nucleic-acid FF (e.g., RNA.OL3/DNA.OL21).
        force_field_ligand : str, default="leaprc.protein.ff19SB"
            LEaP ``source`` line for the ligand/protein/small-molecule FF.
        """
        self.build_string = f"""
                            source {force_field_aptamer}
                            source {force_field_ligand}
                            """
        self.prmtop: app.AmberPrmtopFile | None = None
        self.inpcrd: app.AmberInpcrdFile | None = None
        self.positions: list[mm.Vec3] | None = None
        self.topology: app.Topology | None = None
        self.chains: list[Chain] = []
        self.system: mm.System | None = None
        self.integrator: mm.Integrator | None = None
        self.simulation: app.Simulation | None = None

    # Chains-------------------------------------

    def add_chain(self, sequence: str, structure: Structure) -> None:
        """
        Append a :class:`Chain` with the given sequence and :class:`Structure`.

        Parameters
        ----------
        sequence : str
            **Alias** sequence (space-separated). May be ``''`` for an empty chain
            that will be populated later via :meth:`Chain.create_sequence`.
        structure : Structure
            Template bank providing residue metadata and LEaP init strings.

        Notes
        -----
        The new chain's ``start`` is set to the current total atom count
        (sum of existing chains' lengths). The chain's internal bookkeeping
        (lengths, offsets) is computed immediately if a non-empty sequence
        is provided. Topology/coordinates are only created upon :meth:`build`.
        """
        if self.chains:
            start = sum(chain.length for chain in self.chains)
            chainID = len(self.chains)
        else:
            start = 0
            chainID = 0
        self.chains.append(
            Chain(self, structure, sequence=sequence, start=start, ID=chainID)
        )

    def add_chain_from_pdb(
        self,
        pdb_path: str,
        force_field_aptamer: str,
        force_field_ligand: str,
        structure=None,
        pdb_name: str = "LIG",
        parameterized: bool = False,
    ) -> None:
        """
        Create a single-residue :class:`Structure` from a PDB (or from
        a pre-parameterized building block), then add it as a :class:`Chain`.

        This is typically used for adding the ligand as a one-residue chain.

        Parameters
        ----------
        pdb_path : str
            Path to the input PDB (or a building block). If ``parameterized=True``,
            it should already carry the coordinates you want LEaP to use.
        force_field_aptamer : str
            LEaP ``source`` line for the aptamer FF (e.g., RNA.OL3/DNA.OL21).
        force_field_ligand : str
            LEaP ``source`` line for the ligand/protein FF (e.g., ff19SB or gaff2).
        structure : None
            Unused; kept for signature compatibility.
        pdb_name : str, default="LIG"
            Name for the residue/template generated in LEaP.
        parameterized : bool, default=False
            If True, skip antechamber/parmchk2 and rely on ``loadpdb`` instead.

        Side Effects
        ------------
        Writes ``<pdb_name>.lib`` (and possibly ``.frcmod``) next to the input file.

        See Also
        --------
        maws.prepare.make_lib : Implementation details of the wrapper.
        """
        length = make_lib(
            pdb_path,
            pdb_name,
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            parameterized=parameterized,
        )
        path = str(Path(pdb_path).resolve().parent)
        structure = Structure([pdb_name], residue_length=[length], residue_path=path)
        self.add_chain(pdb_name, structure)

    # LEaP + cache

    def _build_cache_key(self) -> str:
        """
        Compute a deterministic cache key from all build-relevant inputs.

        Returns
        -------
        str
            SHA1 hex digest over:
              - normalized ``build_string`` (whitespace collapsed),
              - each chain's LEaP ``init_string`` with whitespace stripped,
              - each chain's canonical sequence string.

        Notes
        -----
        The cache is stored under ``.maws_cache/`` with filenames
        ``<key>.prmtop`` and ``<key>.inpcrd``.
        """
        payload = {
            "build": " ".join(self.build_string.split()),
            "inits": [
                ("".join(ch.structure.init_string.split())) for ch in self.chains
            ],
            "seqs": [ch.sequence for ch in self.chains],
        }
        return hashlib.sha1(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def build(self, target_path: str = "", file_name: str = "out") -> None:
        """
        Materialize current chains into AMBER ``.prmtop``/``.inpcrd`` using LEaP.
        Uses a content-addressed cache to avoid repeated LEaP runs.

        Parameters
        ----------
        target_path : str, default=""
            Path prefix where LEaP would write outputs **if** the cache is missed.
            Cache files are always under ``.maws_cache/``.
        file_name : str, default="out"
            Base name for outputs when writing (ignored during cache reuse).

        Raises
        ------
        ValueError
            If there are no chains to build.
        RuntimeError
            If LEaP does not produce the expected outputs.

        Notes
        -----
        Build steps:
          1. Construct the LEaP script:
             - prepend ``build_string``,
             - append each chain's ``structure.init_string``,
             - for each chain with a non-empty canonical ``sequence``, emit
               ``CHAIN{i} = sequence { ... }``,
             - ``UNION = combine { CHAIN0 CHAIN1 ... }``,
             - ``saveamberparm UNION <out>.prmtop <out>.inpcrd``.
          2. Compute cache key via :meth:`_build_cache_key`.
          3. If cache miss:
             - write the input script to ``<target_path>/<file_name>.in``,
             - run ``tleap -f`` on it,
             - move results into ``.maws_cache/`` under the cache key.
          4. Load cached results into OpenMM and initialize ``system``,
             ``integrator``, and ``simulation``.
        """
        if not self.chains:
            raise ValueError("Empty Complex! CANNOT build!")

        # Assemble LEaP input
        build_string_base = self.build_string  # keep exact original (whitespace)
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

        # Cache lookup
        cache_dir = Path(".maws_cache")
        cache_dir.mkdir(exist_ok=True)
        key = self._build_cache_key()
        cache_prm = cache_dir / f"{key}.prmtop"
        cache_crd = cache_dir / f"{key}.inpcrd"

        if not (cache_prm.exists() and cache_crd.exists()):
            # Write input near intended outputs for debugging
            in_file = Path(f"{target_path}{file_name}.in")
            in_file.write_text(leap_input)
            # Run tleap directly (no conda run)
            run([find_exe("tleap"), "-f", str(in_file)])

            produced_prm = Path(f"{out_prefix}.prmtop")
            produced_crd = Path(f"{out_prefix}.inpcrd")
            if not (produced_prm.exists() and produced_crd.exists()):
                raise RuntimeError(
                    "LEaP did not produce expected .prmtop/.inpcrd outputs."
                )

            produced_prm.replace(cache_prm)
            produced_crd.replace(cache_crd)

        # Load cached (or newly produced) artifacts
        self.build_string = build_string_base  # restore verbatim preamble
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

    # ------------------------------------------------------------ Geometry ops

    def rebuild(
        self,
        target_path: str = "",
        file_name: str = "out",
        exclusion: list[Chain] | None = None,
    ) -> None:
        """
        Rebuild AMBER artifacts after sequence edits, attempting to preserve
        coordinates for atoms **outside** modified regions.

        Parameters
        ----------
        target_path : str, default=""
            Passed through to :meth:`build`. Only used when the cache misses.
        file_name : str, default="out"
            Passed through to :meth:`build` when writing.
        exclusion : list, default=[]
            Optional list of chains to *exclude* from coordinate mapping. Their
            coordinates will come directly from the fresh build.

        Notes
        -----
        The algorithm:
          - Save old coordinates (``old_positions``).
          - :meth:`build` the updated topology/coordinates (cached).
          - For each chain:
            * If residues were **prepended**:
              align the leading segment using the old/new connection vectors,
              rotate into place around the new bond, fix the bond length, and
              splice back into :attr:`positions`.
            * If residues were **appended**:
              symmetric operation on the trailing segment.
            * If neither:
              splice the chain's *whole* old block back into the new coordinate
              array (preserving internal coordinates).

        This keeps most atoms unmoved and only adjusts the junctions introduced
        by prepend/append operations (driven by ``chain.prepend_history`` or
        ``chain.append_history``).
        """
        exclusion = [] if exclusion is None else list(exclusion)
        old_positions = self.positions[:]
        self.build(target_path=target_path, file_name=file_name)

        for _index, chain in enumerate(self.chains):
            if chain in exclusion:
                continue

            pre_positions = self.positions[chain.start : chain.start_history]
            chain_positions = old_positions[
                chain.start : chain.start + chain.length_history
            ]
            post_positions = self.positions[
                chain.start_history + chain.length_history : chain.start + chain.length
            ]

            # ---- handle prepended atoms ----
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
                s = np.sin(phi_2)
                c = np.cos(phi_2)
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

            # ---- handle appended atoms ----
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

            # ---- no boundary edits: reuse old internal coordinates ----
            if not (chain.append_history or chain.prepend_history):
                self.positions = (
                    self.positions[: chain.start]
                    + old_positions[
                        chain.start_history : chain.start_history + chain.length_history
                    ]
                    + self.positions[chain.start + chain.length :]
                )

    def rotate_element(self, element, angle: float, reverse: bool = False) -> None:
        """
        Rotate a contiguous **global** element of atoms by ``angle`` radians.

        Parameters
        ----------
        element : list[int]
            Triple ``[start, bond, end_exclusive]`` in **global** atom indices.
        angle : float
            Rotation angle in radians.
        reverse : bool, default=False
            Passed through to :meth:`rotate_global` to rotate the complementary
            segment when applicable.

        Raises
        ------
        ValueError
            If :attr:`positions` is not initialized (call :meth:`build` first).

        Notes
        -----
        This convenience method computes the axis as
        ``positions[bond] - positions[start]`` and forwards to
        :meth:`rotate_global` for the actual rotation.
        """
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
    ) -> None:
        """
        Core rotation kernel (Rodrigues-style matrix) for either the whole chain
        (``glob=True``) or a sub-element (``glob=False``).

        Parameters
        ----------
        element : list[int]
            ``[start, bond, end_exclusive]`` in **global** indices. If ``reverse``
            is True, the pivot is taken from ``element[2]`` (end) instead of the
            start/bond pair for the purpose of the pre/post shifts.
        axis : array-like
            Rotation axis (direction is used after normalization).
        angle : float
            Angle in radians.
        reverse : bool, default=False
            If True, rotate the **complement** of the selected range relative to
            the pivot.
        glob : bool, default=True
            If True, rotate from ``element[0]`` to ``element[2]``; if False, start
            from ``element[1]`` (used by :meth:`rotate_element`).

        Raises
        ------
        ValueError
            If :attr:`positions` is not initialized.
        """
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
        s = np.sin(phi_2)
        c = np.cos(phi_2)
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

    def translate_global(self, element, shift) -> None:
        """
        Translate a global element by a displacement vector.

        Parameters
        ----------
        element : list[int]
            ``[start, bond, end_exclusive]`` in **global** indices. The ``bond``
            entry is ignored for translation; only the range ``[start:end)`` is used.
        shift : array-like or openmm.unit.Quantity
            Displacement vector. Quantity units should be compatible with Å.

        Raises
        ------
        ValueError
            If :attr:`positions` is not initialized.
        """
        if not self.positions:
            raise ValueError(
                "This Complex contains no positions! You CANNOT translate!"
            )
        vec_shift = shift
        pos = self.positions[:]
        for j in range(element[0], element[2]):
            pos[j] += vec_shift
        self.positions = pos[:]

    # ---------------------------------------------------------- Energetics/MD

    def get_energy(self) -> tuple[float, list[mm.Vec3]]:
        """
        Compute potential energy for the current coordinates.

        Returns
        -------
        (float, list[mm.Vec3])
            Tuple ``(potential_energy_kJ_per_mol, positions)``.
        """
        self.simulation.context.setPositions(self.positions)
        state = self.simulation.context.getState(
            getPositions=True, getEnergy=True, groups=1
        )
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        return free_E, self.positions

    def minimize(self, max_iterations: int = 100) -> float:
        """
        Local energy minimization.

        Parameters
        ----------
        max_iterations : int, default=100
            Maximum number of minimization iterations.

        Returns
        -------
        float
            Final potential energy (kJ/mol).
        """
        self.simulation.context.setPositions(self.positions)
        self.simulation.minimizeEnergy(maxIterations=max_iterations)
        state = self.simulation.context.getState(
            getPositions=True, getEnergy=True, groups=1
        )
        self.positions = state.getPositions()
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        return free_E

    def step(self, number_of_steps: int) -> tuple[float, list[mm.Vec3]]:
        """
        Advance MD by a number of steps.

        Parameters
        ----------
        number_of_steps : int
            Number of integrator steps to execute.

        Returns
        -------
        (float, list[mm.Vec3])
            Same as :meth:`get_energy`.
        """
        self.simulation.step(number_of_steps)
        self.positions = self.simulation.context.getPositions()
        return self.get_energy()

    def rigid_minimize(
        self, max_iterations: int = 100, max_step_iterations: int = 100
    ) -> None:
        """
        Experimental search: random residue torsions followed by local minimization.

        Parameters
        ----------
        max_iterations : int, default=100
            Outer loop iterations over residues.
        max_step_iterations : int, default=100
            Random torsion proposals per residue before accepting an improvement.

        Notes
        -----
        For each residue, propose random torsions; if a proposal improves the
        energy, keep it, otherwise revert the coordinates. This is a greedy
        heuristic and not guaranteed to converge to a global minimum.
        """
        energy = None
        for _i in range(max_iterations):
            for chain in self.chains:
                for idx, residue in enumerate(chain.sequence_array):
                    for _j in range(max_step_iterations):
                        positions = self.positions[:]
                        rots = chain.structure.rotating_elements[residue]
                        if rots == [None]:  # skip residues with no torsions
                            continue
                        n_rots = len(rots)

                        chain.rotate_in_residue(
                            idx,
                            int(
                                np.random.randint(n_rots)
                            ),  # or: np.random.choice(n_rots)
                            np.random.uniform(-np.pi, np.pi),
                        )

                        free_E = self.get_energy()[0]
                        if free_E < energy or energy is None:
                            energy = free_E
                            self.positions = positions[:]

    def pert_min(self, size: float = 1e-1, iterations: int = 50) -> None:
        """
        Chain-wriggling heuristic: apply small random kicks, then minimize.

        Parameters
        ----------
        size : float, default=1e-1
            Uniform kick magnitude in Å for each coordinate component.
        iterations : int, default=50
            Number of (kick → minimize) cycles to perform.
        """
        for _repeat in range(iterations):
            for i in range(len(self.positions)):
                self.positions[i] += np.random.uniform(-size, size, 3) * unit.angstrom
            self.minimize()

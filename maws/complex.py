"""
maws.complex
============

High-level container that composes one or more Chain objects and materializes
them into AMBER topology/coordinates using AmberTools LEaP, then into OpenMM
objects for simulation.

Builds are content-addressed and cached under ``.maws_cache/``.
"""

from __future__ import annotations

import hashlib
import json
import os
import warnings
from collections import deque
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit
from openmm.unit import Quantity

from maws.chain import Chain
from maws.helpers import angle as ang
from maws.helpers import nostrom
from maws.prepare import make_lib
from maws.structure import Structure
from maws.tools import find_exe, run


class Complex:
    """
    Container for Chain objects that builds AMBER topology and runs OpenMM simulations.

    Attributes
    ----------
    chains : list[Chain]
        Chains in insertion order.
    positions : list[openmm.Vec3] | None
        Flat coordinate array over all chains.
    topology : openmm.app.Topology | None
        OpenMM topology.
    simulation : openmm.app.Simulation | None
        OpenMM Simulation object.
    """

    def __init__(
        self,
        force_field_aptamer: str = "leaprc.RNA.OL3",
        force_field_ligand: str = "leaprc.protein.ff19SB",
        salt_conc: float = 0.15,
    ):
        """
        Parameters
        ----------
        force_field_aptamer : str, default="leaprc.RNA.OL3"
            LEaP ``source`` line for the nucleic-acid FF (e.g., RNA.OL3/DNA.OL21).
        force_field_ligand : str, default="leaprc.protein.ff19SB"
            LEaP ``source`` line for the ligand/protein/small-molecule FF.
        salt_conc : float, default=0.15
            Monovalent salt concentration (mol/L) for the GB implicit-solvent
            Debye-Hückel screening term. The default (~physiological, 150 mM)
            screens the highly charged backbone; ``0.0`` reproduces the old
            unscreened behavior. Note: GB screening is monovalent only and does
            not model divalent ions such as Mg²⁺.
        """
        self.build_string = f"""
                            source {force_field_aptamer}
                            source {force_field_ligand}
                            """
        self.salt_conc = salt_conc
        self.prmtop: app.AmberPrmtopFile | None = None
        self.inpcrd: app.AmberInpcrdFile | None = None
        self.positions: list[mm.Vec3] | None = None
        self.topology: app.Topology | None = None
        self.chains: list[Chain] = []
        self.system: mm.System | None = None
        self.integrator: mm.Integrator | None = None
        self.simulation: app.Simulation | None = None
        self._bond_adjacency: dict[int, set[int]] | None = None
        self._moving_set_cache: dict[tuple[int, int], list[int]] = {}

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

    # ---- Chain access convenience methods ----

    def get_chain(self, index: int) -> Chain:
        """
        Return chain by index with bounds checking.

        Parameters
        ----------
        index : int
            Chain index (0-based). Negative indices work (Python-style).

        Returns
        -------
        Chain
            The chain at the given index.

        Raises
        ------
        IndexError
            If index is out of bounds.
        """
        if not self.chains:
            raise IndexError("Complex has no chains")
        return self.chains[index]

    def aptamer_chain(self) -> Chain:
        """
        Convenience: return chain[0], typically the aptamer.

        Returns
        -------
        Chain
            The first chain, conventionally the aptamer in MAWS workflows.

        Raises
        ------
        IndexError
            If no chains exist.

        Notes
        -----
        This is syntactic sugar for ``complex.chains[0]``. In MAWS, the
        aptamer is always added first via :meth:`add_chain`.
        """
        return self.get_chain(0)

    def ligand_chain(self) -> Chain:
        """
        Convenience: return chain[1], typically the ligand.

        Returns
        -------
        Chain
            The second chain, conventionally the ligand in MAWS workflows.

        Raises
        ------
        IndexError
            If fewer than two chains exist.

        Notes
        -----
        This is syntactic sugar for ``complex.chains[1]``. In MAWS, the
        ligand is typically added second via :meth:`add_chain_from_pdb`.
        """
        return self.get_chain(1)

    # ------------------ Platform selection (CPU/GPU) ------------------------

    def _create_simulation(self) -> app.Simulation:
        """
        Create OpenMM Simulation, preferring GPU if available.

        Order of precedence:
        1. MAWS_OPENMM_PLATFORM env var (if set)
        2. CUDA (if available)
        3. OpenCL (if available)
        4. OpenMM default

        Returns
        -------
        app.Simulation
            Configured simulation object.
        """

        # Load plugins so GPU platforms are visible
        mm.Platform.loadPluginsFromDirectory(mm.Platform.getDefaultPluginsDirectory())

        # Check for user override via env var
        env_platform = os.getenv("MAWS_OPENMM_PLATFORM")
        if env_platform:
            try:
                platform = mm.Platform.getPlatformByName(env_platform)
                return app.Simulation(
                    self.topology, self.system, self.integrator, platform
                )
            except Exception:
                pass  # Fall through to auto-detection

        # Try GPU platforms in order
        for name in ("CUDA", "OpenCL"):
            try:
                platform = mm.Platform.getPlatformByName(name)
                return app.Simulation(
                    self.topology, self.system, self.integrator, platform
                )
            except Exception as e:
                warnings.warn(
                    f"Platform '{name}' unavailable: {e}", RuntimeWarning, stacklevel=2
                )
                continue

        # Let OpenMM choose (usually CPU)
        return app.Simulation(self.topology, self.system, self.integrator)

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
              - each chain's canonical sequence string,
              - the *contents* of each chain's ``.lib``/``.frcmod`` files.

        Notes
        -----
        The cache is stored under ``.maws_cache/`` with filenames
        ``<key>.prmtop`` and ``<key>.inpcrd``.

        The library file contents are included because ``make_lib`` hardcodes the
        ligand residue name ``LIG``: every protein writes ``LIG.lib`` to the same
        path and yields the same canonical sequence ``"LIG"``. Hashing only the
        ``init_string`` (a ``loadoff <path>`` line) would therefore collide across
        *different* proteins and silently serve the first one's cached topology.
        """
        payload = {
            "build": " ".join(self.build_string.split()),
            "inits": [
                ("".join(ch.structure.init_string.split())) for ch in self.chains
            ],
            "seqs": [ch.sequence for ch in self.chains],
            "libs": [self._chain_lib_digests(ch) for ch in self.chains],
        }
        return hashlib.sha1(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    @staticmethod
    def _chain_lib_digests(chain: Chain) -> list[list[str]]:
        """
        Return content digests of the ``.lib``/``.frcmod`` files a chain loads.

        Parameters
        ----------
        chain : Chain
            Chain whose :class:`Structure` references LEaP resource files.

        Returns
        -------
        list[list[str]]
            One ``[name, lib_digest, frcmod_digest]`` row per residue template.
            A missing file contributes the sentinel ``""``. Returns an empty list
            when the structure has no ``residue_path`` (no files to load).

        Notes
        -----
        Used by :meth:`_build_cache_key` so that two ligands sharing the same
        ``LIG.lib`` path but differing in contents produce different cache keys.
        """
        structure = chain.structure
        residue_path = getattr(structure, "residue_path", None)
        if residue_path is None:
            return []

        base = Path(residue_path) if residue_path != "" else Path(".")

        def digest(path: Path) -> str:
            if not path.exists():
                return ""
            return hashlib.sha1(path.read_bytes()).hexdigest()

        return [
            [name, digest(base / f"{name}.lib"), digest(base / f"{name}.frcmod")]
            for name in structure.residue_names
        ]

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
        # A different topology means different connectivity, so anything
        # derived from the previous one is stale.
        self._bond_adjacency = None
        self._moving_set_cache = {}
        self.integrator = mm.LangevinIntegrator(
            300.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
        )
        self.system = self._make_system()
        self.simulation = self._create_simulation()

    def _make_system(self) -> mm.System:
        """
        Create the OpenMM ``System`` from the loaded prmtop.

        Isolated from :meth:`build` so the implicit-solvent settings (notably
        the salt-screening term) can be exercised without running LEaP.
        """
        return self.prmtop.createSystem(
            nonbondedCutoff=5 * unit.angstrom,
            nonbondedMethod=app.NoCutoff,
            constraints=None,
            implicitSolvent=app.OBC1,
            implicitSolventSaltConc=self.salt_conc * unit.molar,
        )

    # ---------------------------------  Geometry ops  ------------------------------

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

    def _adjacency(self) -> dict[int, set[int]]:
        """Return each atom index mapped to the indices bonded to it.

        Read from the topology that :meth:`build` loads, and kept until the
        next build replaces it.

        Returns
        -------
        dict[int, set[int]]
            One entry per atom in the Complex.

        Raises
        ------
        ValueError
            If the Complex has no topology, because :meth:`build` has not run.
        """
        if self._bond_adjacency is None:
            if self.topology is None:
                raise ValueError(
                    "This Complex has no topology; call build() before rotating."
                )
            adjacency: dict[int, set[int]] = {
                atom.index: set() for atom in self.topology.atoms()
            }
            for first, second in self.topology.bonds():
                adjacency[first.index].add(second.index)
                adjacency[second.index].add(first.index)
            self._bond_adjacency = adjacency
        return self._bond_adjacency

    def moving_set(self, fixed: int, pivot: int) -> list[int]:
        """moving_set(fixed, pivot) -> list[int]

        Return the atoms that move when one bond is turned.

        Deleting a single bond splits a molecule into two connected parts.
        Turning that bond moves one part and leaves the other where it is.
        This returns the part containing `pivot`, found by following bonds
        outward from it and never crossing the `fixed`-`pivot` bond.

        The answer comes from which atoms are bonded, so it does not depend
        on the order the atoms happen to be stored in.

        Parameters
        ----------
        fixed : int
            Index of the bonded atom on the side that stays put.
        pivot : int
            Index of the bonded atom on the side that moves. It is included
            in the result, but it lies on the rotation axis and so keeps its
            position.

        Returns
        -------
        list[int]
            Atom indices, ascending.

        Raises
        ------
        ValueError
            If the Complex has no topology, or if the two atoms named are not
            bonded to each other.

        See Also
        --------
        rotate_element : Turns a bond, using this to decide what moves.

        Notes
        -----
        Answers are cached. Connectivity is fixed between builds, and this is
        consulted once per conformation sampled.

        Examples
        --------
        The 14 atoms of a guanine base move when its C1'-N9 bond is turned.

        >>> complex.moving_set(c1_prime_index, n9_index)  # doctest: +SKIP
        [44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57]
        """
        cached = self._moving_set_cache.get((fixed, pivot))
        if cached is not None:
            return cached

        adjacency = self._adjacency()
        if pivot not in adjacency.get(fixed, ()):
            raise ValueError(
                f"Atoms {fixed} and {pivot} are not bonded, so there is no "
                f"bond between them to turn."
            )

        reached = {pivot}
        queue = deque([pivot])
        while queue:
            for neighbour in adjacency[queue.popleft()]:
                if neighbour != fixed and neighbour not in reached:
                    reached.add(neighbour)
                    queue.append(neighbour)

        result = sorted(reached)
        self._moving_set_cache[(fixed, pivot)] = result
        return result

    def rotate_element(
        self, element: Sequence[int], angle: float, reverse: bool = False
    ) -> None:
        """Turn one bond by `angle` radians.

        A torsion is a turn of part of a molecule about one of its own bonds.
        It changes the shape while leaving every bond length and bond angle
        as it was.

        Parameters
        ----------
        element : sequence of int
            ``[first, second]`` atom indices naming the bond to turn, counted
            across the whole Complex. A third entry is accepted and ignored;
            which atoms move is read from the topology.
        angle : float
            How far to turn, in radians.
        reverse : bool, default=False
            Turn the part of the molecule joined to `first` rather than the
            part joined to `second`. Both named atoms sit on the rotation
            axis, so either choice leaves all bond lengths unchanged.

        Raises
        ------
        ValueError
            If :attr:`positions` is unset, if the Complex has no topology, or
            if the two atoms named are not bonded to each other.

        See Also
        --------
        moving_set : Decides which atoms this moves.
        rotate_global : Turns a whole chain about an axis of your choosing.

        Examples
        --------
        >>> before = complex.positions[:]  # doctest: +SKIP
        >>> complex.rotate_element([c1_prime, n9], 1.0)  # doctest: +SKIP
        >>> complex.positions == before  # doctest: +SKIP
        False
        """
        if not self.positions:
            raise ValueError("This Complex contains no positions! You CANNOT rotate!")

        first, second = element[0], element[1]
        fixed, pivot = (second, first) if reverse else (first, second)
        positions = self.positions
        self._rotate_indices(
            self.moving_set(fixed, pivot),
            positions[second] - positions[first],
            angle,
            pivot,
        )

    def _rotate_indices(
        self,
        indices: Sequence[int],
        axis: Quantity,
        angle: float,
        centre: int,
    ) -> None:
        r"""Turn the named atoms about an axis through one of them.

        Parameters
        ----------
        indices : sequence of int
            Atom indices to move, counted across the whole Complex.
        axis : openmm.unit.Quantity
            Direction of the rotation axis, a 3-vector in ångström. Only the
            direction is used; the length is discarded.
        angle : float
            How far to turn, in radians.
        centre : int
            Index of the atom the axis passes through. It stays where it is,
            and so does anything else lying on the axis.

        Notes
        -----
        Rotation matrix built from the half-angle :math:`\phi/2`, with
        :math:`s = \sin(\phi/2)`, :math:`c = \cos(\phi/2)` and the axis
        written as a unit vector :math:`(x, y, z)`:

        .. math::
            R_{00} = 2 (x^2 - 1) s^2 + 1
            \qquad
            R_{01} = 2 x y s^2 - 2 z c s

        with the other entries following by cyclic permutation. Positions are
        multiplied as :math:`v R`, which turns by :math:`-\phi`. Nothing in
        this package depends on that sign, since angles are drawn from a range
        symmetric about zero.

        :attr:`positions` is replaced.
        """
        x, y, z = np.asarray(nostrom(axis)) / np.linalg.norm(np.asarray(nostrom(axis)))
        s = np.sin(angle / 2.0)
        c = np.cos(angle / 2.0)
        rot = np.array(
            [
                [
                    2 * (x**2 - 1) * s**2 + 1,
                    2 * x * y * s**2 - 2 * z * c * s,
                    2 * x * z * s**2 + 2 * y * c * s,
                ],
                [
                    2 * x * y * s**2 + 2 * z * c * s,
                    2 * (y**2 - 1) * s**2 + 1,
                    2 * z * y * s**2 - 2 * x * c * s,
                ],
                [
                    2 * x * z * s**2 - 2 * y * c * s,
                    2 * z * y * s**2 + 2 * x * c * s,
                    2 * (z**2 - 1) * s**2 + 1,
                ],
            ]
        )

        pos = self.positions[:]
        to_origin = mm.Vec3(0, 0, 0) * unit.angstroms - pos[centre]
        for j in indices:
            shifted = (pos[j] + to_origin).value_in_unit(unit.angstrom)
            turned = np.dot(np.array(shifted), rot)
            pos[j] = mm.Vec3(*turned) * unit.angstrom - to_origin
        self.positions = pos[:]

    def rotate_global(
        self,
        element: Sequence[int],
        axis: Quantity,
        angle: float,
        reverse: bool = False,
        glob: bool = True,
    ) -> None:
        """Turn a contiguous run of atoms about an axis of your choosing.

        Unlike :meth:`rotate_element`, the axis here is given rather than
        taken from a bond, so this can reorient a whole chain in space.

        Parameters
        ----------
        element : sequence of int
            ``[start, bond, end]`` atom indices, counted across the whole
            Complex. ``end`` is exclusive.
        axis : openmm.unit.Quantity
            Direction of the rotation axis, a 3-vector in ångström. Attach
            units with :func:`maws.helpers.angstrom` if you have a plain
            array; a unitless one is rejected.
        angle : float
            How far to turn, in radians.
        reverse : bool, default=False
            Put the axis through ``end`` instead of through the first atom of
            the range.
        glob : bool, default=True
            Turn ``[start, end)``. Set ``False`` to turn ``[bond, end)``.

        Raises
        ------
        ValueError
            If :attr:`positions` is not initialized.

        See Also
        --------
        rotate_element : Turns a single bond, working out what moves.

        Examples
        --------
        Point a chain along a different axis, turning about its first atom.

        >>> from maws.helpers import angstrom
        >>> complex.rotate_global(  # doctest: +SKIP
        ...     chain.element, angstrom([0.0, 0.0, 1.0]), 3.14159
        ... )
        """
        if not self.positions:
            raise ValueError("This Complex contains no positions! You CANNOT rotate!")

        first = element[0] if glob else element[1]
        centre = element[2] if reverse else first
        self._rotate_indices(range(first, element[2]), axis, angle, centre)

    def translate_global(self, element: Sequence[int], shift: Quantity) -> None:
        """
        Translate a global element by a displacement vector.

        Parameters
        ----------
        element : Sequence[int]
            ``[start, bond, end_exclusive]`` in **global** indices.
            The ``bond`` entry is ignored for translation; only ``[start:end)`` is used.
        shift : openmm.unit.Quantity
            Displacement vector **with length units** (Å). For example:
            ``np.array([5.0, 0.0, 0.0]) * unit.angstrom``.
            Unitless arrays are **not** accepted.

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

    # ------------------ Energetics/MD --------------------------------------

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
        # ensure the Context starts from our current coordinates
        self.simulation.context.setPositions(self.positions)

        # integrate
        self.simulation.step(number_of_steps)

        # read back positions and energy from the Context
        state = self.simulation.context.getState(
            getPositions=True, getEnergy=True, groups=1
        )
        self.positions = state.getPositions()
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        return free_E, self.positions

    def _freeze_particles(self, indices: Sequence[int]) -> dict[int, Quantity]:
        """
        Temporarily set the mass of ``indices`` to zero so they cannot move.

        OpenMM's :class:`~openmm.LocalEnergyMinimizer` (and the integrators)
        hold zero-mass particles fixed, so this is how a rigid region is
        expressed. Zero mass affects dynamics only - the potential energy of a
        configuration is unchanged.

        Parameters
        ----------
        indices : Sequence[int]
            Global atom indices to immobilize.

        Returns
        -------
        dict[int, openmm.unit.Quantity]
            The original masses, for :meth:`_restore_particle_masses`.

        Notes
        -----
        Either every index is frozen and the caller receives the masses needed
        to undo it, or the System is left exactly as it was found. There is no
        outcome in which particles stay massless with no record of their
        original masses.
        """
        if not indices or self.system is None:
            return {}

        saved = {i: self.system.getParticleMass(i) for i in indices}
        try:
            for i in indices:
                self.system.setParticleMass(i, 0.0)
            if self.simulation is not None:
                self.simulation.context.reinitialize(preserveState=True)
        except BaseException:
            # The saved masses are lost with the propagating exception, so the
            # caller could not undo a partial freeze even if it wanted to.
            self._restore_particle_masses(saved)
            raise
        return saved

    def _restore_particle_masses(self, saved: dict[int, Quantity]) -> None:
        """
        Undo :meth:`_freeze_particles`, restoring the original masses.

        Parameters
        ----------
        saved : dict[int, openmm.unit.Quantity]
            Mapping returned by :meth:`_freeze_particles`.
        """
        if not saved or self.system is None:
            return

        for i, mass in saved.items():
            self.system.setParticleMass(i, mass)
        if self.simulation is not None:
            self.simulation.context.reinitialize(preserveState=True)

    def pert_min(
        self,
        size: float = 1e-1,
        iterations: int = 50,
        atoms: Sequence[int] | None = None,
    ) -> None:
        """
        Chain-wriggling heuristic: apply small random kicks, then minimize.

        Parameters
        ----------
        size : float, default=1e-1
            Uniform kick magnitude in Å for each coordinate component.
        iterations : int, default=50
            Number of (kick → minimize) cycles to perform.
        atoms : Sequence[int] or None, default=None
            Global indices of the atoms allowed to move. Atoms outside this set
            are neither kicked nor relaxed: their masses are set to zero for the
            duration of the call, which pins them during minimization. ``None``
            lets the **whole** complex move.

        Raises
        ------
        IndexError
            If any entry of ``atoms`` falls outside ``[0, len(self.positions))``.

        Notes
        -----
        In a MAWS run the docking target is rigid and only the aptamer chain may
        move, so callers should pass that chain's atom range. Leaving ``atoms``
        unset perturbs the target as well, which both randomizes the input
        structure and lets the distortion accumulate across growth steps via the
        carried-forward coordinates.
        """
        n_atoms = len(self.positions)
        if atoms is None:
            mobile: Sequence[int] = range(n_atoms)
            immobile: Sequence[int] = []
        else:
            mobile = sorted(set(atoms))
            if mobile and not 0 <= mobile[0] <= mobile[-1] < n_atoms:
                # A negative index would otherwise kick an atom at the far end
                # of the array while that same atom stays in ``immobile`` and
                # is frozen - moved and pinned at once, with no error raised.
                raise IndexError(
                    f"atoms must be global indices in [0, {n_atoms}), got "
                    f"[{mobile[0]}, {mobile[-1]}]"
                )
            immobile = sorted(set(range(n_atoms)).difference(mobile))

        saved_masses = self._freeze_particles(immobile)
        try:
            for _repeat in range(iterations):
                for i in mobile:
                    self.positions[i] += (
                        np.random.uniform(-size, size, 3) * unit.angstrom
                    )
                self.minimize()
        finally:
            self._restore_particle_masses(saved_masses)

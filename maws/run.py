"""
maws.run
========

Python API for the MAWS aptamer design algorithm.

This module provides a programmatic interface to run MAWS without using the CLI.
Use :class:`MAWSConfig` to configure the run and :func:`run_maws` to execute it.

Examples
--------
>>> from maws.run import MAWSConfig, run_maws
>>> config = MAWSConfig(
...     name="my_aptamer",
...     pdb_path="data/ligand.pdb",
...     num_nucleotides=15,
...     aptamer_type="RNA",
...     molecule_type="protein",
... )
>>> result = run_maws(config)  # doctest: +SKIP
>>> print(result.sequence)  # doctest: +SKIP
"""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass
from datetime import datetime
from typing import Literal

import numpy as np
from openmm import app

import maws.space as space
from maws.complex import Complex
from maws.dna_structure import load_dna_structure
from maws.helpers import center_of_mass, nostrom
from maws.pdb_cleaner import clean_one_file
from maws.rna_structure import load_rna_structure
from maws.routines import S

AptamerType = Literal["RNA", "DNA"]
MoleculeType = Literal["protein", "organic", "lipid"]


@dataclass
class MAWSConfig:
    """
    Configuration for a MAWS aptamer design run.

    Parameters
    ----------
    pdb_path : str
        Path to the ligand PDB file.
    num_nucleotides : int, default=15
        Number of nucleotides in the aptamer.
    name : str, default="MAWS_aptamer"
        Job name (used for output files).
    aptamer_type : {"RNA", "DNA"}, default="RNA"
        Type of aptamer to design.
    molecule_type : {"protein", "organic", "lipid"}, default="protein"
        Type of ligand molecule.
    beta : float, default=0.01
        Inverse temperature for Boltzmann sampling.
    first_chunk_size : int, default=5000
        Number of samples in the first MAWS step.
    second_chunk_size : int, default=5000
        Number of samples in subsequent steps.
    clean_pdb : bool, default=False
        Whether to clean the input PDB before processing.
    keep_chains : str, default="all"
        Chain policy for cleaner: "all", "one", or comma-separated list.
    remove_h : bool, default=False
        Whether to remove hydrogens during cleaning.
    drop_hetatm : bool, default=False
        Whether to drop HETATM records during cleaning.
    verbose : bool, default=True
        Whether to log progress to console.

    Examples
    --------
    >>> config = MAWSConfig(
    ...     pdb_path="data/ligand.pdb",
    ...     num_nucleotides=10,
    ...     aptamer_type="DNA",
    ... )
    >>> config.aptamer_type
    'DNA'
    """

    pdb_path: str
    num_nucleotides: int
    name: str = "MAWS_aptamer"
    aptamer_type: AptamerType = "RNA"
    molecule_type: MoleculeType = "protein"
    beta: float = 0.01
    first_chunk_size: int = 5000
    second_chunk_size: int = 5000
    clean_pdb: bool = False
    keep_chains: str = "all"
    remove_h: bool = False
    drop_hetatm: bool = False
    verbose: bool = True


@dataclass
class MAWSResult:
    """
    Result of a MAWS aptamer design run.

    Attributes
    ----------
    sequence : str
        Best aptamer sequence found.
    energy : float
        Final energy of the best configuration.
    entropy : float
        Final entropy score.
    complex : Complex
        OpenMM Complex object with final structure.
    config : MAWSConfig
        Configuration used for the run.
    pdb_path : str
        Path to the result PDB file.

    Examples
    --------
    >>> result.sequence
    'GAUCC...'
    >>> result.save_pdb("output.pdb")
    """

    sequence: str
    energy: float
    entropy: float
    complex: Complex
    config: MAWSConfig
    pdb_path: str = ""

    def save_pdb(self, path: str) -> str:
        """
        Save the result structure to a PDB file.

        Parameters
        ----------
        path : str
            Output file path.

        Returns
        -------
        str
            Path to the saved file.
        """
        with open(path, "w") as f:
            app.PDBFile.writeModel(
                self.complex.topology, self.complex.positions, file=f
            )
        return path


def _setup_logger(name: str, verbose: bool) -> logging.Logger:
    """Set up logger for the run."""
    logger = logging.getLogger(f"maws.{name}")
    logger.setLevel(logging.INFO)
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler
    file_handler = logging.FileHandler(f"{name}_output.log", mode="w")
    file_handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(file_handler)

    # Console handler (if verbose)
    if verbose:
        console = logging.StreamHandler()
        console.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(console)

    return logger


def _resolve_pdb_path(
    pdb_path: str,
    molecule_type: str,
    *,
    clean_pdb: bool,
    keep_chains: str,
    remove_h: bool,
    drop_hetatm: bool,
    logger: logging.Logger | None = None,
) -> tuple[str, str]:
    """Decide the PDB path to use (possibly cleaned)."""
    original = pdb_path
    if not clean_pdb:
        return pdb_path, original

    if molecule_type != "protein":
        if logger:
            logger.warning(
                "--clean-pdb is intended for protein inputs; skipping for non-protein."
            )
        return pdb_path, original

    try:
        new_path = clean_one_file(
            pdb_path,
            keep=keep_chains,
            remove_h=remove_h,
            drop_hetatm=drop_hetatm,
        )
        if logger:
            logger.info("Cleaned PDB written to: %s", new_path)
        return new_path, original
    except Exception as e:
        if logger:
            logger.warning("PDB cleaning failed; using original. Reason: %s", e)
        return pdb_path, original


def run_maws(config: MAWSConfig) -> MAWSResult:
    """
    Run the MAWS aptamer design algorithm.

    Parameters
    ----------
    config : MAWSConfig
        Configuration object specifying run parameters.

    Returns
    -------
    MAWSResult
        Object containing the best aptamer sequence and structure.

    Examples
    --------
    >>> from maws.run import MAWSConfig, run_maws
    >>> config = MAWSConfig(
    ...     pdb_path="data/ligand.pdb",
    ...     num_nucleotides=5,
    ...     first_chunk_size=10,
    ...     second_chunk_size=10,
    ... )
    >>> result = run_maws(config)  # doctest: +SKIP
    >>> len(result.sequence) == 5 + 1  # doctest: +SKIP
    True
    """
    logger = _setup_logger(config.name, config.verbose)

    # Version info
    VERSION = "1.0"
    RELEASE_DATE = "2026"
    N_ELEMENTS = 4  # rotatable backbone torsions per residue

    logger.info("MAWS - Making Aptamers With Software")
    logger.info("Active version: %s (released: %s)", VERSION, RELEASE_DATE)
    logger.info("Type of aptamer: %s", config.aptamer_type)
    logger.info("Type of ligand molecule: %s", config.molecule_type)
    logger.info("Job: %s", config.name)

    # Resolve PDB path
    pdb_path, original_pdb = _resolve_pdb_path(
        config.pdb_path,
        config.molecule_type,
        clean_pdb=config.clean_pdb,
        keep_chains=config.keep_chains,
        remove_h=config.remove_h,
        drop_hetatm=config.drop_hetatm,
        logger=logger,
    )
    logger.info("Input file: %s", original_pdb)

    # Choose force fields
    if config.aptamer_type == "RNA":
        molecule = load_rna_structure()
        nt_list = "GAUC"
        force_field_aptamer = "leaprc.RNA.OL3"
    else:
        molecule = load_dna_structure()
        nt_list = "GATC"
        force_field_aptamer = "leaprc.DNA.OL21"

    if config.molecule_type == "protein":
        force_field_ligand = "leaprc.protein.ff19SB"
        parameterized = True
    elif config.molecule_type == "organic":
        force_field_ligand = "leaprc.gaff2"
        parameterized = False
    else:
        force_field_ligand = "leaprc.lipid21"
        parameterized = False

    logger.info("Force field (aptamer): %s", force_field_aptamer)
    logger.info("Force field (ligand): %s", force_field_ligand)
    logger.info("Start time: %s", datetime.now())

    # Build template complex
    cpx = Complex(
        force_field_aptamer=force_field_aptamer,
        force_field_ligand=force_field_ligand,
    )
    cpx.add_chain("", molecule)  # empty aptamer chain
    cpx.add_chain_from_pdb(
        pdb_path=pdb_path,
        force_field_aptamer=force_field_aptamer,
        force_field_ligand=force_field_ligand,
        parameterized=parameterized,
    )

    # Build ligand-only complex for COM calculation
    c = Complex(
        force_field_aptamer=force_field_aptamer,
        force_field_ligand=force_field_ligand,
    )
    c.add_chain_from_pdb(
        pdb_path=pdb_path,
        force_field_aptamer=force_field_aptamer,
        force_field_ligand=force_field_ligand,
        parameterized=parameterized,
    )
    c.build()

    # Sampling spaces
    cube = space.Cube(20.0, center_of_mass(np.asarray(nostrom(c.positions))))
    rotations = space.NAngles(N_ELEMENTS)

    # Tracking
    best_entropy = None
    best_sequence = None
    best_positions = None

    logger.info("Initialized; starting step 1.")

    # ---- Step 1: choose first nucleotide ----
    for ntide in nt_list:
        logger.info("Step 1: sampling for '%s'", ntide)
        energies = []
        free_E = None
        position = None

        cx = copy.deepcopy(cpx)
        aptamer = cx.aptamer_chain()
        aptamer.create_sequence(ntide)
        cx.build()

        for _ in range(config.first_chunk_size):
            cx.positions = cx.positions[:]
            sample = cube.generator()
            cx.translate_global(aptamer.element, sample[:3])
            cx.rotate_global(aptamer.element, sample[3:6], sample[6])

            for i in range(N_ELEMENTS):
                rotation_sample = rotations.generator()
                aptamer.rotate_in_residue(0, i, rotation_sample[i])

            cx.minimize()
            E = cx.get_energy()
            energies.append(E)
            if free_E is None or E < free_E:
                free_E = E
                position = cx.positions[:]

            cx.rebuild()

        entropy = float(S(energies, config.beta))
        if best_entropy is None or entropy > best_entropy:
            best_entropy = entropy
            best_sequence = ntide
            best_positions = position[:]

    logger.info("Step 1 completed. Selected: %s", best_sequence)

    # ---- Subsequent steps ----
    for step_num in range(config.num_nucleotides):
        logger.info("Step %d: extending sequence", step_num + 2)
        best_old_sequence = best_sequence
        best_old_positions = best_positions[:]

        for ntide in nt_list:
            energies = []
            free_E = None
            position = None

            cx = copy.deepcopy(cpx)
            aptamer = cx.aptamer_chain()
            aptamer.create_sequence(best_old_sequence)
            cx.build()
            cx.positions = best_old_positions[:]

            aptamer.append_sequence(ntide)
            cx.rebuild()

            for _ in range(config.second_chunk_size):
                cx.positions = cx.positions[:]
                for i in range(N_ELEMENTS):
                    rotation_sample = rotations.generator()
                    aptamer.rotate_in_historic_residue(
                        len(aptamer.sequence_array) - 1, i, rotation_sample[i]
                    )

                cx.minimize()
                E = cx.get_energy()
                energies.append(E)
                if free_E is None or E < free_E:
                    free_E = E
                    position = cx.positions[:]

                cx.rebuild()

            entropy = float(S(energies, config.beta))
            if entropy > best_entropy:
                best_entropy = entropy
                best_sequence = " ".join(aptamer.alias_sequence_array)
                best_positions = position[:]

        logger.info("Step %d completed. Sequence: %s", step_num + 2, best_sequence)

    # Final result
    result_complex = copy.deepcopy(cpx)
    aptamer = result_complex.aptamer_chain()
    aptamer.create_sequence(best_sequence)
    result_complex.build()
    result_complex.positions = best_positions[:]

    result_pdb = f"{config.name}_RESULT.pdb"
    with open(result_pdb, "w") as f:
        app.PDBFile.writeModel(
            result_complex.topology, result_complex.positions, file=f
        )

    logger.info("Completed. Final sequence: %s", best_sequence)
    logger.info("Result saved to: %s", result_pdb)

    return MAWSResult(
        sequence=best_sequence,
        energy=float(S([cx.get_energy()], config.beta)),
        entropy=best_entropy,
        complex=result_complex,
        config=config,
        pdb_path=result_pdb,
    )

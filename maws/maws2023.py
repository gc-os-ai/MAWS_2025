#!/usr/bin/env python3
# MAWS is part of the sharksome software suite
# This software is published under MIT license
# COPYRIGHT 2017 Michael Jendrusch
# Authors: Michael Jendrusch, Stefan Holderbach
#
# Modifications by DTU Biobuilders (2023) and additional refactor (2025)

import argparse
import copy
import logging
from datetime import datetime

import numpy as np
from openmm import app, unit

import maws.space as space
from maws.complex import Complex
from maws.dna_structure import load_dna_structure
from maws.helpers import nostrom
from maws.kernels import center_of_mass
from maws.pdb_cleaner import clean_one_file
from maws.rna_structure import load_rna_structure
from maws.routines import S

# VERSION = "2.1" # Original Authoras To-do: cite in readme
# RELEASE_DATE = "2017" # Original Authors To-do: cite in readme
VERSION = "2.2"  # Siddharth
RELEASE_DATE = "2025"  # Siddharth
METHOD = "Kullback-Leibler"


def _resolve_pdb_path(
    pdb_path: str,
    molecule_type: str,
    *,
    clean_pdb: bool,
    keep_chains: str,
    remove_h: bool,
    drop_hetatm: bool,
    logger=None,
) -> tuple[str, str]:
    """
    Decide the PDB path to use (possibly cleaned). Returns (final_path, original_path).

    - Cleans only if `clean_pdb=True` AND `molecule_type == "protein"`.
    - Logs to `logger(msg: str)` if provided (same file handle as your run log).
    """
    original = pdb_path
    if not clean_pdb:
        return pdb_path, original

    if molecule_type != "protein":
        if logger:
            logger(
                "[WARN] --clean-pdb is intended for protein inputs;"
                "skipping for non-protein types.\n"
            )
        return pdb_path, original

    try:
        new_path = clean_one_file(
            pdb_path,
            keep=keep_chains,  # "all" | "one" | "A,B"
            remove_h=remove_h,  # True/False
            drop_hetatm=drop_hetatm,  # True/False
        )
        if logger:
            logger(f"[CLEAN] Cleaned PDB written to: {new_path}\n")
        return new_path, original
    except Exception as e:
        if logger:
            logger(f"[WARN] PDB cleaning failed; using original file. Reason: {e}\n")
        return pdb_path, original


def parse_args():
    parser = argparse.ArgumentParser(description="MAWS - Making Aptamers With Software")
    parser.add_argument(
        "-n", "--name", type=str, default="MAWS_aptamer", help="Job name."
    )
    parser.add_argument(
        "-nt",
        "--ntides",
        type=int,
        default=15,
        help="Number of nucleotides in the aptamer.",
    )
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        default="./data/pfoa.pdb",
        help="Path to your ligand PDB file.",
    )
    parser.add_argument(
        "-ta",
        "--aptamertype",
        type=str,
        default="RNA",
        choices=["RNA", "DNA"],
        help="Type of aptamer (DNA or RNA).",
    )
    parser.add_argument(
        "-tm",
        "--moleculetype",
        type=str,
        default="protein",
        choices=["protein", "organic", "lipid"],
        help="Type of ligand molecule.",
    )
    parser.add_argument(
        "-b", "--beta", type=float, default=0.01, help="Inverse temperature."
    )
    parser.add_argument(
        "-c1",
        "--firstchunksize",
        type=int,
        default=5000,
        help="Number of samples in the first MAWS step.",
    )
    parser.add_argument(
        "-c2",
        "--secondchunksize",
        type=int,
        default=5000,
        help="Number of samples in all subsequent MAWS steps.",
    )
    parser.add_argument(
        "--clean-pdb",
        action="store_true",
        help="Clean the input PDB before LEaP (use for protein-type inputs).",
    )
    parser.add_argument(
        "--keep-chains",
        type=str,
        default="all",
        help='Chain policy for cleaner: "all", "one", or comma list like "A,B".',
    )
    parser.add_argument(
        "--remove-h", action="store_true", help="Cleaner: remove hydrogens."
    )
    parser.add_argument(
        "--drop-hetatm",
        action="store_true",
        help="Cleaner: drop all HETATM records (NOT recommended for small molecules).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Params
    JOB_NAME = args.name
    BETA = args.beta
    FIRST_CHUNK_SIZE = args.firstchunksize
    SECOND_CHUNK_SIZE = args.secondchunksize
    N_NTIDES = args.ntides
    PDB_PATH = args.path
    APTAMER_TYPE = args.aptamertype
    MOLECULE_TYPE = args.moleculetype
    N_ELEMENTS = 4  # rotatable backbone torsions per residue

    # ---------------- Logging configuration ----------------
    logger = logging.getLogger("maws")
    logger.setLevel(logging.INFO)

    # Avoid duplicate handlers if main() is called multiple times
    if logger.hasHandlers():
        logger.handlers.clear()

    file_handler = logging.FileHandler(f"{JOB_NAME}_output.log", mode="w")
    file_handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(file_handler)

    # Logs: only entropy + step cache are opened as files now
    with (
        open(f"{JOB_NAME}_entropy.log", "w") as entropy_log,
        open(f"{JOB_NAME}_step_cache.pdb", "w") as step,
    ):
        # Header
        output.write("MAWS - Making Aptamers With Software\n")
        output.write(f"Active version: {VERSION} (released:_{RELEASE_DATE})\n")
        output.write(f"Computational method: {METHOD}\n")
        output.write(f"Type of aptamer: {APTAMER_TYPE}\n")
        output.write(f"Type of ligand molecule: {MOLECULE_TYPE}\n")
        output.write(f"Job: {JOB_NAME}\n")
        # Resolve (and optionally clean) the PDB before any LEaP calls
        PDB_PATH, ORIGINAL_PDB_PATH = _resolve_pdb_path(
            PDB_PATH,
            MOLECULE_TYPE,
            clean_pdb=args.clean_pdb,
            keep_chains=args.keep_chains,
            remove_h=args.remove_h,
            drop_hetatm=args.drop_hetatm,
            logger=output.write,
        )
        output.write(f"Input file (original): {ORIGINAL_PDB_PATH}\n")
        if PDB_PATH != ORIGINAL_PDB_PATH:
            output.write(f"Input file (cleaned):  {PDB_PATH}\n")
        output.write(f"Sample number in initial step: {FIRST_CHUNK_SIZE}\n")
        output.write(f"Sample number per further steps: {SECOND_CHUNK_SIZE}\n")
        output.write(
            f"Number of further steps: {N_NTIDES} (sequence length = {N_NTIDES + 1})\n"
        )
        logger.info("Input file (original): %s", ORIGINAL_PDB_PATH)
        if PDB_PATH != ORIGINAL_PDB_PATH:
            logger.info("Input file (cleaned): %s", PDB_PATH)
        logger.info("Sample number in initial step: %d", FIRST_CHUNK_SIZE)
        logger.info("Sample number per further steps: %d", SECOND_CHUNK_SIZE)
        logger.info(
            "Number of further steps: %d (sequence length = %d)",
            N_NTIDES,
            N_NTIDES + 1,
        )
        logger.info("Value of beta: %s", BETA)
        logger.info("Start time: %s", datetime.now())

        # Choose aptamer FF and residue
        if APTAMER_TYPE == "RNA":
            molecule = load_rna_structure()
            nt_list = "GAUC"
            force_field_aptamer = "leaprc.RNA.OL3"
        else:  # DNA
            molecule = load_dna_structure()
            nt_list = "GATC"
            force_field_aptamer = "leaprc.DNA.OL21"
        logger.info("Force field selected for the aptamer: %s", force_field_aptamer)

        # Choose ligand FF
        if MOLECULE_TYPE == "protein":
            force_field_ligand = "leaprc.protein.ff19SB"
            parameterized = (
                True  # proteins are pre-parameterized in ff19SB → skip antechamber
            )
        elif MOLECULE_TYPE == "organic":
            force_field_ligand = "leaprc.gaff2"
            parameterized = False  # small molecules need antechamber/parmchk2
        else:
            force_field_ligand = "leaprc.lipid21"
            parameterized = False  # standard lipids in lipid21 → skip antechamber
        output.write(
            f"Force field selected for the ligand molecule: {force_field_ligand}\n"
        )

        # Complex template with an empty aptamer chain + the ligand from PDB
        cpx = Complex(
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
        )
        cpx.add_chain("", molecule)  # empty aptamer chain (sequence added later)
        cpx.add_chain_from_pdb(
            pdb_path=PDB_PATH,
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            parameterized=parameterized,
        )

        # Build a separate complex with just the ligand to compute COM for sampling
        c = Complex(
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
        )
        c.add_chain_from_pdb(
            pdb_path=PDB_PATH,
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            parameterized=parameterized,
        )
        c.build()

        # Sampling spaces: cube of width 20 Å around ligand COM
        cube = space.Cube(
            20.0,
            center_of_mass(np.asarray(nostrom(c.positions))),
        )
        rotations = space.NAngles(N_ELEMENTS)

        # Tracking best candidate
        best_entropy = None
        best_sequence = None
        best_positions = None
        # best_ntide = None
        best_topology = None

        logger.info("Initialized successfully; starting step 1.")

        # ---- Step 1: choose the first nucleotide ----------------------------
        for ntide in nt_list:
            logger.info("Step 1: starting initial sampling for '%s'", ntide)
            energies = []
            free_E = None
            position = None

            # Clone the template complex
            cx = copy.deepcopy(cpx)
            aptamer = cx.aptamer_chain()

            # Seed sequence and build (LEaP build will hit cache after first time)
            aptamer.create_sequence(ntide)
            cx.build()

            # Remember initial positions
            positions0 = cx.positions[:]

            # Sample orientations/rotations
            for _ in range(FIRST_CHUNK_SIZE):
                orientation = cube.generator()
                rotation = rotations.generator()

                aptamer.translate_global(orientation[0:3] * unit.angstrom)
                aptamer.rotate_global(
                    orientation[3:-1] * unit.angstrom, orientation[-1]
                )

                for j in range(N_ELEMENTS):
                    aptamer.rotate_in_residue(0, j, rotation[j])

                energy = cx.get_energy()[0]
                if free_E is None or energy < free_E:
                    free_E = energy
                    position = cx.positions[:]
                energies.append(energy)

                # Reset for next sample
                cx.positions = positions0[:]

            entropy = S(energies, beta=BETA)

            # # Outputs
            # with open(f"{JOB_NAME}_1_{ntide}.pdb", "w") as pdblog:
            #     app.PDBFile.writeModel(
            #         copy.deepcopy(cx.topology),
            #         position[:],
            #         file=pdblog,
            #         modelIndex=1,
            #     )

            logger.info(
                "Step 1: finished candidate '%s' (entropy=%s, best_E=%s)",
                aptamer.alias_sequence,
                entropy,
                free_E,
            )

            entropy_log.write(
                f"SEQUENCE: {aptamer.alias_sequence} "
                f"ENTROPY: {entropy} ENERGY: {free_E}\n"
            )

            if best_entropy is None or entropy < best_entropy:
                best_entropy = entropy
                best_sequence = ntide
                best_positions = position[:]
                best_topology = copy.deepcopy(cx.topology)

        # Cache best-of-step to file
        app.PDBFile.writeModel(best_topology, best_positions, file=step, modelIndex=1)
        logger.info("Completed first step. Selected nucleotide: %s", best_sequence)
        logger.info("Starting further steps to append %d nucleotides", N_NTIDES)

        # ---- Steps 2..N: grow sequence --------------------------------------
        for i in range(1, N_NTIDES):
            best_old_sequence = best_sequence
            best_old_positions = best_positions[:]
            best_entropy = None

            logger.info(
                "Step %d: starting with current best sequence %s",
                i + 1,
                best_old_sequence,
            )

            for ntide in nt_list:
                for append in [True, False]:
                    energies = []
                    free_E = None
                    position = None

                    cx = copy.deepcopy(cpx)
                    aptamer = cx.aptamer_chain()
                    aptamer.create_sequence(best_old_sequence)

                    cx.build()  # cached
                    cx.positions = best_old_positions[:]

                    if append:
                        aptamer.append_sequence(ntide)
                    else:
                        aptamer.prepend_sequence(ntide)

                    cx.rebuild()  # cached build + coordinate mapping
                    cx.pert_min(size=0.5)  # light shake to find nearby minima

                    positions0 = cx.positions[:]

                    for _ in range(SECOND_CHUNK_SIZE):
                        rotation = rotations.generator()

                        # Forward rotations on the new residue’s internal bonds
                        for j in range(N_ELEMENTS - 1):
                            if append:
                                aptamer.rotate_in_residue(-1, j, rotation[j])
                            else:
                                aptamer.rotate_in_residue(
                                    0, j, rotation[j], reverse=True
                                )

                        # Backward rotation (C3'-O3')
                        if append:
                            aptamer.rotate_in_residue(-2, 3, rotation[3])
                        else:
                            aptamer.rotate_in_residue(0, 3, rotation[3], reverse=True)

                        energy = cx.get_energy()[0]
                        if free_E is None or energy < free_E:
                            free_E = energy
                            position = cx.positions[:]
                        energies.append(energy)

                        cx.positions = positions0[:]

                    entropy = S(energies, beta=BETA)

                    logger.info(
                        "Step %d: finished candidate sequence %s "
                        "(%s '%s') (entropy=%s, best_E=%s)",
                        i + 1,
                        aptamer.alias_sequence,
                        "append" if append else "prepend",
                        ntide,
                        entropy,
                        free_E,
                    )

                    entropy_log.write(
                        f"SEQUENCE: {aptamer.alias_sequence} "
                        f"ENTROPY: {entropy} ENERGY: {free_E}\n"
                    )

                    if best_entropy is None or entropy < best_entropy:
                        best_entropy = entropy
                        best_positions = position[:]
                        best_sequence = aptamer.alias_sequence
                        best_topology = copy.deepcopy(cx.topology)

                        logger.info(
                            "Step %d: new best sequence %s "
                            "(added '%s' on %s end; entropy=%s, best_E=%s)",
                            i + 1,
                            best_sequence,
                            ntide,
                            "3'" if append else "5'",
                            best_entropy,
                            free_E,
                        )

            app.PDBFile.writeModel(
                best_topology, best_positions, file=step, modelIndex=1
            )
            logger.info(
                "Completed step %d. Selected sequence: %s",
                i + 1,
                best_sequence,
            )

        # ---- Final render ----------------------------------------------------
        result_complex = copy.deepcopy(cpx)
        aptamer = result_complex.aptamer_chain()
        aptamer.create_sequence(best_sequence)
        result_complex.build()  # cached
        result_complex.positions = best_positions[:]

        result_pdb_name = f"{JOB_NAME}_RESULT.pdb"
        with open(result_pdb_name, "w") as pdb_result:
            app.PDBFile.writeModel(
                result_complex.topology, result_complex.positions, file=pdb_result
            )

        logger.info(
            "Run completed. Final sequence: %s (result written to %s)",
            best_sequence,
            result_pdb_name,
        )


if __name__ == "__main__":
    main()

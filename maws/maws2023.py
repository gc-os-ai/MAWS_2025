#!/usr/bin/env python3
# MAWS is part of the sharksome software suite
# This software is published under MIT license
# COPYRIGHT 2017 Michael Jendrusch
# Authors: Michael Jendrusch, Stefan Holderbach
#
# Modifications by DTU Biobuilders (2023) and additional refactor (2025)

import argparse
import copy
from datetime import datetime

import numpy as np
from openmm import app, unit

import maws.space as space
from maws.complex import Complex
from maws.dna_structure import load_dna_structure
from maws.helpers import nostrom
from maws.kernels import center_of_mass
from maws.rna_structure import load_rna_structure
from maws.routines import S

VERSION = "2.1"
RELEASE_DATE = "2017"
METHOD = "Kullback-Leibler"


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

    # Logs
    with (
        open(f"{JOB_NAME}_output.log", "w") as output,
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
        output.write(f"Input file: {PDB_PATH}\n")
        output.write(f"Sample number in initial step: {FIRST_CHUNK_SIZE}\n")
        output.write(f"Sample number per further steps: {SECOND_CHUNK_SIZE}\n")
        output.write(
            f"Number of further steps: {N_NTIDES} (sequence length = {N_NTIDES + 1})\n"
        )
        output.write(f"Value of beta: {BETA}\n")
        output.write(f"Start time: {datetime.now()}\n")

        # Choose aptamer FF and residue
        if APTAMER_TYPE == "RNA":
            molecule = load_rna_structure()
            nt_list = "GAUC"
            force_field_aptamer = "leaprc.RNA.OL3"
        else:  # DNA
            molecule = load_dna_structure()
            nt_list = "GATC"
            force_field_aptamer = "leaprc.DNA.OL21"
        output.write(f"Force field selected for the aptamer: {force_field_aptamer}\n")

        # Choose ligand FF
        if MOLECULE_TYPE == "protein":
            force_field_ligand = "leaprc.protein.ff19SB"
        elif MOLECULE_TYPE == "organic":
            force_field_ligand = "leaprc.gaff2"
        else:
            force_field_ligand = "leaprc.lipid21"
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
            parameterized=False,
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
            parameterized=False,
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
        best_ntide = None
        best_topology = None

        output.write("Initialized successfully!\n")

        # ---- Step 1: choose the first nucleotide ----------------------------
        for ntide in nt_list:
            output.write(f"{datetime.now()}: starting initial step for '{ntide}'\n")
            energies = []
            free_E = None
            position = None

            # Clone the template complex
            cx = copy.deepcopy(cpx)
            aptamer = cx.chains[0]

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

            # Outputs
            with open(f"{JOB_NAME}_1_{ntide}.pdb", "w") as pdblog:
                app.PDBFile.writeModel(
                    copy.deepcopy(cx.topology),
                    position[:],
                    file=pdblog,
                    modelIndex=1,
                )

            entropy_log.write(
                f"SEQUENCE: {aptamer.alias_sequence} "
                f"ENTROPY: {entropy} ENERGY: {free_E}\n"
            )

            if best_entropy is None or entropy < best_entropy:
                best_entropy = entropy
                best_sequence = ntide
                best_ntide = ntide
                best_positions = position[:]
                best_topology = copy.deepcopy(cx.topology)

        # Cache best-of-step to file
        app.PDBFile.writeModel(best_topology, best_positions, file=step, modelIndex=1)
        with open(f"{JOB_NAME}_best_1_{best_ntide}.pdb", "w") as pdblog:
            app.PDBFile.writeModel(
                best_topology, best_positions, file=pdblog, modelIndex=1
            )

        output.write(
            f"{datetime.now()}: Completed first step. "
            f"Selected nucleotide: {best_sequence}\n"
        )
        output.write(
            f"{datetime.now()}: Starting further steps to append "
            f"{N_NTIDES} nucleotides\n"
        )

        # ---- Steps 2..N: grow sequence --------------------------------------
        for i in range(1, N_NTIDES):
            best_old_sequence = best_sequence
            best_old_positions = best_positions[:]
            best_entropy = None

            for ntide in nt_list:
                for append in [True, False]:
                    energies = []
                    free_E = None
                    position = None

                    cx = copy.deepcopy(cpx)
                    aptamer = cx.chains[0]
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

                    with open(f"{JOB_NAME}_{i + 1}_{ntide}.pdb", "w") as pdblog:
                        app.PDBFile.writeModel(
                            copy.deepcopy(cx.topology),
                            position[:],
                            file=pdblog,
                            modelIndex=1,
                        )

                    entropy_log.write(
                        f"SEQUENCE: {aptamer.alias_sequence} "
                        f"ENTROPY: {entropy} ENERGY: {free_E}\n"
                    )

                    if best_entropy is None or entropy < best_entropy:
                        best_entropy = entropy
                        best_positions = position[:]
                        best_ntide = ntide
                        best_sequence = aptamer.alias_sequence
                        best_topology = copy.deepcopy(cx.topology)

            app.PDBFile.writeModel(
                best_topology, best_positions, file=step, modelIndex=1
            )
            output.write(
                f"{datetime.now()}: Completed step {i + 1}. "
                f"Selected sequence: {best_sequence}\n"
            )
            with open(f"{JOB_NAME}_best_{i + 1}_{best_ntide}.pdb", "w") as pdblog:
                app.PDBFile.writeModel(
                    best_topology, best_positions, file=pdblog, modelIndex=1
                )

        # ---- Final render ----------------------------------------------------
        result_complex = copy.deepcopy(cpx)
        aptamer = result_complex.chains[0]
        aptamer.create_sequence(best_sequence)
        result_complex.build()  # cached
        result_complex.positions = best_positions[:]

        with open(f"{JOB_NAME}_RESULT.pdb", "w") as pdb_result:
            app.PDBFile.writeModel(
                result_complex.topology, result_complex.positions, file=pdb_result
            )

        output.write(f"{datetime.now()}: Run completed. Thank you for using MAWS!\n\n")
        output.write(f"Final sequence: {best_sequence}\n")


if __name__ == "__main__":
    main()

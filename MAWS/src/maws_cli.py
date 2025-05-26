#!/usr/bin/env python3
# MAWS is part of the sharksome software suite
# This software is published under the MIT licence
# Copyright © 2017 Michael Jendrusch & Stefan Holderbach

# Modern refactor 2025 – siddharth
#
# VERSION information

from __future__ import annotations


VERSION       = "2.1"          # bump whenever code changes!
RELEASE_DATE  = "2025-05-26"
METHOD        = "Kullback-Leibler"

# --------------------------------------------------------------------------- #
# Imports (modernised)
# --------------------------------------------------------------------------- #


import argparse
import copy
import logging
from datetime import datetime
from pathlib import Path
from typing import List

import numpy as np
import openmm as mm
from openmm import app, unit

from MAWS.src.LoadFrom import xml_structure as XMLStructure
from MAWS.src.helpers   import nostrom,from_angstrom
from MAWS.src.Complex   import Complex
from MAWS.src.Routines  import S                                 # entropy
from MAWS.src.Kernels   import centerOfMass
from MAWS.src.Spaces    import Cube, NAngles

# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="MAWS – Making Aptamers With Software")
    p.add_argument("-n", "--name", default="MAWS_job", help="job name prefix")
    p.add_argument("-b", "--beta", type=float, default=0.01, help="inverse temperature β")
    p.add_argument("-c1", "--firstchunksize", type=int, default=5000, help="# samples in first step")
    p.add_argument("-c2", "--secondchunksize", type=int, default=5000, help="# samples in later steps")
    p.add_argument("-t", "--ntides", type=int, default=60, help="# further nucleotides to add")
    p.add_argument("-p", "--pdb", required=True, help="ligand PDB file")
    p.add_argument("-a", "--aptamer-type", choices=("DNA", "RNA"), default="DNA", help="residue library")
    p.add_argument("--log", default="INFO", help="logging level (DEBUG, INFO, …)")
    return p.parse_args()


# --------------------------------------------------------------------------- #
# Helper to emit progress into log & console
# --------------------------------------------------------------------------- #
def _setup_logging(level: str, job: str) -> logging.Logger:
    logging.basicConfig(
        level=level.upper(),
        format="%(asctime)s  %(levelname)-8s %(message)s",
        handlers=[
            logging.FileHandler(f"{job}_output.log", mode="w"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger("maws")


# --------------------------------------------------------------------------- #
# MAIN
# --------------------------------------------------------------------------- #
def main() -> None:
    args = _parse_cli()
    log  = _setup_logging(args.log, args.name)

    log.info("MAWS %s (%s) – %s", VERSION, RELEASE_DATE, METHOD)
    log.info("Job '%s' started with PDB '%s'", args.name, args.pdb)

    # ------------------------------------------------------------------ #
    # Residue library & nucleotide list
    # ------------------------------------------------------------------ #
    # xml_file = Path(__file__).with_suffix(".xml").with_name(f"{args.aptamer_type}.xml")
    xml_file = (Path(__file__).resolve().parent.parent
                / "xml_files"
                / f"{args.aptamer_type}.xml")
    struct   = XMLStructure(xml_file)
    # nucleotides = ["G", "A", "U", "C"] if args.aptamer_type == "RNA" else ["DGN", "DAN", "DTN", "DCN"]
    nucleotides = ["G", "A", "U", "C"] if args.aptamer_type == "RNA" else ["DG", "DA", "DT", "DC"]
    # choose both leaprcs for DNA (or RNA)
    ff = (
        "leaprc.protein.ff14SB\nsource leaprc.DNA.bsc1"
        if args.aptamer_type == "DNA"
        else "leaprc.protein.ff14SB\nsource leaprc.RNA.OL3"
    )
    # ------------------------------------------------------------------ #
    # Prototype Complex (aptamer chain + ligand)
    # ------------------------------------------------------------------ #
    # proto = Complex("leaprc.protein.ff14SB")
    proto = Complex(ff)
    proto.add_chain("", struct)                                    # empty aptamer placeholder
    proto.add_chain_from_pdb(args.pdb, protein_name=args.name, parameterized=True)
    print("Structure library contains units: %s", struct.residue_names)


    # Build helper complex 'c' to compute CoM for sampling cube
    # c = Complex("leaprc.protein.ff14SB")
    c = Complex(ff)
    c.add_chain_from_pdb(args.pdb, protein_name=args.name, parameterized=True)
    c.build(args.name)
    # cube_center = center_of_mass(np.asarray(nostrom(c.positions)))
    positions_np = np.array([from_angstrom(p) for p in c.positions])
    cube_center  = centerOfMass(positions_np)
    cube        = Cube(20.0, cube_center)
    rotations   = NAngles(4)

    # ------------------------------------------------------------------ #
    # FILE handles for entropy & step cache (retain original behaviour)
    # ------------------------------------------------------------------ #
    entropy_log = open(f"{args.name}_entropy.log", "w")
    step_cache  = open(f"{args.name}_step_cache.pdb", "w")

    # ------------------------------------------------------------------ #
    # FIRST STEP – choose best mono-nucleotide
    # ------------------------------------------------------------------ #
    best_entropy: float | None               = None
    best_sequence: str | None                = None
    best_positions: List[unit.Quantity] | None = None
    best_topology = None

    for nt in nucleotides:
        log.info("Initial sampling for nucleotide %s", nt)
        energies = []
        free_E   = None
        position = None

        cplx   = copy.deepcopy(proto)
        chain  = cplx.chains[0]
        chain.create_sequence(nt)
        print("Building chain with sequence tokens: %s", chain.sequence_tokens)
        cplx.build(args.name)
        pos0 = cplx.positions[:]

        for k in range(args.firstchunksize):
            samp  = cube.generator()
            theta = rotations.generator()

            chain.translate_global(unit.Quantity(samp[:3], unit.angstrom))
            chain.rotate_global(unit.Quantity(samp[3:6], unit.angstrom), samp[6])
            for j in range(4):
                chain.rotate_in_residue(0, j, theta[j])

            e = cplx.get_energy()[0]
            energies.append(e)
            if free_E is None or e < free_E:
                free_E, position = e, cplx.positions[:]
                log.debug("FIRST %s chunk %d: energy %.2f kJ", nt, k + 1, e)

            cplx.positions = pos0[:]

        H = S(energies, beta=args.beta)
        entropy_log.write(f"SEQUENCE {nt}  ENTROPY {H:.6f}  ENERGY {free_E:.2f}\n")

        if best_entropy is None or H < best_entropy:
            best_entropy, best_sequence, best_positions = H, nt, position[:]
            best_topology = copy.deepcopy(cplx.topology)

    assert best_sequence and best_positions
    log.info("Best mono-nucleotide: %s (entropy %.4f)", best_sequence, best_entropy)

    # ------------------------------------------------------------------ #
    # FURTHER STEPS – grow chain to N_NTIDES
    # ------------------------------------------------------------------ #
    for step_idx in range(args.ntides):
        log.info("Iteration %d / %d – current best seq '%s'", step_idx + 2, args.ntides + 1, best_sequence)
        best_old_seq = best_sequence
        best_old_pos = best_positions[:]
        best_entropy = None

        for nt in nucleotides:
            for append in (True, False):
                cplx = copy.deepcopy(proto)
                chain = cplx.chains[0]
                chain.create_sequence(best_old_seq)
                cplx.build(args.name)
                cplx.positions = best_old_pos[:]

                if append:
                    chain.append_sequence(nt)
                else:
                    chain.prepend_sequence(nt)

                cplx.rebuild(args.name)
                try:
                    cplx.pert_min(size=0.5)
                except Exception:
                    pass

                pos0 = cplx.positions[:]
                energies = []
                free_E, position = None, None

                for k in range(args.secondchunksize):
                    theta = rotations.generator()
                    # rotate new residue’s own rotors
                    for j in range(3):
                        if append:
                            chain.rotate_in_residue(-1, j, theta[j])
                        else:
                            chain.rotate_in_residue(0, j, theta[j], reverse=True)
                    # rotate backbone C3’O3’ / previous residue
                    if append:
                        chain.rotate_in_residue(-2, 3, theta[3])
                    else:
                        chain.rotate_in_residue(0, 3, theta[3], reverse=True)

                    e = cplx.get_energy()[0]
                    energies.append(e)
                    if free_E is None or e < free_E:
                        free_E, position = e, cplx.positions[:]
                        log.debug("STEP %d  %s %s chunk %d: E=%.2f",
                                  step_idx + 2, "APPEND" if append else "PREPEND", nt, k + 1, e)
                    cplx.positions = pos0[:]

                H = S(energies, beta=args.beta)
                entropy_log.write(f"SEQUENCE {chain.alias_sequence}  ENTROPY {H:.6f}  ENERGY {free_E:.2f}\n")

                if best_entropy is None or H < best_entropy:
                    best_entropy = H
                    best_sequence = chain.alias_sequence
                    best_positions = position[:]
                    best_topology = copy.deepcopy(cplx.topology)

        log.info("Completed step %d; best sequence now '%s'", step_idx + 2, best_sequence)

    # ------------------------------------------------------------------ #
    # WRITE final model
    # ------------------------------------------------------------------ #
    result = copy.deepcopy(proto)
    result.chains[0].create_sequence(best_sequence)
    result.build(file_name=args.name)
    result.positions = best_positions[:]

    with open(f"{args.name}_RESULT.pdb", "w") as fh:
        app.PDBFile.writeModel(result.topology, result.positions, file=fh)
    log.info("Run completed – final sequence: %s", best_sequence)

    # tidy
    entropy_log.close()
    step_cache.close()


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    main()

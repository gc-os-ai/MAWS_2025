from __future__ import annotations

import copy
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from openmm import app, unit

import maws.space as space
from maws.complex import Complex
from maws.dna_structure import load_dna_structure
from maws.pdb_cleaner import resolve_pdb_path
from maws.rna_structure import load_rna_structure
from maws.routines import S

AptamerType = Literal["RNA", "DNA"]
MoleculeType = Literal["protein", "organic", "lipid"]
PDBInput = str | Path


@dataclass(frozen=True)
class MawsResult:
    """
    Result of a MAWS aptamer design run.

    Attributes
    ----------
    sequence : str
        Best aptamer sequence found.
    energy : float
        Energy of the final best configuration.
    entropy : float
        Entropy score used for selection.
    pdb_path : str
        Path to the saved result PDB file produced by the run.
    """

    sequence: str
    energy: float
    entropy: float
    pdb_path: str | None = None


class MawsRunner:
    def __init__(
        self,
        *,
        num_nucleotides: int,
        aptamer_type: AptamerType,
        molecule_type: MoleculeType,
        beta: float = 0.01,
        first_chunck_size: int = 5000,
        second_chunck_size: int = 5000,
        clean_pdb: bool = False,
        keep_chains: str = "all",
        remove_h: bool = False,
        drop_hetatm: bool = False,
        verbose: bool = False,
        shape: Literal["cube", "sphere", "shell"] = "shell",
        reach: float = 10.0,
        probe: float = 1.4,
    ) -> None:
        if num_nucleotides <= 0:
            raise ValueError("num_nucleotides couldn't be less than 0")
        if first_chunck_size <= 0 or second_chunck_size <= 0:
            raise ValueError("Chunck size must be greater than 0")
        if reach < 0:
            raise ValueError(f"reach must be >= 0, got {reach}")
        if probe < 0:
            raise ValueError(f"probe must be >= 0, got {probe}")

        self.num_nucleotides = num_nucleotides
        self.aptamer_type = aptamer_type
        self.molecule_type = molecule_type
        self.beta = beta
        self.first_chunk_size = first_chunck_size
        self.second_chunk_size = second_chunck_size
        self.clean_pdb = clean_pdb
        self.keep_chains = keep_chains
        self.remove_h = remove_h
        self.drop_hetatm = drop_hetatm
        self.verbose = verbose
        self.shape = shape
        self.reach = reach
        self.probe = probe

    def run(
        self,
        *,
        pdb: PDBInput,
        name: str = "MAWS_aptamer",
        output_pdb: str | Path | None = None,
    ) -> MawsResult:
        """
        Run the MAWS algorithm.

        Parameters
        ----------
        pdb : str | Path
            Input ligand PDB file path.
        name : str
            Run name used only for log context and artifact naming.
        output_pdb : str | Path | None
            If provided:
              - if it's an existing directory -> writes `{name}_RESULT.pdb` inside it
              - otherwise treated as the exact output file path (parent dirs created)
            If None, no PDB is written.
        logger : logging.Logger | None
            Logger to use. If None, uses module logger.

        Returns
        -------
        MAWSResult
        """
        N_BACKBONE_TORSIONS = (
            4  # MAWS rotates 4 backbone torsions per residue in this implementation
        )
        log = logging.getLogger(__name__)

        if self.verbose:
            log.info("MAWS run started: name=%s", name)
        log.debug(
            "Config: num_nucleotides=%d aptamer_type=%s molecule_type=%s beta=%s "
            "c1=%d c2=%d clean_pdb=%s keep_chains=%s remove_h=%s drop_hetatm=%s",
            self.num_nucleotides,
            self.aptamer_type,
            self.molecule_type,
            self.beta,
            self.first_chunk_size,
            self.second_chunk_size,
            self.clean_pdb,
            self.keep_chains,
            self.remove_h,
            self.drop_hetatm,
        )

        # Resolve (and optionally clean) the PDB path before LEaP calls
        pdb_path, original_pdb_path = resolve_pdb_path(
            str(pdb),
            self.molecule_type,
            clean_pdb=self.clean_pdb,
            keep_chains=self.keep_chains,
            remove_h=self.remove_h,
            drop_hetatm=self.drop_hetatm,
            logger=log,
        )
        log.debug("Input PDB original=%s final=%s", original_pdb_path, pdb_path)

        # Choose aptamer FF and residue template
        if self.aptamer_type == "RNA":
            molecule = load_rna_structure()
            nt_list = "GAUC"
            force_field_aptamer = "leaprc.RNA.OL3"
        else:  # DNA
            molecule = load_dna_structure()
            nt_list = "GATC"
            force_field_aptamer = "leaprc.DNA.OL21"

        # Choose ligand FF
        if self.molecule_type == "protein":
            force_field_ligand = "leaprc.protein.ff19SB"
            parameterized = True
        elif self.molecule_type == "organic":
            force_field_ligand = "leaprc.gaff2"
            parameterized = False
        else:
            force_field_ligand = "leaprc.lipid21"
            parameterized = False

        log.debug(
            "Forcefields: aptamer=%s ligand=%s parameterized=%s nt_list=%s",
            force_field_aptamer,
            force_field_ligand,
            parameterized,
            nt_list,
        )

        # Template complex with empty aptamer chain + ligand from PDB
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

        # Ligand-only complex for COM sampling center
        ligand_only = Complex(
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
        )
        ligand_only.add_chain_from_pdb(
            pdb_path=pdb_path,
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            parameterized=parameterized,
        )
        ligand_only.build()

        sampler = space.make_sampler(
            self.shape,
            ligand_only,
            reach=self.reach,
            probe=self.probe,
        )
        rotations = space.NAngles(N_BACKBONE_TORSIONS)

        # Track best candidate across steps
        best_entropy = None
        best_energy = None
        best_sequence = None
        best_positions = None

        if self.verbose:
            log.info("MAWS step 1: selecting first nucleotide")
        else:
            log.debug("Step 1 start")

        # ---- Step 1: choose first nucleotide ----
        for ntide in nt_list:
            energies = []
            free_E = None
            position = None

            cx = copy.deepcopy(cpx)
            aptamer = cx.aptamer_chain()
            aptamer.create_sequence(ntide)
            cx.build()

            positions0 = cx.positions[:]

            for _ in range(self.first_chunk_size):
                pose = sampler.generator()
                rotation = rotations.generator()

                cx.translate_global(aptamer.element, pose.position * unit.angstrom)
                cx.rotate_global(aptamer.element, pose.axis * unit.angstrom, pose.angle)

                for j in range(N_BACKBONE_TORSIONS):
                    aptamer.rotate_in_residue(0, j, rotation[j])

                energy = cx.get_energy()[0]
                if free_E is None or energy < free_E:
                    free_E = energy
                    position = cx.positions[:]
                energies.append(energy)

                cx.positions = positions0[:]

            entropy = S(energies, beta=self.beta)
            log.debug(
                "Step1 candidate=%s entropy=%s best_E=%s",
                aptamer.alias_sequence,
                entropy,
                free_E,
            )

            if best_entropy is None or entropy < best_entropy:
                best_entropy = entropy
                best_energy = free_E
                best_sequence = ntide
                best_positions = position[:]

        log.debug(
            "After step1 best_sequence=%s best_entropy=%s", best_sequence, best_entropy
        )

        # ---- Steps 2..N: grow sequence (append or prepend) ----
        if self.verbose:
            log.info("MAWS steps 2..N: growing sequence")
        for i in range(1, self.num_nucleotides):
            best_old_sequence = best_sequence
            best_old_positions = best_positions[:]

            # per-step selection (same as CLI)
            best_entropy = None
            best_energy = None

            log.debug("Step%d starting best_old_sequence=%s", i + 1, best_old_sequence)

            for ntide in nt_list:
                for append in (True, False):
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

                    cx.rebuild()
                    cx.pert_min(size=0.5)

                    positions0 = cx.positions[:]

                    for _ in range(self.second_chunk_size):
                        rotation = rotations.generator()

                        # Forward rotations on the new residue’s internal bonds
                        for j in range(N_BACKBONE_TORSIONS - 1):
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

                    entropy = S(energies, beta=self.beta)
                    log.debug(
                        "Step%d candidate=%s (%s %s) entropy=%s best_E=%s",
                        i + 1,
                        aptamer.alias_sequence,
                        "append" if append else "prepend",
                        ntide,
                        entropy,
                        free_E,
                    )

                    if best_entropy is None or entropy < best_entropy:
                        best_entropy = entropy
                        best_energy = free_E
                        best_sequence = aptamer.alias_sequence
                        best_positions = position[:]

        # Optional final PDB artifact
        written_pdb = None
        if output_pdb is not None:
            out = Path(output_pdb)
            if out.exists() and out.is_dir():
                out = out / f"{name}_RESULT.pdb"
            else:
                out.parent.mkdir(parents=True, exist_ok=True)

            result_complex = copy.deepcopy(cpx)
            aptamer = result_complex.aptamer_chain()
            aptamer.create_sequence(best_sequence)
            result_complex.build()
            result_complex.positions = best_positions[:]

            with open(out, "w") as f:
                app.PDBFile.writeModel(
                    result_complex.topology,
                    result_complex.positions,
                    file=f,
                )
            written_pdb = str(out)
            log.debug("Wrote PDB artifact: %s", written_pdb)

        if self.verbose:
            log.info("MAWS run finished: name=%s sequence=%s", name, best_sequence)
        else:
            log.debug("MAWS run finished: name=%s sequence=%s", name, best_sequence)

        return MawsResult(
            sequence=str(best_sequence),
            energy=float(best_energy) if best_energy is not None else float("nan"),
            entropy=float(best_entropy) if best_entropy is not None else float("nan"),
            pdb_path=written_pdb,
        )

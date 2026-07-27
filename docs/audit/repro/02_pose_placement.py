"""Reproduce findings B1, B2 and B3: the sampled pose is a displacement, and the
SAS filter does not prevent aptamer-target clashes.

Requires a prebuilt ligand library. The repo ships one at data/LIG.lib (1408 atoms,
built from data/1BRQ_cleaned.pdb), which this script uses so that no antechamber
run is needed.

    conda activate maws
    python docs/audit/repro/02_pose_placement.py

Runtime: ~1 minute, dominated by two tleap builds. Artifacts go to a temp dir.
"""

import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
from openmm import unit
from scipy.spatial import KDTree

import maws.space as space
from maws.complex import Complex
from maws.helpers import nostrom
from maws.rna_structure import load_rna_structure
from maws.structure import Structure

REPO = Path(__file__).resolve().parents[3]
LIG_LIB = REPO / "data" / "LIG.lib"
LIG_ATOMS = 1408
N_POSES = 300


def main():
    warnings.filterwarnings("ignore")
    # maws seeds nothing (finding B5); seed the global state so this script is
    # reproducible. Without this the clash fractions below move by a few points
    # between runs - which is itself worth noticing.
    np.random.seed(0)
    if not LIG_LIB.exists():
        raise SystemExit(f"missing {LIG_LIB}; build it first with a protein-type run")

    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        shutil.copy(LIG_LIB, d / "LIG.lib")
        cwd = Path.cwd()
        import os

        os.chdir(d)  # Complex.build writes out.in / .maws_cache into cwd
        try:
            lig = Structure(["LIG"], residue_length=[LIG_ATOMS], residue_path=".")

            cpx = Complex("leaprc.RNA.OL3", "leaprc.protein.ff19SB")
            cpx.add_chain("", load_rna_structure())
            cpx.add_chain("LIG", lig)
            apt = cpx.aptamer_chain()
            apt.create_sequence("G")
            cpx.build()
            a0, a2 = apt.element[0], apt.element[2]
            pristine = cpx.positions[:]

            ligand_only = Complex("leaprc.RNA.OL3", "leaprc.protein.ff19SB")
            ligand_only.add_chain("LIG", lig)
            ligand_only.build()

            dims = space.compute_envelope_dims(ligand_only, reach=10.0)
            sampler = space.make_sampler(ligand_only, reach=10.0, probe=1.4)

            pos = np.array(nostrom(pristine))
            print("=" * 78)
            print("B1: the sampled pose is applied as a DISPLACEMENT, not a target")
            print("=" * 78)
            print(
                f"  aptamer centroid as LEaP built it : "
                f"{pos[a0:a2].mean(axis=0).round(2)}"
            )
            print(
                f"  target centroid                   : "
                f"{pos[a2:].mean(axis=0).round(2)}"
            )
            print(
                f"  envelope: centre {np.round(dims['centre'], 2)} "
                f"radius {dims['radius']:.2f} A"
            )

            target = pos[a2:]
            tree = KDTree(target)
            offs, mind = [], []
            for _ in range(N_POSES):
                cpx.positions = pristine[:]
                p = sampler.generator()
                cpx.translate_global(apt.element, p.position * unit.angstrom)
                cpx.rotate_global(apt.element, p.axis * unit.angstrom, p.angle)
                a = np.array(nostrom(cpx.positions))[a0:a2]
                offs.append(np.linalg.norm(a.mean(axis=0) - p.position))
                mind.append(tree.query(a, k=1)[0].min())
            offs, mind = np.array(offs), np.array(mind)

            print(
                f"\n  |aptamer centroid after move - sampled point|, {N_POSES} poses:"
            )
            print(
                f"     mean {offs.mean():.2f}  min {offs.min():.2f}  "
                f"max {offs.max():.2f} A"
            )
            print("  -> should be ~0 if the pose were used as a target position.")
            print("     (B3 adds to this: rotate_global pivots on the aptamer's FIRST")
            print(
                "      atom, not its centroid, so the rotation displaces it further.)"
            )

            print()
            print("=" * 78)
            print("B2: poses the SAS filter accepts as 'clear' still clash")
            print("=" * 78)
            print("  true min aptamer-atom -> target-atom distance:")
            print(f"     mean {mind.mean():.2f} A   min {mind.min():.2f} A")
            for thr in (1.0, 2.0, 2.6, 3.0):
                print(
                    f"     fraction of ACCEPTED poses with a contact < {thr:.1f} A : "
                    f"{(mind < thr).mean():.3f}"
                )
            print("  -> Excluder.is_clear() tests the sampled POINT, not the ~30 atoms")
            print("     that actually get placed there.")
        finally:
            os.chdir(cwd)


if __name__ == "__main__":
    main()

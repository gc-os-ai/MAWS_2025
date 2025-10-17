# examples/complex_aptamer.py
#
# Linear walkthrough of maws.complex.Complex:
# - Build a complex (empty DNA chain + ligand from PDB)
# - Place chain (translate/rotate)
# - Probe a torsion
# - Append a base and rebuild (coordinate splicing)
# - Minimize and do a few MD steps
#
# Requirements:
#   - AmberTools on PATH (tleap, antechamber, parmchk2) for add_chain_from_pdb()
#   - OpenMM installed

from __future__ import annotations

import copy
import math
from pathlib import Path

import numpy as np
from openmm import app, unit

from maws.complex import Complex
from maws.dna_structure import load_dna_structure  # switch to RNA if you like

# ---- config -----------------------------------------------------------------
DATA_PDB = Path(__file__).resolve().parent.parent / "data" / "pfoa.pdb"
FF_APT = "leaprc.DNA.OL21"  # DNA aptamer
FF_LIG = "leaprc.gaff2"  # organic ligand

# ---- 1) make complex: empty aptamer chain + ligand from PDB ------------------
dna = load_dna_structure()  # Structure defines lengths/aliases/torsions
cpx = Complex(FF_APT, FF_LIG)  # sets build preamble (LEaP sources)

cpx.add_chain("", dna)  # Chain[0] = aptamer (empty for now)
cpx.add_chain_from_pdb(  # Chain[1] = ligand (1 residue)
    pdb_path=str(DATA_PDB),
    force_field_aptamer=FF_APT,
    force_field_ligand=FF_LIG,
    pdb_name="LIG",
    parameterized=False,
)

apt = cpx.chains[0]

# ---- 2) seed a 1-nt sequence and build (hits cache on repeats) --------------
apt.create_sequence("G")  # alias → canonical via Structure.translate
cpx.build()
E, _ = cpx.get_energy()
print(f"[build] seq={apt.alias_sequence}  E={E:.2f} kJ/mol")

with open("01_initial_G.pdb", "w") as fh:
    app.PDBFile.writeModel(cpx.topology, cpx.positions, file=fh, modelIndex=1)

# ---- 3) place the aptamer near the ligand (translate + rotate) --------------
apt.translate_global(np.array([5.0, 0.0, 0.0]) * unit.angstrom)  # +5 Å along x
apt.rotate_global(
    np.array([0.0, 1.0, 0.0]) * unit.angstrom, math.radians(30)
)  # 30° about y
E, _ = cpx.get_energy()
print(f"[place] +5 Å (x), rot 30° (y)  E={E:.2f} kJ/mol")

with open("02_placed.pdb", "w") as fh:
    app.PDBFile.writeModel(cpx.topology, cpx.positions, file=fh, modelIndex=1)

# ---- 4) quick torsion probe: residue 0, torsion index 0 ---------------------
pos0 = copy.deepcopy(cpx.positions)
for theta in (0, 60, 120, 180):
    cpx.positions = copy.deepcopy(pos0)
    apt.rotate_in_residue(0, 0, math.radians(theta))  # (res_idx=0, torsion=0)
    E, _ = cpx.get_energy()
    print(f"[torsion] residue 0, torsion 0 @ {theta:>3}°  ->  E={E:.2f} kJ/mol")
cpx.positions = pos0  # restore

# ---- 5) grow sequence and rebuild (coordinate splicing around the edit) -----
apt.append_sequence("A")  # alias base, e.g. "A"
cpx.rebuild()  # cached AMBER + splice old/new coords at the junction
E, _ = cpx.get_energy()
print(f"[rebuild] appended 'A' -> seq={apt.alias_sequence}  E={E:.2f} kJ/mol")

with open("03_after_append_A.pdb", "w") as fh:
    app.PDBFile.writeModel(cpx.topology, cpx.positions, file=fh, modelIndex=1)

# ---- 6) minimize and a tiny MD nudge ----------------------------------------
Emin = cpx.minimize(max_iterations=700)
print(f"[minimize] 700 iters  ->  E={Emin:.2f} kJ/mol")

cpx.step(1000)  # very short MD
E, _ = cpx.get_energy()
print(f"[md] 1000 steps  ->  E={E:.2f} kJ/mol")

with open("04_relaxed.pdb", "w") as fh:
    app.PDBFile.writeModel(cpx.topology, cpx.positions, file=fh, modelIndex=1)

# ---- 7) (optional) small random kicks + minimize ----------------------------
cpx.pert_min(size=0.2, iterations=5)
E, _ = cpx.get_energy()
print(f"[pert_min] size=0.2 Å, cycles=5  ->  E={E:.2f} kJ/mol")

with open("05_pert_min.pdb", "w") as fh:
    app.PDBFile.writeModel(cpx.topology, cpx.positions, file=fh, modelIndex=1)

print("\nDone. Cache is under .maws_cache/ (reused automatically).")

# examples/basic_complex.py

from maws.complex import Complex
from maws.dna_structure import load_dna_structure

dna = load_dna_structure()  # Structure
pdb_path = "../MAWS_2025/data/pfoa.pdb"

force_field_aptamer = "leaprc.DNA.OL21"
force_field_ligand = "leaprc.gaff2"

cpx = Complex(
    force_field_aptamer=force_field_aptamer,
    force_field_ligand=force_field_ligand,
)

# Start with an empty DNA chain (sequence added later)
cpx.add_chain("", dna)

# Add the ligand as a one-residue chain from PDB
# (writes LIG.lib / LIG.frcmod next to the PDB)
cpx.add_chain_from_pdb(
    pdb_path=pdb_path,
    force_field_aptamer=force_field_aptamer,
    force_field_ligand=force_field_ligand,
    pdb_name="LIG",
    parameterized=False,
)

cpx.build()

E, _ = cpx.get_energy()
print("Potential energy (kJ/mol):", E)

"""Reproduce findings E1, E2 and E3: the PDB cleaner silently deletes and reorders.

    conda activate maws
    python docs/audit/repro/04_pdb_cleaner.py

Runtime: seconds. Writes into a temp dir; data/ is not modified.
"""

import collections
import logging
import shutil
import tempfile
from pathlib import Path

from maws.pdb_cleaner import resolve_pdb_path

REPO = Path(__file__).resolve().parents[3]
LOG = logging.getLogger("repro")


def summarize(path, tag):
    lines = Path(path).read_text().splitlines()
    atoms = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    ters = [ln for ln in lines if ln.startswith("TER")]
    het = collections.Counter(
        ln[17:20].strip() for ln in atoms if ln.startswith("HETATM")
    )
    chains = collections.Counter(ln[21] for ln in atoms)
    ins = {(ln[21], ln[22:27].strip()) for ln in atoms if ln[26] != " "}
    print(f"  [{tag}] atoms={len(atoms)}  TER={len(ters)}  chains={dict(chains)}")
    print(f"           HETATM residues={dict(het)}  insertion-coded atoms={len(ins)}")


def atom_line(serial, name, res, chain, seq, alt=" ", occ=1.00, el="C"):
    return (
        f"ATOM  {serial:5d} {name:<4s}{alt}{res} {chain}{seq:4d}    "
        f"{serial:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{0.0:6.2f}          {el:>2s}\n"
    )


def main():
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)

        print("=" * 78)
        print("E2 + E3: real target, data/1HAO.pdb (thrombin + 15mer + inhibitor)")
        print("=" * 78)
        src = REPO / "data" / "1HAO.pdb"
        shutil.copy(src, d / "1HAO.pdb")
        summarize(src, "ORIGINAL")
        out, _ = resolve_pdb_path(
            str(d / "1HAO.pdb"),
            "protein",
            clean_pdb=True,
            keep_chains="all",
            remove_h=False,
            drop_hetatm=False,
            logger=LOG,
        )
        summarize(out, "CLEANED ")
        print("  -> E2: insertion-coded residues deleted, including a TER record, so")
        print("         3 chains become 2 and chain L is now bonded to chain H.")
        print("  -> E3: the 0G6 inhibitor and all waters are gone despite")
        print("         drop_hetatm=False, because they follow the last TER.")

        print()
        print("=" * 78)
        print("E1: a single altLoc re-sorts the whole file")
        print("=" * 78)
        rows, serial = [], 1
        for seq, res in [(8, "ALA"), (9, "GLY"), (10, "SER"), (11, "VAL")]:
            for nm, el in [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]:
                if seq == 10:  # one altLoc pair
                    rows.append(atom_line(serial, nm, res, "A", seq, "A", 0.40, el))
                    serial += 1
                    rows.append(atom_line(serial, nm, res, "A", seq, "B", 0.60, el))
                    serial += 1
                else:
                    rows.append(atom_line(serial, nm, res, "A", seq, el=el))
                    serial += 1
        rows.append(f"TER   {serial:5d}      VAL A  11\nEND\n")
        (d / "alt.pdb").write_text("".join(rows))

        order_in = [ln[22:26].strip() for ln in rows if ln.startswith("ATOM")][::4]
        print(f"  INPUT residue order : {order_in}")
        out, _ = resolve_pdb_path(
            str(d / "alt.pdb"),
            "protein",
            clean_pdb=True,
            keep_chains="all",
            remove_h=False,
            drop_hetatm=False,
            logger=LOG,
        )
        print("  CLEANED file:")
        for ln in Path(out).read_text().splitlines():
            print(f"    {ln}")
        print("  -> residues emitted 10, 11, 8, 9; TER placed before its own residue;")
        print("     atoms alphabetised within each residue (C, CA, N, O).")
        print("     tleap builds connectivity from file order.")


if __name__ == "__main__":
    main()

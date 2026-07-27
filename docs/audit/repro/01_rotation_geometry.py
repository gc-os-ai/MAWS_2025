"""Reproduce findings A1, A2 and A4: the rotation kernel breaks covalent geometry.

Builds a real RNA G5-A3 dinucleotide with tleap, then drives MAWS's own
Chain.rotate_in_residue over it and measures (a) how many bonds change length and
(b) how much rotation is actually applied versus how much was requested.

    conda activate maws
    python docs/audit/repro/01_rotation_geometry.py

Runtime: a few seconds. Writes scratch files to a temp dir.
"""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from openmm import unit
from openmm.app import AmberPrmtopFile, PDBFile

from maws import complex as maws_complex
from maws.complex import Complex
from maws.rna_structure import load_rna_structure

LEAP = """source leaprc.RNA.OL3
source leaprc.DNA.OL21
r1 = sequence {{ G5 A3 }}
savepdb r1 {d}/rna_GA.pdb
saveamberparm r1 {d}/rna_GA.prmtop {d}/rna_GA.inpcrd
d1 = sequence {{ DG5 DA3 }}
savepdb d1 {d}/dna_GA.pdb
g1 = sequence {{ GN }}
savepdb g1 {d}/rna_GN.pdb
dgn = sequence {{ DGN }}
savepdb dgn {d}/dna_GN.pdb
quit
"""


def build_reference(d: Path):
    (d / "leap.in").write_text(LEAP.format(d=d))
    subprocess.run(
        ["tleap", "-f", str(d / "leap.in")],
        cwd=d,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def dihedral(positions, a, b, c, e):
    p = np.array(positions.value_in_unit(unit.angstrom))
    b0, b1, b2 = p[a] - p[b], p[c] - p[b], p[e] - p[c]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def main():
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        build_reference(d)

        pdb = PDBFile(str(d / "rna_GA.pdb"))
        prm = AmberPrmtopFile(str(d / "rna_GA.prmtop"))
        bonds = [(b[0].index, b[1].index) for b in prm.topology.bonds()]

        def fresh():
            c = Complex()
            c.add_chain("G A", load_rna_structure())
            c.positions = pdb.positions[:]
            return c

        def bond_lengths(pos):
            a = np.array(pos.value_in_unit(unit.angstrom))
            return np.array([np.linalg.norm(a[i] - a[j]) for i, j in bonds])

        ref = bond_lengths(pdb.positions)
        theta = 0.7

        print("=" * 78)
        print("A4: atom ordering as LEaP actually builds it")
        print("=" * 78)
        for f in ("rna_GA.pdb", "dna_GA.pdb", "rna_GN.pdb", "dna_GN.pdb"):
            for r in PDBFile(str(d / f)).topology.residues():
                names = [a.name for a in r.atoms()]
                print(
                    f"  {f:12s} {r.name:4s} n={len(names):2d}  "
                    f"idx-1={names[-1]:6s} idx-2={names[-2]:6s} "
                    f"idx-4={names[-4]:6s} idx-7={names[-7]:6s}"
                )
        print("  -> DNA XN/X3 use rotation spec (-4, -2) = H2' -> O3', not C3' -> O3'.")

        print()
        print("=" * 78)
        print("A1: how many times is the kernel called, over which atom range?")
        print("=" * 78)
        calls = []
        original = maws_complex.Complex.rotate_global

        def spy(self, element, axis, angle, reverse=False, glob=True):
            calls.append((list(element), glob))
            return original(self, element, axis, angle, reverse=reverse, glob=glob)

        maws_complex.Complex.rotate_global = spy
        for j, label in enumerate(["j=0 5'-OH", "j=1 beta", "j=2 chi", "j=3 epsilon"]):
            calls.clear()
            fresh().aptamer_chain().rotate_in_residue(0, j, theta)
            ranges = " ".join(f"[{e[0] if g else e[1]}:{e[2]})" for e, g in calls)
            print(f"  {label:12s} -> {len(calls)} calls, ranges: {ranges}")
        maws_complex.Complex.rotate_global = original
        print("  -> 3 calls per torsion; negative starts wrap to the end of the array.")

        print()
        print("=" * 78)
        print("A1/A2: covalent bonds after ONE requested torsion (theta = 0.7 rad)")
        print("=" * 78)
        print(f"  reference bond lengths span {ref.min():.3f} - {ref.max():.3f} A")

        def report(tag, c):
            d_ = np.abs(bond_lengths(c.positions) - ref)
            print(
                f"  {tag:50s} max|dL|={d_.max():7.3f} A  "
                f"broken={int((d_ > 0.05).sum()):2d}/{len(bonds)}"
            )

        print("\n  -- forward (append-style) --")
        for res, j, name in [
            (1, 0, "alpha"),
            (1, 1, "beta"),
            (1, 2, "chi"),
            (0, 3, "epsilon"),
        ]:
            c = fresh()
            c.aptamer_chain().rotate_in_residue(res, j, theta)
            report(f"rotate_in_residue({res}, j={j})  {name}", c)

        print("\n  -- reverse=True (every prepend move) --")
        for j in range(4):
            c = fresh()
            c.aptamer_chain().rotate_in_residue(0, j, theta, reverse=True)
            report(f"rotate_in_residue(0, j={j}, reverse=True)", c)

        print()
        print("=" * 78)
        print("A1: net angle applied vs requested, on a torsion that does not break")
        print("=" * 78)
        for req in (0.3, 0.7, 1.5):
            c = fresh()
            before = dihedral(c.positions, 0, 1, 2, 5)
            c.aptamer_chain().rotate_in_residue(0, 1, req)
            after = dihedral(c.positions, 0, 1, 2, 5)
            delta = (after - before + 180) % 360 - 180
            wrapped = "  (3*theta wrapped past 180 deg)" if abs(req * 3) > np.pi else ""
            print(
                f"  requested {np.degrees(req):7.2f} deg -> applied {delta:8.2f} deg"
                f"   ratio {delta / np.degrees(req):+.3f}{wrapped}"
            )
        print("  -> ratio -3.000: the torsion is applied three times. The sign is the")
        print(
            "     kernel's transposed-matrix convention (v.R, i.e. rotation by -theta)."
        )


if __name__ == "__main__":
    main()

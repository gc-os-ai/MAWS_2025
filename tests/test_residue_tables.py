"""The RNA/DNA residue tables must describe the molecules LEaP actually builds.

rna_structure.py and dna_structure.py hardcode atom counts and torsion
definitions as bare integer offsets. Nothing in the library checks them
against a real structure, so a transcription slip stays invisible until it
shows up as bad geometry. See issue #47 (finding A4).
"""

import subprocess
import tempfile
from pathlib import Path

import pytest
from openmm.app import AmberPrmtopFile

try:
    from maws.dna_structure import load_dna_structure
    from maws.rna_structure import load_rna_structure
    from maws.tools import find_exe

    try:
        find_exe("tleap")
        HAS_AMBERTOOLS = True
    except Exception:
        HAS_AMBERTOOLS = False
    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False
    HAS_AMBERTOOLS = False

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not HAS_DEPS, reason="OpenMM/maws not available"),
    pytest.mark.skipif(not HAS_AMBERTOOLS, reason="AmberTools (tleap) not available"),
]


def _leap_units():
    """(sequence, unit name, residue names in the order we asked for them).

    A trinucleotide gives the 5', internal and 3' forms; a lone residue gives
    the N form. Between them that is all 16 templates per polymer.
    """
    units = []
    for prefix, bases in (("", "GAUC"), ("D", "GATC")):
        for base in bases:
            p = f"{prefix}{base}"
            units.append((f"{p}5 {p} {p}3", f"tri{p}", [f"{p}5", p, f"{p}3"]))
            units.append((f"{p}N", f"mono{p}", [f"{p}N"]))
    return units


@pytest.fixture(scope="session")
def residue_templates():
    """resname -> (atom names, intra-residue bonds as residue-local pairs).

    Residues are keyed by the name we asked LEaP for, not the name it reports
    back: Amber labels both G5 and G3 as plain "G", so trusting the label
    would silently collapse three templates into one.
    """
    units = _leap_units()
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        lines = ["source leaprc.RNA.OL3", "source leaprc.DNA.OL21"]
        for sequence, unit_name, _ in units:
            lines += [
                f"{unit_name} = sequence {{ {sequence} }}",
                f"saveamberparm {unit_name} {d}/{unit_name}.prmtop "
                f"{d}/{unit_name}.inpcrd",
            ]
        lines.append("quit")
        (d / "leap.in").write_text("\n".join(lines))
        subprocess.run(
            [find_exe("tleap"), "-f", str(d / "leap.in")],
            cwd=d,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        templates = {}
        for _sequence, unit_name, wanted in units:
            topology = AmberPrmtopFile(str(d / f"{unit_name}.prmtop")).topology
            bonds = {(b[0].index, b[1].index) for b in topology.bonds()}
            for residue, resname in zip(topology.residues(), wanted, strict=True):
                indices = [a.index for a in residue.atoms()]
                first = min(indices)
                templates[resname] = (
                    [a.name for a in residue.atoms()],
                    {
                        (min(i, j) - first, max(i, j) - first)
                        for i, j in bonds
                        if i in indices and j in indices
                    },
                )
        return templates


def _all_residues():
    """(loader, polymer tag, residue name) for every template in both tables."""
    if not HAS_DEPS:
        return []
    cases = []
    for loader, tag in ((load_rna_structure, "RNA"), (load_dna_structure, "DNA")):
        cases += [(loader, tag, name) for name in loader().residue_names]
    return cases


RESIDUES = _all_residues()
RESIDUE_IDS = [f"{tag}-{name}" for _loader, tag, name in RESIDUES]


@pytest.mark.parametrize(("loader", "tag", "resname"), RESIDUES, ids=RESIDUE_IDS)
def test_torsion_axes_are_covalent_bonds(residue_templates, loader, tag, resname):
    """Every torsion must turn about a real bond.

    Rotating about a pair of atoms that are not bonded is not a torsion at
    all - it swings a fragment about an arbitrary line through the molecule.
    """
    structure = loader()
    atom_names, bonds = residue_templates[resname]
    length = structure.residue_length[resname]

    failures = []
    for index, spec in enumerate(structure.rotating_elements[resname]):
        start = spec[0] + length if spec[0] < 0 else spec[0]
        bond = spec[1] + length if spec[1] < 0 else spec[1]
        if (min(start, bond), max(start, bond)) not in bonds:
            failures.append(
                f"j={index} spec {spec} -> atoms {start},{bond} = "
                f"{atom_names[start]}-{atom_names[bond]}"
            )

    assert not failures, f"{tag} {resname}: axis is not a covalent bond; " + "; ".join(
        failures
    )


@pytest.mark.parametrize(("loader", "tag", "resname"), RESIDUES, ids=RESIDUE_IDS)
def test_residue_length_matches_the_built_template(
    residue_templates, loader, tag, resname
):
    """RESIDUE_LENGTH is what every negative index is resolved against.

    If it drifts from the real atom count, every negative offset in the
    torsion and connectivity tables silently points at the wrong atom.
    """
    atom_names, _bonds = residue_templates[resname]
    declared = loader().residue_length[resname]
    assert declared == len(atom_names), (
        f"{tag} {resname}: table says {declared} atoms, LEaP builds {len(atom_names)}"
    )

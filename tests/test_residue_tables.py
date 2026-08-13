"""Checks the residue tables against the molecules AmberTools actually builds.

A *residue* is one letter of an aptamer: G, A, U or C for RNA, and the same
with a leading D for DNA. Each comes in four forms depending on where it sits
in the strand, which is why there are sixteen per polymer. A *torsion* is a
turn of part of a molecule about one of its own bonds.

:mod:`maws.rna_structure` and :mod:`maws.dna_structure` describe every residue
as bare integers: how many atoms it has, and which atom positions form each
torsion. Those integers are written by hand, and a wrong one still runs. These
tests build all thirty-two residues with LEaP, the AmberTools program that
turns a sequence into a molecule, and compare the tables against the result.
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
    """Return the LEaP units needed to cover every residue form.

    A strand of three residues yields the start form, the middle form and the
    end form; a strand of one yields the standalone form. Four bases times two
    polymers times those two strands reaches all thirty-two.

    Returns
    -------
    list of tuple
        ``(sequence, unit_name, resnames)`` per unit, where `sequence` is the
        LEaP sequence string, `unit_name` names the files it writes, and
        `resnames` lists the residues in the order they were requested.
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
    """Build every residue once and return its atoms and internal bonds.

    Returns
    -------
    dict
        Residue name to ``(atom_names, bonds)``. `atom_names` is in file
        order; `bonds` holds index pairs counted from the start of that
        residue, covering only bonds with both ends inside it.

    Notes
    -----
    Residues are keyed by the name requested from LEaP rather than the name it
    reports back. Amber labels both the start form and the end form of G as
    plain ``"G"``, so keying on the reported name merges three residues into
    one.
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
    """Return one parametrise case per residue across both tables.

    Returns
    -------
    list of tuple
        ``(loader, tag, resname)``, where `loader` builds the
        :class:`~maws.structure.Structure` holding that residue and `tag` is
        ``"RNA"`` or ``"DNA"``.
    """
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
    """Every torsion turns about two atoms that are bonded to each other."""
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
    """Every declared residue length matches the atom count LEaP builds."""
    atom_names, _bonds = residue_templates[resname]
    declared = loader().residue_length[resname]
    assert declared == len(atom_names), (
        f"{tag} {resname}: table says {declared} atoms, LEaP builds {len(atom_names)}"
    )

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
from collections import defaultdict, deque
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


@pytest.fixture(scope="session")
def assembled_chains():
    """Build every LEaP unit once and return it as a whole strand.

    Returns
    -------
    dict
        Unit name to ``(resnames, offsets, n_atoms, bonds)``. `offsets` gives
        each residue's first atom index within the strand, and `bonds` holds
        index pairs counted across the whole strand.

    See Also
    --------
    residue_templates : The same molecules split into single residues.

    Notes
    -----
    Whole strands are needed because a torsion whose table entry ends in
    ``None`` moves every atom to the end of the strand, not to the end of the
    residue it belongs to.
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

        chains = {}
        for _sequence, unit_name, wanted in units:
            topology = AmberPrmtopFile(str(d / f"{unit_name}.prmtop")).topology
            offsets, tally = [], 0
            for residue in topology.residues():
                offsets.append(tally)
                tally += sum(1 for _ in residue.atoms())
            chains[unit_name] = (
                wanted,
                offsets,
                topology.getNumAtoms(),
                {(b[0].index, b[1].index) for b in topology.bonds()},
            )
        return chains


def moving_set_from_bonds(bonds, fixed, pivot):
    """Return the atoms joined to `pivot` once the `fixed`-`pivot` bond is cut.

    Written out here rather than calling
    :meth:`maws.complex.Complex.moving_set`, which computes the same thing: a
    test that checks a function against itself passes whatever the function
    does.

    Parameters
    ----------
    bonds : iterable of tuple of int
        Index pairs, one per bond.
    fixed : int
        Index of the bonded atom whose side is discarded.
    pivot : int
        Index of the bonded atom whose side is returned. Included in the
        result.

    Returns
    -------
    set of int
        Atom indices reachable from `pivot` without crossing the cut bond.

    Examples
    --------
    A four-atom chain cut in the middle splits into two pairs.

    >>> sorted(moving_set_from_bonds([(0, 1), (1, 2), (2, 3)], 1, 2))
    [2, 3]
    """
    adjacency = defaultdict(set)
    for a, b in bonds:
        adjacency[a].add(b)
        adjacency[b].add(a)

    seen, queue, reached = {fixed}, deque([pivot]), {pivot}
    while queue:
        for neighbour in adjacency[queue.popleft()]:
            if neighbour not in seen and neighbour not in reached:
                reached.add(neighbour)
                queue.append(neighbour)
    return reached


UNIT_NAMES = [name for _seq, name, _res in _leap_units()] if HAS_DEPS else []


@pytest.mark.parametrize("unit_name", UNIT_NAMES)
def test_table_ranges_match_the_bond_graph(assembled_chains, unit_name):
    """Every hand-written atom range names the atoms the bonds say it should."""
    resnames, offsets, n_atoms, bonds = assembled_chains[unit_name]
    is_dna = resnames[0].startswith("D")
    structure = (load_dna_structure if is_dna else load_rna_structure)()

    mismatches = []
    for residue_index, resname in enumerate(resnames):
        offset = offsets[residue_index]
        length = structure.residue_length[resname]
        for j, spec in enumerate(structure.rotating_elements[resname]):
            start = (spec[0] + length if spec[0] < 0 else spec[0]) + offset
            pivot = (spec[1] + length if spec[1] < 0 else spec[1]) + offset
            end = (
                n_atoms
                if spec[2] is None
                else (spec[2] + length if spec[2] < 0 else spec[2]) + offset
            )

            from_table = set(range(pivot, end))
            from_graph = moving_set_from_bonds(bonds, start, pivot)
            if from_table != from_graph:
                mismatches.append(
                    f"{resname} j={j} spec {spec}: table {len(from_table)} atoms, "
                    f"graph {len(from_graph)}, "
                    f"differ by {sorted(from_table ^ from_graph)}"
                )

    assert not mismatches, f"{unit_name}: " + "; ".join(mismatches)


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

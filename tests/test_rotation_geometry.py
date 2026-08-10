"""Geometric invariants for the rotation kernel.

A torsion rotation is allowed to change dihedral angles. It is never
allowed to change a covalent bond length. See issue #47.
"""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest
from openmm import unit
from openmm.app import AmberPrmtopFile, PDBFile

try:
    from maws.complex import Complex
    from maws.rna_structure import load_rna_structure
    from maws.tools import find_exe

    try:
        find_exe("tleap")
        HAS_AMBERTOOLS = True
    except Exception:
        HAS_AMBERTOOLS = False
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
    HAS_AMBERTOOLS = False

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available"),
    pytest.mark.skipif(not HAS_AMBERTOOLS, reason="AmberTools (tleap) not available"),
]

LEAP_INPUT = """source leaprc.RNA.OL3
r1 = sequence {{ G5 A3 }}
savepdb r1 {d}/ref.pdb
saveamberparm r1 {d}/ref.prmtop {d}/ref.inpcrd
quit
"""

# How much is a bond allowed to change length when we rotate a torsion?
#
# In theory: not at all. Rotating swings atoms around an axis. It should
# never stretch or squash a chemical bond. But we cannot test for exactly
# zero, because computers round every arithmetic result a tiny bit.
#
# So we measured both ends on real structures:
#
#   a rotation that WORKS  changes bonds by 0.000000000000002 A  <- rounding
#   the bugs in #47        changed bonds by 0.054 A to 8.484 A   <- real damage
#
# Any threshold between those two numbers does the job. We picked 0.001 A:
#   - large enough that rounding noise never trips it (500 million x bigger)
#   - small enough that even the mildest real bug did trip it (54 x smaller)
#
# If this ever looks too strict, re-measure both ends before changing it.
BOND_TOLERANCE_ANGSTROM = 1e-3


def dihedral(positions, a, b, c, d):
    """Torsion angle a-b-c-d in radians, about the b-c bond.

    Do not "simplify" this to maws.helpers. helpers.angle is unsigned, so
    it cannot see a rotation's direction at all, and directed_angle shares
    a module - and a sign convention - with the kernel under test.
    """
    xyz = np.asarray(positions.value_in_unit(unit.angstrom), dtype=float)
    arm_in = xyz[a] - xyz[b]
    axis = xyz[c] - xyz[b]
    arm_out = xyz[d] - xyz[c]

    axis = axis / np.linalg.norm(axis)
    # Keep only what is perpendicular to the axis - that is what swings.
    v = arm_in - np.dot(arm_in, axis) * axis
    w = arm_out - np.dot(arm_out, axis) * axis
    return float(np.arctan2(np.dot(np.cross(axis, v), w), np.dot(v, w)))


def angle_difference(after, before):
    """after - before, wrapped into (-pi, pi]."""
    return (after - before + np.pi) % (2 * np.pi) - np.pi


def bond_lengths(positions, bonds):
    """Length of every covalent bond, in angstrom.

    OpenMM stores nanometres; forgetting the conversion makes every length
    10x too small and the tolerance meaningless.
    """
    xyz = np.asarray(positions.value_in_unit(unit.angstrom), dtype=float)
    left = xyz[[i for i, _ in bonds]]
    right = xyz[[j for _, j in bonds]]
    return np.linalg.norm(left - right, axis=1)


@pytest.fixture(scope="session")
def reference_dinucleotide():
    """LEaP-built RNA G5-A3. Returns (positions, bonds).

    The bond list comes from the prmtop, i.e. from AmberTools - never
    hand-written and never derived from the code under test.
    """
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "leap.in").write_text(LEAP_INPUT.format(d=d))
        subprocess.run(
            [find_exe("tleap"), "-f", str(d / "leap.in")],
            cwd=d,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        pdb = PDBFile(str(d / "ref.pdb"))
        prm = AmberPrmtopFile(str(d / "ref.prmtop"))
        bonds = [(b[0].index, b[1].index) for b in prm.topology.bonds()]
        return pdb.positions, bonds


@pytest.fixture
def dinucleotide_complex(reference_dinucleotide):
    """A fresh Complex per test, on the reference coordinates.

    positions[:] copies: the session fixture is shared, so in-place
    mutation would leak between tests.
    """
    positions, _bonds = reference_dinucleotide
    cpx = Complex()
    cpx.add_chain("G A", load_rna_structure())
    cpx.positions = positions[:]
    return cpx


def assert_bonds_preserved(before, after, bonds, what):
    """Fail with the worst offending bond named, not just a boolean."""
    deviation = np.abs(after - before)
    worst = int(deviation.argmax())
    n_broken = int((deviation > BOND_TOLERANCE_ANGSTROM).sum())
    assert deviation[worst] < BOND_TOLERANCE_ANGSTROM, (
        f"{what} changed bond {bonds[worst]} by {deviation[worst]:.3f} A "
        f"({before[worst]:.3f} -> {after[worst]:.3f} A); "
        f"{n_broken} of {len(bonds)} bonds broken"
    )


# The four torsions run.py drives when it APPENDS a nucleotide: j=0,1,2 on
# the new 3' residue, j=3 on the one before it (run.py:315-327).
APPEND_TORSIONS = [
    (1, 0, "alpha  P-O5'"),
    (1, 1, "beta   O5'-C5'"),
    (1, 2, "chi    C1'-N9"),
    (0, 3, "epsilon C3'-O3'"),
]

# The same four when it PREPENDS: all on the new 5' residue, reverse=True.
PREPEND_TORSIONS = [(0, j, APPEND_TORSIONS[j][2]) for j in range(4)]


@pytest.mark.parametrize(
    ("residue", "torsion", "label"), APPEND_TORSIONS, ids=lambda v: str(v)
)
def test_append_torsion_preserves_bond_lengths(
    reference_dinucleotide, dinucleotide_complex, residue, torsion, label
):
    """No torsion on the append path may change a bond length."""
    _positions, bonds = reference_dinucleotide
    before = bond_lengths(dinucleotide_complex.positions, bonds)

    dinucleotide_complex.aptamer_chain().rotate_in_residue(residue, torsion, 0.7)

    after = bond_lengths(dinucleotide_complex.positions, bonds)
    assert_bonds_preserved(before, after, bonds, f"append {label}")


@pytest.mark.parametrize(
    ("residue", "torsion", "label"), PREPEND_TORSIONS, ids=lambda v: str(v)
)
def test_prepend_torsion_preserves_bond_lengths(
    reference_dinucleotide, dinucleotide_complex, residue, torsion, label
):
    """Same invariant on the prepend path, which uses reverse=True."""
    _positions, bonds = reference_dinucleotide
    before = bond_lengths(dinucleotide_complex.positions, bonds)

    dinucleotide_complex.aptamer_chain().rotate_in_residue(
        residue, torsion, 0.7, reverse=True
    )

    after = bond_lengths(dinucleotide_complex.positions, bonds)
    assert_bonds_preserved(before, after, bonds, f"prepend {label}")


def moved_atom_count(before, after):
    """How many atoms actually changed position."""
    displacement = np.linalg.norm(
        np.asarray(after.value_in_unit(unit.angstrom), dtype=float)
        - np.asarray(before.value_in_unit(unit.angstrom), dtype=float),
        axis=1,
    )
    return int((displacement > 1e-6).sum())


# Every rotation the search relies on, minus prepend alpha - see the test
# below for why that one is necessarily inert.
ROTATIONS_THAT_MUST_MOVE = [
    (1, 0, False, "append alpha"),
    (1, 1, False, "append beta"),
    (1, 2, False, "append chi"),
    (0, 3, False, "append epsilon"),
    (0, 1, True, "prepend beta"),
    (0, 2, True, "prepend chi"),
    (0, 3, True, "prepend epsilon"),
]


@pytest.mark.parametrize(
    ("residue", "torsion", "reverse", "label"),
    ROTATIONS_THAT_MUST_MOVE,
    ids=lambda v: str(v),
)
def test_rotation_actually_moves_atoms(
    dinucleotide_complex, residue, torsion, reverse, label
):
    """A rotation must rotate something.

    Without this, an implementation that silently did nothing would satisfy
    every bond-length assertion in this file.
    """
    before = dinucleotide_complex.positions[:]

    dinucleotide_complex.aptamer_chain().rotate_in_residue(
        residue, torsion, 0.7, reverse=reverse
    )

    n_moved = moved_atom_count(before, dinucleotide_complex.positions)
    assert n_moved > 0, f"{label} moved no atoms at all"


def test_prepend_alpha_is_inert_by_geometry(dinucleotide_complex):
    """Prepend alpha moves nothing, and cannot.

    It turns the fragment upstream of the HO5'-O5' bond. On a 5' terminus
    the only atom upstream is HO5' itself, which lies on the rotation axis.
    So the prepend path has three usable torsions, not four (audit A3).

    Pinned here so the day it starts moving atoms, someone has to explain
    why - that would mean the range or the pivot has shifted.
    """
    before = dinucleotide_complex.positions[:]

    dinucleotide_complex.aptamer_chain().rotate_in_residue(0, 0, 0.7, reverse=True)

    assert moved_atom_count(before, dinucleotide_complex.positions) == 0


def test_torsion_applies_the_requested_angle(dinucleotide_complex):
    """Asking for a 0.7 rad turn must turn the bond by 0.7 rad.

    Magnitude, not sign: the kernel rotates by -theta, nothing depends on
    that, so asserting the sign would fire on a harmless convention change
    while still missing a wrong amount.
    """
    requested = 0.7
    # HO5'-O5'-C5'-C4' of residue 0: the torsion that j=1 drives.
    before = dihedral(dinucleotide_complex.positions, 0, 1, 2, 5)

    dinucleotide_complex.aptamer_chain().rotate_in_residue(0, 1, requested)

    after = dihedral(dinucleotide_complex.positions, 0, 1, 2, 5)
    applied = angle_difference(after, before)
    assert abs(applied) == pytest.approx(requested, abs=1e-6), (
        f"asked for {requested:.3f} rad, applied {applied:+.4f} rad "
        f"({applied / requested:+.3f}x the requested angle)"
    )

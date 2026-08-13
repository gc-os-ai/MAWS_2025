"""
Turning bonds on a real structure, and checking the molecule survives it.

The rest of the suite works against a stand-in builder that lays atoms out on a
grid, which is enough for spans, counts and bookkeeping but says nothing about
chemistry. Here the strand is built by AmberTools, so the atoms are where a
nucleic acid's atoms go and the list of which are bonded to which is the real
one, read out of what the build produced.

Every turnable bond the search would use is turned, and every bond length has
to come back unchanged. A turn that moves one atom of a bond and not the other
leaves a stretched bond, which still produces an energy and a score —
meaningless ones — and which nothing else in the package would notice.

Everything here is marked ``integration`` and skipped unless you ask for it::

    pytest -m integration
"""

from __future__ import annotations

import shutil

import numpy as np
import pytest
from bond_checks import assert_bonds_survive, bonds_from_topology

from maws.build import LeapBuilder
from maws.forcefield import ForceField
from maws.libraries import dna, rna
from maws.topology import Assembly

pytestmark = pytest.mark.integration

HAS_LEAP = shutil.which("tleap") is not None
"""Whether the structure-building program is installed."""

needs_leap = pytest.mark.skipif(HAS_LEAP is False, reason="tleap is not on PATH")

SEQUENCE = "G A U C"
"""Four different nucleotides, so every residue table in use gets exercised."""


@pytest.fixture(scope="module")
def strand(tmp_path_factory):
    """A four-nucleotide strand built by AmberTools.

    Built once for the whole file, because running ``tleap`` takes seconds and
    nothing here changes the structure — every test starts from these same
    positions and produces new ones. The cache is a fresh directory, so no
    earlier run can supply the answer.
    """
    builder = LeapBuilder(cache_dir=tmp_path_factory.mktemp("cache"))
    return builder.build(
        Assembly().with_aptamer(rna(), SEQUENCE),
        ForceField.for_target("RNA", "protein"),
    )


@pytest.fixture(scope="module")
def bonds(strand):
    """Every bond of that strand, as pairs of atom numbers."""
    return bonds_from_topology(strand.amber.topology)


@needs_leap
class TestTheBondListIsUsable:
    """The checks below say nothing unless the bond list is the real one."""

    def test_a_nucleic_acid_has_about_one_bond_per_atom(self, strand, bonds):
        """A chain molecule has roughly as many bonds as atoms.

        Far fewer would mean the list came back nearly empty and every check
        below is passing on nothing.
        """
        assert 0.8 * strand.n_atoms < len(bonds) < 1.5 * strand.n_atoms

    def test_every_bond_is_between_two_atoms_of_the_strand(self, strand, bonds):
        """The bond list and the positions must be numbered the same way.

        Both come from the same build, so an atom number means the same atom in
        each. If that ever stopped being true the lengths measured below would
        be between unrelated atoms and would look broken at random.
        """
        assert bonds.min() >= 0
        assert bonds.max() < strand.n_atoms

    def test_no_bond_is_longer_than_two_angstrom(self, strand, bonds):
        """A covalent bond between two of these elements is 0.9 to 1.8 Å.

        This is what says the numbering really does line up: if it did not,
        some of these distances would come out at the width of the molecule.
        """
        from bond_checks import bond_lengths

        assert bond_lengths(strand.pose.xyz, bonds).max() < 2.0


@needs_leap
class TestTheTurnableBondsAreBonds:
    """Every bond the residue tables declare has to be a bond in the molecule.

    The tables give atom numbers counted from the start of a residue, worked
    out by reading the order LEaP lays the atoms out in. A number that is off
    by a few still names two atoms, and turning about the line between them
    still produces a structure and an energy — it simply is not a torsion. The
    bond angles around it are wrenched instead, which no other check notices,
    because nothing is stretched: everything that moves stays the same distance
    from the line it turns about.

    This is the check that says the numbers in the tables are the right
    numbers, and it can only be made against a real build.
    """

    @pytest.mark.parametrize(
        ("library", "sequence", "aptamer"),
        [
            (rna(), "G A U C", "RNA"),
            (rna(), "G", "RNA"),
            (dna(), "G A T C", "DNA"),
            (dna(), "G", "DNA"),
        ],
        ids=["rna-four", "rna-one", "dna-four", "dna-one"],
    )
    def test_every_declared_torsion_names_two_bonded_atoms(
        self, tmp_path, library, sequence, aptamer
    ):
        """Across both nucleic acids, and both the four end-of-strand forms.

        A one-nucleotide strand uses the standalone forms of the residues and a
        four-nucleotide one uses the 5' form, the middle form and the 3' form,
        so between them the four rows each residue has in the tables are all
        covered.
        """
        system = LeapBuilder(cache_dir=tmp_path).build(
            Assembly().with_aptamer(library, sequence),
            ForceField.for_target(aptamer, "protein"),
        )
        bonds = {tuple(pair) for pair in bonds_from_topology(system.amber.topology)}
        chain = system.chain("aptamer")
        for index in range(chain.n_residues):
            residue = chain.residue(index)
            for which in range(residue.n_torsions):
                torsion = residue.torsion(which, "3prime")
                pair = (
                    min(torsion.pivot, torsion.bond),
                    max(torsion.pivot, torsion.bond),
                )
                gap = float(
                    np.linalg.norm(
                        system.pose.xyz[torsion.pivot] - system.pose.xyz[torsion.bond]
                    )
                )
                assert pair in bonds, (
                    f"{residue.name} torsion {which} names atoms "
                    f"{pair[0]}-{pair[1]}, which are {gap:.2f} Å apart and not "
                    f"bonded to each other"
                )


@needs_leap
class TestTurningOneBond:
    """One turn at a time, over every bond of every residue."""

    @pytest.mark.parametrize("angle", [0.3, 1.0, 2.5])
    @pytest.mark.parametrize("direction", ["3prime", "5prime"])
    def test_every_turn_keeps_every_bond_length(self, strand, bonds, angle, direction):
        """A turn is a rigid motion of the atoms it moves, so no bond changes.

        Both directions are checked. They turn different parts of the strand
        about the same bond, and a 5' form whose moving atoms are not the exact
        complement of the 3' form's would move one atom of some bond and not
        the other.
        """
        chain = strand.chain("aptamer")
        for residue_index in range(chain.n_residues):
            for torsion in chain.residue(residue_index).torsions(direction):
                turned = strand.pose.rotate(torsion, angle)
                assert_bonds_survive(
                    strand.pose,
                    turned,
                    f"turning residue {residue_index} bond "
                    f"{torsion.pivot}-{torsion.bond} by {angle} rad towards "
                    f"the {direction} end",
                    pairs=bonds,
                )


@needs_leap
class TestTurningTheWholeSetOfBonds:
    """Many turns in a row, the way the search actually shapes a strand."""

    def test_a_sampled_shape_keeps_every_bond(self, strand, bonds, rng):
        """The shapes the search tries must all still be molecules.

        One turn at a time is the check above; this is the same thing done as
        the search does it, several turns applied one after another from one
        draw of random angles. Each turn carries the ones after it, so an error
        in how the moving atoms are worked out compounds here rather than
        cancelling.
        """
        torsions = strand.chain("aptamer").residue(-1).torsions("3prime")
        for _ in range(20):
            angles = rng.uniform(0.0, 2.0 * np.pi, len(torsions))
            shaped = strand.pose.rotate_all(torsions, angles)
            assert_bonds_survive(strand.pose, shaped, "a sampled shape", pairs=bonds)


@needs_leap
class TestEveryBondOfferedDoesSomething:
    """A bond that moves nothing is a wasted degree of freedom.

    The search spends one of its random angles on each bond it is handed, so
    one that cannot change the shape quietly narrows the set of shapes tried.
    Whether a bond can change the shape depends on where the atoms are: if
    every atom it moves happens to sit on the line through the two atoms of the
    bond, turning it does nothing at all. That is a question about real
    positions, which is why it is asked here and not in the unit tier.
    """

    @pytest.mark.parametrize("direction", ["3prime", "5prime"])
    def test_turning_any_offered_bond_changes_the_shape(self, strand, direction):
        """Some atom must end up somewhere else."""
        chain = strand.chain("aptamer")
        for residue_index in range(chain.n_residues):
            for torsion in chain.residue(residue_index).torsions(direction):
                turned = strand.pose.rotate(torsion, 1.0)
                moved = float(np.abs(turned.xyz - strand.pose.xyz).max())
                assert moved > 1e-6, (
                    f"residue {residue_index} bond {torsion.pivot}-"
                    f"{torsion.bond} read towards the {direction} end moves no "
                    f"atom off its own axis"
                )

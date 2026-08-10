"""
Tests that need AmberTools and OpenMM installed.

Everything here is marked ``integration`` and skipped unless you ask for it::

    pytest -m integration

These are the only tests that run ``tleap``, build a real structure, or compute
a real energy. They are slow — seconds to minutes each — which is why the rest
of the suite works against stand-ins instead. What they are for is checking the
joins the stand-ins cannot: that a structure built by AmberTools has the atom
counts the residue tables predict, that its parameters load into OpenMM, and
that a real energy comes back as a finite number.
"""

from __future__ import annotations

import shutil

import numpy as np
import pytest

from maws.build import LeapBuilder
from maws.forcefield import ForceField
from maws.libraries import rna
from maws.topology import Assembly
from maws.values import NucleotideSequence

pytestmark = pytest.mark.integration

HAS_LEAP = shutil.which("tleap") is not None
"""Whether the structure-building program is installed."""

needs_leap = pytest.mark.skipif(HAS_LEAP is False, reason="tleap is not on PATH")


@pytest.fixture
def leap_builder(tmp_path) -> LeapBuilder:
    """A builder writing into a fresh cache, so no earlier run interferes."""
    return LeapBuilder(cache_dir=tmp_path / "cache")


@pytest.fixture
def forcefield() -> ForceField:
    """The usual RNA-against-protein setup."""
    return ForceField.for_target("RNA", "protein")


@needs_leap
class TestRealBuild:
    """Building an actual structure with AmberTools."""

    def test_a_built_strand_has_the_atom_count_the_tables_predict(
        self, leap_builder, forcefield
    ):
        """The residue tables and the builder must agree, or every span is wrong.

        This is the one check the stand-in builder cannot make: it takes its
        atom counts from the same tables, so it can never disagree with them.
        """
        assembly = Assembly().with_aptamer(rna(), "G A U")
        system = leap_builder.build(assembly, forcefield)
        expected = rna().atom_count(NucleotideSequence.parse("G A U").canonical(rna()))
        assert system.n_atoms == expected

    def test_the_built_positions_are_finite(self, leap_builder, forcefield):
        """A structure with a stray infinity would poison every later step."""
        system = leap_builder.build(Assembly().with_aptamer(rna(), "G A"), forcefield)
        assert np.isfinite(system.pose.xyz).all()

    def test_no_two_atoms_share_a_position(self, leap_builder, forcefield):
        """Two atoms in one place would leave a rotation axis undefined."""
        system = leap_builder.build(Assembly().with_aptamer(rna(), "G A"), forcefield)
        distances = np.linalg.norm(
            system.pose.xyz[:, None, :] - system.pose.xyz[None, :, :], axis=-1
        )
        np.fill_diagonal(distances, np.inf)
        assert distances.min() > 0.1

    def test_building_the_same_strand_twice_reuses_the_stored_result(
        self, leap_builder, forcefield
    ):
        """The second build is served from the cache and matches the first."""
        assembly = Assembly().with_aptamer(rna(), "G A")
        first = leap_builder.build(assembly, forcefield)
        second = leap_builder.build(assembly, forcefield)
        np.testing.assert_array_equal(first.pose.xyz, second.pose.xyz)

    def test_a_built_structure_carries_its_amber_files(self, leap_builder, forcefield):
        """Those files are what lets the structure be scored with real physics."""
        system = leap_builder.build(Assembly().with_aptamer(rna(), "G"), forcefield)
        assert system.amber is not None
        assert system.amber.prmtop_path.exists()

    def test_element_symbols_come_from_the_built_structure(
        self, leap_builder, forcefield
    ):
        """A nucleic acid is mostly carbon, oxygen, nitrogen and hydrogen."""
        system = leap_builder.build(Assembly().with_aptamer(rna(), "G"), forcefield)
        assert set(system.elements) <= {"C", "N", "O", "P", "H"}


@needs_leap
class TestRealEnergy:
    """Scoring a real structure with OpenMM."""

    def test_a_real_energy_is_a_finite_number(self, leap_builder, forcefield):
        """If this fails, the parameters did not load and nothing else can work."""
        system = leap_builder.build(Assembly().with_aptamer(rna(), "G A"), forcefield)
        model = system.energy_model(platform="CPU")
        assert np.isfinite(model.evaluate(system.pose))

    def test_scoring_the_same_positions_twice_gives_the_same_number(
        self, leap_builder, forcefield
    ):
        """The scorer is a function of the positions and nothing else.

        Not bit-for-bit: the compute backend sums the same terms in whatever
        order its threads finish in, so the last few digits move between calls.
        """
        system = leap_builder.build(Assembly().with_aptamer(rna(), "G"), forcefield)
        model = system.energy_model(platform="CPU")
        assert model.evaluate(system.pose) == pytest.approx(
            model.evaluate(system.pose), rel=1e-6
        )

    def test_settling_a_structure_does_not_raise_its_energy(
        self, leap_builder, forcefield
    ):
        """Settling rolls downhill, so it can only improve or hold steady."""
        system = leap_builder.build(Assembly().with_aptamer(rna(), "G A"), forcefield)
        model = system.energy_model(platform="CPU")
        before = model.evaluate(system.pose)
        after = model.minimize(system.pose, max_iterations=20)
        assert after.energy <= before + 1e-6

    def test_salt_concentration_changes_the_energy(self, leap_builder):
        """Salt damps the pull between charges, so it must alter the number.

        A nucleic acid backbone carries a charge on every residue, so this is
        one of the larger terms in the score rather than a refinement.
        """
        assembly = Assembly().with_aptamer(rna(), "G A")
        unscreened = ForceField.for_target("RNA", "protein", salt_conc=0.0)
        screened = ForceField.for_target("RNA", "protein", salt_conc=0.5)
        system = leap_builder.build(assembly, unscreened)

        from maws.energy import OpenMMEnergy

        without = OpenMMEnergy.from_prmtop(
            system.amber.prmtop_path, unscreened, platform="CPU"
        ).evaluate(system.pose)
        with_salt = OpenMMEnergy.from_prmtop(
            system.amber.prmtop_path, screened, platform="CPU"
        ).evaluate(system.pose)
        assert without != pytest.approx(with_salt)


@needs_leap
class TestRealGrowth:
    """Growing a real strand, which rebuilds it each time."""

    def test_growing_keeps_every_atom_that_already_existed(
        self, leap_builder, forcefield
    ):
        """The shape found so far is the whole result of the search to that point.

        Rebuilding must therefore put the existing atoms back exactly, not
        merely close by.
        """
        from maws.regrow import grow_chain

        system = leap_builder.build(Assembly().with_aptamer(rna(), "G"), forcefield)
        grown = grow_chain(
            system,
            system.pose,
            role="aptamer",
            token="A",
            direction="3prime",
            builder=leap_builder,
        )
        kept = len(system.chain("aptamer").span)
        np.testing.assert_array_equal(grown.pose.xyz[:kept], system.pose.xyz[:kept])

    def test_the_new_residue_is_joined_at_a_sensible_distance(
        self, leap_builder, forcefield
    ):
        """A new residue placed far from the strand would be nonsense.

        The join is a chemical bond, so the nearest approach between the new
        residue and the rest of the strand should be bond-length scale.
        """
        from maws.regrow import grow_chain

        system = leap_builder.build(Assembly().with_aptamer(rna(), "G"), forcefield)
        grown = grow_chain(
            system,
            system.pose,
            role="aptamer",
            token="A",
            direction="3prime",
            builder=leap_builder,
        )
        kept = len(system.chain("aptamer").span)
        old_part = grown.pose.xyz[:kept]
        new_part = grown.pose.xyz[kept:]
        gaps = np.linalg.norm(old_part[:, None, :] - new_part[None, :, :], axis=-1)
        assert gaps.min() < 2.5

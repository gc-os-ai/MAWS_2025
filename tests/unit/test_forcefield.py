"""
Tests for :mod:`maws.forcefield`.

A force field value is three settings that have to travel together: which
parameters describe the aptamer, which describe the target, and how much salt
the surrounding water holds. The tests here pin down the pairings, the four
nucleotides each pairing allows, and the fact that none of it can be changed
once the value exists.
"""

from __future__ import annotations

import dataclasses

import pytest

from maws.errors import ConfigurationError
from maws.forcefield import ForceField


class TestForTarget:
    """Choosing a pair of parameter sets by saying what the two molecules are."""

    @pytest.mark.parametrize(
        ("aptamer", "molecule", "aptamer_source", "ligand_source"),
        [
            ("RNA", "protein", "leaprc.RNA.OL3", "leaprc.protein.ff19SB"),
            ("RNA", "organic", "leaprc.RNA.OL3", "leaprc.gaff2"),
            ("RNA", "lipid", "leaprc.RNA.OL3", "leaprc.lipid21"),
            ("DNA", "protein", "leaprc.DNA.OL21", "leaprc.protein.ff19SB"),
            ("DNA", "organic", "leaprc.DNA.OL21", "leaprc.gaff2"),
            ("DNA", "lipid", "leaprc.DNA.OL21", "leaprc.lipid21"),
        ],
    )
    def test_each_pairing_names_both_parameter_sets(
        self, aptamer, molecule, aptamer_source, ligand_source
    ):
        """The two molecules are chosen separately and neither disturbs the other."""
        forcefield = ForceField.for_target(aptamer, molecule)
        assert forcefield.aptamer_source == aptamer_source
        assert forcefield.ligand_source == ligand_source

    @pytest.mark.parametrize(
        ("molecule", "expected"),
        [("protein", True), ("organic", False), ("lipid", False)],
    )
    def test_only_a_protein_target_is_ready_to_build_straight_away(
        self, molecule, expected
    ):
        """A protein target needs no parameters worked out for it.

        The protein parameter set already covers proteins, so the
        parameter-fitting programs do not have to run. The other two kinds
        describe families too broad to cover every member, so for those the
        parameters for that particular molecule are worked out first.
        """
        assert ForceField.for_target("RNA", molecule).parameterized is expected

    def test_the_salt_concentration_is_carried_through_unchanged(self):
        """The concentration asked for is the one the value ends up holding."""
        assert ForceField.for_target("RNA", "protein", salt_conc=0.5).salt_conc == 0.5

    def test_the_salt_concentration_defaults_to_roughly_that_of_blood(self):
        """A run that says nothing about salt gets 0.15 mol/L."""
        assert ForceField.for_target("RNA", "protein").salt_conc == 0.15


class TestForTargetRejections:
    """What happens when a run names a molecule nobody has parameters for."""

    def test_an_unknown_nucleic_acid_is_rejected(self):
        """Only RNA and DNA can be grown, so anything else stops the run."""
        with pytest.raises(ConfigurationError, match="aptamer must be one of"):
            ForceField.for_target("PNA", "protein")

    def test_an_unknown_nucleic_acid_lists_the_ones_that_work(self):
        """The message saves a trip to the source to find the spelling."""
        with pytest.raises(ConfigurationError) as caught:
            ForceField.for_target("PNA", "protein")
        assert "'DNA', 'RNA'" in str(caught.value)

    def test_an_unknown_target_kind_is_rejected(self):
        """A target has to be one of the three kinds with a parameter set."""
        with pytest.raises(ConfigurationError, match="molecule must be one of"):
            ForceField.for_target("RNA", "metal")

    def test_an_unknown_target_kind_lists_the_ones_that_work(self):
        """The message names all three choices."""
        with pytest.raises(ConfigurationError) as caught:
            ForceField.for_target("RNA", "metal")
        assert "'lipid', 'organic', 'protein'" in str(caught.value)

    def test_a_negative_salt_concentration_is_rejected(self):
        """There is no such thing as less than no salt."""
        with pytest.raises(ConfigurationError, match="cannot be negative"):
            ForceField.for_target("RNA", "protein", salt_conc=-0.1)

    def test_no_salt_at_all_is_allowed(self):
        """Zero is a real setting: it removes the damping between charges."""
        assert ForceField.for_target("RNA", "protein", salt_conc=0.0).salt_conc == 0.0


class TestAlphabet:
    """The four nucleotides a search may grow the aptamer with."""

    def test_rna_is_built_from_g_a_u_and_c(self):
        """RNA has uracil where DNA has thymine."""
        assert ForceField.for_target("RNA", "protein").alphabet == "GAUC"

    def test_dna_is_built_from_g_a_t_and_c(self):
        """DNA has thymine where RNA has uracil."""
        assert ForceField.for_target("DNA", "protein").alphabet == "GATC"

    @pytest.mark.parametrize(
        ("aptamer_source", "expected"),
        [("leaprc.RNA.OL3", "GAUC"), ("leaprc.DNA.OL21", "GATC")],
    )
    def test_the_alphabet_comes_from_the_parameter_set_name(
        self, aptamer_source, expected
    ):
        """The name of the loaded parameters is what settles the alphabet."""
        forcefield = ForceField(aptamer_source, "leaprc.gaff2")
        assert forcefield.alphabet == expected

    def test_the_alphabet_is_not_a_field_of_its_own(self):
        """Worked out on request, so it cannot disagree with the parameters."""
        names = {field.name for field in dataclasses.fields(ForceField)}
        assert "alphabet" not in names


class TestLeapPreamble:
    """The lines that open the input script for the structure-building program."""

    def test_both_parameter_sets_are_loaded(self):
        """A script has to load every parameter set before naming a molecule."""
        assert ForceField.for_target("RNA", "protein").leap_preamble() == [
            "source leaprc.RNA.OL3",
            "source leaprc.protein.ff19SB",
        ]

    def test_the_aptamer_parameters_come_first(self):
        """The order is fixed, so two identical runs write identical scripts."""
        lines = ForceField.for_target("DNA", "organic").leap_preamble()
        assert lines[0] == "source leaprc.DNA.OL21"
        assert lines[1] == "source leaprc.gaff2"


class TestImmutability:
    """A force field value settled once and used at every step of a run."""

    @pytest.mark.parametrize(
        ("field", "value"),
        [
            ("aptamer_source", "leaprc.DNA.OL21"),
            ("ligand_source", "leaprc.gaff2"),
            ("salt_conc", 1.0),
            ("parameterized", False),
        ],
    )
    def test_no_field_can_be_reassigned(self, field, value):
        """One run cannot build under one setup and then score under another."""
        forcefield = ForceField.for_target("RNA", "protein")
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(forcefield, field, value)

    def test_two_setups_asked_for_the_same_way_are_equal(self):
        """Being a plain value makes a force field safe to use as a cache key."""
        assert ForceField.for_target("RNA", "protein") == ForceField.for_target(
            "RNA", "protein"
        )

"""
Tests for :mod:`maws.values`.

These are the types every other part of MAWS is written in terms of, so the
tests here are deliberately concrete: small numbers, hand-worked answers, and
one behaviour per test.
"""

from __future__ import annotations

import pytest

from maws.errors import ConfigurationError
from maws.libraries import dna, rna
from maws.values import (
    AliasSet,
    AtomRange,
    BackboneTorsion,
    Connection,
    LocalTorsion,
    NucleotideSequence,
    ResidueLibrary,
    ResidueTemplate,
    Torsion,
)


class TestAtomRange:
    """A run of atom indices, from start up to but not including stop."""

    def test_length_counts_the_atoms_in_the_run(self):
        """A range holds stop minus start atoms."""
        assert len(AtomRange(10, 25)) == 15

    def test_a_range_that_starts_where_it_stops_is_empty(self):
        """An empty range is allowed, and holds nothing."""
        assert len(AtomRange(7, 7)) == 0

    def test_membership_excludes_the_stop_index(self):
        """The stop index is one past the end, so it is not in the range."""
        span = AtomRange(3, 6)
        assert 3 in span
        assert 5 in span
        assert 6 not in span

    def test_shifting_moves_both_ends_by_the_same_amount(self):
        """Shifting relocates a run without changing how many atoms it holds."""
        assert AtomRange(10, 25).shifted(100) == AtomRange(110, 125)

    def test_shifting_returns_a_new_range(self):
        """The range shifted is left as it was."""
        original = AtomRange(0, 5)
        original.shifted(10)
        assert original == AtomRange(0, 5)

    def test_as_slice_selects_the_same_atoms(self):
        """The slice covers exactly the indices the range lists."""
        assert list(range(20))[AtomRange(4, 8).as_slice()] == [4, 5, 6, 7]

    def test_a_range_cannot_end_before_it_begins(self):
        """An inverted range is rejected where it is written, not later."""
        with pytest.raises(ConfigurationError, match="cannot end before"):
            AtomRange(9, 4)

    def test_atom_indices_are_never_negative(self):
        """Counting backwards is resolved when a library loads, not here."""
        with pytest.raises(ConfigurationError, match="never negative"):
            AtomRange(-1, 5)


class TestTorsion:
    """A bond that can be turned, in whole-design atom numbers."""

    def test_shifting_moves_every_index_together(self):
        """A torsion keeps its shape when the residue holding it moves."""
        torsion = Torsion(pivot=0, bond=3, moving=AtomRange(3, 34))
        assert torsion.shifted(1000) == Torsion(
            pivot=1000, bond=1003, moving=AtomRange(1003, 1034)
        )


class TestTorsionTemplates:
    """The two kinds of recipe a residue can declare for a turnable bond."""

    def test_a_backbone_bond_moves_everything_to_the_end_of_the_chain(self):
        """A backbone bond swings every atom after it, however far that is."""
        placed = BackboneTorsion(pivot=0, bond=3).placed(100, AtomRange(100, 300))
        assert placed == Torsion(pivot=100, bond=103, moving=AtomRange(103, 300))

    def test_a_local_bond_moves_only_the_atoms_it_names(self):
        """A local bond leaves the rest of the strand alone."""
        placed = LocalTorsion(pivot=8, bond=10, stop=25).placed(
            100, AtomRange(100, 300)
        )
        assert placed == Torsion(pivot=108, bond=110, moving=AtomRange(110, 125))

    def test_the_five_prime_form_turns_the_other_side_of_the_same_bond(self):
        """Turning the 5' side moves everything before the bond instead."""
        reversed_form = BackboneTorsion(pivot=0, bond=3).reversed(
            100, AtomRange(0, 300)
        )
        assert reversed_form.moving == AtomRange(0, 103)

    def test_the_five_prime_form_runs_its_axis_the_other_way(self):
        """The two forms share a bond but point their axis oppositely."""
        recipe = BackboneTorsion(pivot=0, bond=3)
        chain = AtomRange(0, 300)
        forward = recipe.placed(0, chain)
        backward = recipe.reversed(0, chain)
        assert (forward.pivot, forward.bond) == (backward.bond, backward.pivot)

    @pytest.mark.parametrize("direction", ["3prime", "5prime"])
    def test_both_forms_agree_on_which_two_atoms_the_bond_is(self, direction):
        """Whichever side moves, it is the same physical bond being turned."""
        recipe = LocalTorsion(pivot=8, bond=10, stop=25)
        chain = AtomRange(0, 300)
        placed = (
            recipe.placed(0, chain)
            if direction == "3prime"
            else recipe.reversed(0, chain)
        )
        assert {placed.pivot, placed.bond} == {8, 10}


class TestResidueTemplate:
    """A description of one kind of residue: its size and its moving parts."""

    def test_a_residue_must_have_atoms(self):
        """A residue with no atoms cannot be built and is rejected."""
        with pytest.raises(ConfigurationError, match="at least one atom"):
            ResidueTemplate(
                name="X",
                n_atoms=0,
                torsions=(),
                head=Connection(0, -1, 1.6),
                tail=Connection(-2, 0, 1.6),
            )

    def test_a_bond_cannot_name_an_atom_the_residue_does_not_have(self):
        """A typo in a residue table is caught as the template is made."""
        with pytest.raises(ConfigurationError, match="only has atoms 0 to 9"):
            ResidueTemplate(
                name="X",
                n_atoms=10,
                torsions=(BackboneTorsion(pivot=0, bond=40),),
                head=Connection(0, -1, 1.6),
                tail=Connection(-2, 0, 1.6),
            )

    def test_a_local_bond_must_actually_move_something(self):
        """A bond whose moving range is empty would do nothing at all."""
        with pytest.raises(ConfigurationError, match="would move no atoms"):
            ResidueTemplate(
                name="X",
                n_atoms=10,
                torsions=(LocalTorsion(pivot=0, bond=5, stop=5),),
                head=Connection(0, -1, 1.6),
                tail=Connection(-2, 0, 1.6),
            )


class TestResidueLibraryFromTables:
    """Building a library out of the flat per-residue tables."""

    def _tables(self, **overrides):
        """Return a minimal set of valid tables, with any part replaced.

        Parameters
        ----------
        **overrides
            Table name to replacement value.

        Returns
        -------
        dict
            Keyword arguments for :meth:`ResidueLibrary.from_tables`.
        """
        tables = {
            "names": ["G", "A"],
            "lengths": [34, 33],
            "rotations": [("G", 0, 3, None), ("A", 0, 3, -7)],
            "connect": [[[0, -1], [-2, 0], 1.6, 1.6]] * 2,
            "aliases": [["G", "GN", "G5", "G", "G3"]],
        }
        tables.update(overrides)
        return tables

    def test_a_row_with_no_end_index_becomes_a_backbone_bond(self):
        """``None`` in the table means "move the rest of the strand"."""
        library = ResidueLibrary.from_tables(**self._tables())
        assert isinstance(library["G"].torsions[0], BackboneTorsion)

    def test_a_row_with_an_end_index_becomes_a_local_bond(self):
        """A number in the table bounds the atoms that move."""
        library = ResidueLibrary.from_tables(**self._tables())
        assert isinstance(library["A"].torsions[0], LocalTorsion)

    def test_counting_backwards_is_resolved_when_the_table_loads(self):
        """``-7`` in a 33-atom residue means atom 26, worked out once."""
        library = ResidueLibrary.from_tables(**self._tables())
        assert library["A"].torsions[0].stop == 26

    def test_tables_of_different_lengths_are_rejected(self):
        """A missing atom count is caught here, not during a build."""
        with pytest.raises(ConfigurationError, match="do not line up"):
            ResidueLibrary.from_tables(**self._tables(lengths=[34]))

    def test_the_error_says_where_the_shorter_table_stopped(self):
        """The message points at the row to go and look at."""
        with pytest.raises(ConfigurationError, match="stops at index 1"):
            ResidueLibrary.from_tables(**self._tables(lengths=[34]))

    def test_a_rotation_row_for_an_unknown_residue_is_rejected(self):
        """A bond cannot be declared for a residue with no atom count."""
        tables = self._tables(rotations=[("Z", 0, 1, None)])
        with pytest.raises(ConfigurationError, match="does not appear"):
            ResidueLibrary.from_tables(**tables)

    def test_a_malformed_alias_row_is_rejected(self):
        """A short alias row means a typo, not a residue to leave as it is."""
        tables = self._tables(aliases=[["G", "GN", "G5"]])
        with pytest.raises(ConfigurationError, match="alias row must be"):
            ResidueLibrary.from_tables(**tables)


class TestResidueLibrary:
    """Looking residues and written nucleotides up in a library."""

    def test_an_unknown_residue_names_the_ones_that_exist(self):
        """The error message saves a trip to the source tables."""
        with pytest.raises(ConfigurationError, match="It knows: "):
            rna()["nonsense"]

    def test_an_unknown_token_names_the_ones_that_exist(self):
        """Same for a written nucleotide."""
        with pytest.raises(ConfigurationError, match="It knows: "):
            rna().alias("Z")

    def test_atom_count_adds_up_a_run_of_residues(self):
        """A chain's size is the sum of its residues' sizes."""
        assert ResidueLibrary.single("LIG", 40).atom_count(["LIG", "LIG"]) == 80

    def test_a_single_residue_library_declares_no_turnable_bonds(self):
        """A target molecule is moved as a whole and never reshaped."""
        assert ResidueLibrary.single("LIG", 40)["LIG"].n_torsions == 0

    def test_a_residue_with_no_alias_row_stands_for_itself(self):
        """A residue nobody aliased keeps its own name everywhere."""
        library = ResidueLibrary.single("LIG", 40)
        assert library.alias("LIG") == AliasSet("LIG", "LIG", "LIG", "LIG")


class TestNucleotideSequence:
    """A strand written as nucleotides, 5' end first."""

    def test_parsing_splits_on_any_run_of_whitespace(self):
        """Two spaces do not produce an empty nucleotide."""
        assert NucleotideSequence.parse("G  A").tokens == ("G", "A")

    def test_parsing_nothing_gives_an_empty_strand(self):
        """An empty strand is what a search starts from."""
        assert NucleotideSequence.parse("") == NucleotideSequence(())

    def test_an_empty_strand_is_falsey(self):
        """``if sequence`` reads as "if there is a strand at all"."""
        assert not NucleotideSequence(())
        assert NucleotideSequence(("G",))

    def test_growing_returns_a_new_strand(self):
        """The strand grown from is left exactly as it was."""
        original = NucleotideSequence.parse("G A")
        grown = original.appended("U")
        assert str(grown) == "G A U"
        assert str(original) == "G A"

    @pytest.mark.parametrize(
        ("direction", "expected"),
        [("3prime", "G A"), ("5prime", "A G")],
    )
    def test_growing_at_either_end(self, direction, expected):
        """Which end grows decides where the new nucleotide is written."""
        grown = NucleotideSequence.parse("G").grown("A", direction)
        assert str(grown) == expected

    def test_a_lone_nucleotide_gets_its_standalone_name(self):
        """A one-residue strand is chemically its own special case."""
        assert NucleotideSequence.parse("G").canonical(rna()) == ("GN",)

    def test_the_ends_of_a_strand_get_end_specific_names(self):
        """The 5' and 3' residues differ chemically from the middle ones."""
        canonical = NucleotideSequence.parse("G A U").canonical(rna())
        assert canonical == ("G5", "A", "U3")

    def test_dna_is_written_with_tokens_that_are_not_residue_names(self):
        """``G`` is a DNA nucleotide but ``DG5`` is what the builder wants."""
        assert NucleotideSequence.parse("G A T").canonical(dna()) == (
            "DG5",
            "DA",
            "DT3",
        )

    def test_an_unknown_nucleotide_is_reported_by_name(self):
        """A typo in a sequence says which letter was not understood."""
        with pytest.raises(ConfigurationError, match="token 'Z'"):
            NucleotideSequence.parse("G Z").canonical(rna())

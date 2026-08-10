"""
Tests for :mod:`maws.topology`.

Two things are pinned down here. First, that describing a design never changes
a description that already exists: the search tries several ways of extending
the same strand and needs the one it started from intact after each attempt.
Second, that a description which cannot answer a question says so plainly
rather than guessing.

The numbers are the real RNA residue sizes, so they can be checked against the
residue tables: G5 has 32 atoms, A3 has 34, and a guanine standing alone as a
whole strand, GN, has 33.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from maws.errors import ConfigurationError
from maws.forcefield import ForceField
from maws.libraries import rna
from maws.topology import (
    Assembly,
    BuiltSystem,
    PdbChain,
    ResidueChain,
    compute_spans,
    default_residue_name,
)
from maws.values import NucleotideSequence, ResidueLibrary


def _pdb_chain(path, role="ligand"):
    """Return a chain that is still only a file on disk.

    Parameters
    ----------
    path : pathlib.Path
        Where the PDB file would be. It is never read.
    role : str, default="ligand"
        The chain's name within its assembly.

    Returns
    -------
    maws.topology.PdbChain
        A chain description with no atom count.
    """
    return PdbChain(role=role, path=path, residue_name="LIG", parameterized=False)


class TestResidueChain:
    """A chain whose nucleotides are known, so its size is known too."""

    def test_atom_count_adds_up_the_residues_of_the_sequence(self):
        """A two-nucleotide strand is 32 atoms of G5 plus 34 of A3."""
        chain = ResidueChain("aptamer", rna(), NucleotideSequence.parse("G A"))
        assert chain.n_atoms == 66

    def test_a_one_nucleotide_strand_uses_the_standalone_residue(self):
        """A lone guanine is a whole strand in itself and has 33 atoms."""
        chain = ResidueChain("aptamer", rna(), NucleotideSequence.parse("G"))
        assert chain.n_atoms == 33

    def test_the_two_ends_of_a_strand_get_end_specific_names(self):
        """The residues at 5' and 3' carry groups a middle residue does not."""
        chain = ResidueChain("aptamer", rna(), NucleotideSequence.parse("G A U"))
        assert chain.canonical == ("G5", "A", "U3")

    def test_replacing_the_sequence_gives_a_chain_of_the_new_size(self):
        """Growing a strand by one nucleotide grows its atom count with it."""
        chain = ResidueChain("aptamer", rna(), NucleotideSequence.parse("G"))
        grown = chain.with_sequence(NucleotideSequence.parse("G A"))
        assert grown.n_atoms == 66

    def test_replacing_the_sequence_leaves_the_original_chain_alone(self):
        """The chain replaced from still describes the strand it always did."""
        chain = ResidueChain("aptamer", rna(), NucleotideSequence.parse("G"))
        chain.with_sequence(NucleotideSequence.parse("G A"))
        assert chain.n_atoms == 33

    def test_replacing_the_sequence_keeps_everything_else(self):
        """Only the nucleotides change; the name and the residue tables stay."""
        chain = ResidueChain("aptamer", rna(), NucleotideSequence.parse("G"))
        grown = chain.with_sequence(NucleotideSequence.parse("G A"))
        assert (grown.role, grown.library) == (chain.role, chain.library)


class TestAssemblyInspection:
    """Asking a description what it contains, without building anything."""

    def test_roles_are_listed_in_the_order_the_chains_were_added(self):
        """The first chain added takes the lowest atom indices of the design."""
        assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(40)
        assert assembly.roles == ("aptamer", "ligand")

    def test_atom_count_is_the_total_across_every_chain(self):
        """A 66-atom strand beside a 40-atom target comes to 106 atoms."""
        assembly = Assembly().with_aptamer(rna(), "G A").with_ligand_stub(40)
        assert assembly.n_atoms() == 106

    def test_an_empty_assembly_has_no_atoms(self):
        """A description with no chains in it is what a run starts from."""
        assert Assembly().n_atoms() == 0

    def test_a_chain_is_found_by_the_name_it_was_given(self):
        """Chains are looked up by name, so a third chain changes nothing."""
        assembly = Assembly().with_aptamer(rna(), "G A").with_ligand_stub(40)
        assert assembly.chain("aptamer").n_atoms == 66

    def test_an_unknown_chain_name_is_rejected(self):
        """Asking for a chain that is not there is a mistake worth reporting."""
        assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(40)
        with pytest.raises(ConfigurationError, match="no chain named 'protein'"):
            assembly.chain("protein")

    def test_an_unknown_chain_name_lists_the_names_that_exist(self):
        """The message shows what could have been asked for instead."""
        assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(40)
        with pytest.raises(ConfigurationError, match="Chains: aptamer, ligand"):
            assembly.chain("protein")


class TestAssemblyAtomCountWithAPdbChain:
    """A target still on disk has no atom count, and nothing invents one."""

    def test_the_total_is_refused_while_a_file_is_unresolved(self, tmp_path):
        """Nobody has counted the atoms in that file, so there is no total."""
        assembly = Assembly().with_chain(_pdb_chain(tmp_path / "target.pdb"))
        with pytest.raises(ConfigurationError, match="atom count is not known yet"):
            assembly.n_atoms()

    def test_the_message_names_the_chain_and_its_file(self, tmp_path):
        """Which of several chains is holding things up, and where it lives."""
        path = tmp_path / "target.pdb"
        assembly = Assembly().with_aptamer(rna(), "G").with_chain(_pdb_chain(path))
        with pytest.raises(ConfigurationError) as caught:
            assembly.n_atoms()
        assert re.search(
            rf"chain 'ligand' comes from {re.escape(str(path))}", str(caught.value)
        )


class TestAssemblyIsNeverChangedInPlace:
    """Every way of describing more still leaves the earlier description whole.

    The search depends on this: it builds one candidate assembly per nucleotide
    it might add, all from the same starting point, and that starting point has
    to survive every one of them.
    """

    def test_adding_a_chain_returns_a_new_assembly(self):
        """The assembly added to keeps the chains it had."""
        start = Assembly().with_aptamer(rna(), "G")
        grown = start.with_ligand_stub(40)
        assert start.roles == ("aptamer",)
        assert grown.roles == ("aptamer", "ligand")

    def test_adding_an_aptamer_returns_a_new_assembly(self):
        """An empty description stays empty after a strand is described."""
        empty = Assembly()
        one = empty.with_aptamer(rna(), "G")
        assert empty.roles == ()
        assert one.roles == ("aptamer",)

    def test_adding_a_stand_in_target_returns_a_new_assembly(self):
        """A stand-in target is added to a copy, never to the original."""
        empty = Assembly()
        with_target = empty.with_ligand_stub(40)
        assert len(empty) == 0
        assert len(with_target) == 1

    def test_adding_a_target_from_a_file_returns_a_new_assembly(self, tmp_path):
        """Naming a PDB file reads nothing and changes nothing."""
        empty = Assembly()
        with_target = empty.with_ligand(
            tmp_path / "target.pdb", ForceField.for_target("RNA", "organic")
        )
        assert empty.roles == ()
        assert with_target.roles == ("ligand",)

    def test_replacing_a_sequence_returns_a_new_assembly(self):
        """Growing the strand leaves the shorter description ready to reuse."""
        start = Assembly().with_aptamer(rna(), "G").with_ligand_stub(40)
        grown = start.with_sequence("aptamer", "G A")
        assert start.n_atoms() == 73
        assert grown.n_atoms() == 106

    @pytest.mark.parametrize(
        "extend",
        [
            lambda a: a.with_chain(
                ResidueChain("extra", rna(), NucleotideSequence(()))
            ),
            lambda a: a.with_aptamer(rna(), "G A U", role="extra"),
            lambda a: a.with_ligand_stub(5, role="extra"),
            lambda a: a.with_sequence("aptamer", "G A U"),
        ],
    )
    def test_every_way_of_describing_more_gives_back_a_different_object(self, extend):
        """A new description is a new object, so the two cannot alias."""
        start = Assembly().with_aptamer(rna(), "G").with_ligand_stub(40)
        assert extend(start) is not start


class TestAssemblyWithChain:
    """Adding a chain to a description."""

    def test_two_chains_cannot_share_a_name(self):
        """Names are how chains are found, so a duplicate would hide one."""
        start = Assembly().with_aptamer(rna(), "G")
        with pytest.raises(ConfigurationError, match="already has a chain"):
            start.with_aptamer(rna(), "G A")

    def test_the_message_names_the_role_that_clashed(self):
        """Which name was taken, so the fix is one edit."""
        start = Assembly().with_ligand_stub(40)
        with pytest.raises(ConfigurationError, match="named 'ligand'"):
            start.with_ligand_stub(20)


class TestAssemblyWithSequence:
    """Replacing the nucleotides of one chain."""

    def test_only_the_named_chain_changes(self):
        """A target beside the strand is the same target afterwards."""
        start = Assembly().with_aptamer(rna(), "G").with_ligand_stub(40)
        grown = start.with_sequence("aptamer", "G A")
        assert grown.chain("ligand") == start.chain("ligand")

    def test_the_named_chain_gets_the_new_nucleotides(self):
        """The strand asked for is the strand the copy describes."""
        start = Assembly().with_aptamer(rna(), "G").with_ligand_stub(40)
        grown = start.with_sequence("aptamer", "G A")
        assert grown.chain("aptamer").canonical == ("G5", "A3")

    def test_the_chain_order_is_kept(self):
        """Atom indices depend on the order, which editing must not disturb."""
        start = Assembly().with_aptamer(rna(), "G").with_ligand_stub(40)
        assert start.with_sequence("aptamer", "G A").roles == ("aptamer", "ligand")

    def test_a_chain_read_from_a_file_has_no_sequence_to_edit(self, tmp_path):
        """A measured molecule is whatever the file says; it is not written."""
        assembly = Assembly().with_chain(_pdb_chain(tmp_path / "target.pdb"))
        with pytest.raises(ConfigurationError, match="no editable sequence"):
            assembly.with_sequence("ligand", "G A")

    def test_an_unknown_chain_name_is_rejected(self):
        """Editing a chain that is not there is caught where it is written."""
        assembly = Assembly().with_aptamer(rna(), "G")
        with pytest.raises(ConfigurationError, match="no chain named 'ligand'"):
            assembly.with_sequence("ligand", "G A")


class TestDefaultResidueName:
    """The name a target molecule is known by while it is being built.

    The name decides what the molecule's parameter file is called. Taking it
    from the file's contents is what stops two different targets writing over
    each other's parameters.
    """

    def test_two_files_holding_different_atoms_get_different_names(self, tmp_path):
        """Two targets in one working directory cannot collide."""
        first = tmp_path / "first.pdb"
        second = tmp_path / "second.pdb"
        first.write_text("ATOM      1  C1  LIG A   1       0.000   0.000   0.000\n")
        second.write_text("ATOM      1  N1  LIG A   1       1.000   0.000   0.000\n")
        assert default_residue_name(first) != default_residue_name(second)

    def test_two_files_holding_the_same_atoms_get_the_same_name(self, tmp_path):
        """One molecule keeps one name, so its parameter file is found again."""
        contents = "ATOM      1  C1  LIG A   1       0.000   0.000   0.000\n"
        first = tmp_path / "here.pdb"
        second = tmp_path / "elsewhere" / "copied.pdb"
        second.parent.mkdir()
        first.write_text(contents)
        second.write_text(contents)
        assert default_residue_name(first) == default_residue_name(second)

    def test_the_same_file_gives_the_same_name_every_time(self, tmp_path):
        """Nothing about the answer depends on when it was asked for."""
        path = tmp_path / "target.pdb"
        path.write_text("ATOM      1  C1  LIG A   1       0.000   0.000   0.000\n")
        assert default_residue_name(path) == default_residue_name(path)

    def test_editing_a_file_changes_the_name_it_yields(self, tmp_path):
        """A target that has been changed is a different molecule to build."""
        path = tmp_path / "target.pdb"
        path.write_text("ATOM      1  C1  LIG A   1       0.000   0.000   0.000\n")
        before = default_residue_name(path)
        path.write_text("ATOM      1  C1  LIG A   1       9.000   0.000   0.000\n")
        assert default_residue_name(path) != before

    def test_a_file_that_is_not_there_yet_still_gets_a_name(self, tmp_path):
        """A run can be described before the target file has been written."""
        assert default_residue_name(tmp_path / "not_written_yet.pdb")

    def test_a_name_is_a_letter_and_six_hex_digits(self, tmp_path):
        """Short enough for the residue-name column of a PDB-style file."""
        name = default_residue_name(tmp_path / "target.pdb")
        assert re.fullmatch(r"L[0-9A-F]{6}", name)


class TestComputeSpans:
    """Working out which stretch of the atom array each chain occupies."""

    def _chains(self):
        """Return a two-chain layout of 33 and then 40 atoms.

        Returns
        -------
        list of maws.topology.ResidueChain
            A lone guanine followed by a 40-atom stand-in target.
        """
        return [
            ResidueChain("aptamer", rna(), NucleotideSequence.parse("G")),
            ResidueChain(
                "ligand",
                ResidueLibrary.single("LIG", 40),
                NucleotideSequence(("LIG",)),
            ),
        ]

    def test_the_first_chain_starts_at_the_beginning_of_the_array(self):
        """Atom numbering starts at zero, not at the first chain's own offset."""
        assert compute_spans(self._chains())["aptamer"].start == 0

    def test_each_chain_starts_where_the_one_before_it_stopped(self):
        """Chains are laid end to end with no gap and no overlap."""
        spans = compute_spans(self._chains())
        assert spans["ligand"].start == spans["aptamer"].stop

    def test_each_span_is_as_long_as_its_chain(self):
        """A 33-atom strand covers 33 rows of the array."""
        spans = compute_spans(self._chains())
        assert len(spans["aptamer"]) == 33
        assert len(spans["ligand"]) == 40

    def test_the_spans_together_cover_every_atom(self):
        """No atom belongs to no chain, and none belongs to two."""
        chains = self._chains()
        spans = compute_spans(chains)
        assert sum(len(span) for span in spans.values()) == sum(
            chain.n_atoms for chain in chains
        )

    def test_no_chains_means_no_spans(self):
        """An empty design has nothing to lay out."""
        assert compute_spans([]) == {}


def _system(assembly, forcefield, n_atoms=None, elements=None, masses=None):
    """Build a system by hand, so a table can be given the wrong length.

    Parameters
    ----------
    assembly : maws.topology.Assembly
        The description to build.
    forcefield : maws.forcefield.ForceField
        The physics the structure is recorded as having been built under.
    n_atoms : int, optional
        How many rows of coordinates to supply. Defaults to what `assembly`
        says it needs.
    elements : sequence of str, optional
        Element symbols. Defaults to one carbon per atom.
    masses : sequence of float, optional
        Masses in daltons. Defaults to carbon's, one per atom.

    Returns
    -------
    maws.topology.BuiltSystem
        The system, if the arguments agree with each other.
    """
    count = assembly.n_atoms() if n_atoms is None else n_atoms
    return BuiltSystem(
        assembly,
        forcefield,
        np.zeros((count, 3)),
        elements=["C"] * count if elements is None else elements,
        masses=[12.011] * count if masses is None else masses,
    )


class TestBuiltSystemConstruction:
    """What a built system insists on before it will exist at all."""

    def test_a_chain_still_on_disk_cannot_be_part_of_a_built_system(
        self, tmp_path, forcefield
    ):
        """Building means every chain's atoms are known; one is not."""
        assembly = Assembly().with_chain(_pdb_chain(tmp_path / "target.pdb"))
        with pytest.raises(ConfigurationError, match="not parameterized yet"):
            _system(assembly, forcefield, n_atoms=10)

    def test_one_element_symbol_short_is_rejected(self, forcefield):
        """A symbol per atom, or nothing reading them can be trusted."""
        assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(10)
        with pytest.raises(ConfigurationError, match="elements has 5 entries"):
            _system(assembly, forcefield, elements=["C"] * 5)

    def test_one_mass_short_is_rejected(self, forcefield):
        """A mass per atom, or a centre of mass comes out somewhere wrong."""
        assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(10)
        with pytest.raises(ConfigurationError, match="masses has 2 entries"):
            _system(assembly, forcefield, masses=[12.011, 12.011])

    def test_the_message_says_how_many_atoms_were_expected(self, forcefield):
        """33 atoms of guanine and 10 of target come to 43."""
        assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(10)
        with pytest.raises(ConfigurationError, match="the system has 43 atoms"):
            _system(assembly, forcefield, elements=["C"] * 5)

    def test_chains_and_coordinates_must_agree_on_the_atom_count(self, forcefield):
        """Coordinates for a design other than this one are caught here."""
        assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(10)
        with pytest.raises(ConfigurationError, match="assembly and the build disagree"):
            _system(assembly, forcefield, n_atoms=50)


class TestBuiltSystemInspection:
    """Reading a built system's layout back out."""

    def test_the_atom_count_is_the_strand_plus_the_target(self, two_residue_system):
        """66 atoms of strand beside a 20-atom stand-in come to 86."""
        assert two_residue_system.n_atoms == 86

    def test_a_span_is_given_for_every_chain(self, two_residue_system):
        """Each chain knows which stretch of the array is its own."""
        assert set(two_residue_system.spans) == {"aptamer", "ligand"}

    def test_the_spans_add_up_to_the_atom_count(self, two_residue_system):
        """Every atom belongs to exactly one chain."""
        total = sum(len(span) for span in two_residue_system.spans.values())
        assert total == two_residue_system.n_atoms

    def test_the_spans_are_laid_out_in_the_order_the_chains_were_added(
        self, two_residue_system
    ):
        """The aptamer was described first, so it takes the lowest indices."""
        spans = two_residue_system.spans
        assert (spans["aptamer"].start, spans["aptamer"].stop) == (0, 66)
        assert (spans["ligand"].start, spans["ligand"].stop) == (66, 86)

    def test_there_is_one_element_symbol_per_atom(self, two_residue_system):
        """Anything judging how much room an atom takes up has a symbol to read."""
        assert len(two_residue_system.elements) == two_residue_system.n_atoms

    def test_there_is_one_mass_per_atom(self, two_residue_system):
        """A centre of mass needs every atom's mass, with none missing."""
        assert len(two_residue_system.masses) == two_residue_system.n_atoms

    def test_a_chain_window_covers_that_chain_and_nothing_else(
        self, two_residue_system
    ):
        """A window onto a chain names the same atoms its span does."""
        view = two_residue_system.chain("ligand")
        assert view.span == two_residue_system.spans["ligand"]

    def test_an_empty_strand_leaves_only_the_target(self, empty_system):
        """A search starts with no nucleotides at all, which is allowed."""
        assert empty_system.n_atoms == 20
        assert len(empty_system.spans["aptamer"]) == 0

    def test_an_unknown_chain_name_is_rejected(self, two_residue_system):
        """A misspelt name is a mistake, not an empty result."""
        with pytest.raises(ConfigurationError, match="no chain named 'protein'"):
            two_residue_system.chain("protein")

    def test_an_unknown_chain_name_lists_the_names_that_exist(self, two_residue_system):
        """The message shows what could have been asked for instead."""
        with pytest.raises(ConfigurationError, match="Chains: aptamer, ligand"):
            two_residue_system.chain("protein")


class TestBuiltSystemEnergyModel:
    """Asking a system for something that can calculate its energy."""

    def test_a_system_of_stand_ins_has_no_energy_model(self, two_residue_system):
        """Atoms on a grid are not a molecule, so no real physics applies."""
        with pytest.raises(ConfigurationError, match="no AMBER topology"):
            two_residue_system.energy_model()

    def test_the_message_points_at_the_stand_in_scorer(self, two_residue_system):
        """There is something that will score a stand-in, and it is named."""
        with pytest.raises(ConfigurationError, match="StubEnergy"):
            two_residue_system.energy_model()

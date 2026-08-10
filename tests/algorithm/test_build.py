"""
Tests for :mod:`maws.build`.

Everything here runs with nothing installed. :class:`~maws.build.FakeBuilder`
puts atoms on a grid, so its whole behaviour can be checked by arithmetic.

:class:`~maws.build.LeapBuilder` needs AmberTools to build anything, so only the
two parts of it that do not are exercised: the name it stores a result under and
the script it hands to ``tleap``. Both are private methods, and they are called
directly here because reaching them through :meth:`~maws.build.LeapBuilder.build`
would mean running the programs this test suite deliberately does without.
"""

from __future__ import annotations

import numpy as np
import pytest

from maws.build import (
    CHAIN_SEPARATION,
    GRID_SPACING,
    FakeBuilder,
    LeapBuilder,
    build,
)
from maws.errors import ConfigurationError
from maws.forcefield import ForceField
from maws.libraries import rna
from maws.topology import (
    DEFAULT_STUB_ELEMENT,
    DEFAULT_STUB_MASS,
    Assembly,
    PdbChain,
    ResidueChain,
)
from maws.values import AtomRange, NucleotideSequence

GUANINE_ATOMS = 33
"""How many atoms a lone guanine has, per the RNA residue table."""

STUB_ATOMS = 20
"""How many atoms the stand-in target used in this file has."""


def two_chain_assembly() -> Assembly:
    """Return one guanine beside a twenty-atom stand-in target.

    Returns
    -------
    maws.topology.Assembly
        A description with an ``aptamer`` chain and a ``ligand`` chain.
    """
    return Assembly().with_aptamer(rna(), "G").with_ligand_stub(STUB_ATOMS)


def write_pdb(path, n_atoms: int, *, record: str = "ATOM") -> None:
    """Write a PDB file holding `n_atoms` atom records and nothing else.

    Parameters
    ----------
    path : pathlib.Path
        Where to write the file.
    n_atoms : int
        How many atom records to write.
    record : {"ATOM", "HETATM"}, default="ATOM"
        Which record name to give every line.
    """
    lines = ["REMARK a stand-in for a measured structure"]
    lines += [
        f"{record:<6}{index + 1:>5}  C   LIG A   1"
        f"{0.0:>12.3f}{0.0:>8.3f}{0.0:>8.3f}  1.00  0.00           C"
        for index in range(n_atoms)
    ]
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


class TestFakeBuilderLayout:
    """Where a grid builder puts the atoms of a two-chain design."""

    def test_the_structure_holds_every_chain_s_atoms(self, builder, forcefield):
        """A guanine of 33 atoms beside a 20-atom target comes to 53."""
        system = builder.build(two_chain_assembly(), forcefield)
        assert system.n_atoms == GUANINE_ATOMS + STUB_ATOMS

    def test_each_chain_owns_a_run_of_the_coordinate_array(self, builder, forcefield):
        """The first chain takes the lowest indices and the second follows it."""
        system = builder.build(two_chain_assembly(), forcefield)
        assert system.chain("aptamer").span == AtomRange(0, GUANINE_ATOMS)
        assert system.chain("ligand").span == AtomRange(
            GUANINE_ATOMS, GUANINE_ATOMS + STUB_ATOMS
        )

    def test_no_two_atoms_share_a_position(self, builder, forcefield):
        """A turnable bond is a line through two atoms, so a repeated position
        would leave that line undefined."""
        system = builder.build(two_chain_assembly(), forcefield)
        distinct = np.unique(system.pose.xyz, axis=0)
        assert len(distinct) == system.n_atoms

    def test_the_two_chains_are_placed_clear_of_each_other(self, builder, forcefield):
        """No atom of one chain sits on top of an atom of the other."""
        system = builder.build(two_chain_assembly(), forcefield)
        first = system.pose.atoms(system.chain("aptamer").span)
        second = system.pose.atoms(system.chain("ligand").span)
        gaps = np.linalg.norm(first[:, None, :] - second[None, :, :], axis=-1)
        assert gaps.min() > 0.0

    def test_the_gap_between_the_chains_is_about_the_separation_asked_for(
        self, builder, forcefield
    ):
        """The two grids start 25 Å apart, and the first is only a few ångström
        wide, so the nearest pair of atoms is a little under that."""
        system = builder.build(two_chain_assembly(), forcefield)
        first = system.pose.atoms(system.chain("aptamer").span)
        second = system.pose.atoms(system.chain("ligand").span)
        gaps = np.linalg.norm(first[:, None, :] - second[None, :, :], axis=-1)
        assert gaps.min() == pytest.approx(CHAIN_SEPARATION, rel=0.25)
        assert gaps.min() < CHAIN_SEPARATION

    def test_wider_spacing_spreads_a_chain_out_proportionally(self, forcefield):
        """Doubling the spacing doubles every coordinate of the first chain,
        which starts at the origin."""
        assembly = Assembly().with_ligand_stub(STUB_ATOMS)
        narrow = FakeBuilder(spacing=GRID_SPACING).build(assembly, forcefield)
        wide = FakeBuilder(spacing=2 * GRID_SPACING).build(assembly, forcefield)
        assert wide.pose.xyz == pytest.approx(2 * narrow.pose.xyz)

    def test_the_separation_decides_where_the_second_chain_starts(self, forcefield):
        """The second chain's grid begins exactly one separation along x."""
        system = FakeBuilder(separation=40.0).build(two_chain_assembly(), forcefield)
        second = system.pose.atoms(system.chain("ligand").span)
        assert second[:, 0].min() == pytest.approx(40.0)

    def test_building_the_same_design_twice_gives_the_same_positions(self, forcefield):
        """There is no randomness in a grid, so two builds agree exactly."""
        once = FakeBuilder().build(two_chain_assembly(), forcefield)
        twice = FakeBuilder().build(two_chain_assembly(), forcefield)
        assert once.pose.xyz == pytest.approx(twice.pose.xyz)

    def test_a_design_with_no_chains_is_rejected(self, builder, forcefield):
        """There is nothing to lay out, so this is caught rather than producing
        an empty structure."""
        with pytest.raises(ConfigurationError, match="no chains"):
            builder.build(Assembly(), forcefield)

    def test_every_atom_is_given_a_symbol_and_a_mass(self, builder, forcefield):
        """Samplers read one element symbol per atom and one mass per atom."""
        system = builder.build(two_chain_assembly(), forcefield)
        assert system.elements == (DEFAULT_STUB_ELEMENT,) * system.n_atoms
        assert system.masses == (DEFAULT_STUB_MASS,) * system.n_atoms


class TestFakeBuilderPrepare:
    """Resolving a chain that is still a file on disk."""

    @pytest.mark.parametrize("record", ["ATOM", "HETATM"])
    def test_both_kinds_of_atom_record_are_counted(
        self, tmp_path, builder, forcefield, record
    ):
        """A measured structure may list its atoms under either record name."""
        path = tmp_path / "target.pdb"
        write_pdb(path, 7, record=record)
        chain = PdbChain(
            role="ligand", path=path, residue_name="LIG", parameterized=True
        )
        assert builder.prepare(chain, forcefield).n_atoms == 7

    def test_lines_that_are_not_atom_records_are_ignored(
        self, tmp_path, builder, forcefield
    ):
        """The remark and end lines of the file are not atoms."""
        path = tmp_path / "target.pdb"
        write_pdb(path, 3)
        chain = PdbChain(
            role="ligand", path=path, residue_name="LIG", parameterized=True
        )
        assert builder.prepare(chain, forcefield).n_atoms == 3

    def test_a_file_with_no_atoms_is_rejected(self, tmp_path, builder, forcefield):
        """A file naming no atoms describes nothing that can be built."""
        path = tmp_path / "empty.pdb"
        path.write_text("REMARK nothing here\nEND\n")
        chain = PdbChain(
            role="ligand", path=path, residue_name="LIG", parameterized=True
        )
        with pytest.raises(ConfigurationError, match="holds no ATOM or HETATM"):
            builder.prepare(chain, forcefield)

    def test_a_file_that_cannot_be_read_is_reported_by_name(
        self, tmp_path, builder, forcefield
    ):
        """The message names the path, which is normally where the mistake is."""
        chain = PdbChain(
            role="ligand",
            path=tmp_path / "not-there.pdb",
            residue_name="LIG",
            parameterized=True,
        )
        with pytest.raises(ConfigurationError, match="cannot read"):
            builder.prepare(chain, forcefield)


class TestBuildFunction:
    """Choosing which builder does the work."""

    def test_the_builder_given_is_the_one_used(self, forcefield):
        """Passing a grid builder builds a grid, with no AmberTools involved."""
        system = build(two_chain_assembly(), forcefield, builder=FakeBuilder())
        assert system.n_atoms == GUANINE_ATOMS + STUB_ATOMS
        assert system.amber is None

    def test_the_force_field_is_carried_onto_the_result(self, forcefield):
        """Whatever scores the structure later gets the physics it was built
        under."""
        system = build(two_chain_assembly(), forcefield, builder=FakeBuilder())
        assert system.forcefield is forcefield


def residue_chains() -> tuple[ResidueChain, ...]:
    """Return an aptamer chain and a stand-in target chain.

    Returns
    -------
    tuple of maws.topology.ResidueChain
        Two chains that both know their residues.
    """
    return tuple(two_chain_assembly().chains)


class TestLeapCacheKey:
    """The name a built structure is stored and looked up under."""

    def test_the_same_design_always_gets_the_same_name(self, forcefield):
        """A structure already on disk is found again instead of being rebuilt."""
        builder = LeapBuilder()
        chains = residue_chains()
        assert builder._cache_key(chains, forcefield) == builder._cache_key(
            chains, forcefield
        )

    def test_a_different_sequence_gets_a_different_name(self, forcefield):
        """Growing the strand by one nucleotide must not reuse the old
        structure."""
        builder = LeapBuilder()
        one = two_chain_assembly()
        two = one.with_sequence("aptamer", "G A")
        assert builder._cache_key(tuple(one.chains), forcefield) != builder._cache_key(
            tuple(two.chains), forcefield
        )

    def test_a_different_force_field_gets_a_different_name(self, forcefield):
        """The same nucleotides built under different physics are different
        structures."""
        builder = LeapBuilder()
        chains = residue_chains()
        other = ForceField.for_target("DNA", "protein")
        assert builder._cache_key(chains, forcefield) != builder._cache_key(
            chains, other
        )

    def test_the_name_is_forty_hexadecimal_characters(self, forcefield):
        """It becomes a file name, so it has to be short and safe to write."""
        key = LeapBuilder()._cache_key(residue_chains(), forcefield)
        assert len(key) == 40
        assert set(key) <= set("0123456789abcdef")


class TestLeapScript:
    """The input script that says what to build."""

    def script(self, chains, forcefield, tmp_path) -> str:
        """Return the script for `chains`, written to files under `tmp_path`.

        Parameters
        ----------
        chains : sequence of maws.topology.ResidueChain
            The chains to build.
        forcefield : maws.forcefield.ForceField
            Which parameters to build with.
        tmp_path : pathlib.Path
            Where the two output files would be written.

        Returns
        -------
        str
            The script.
        """
        return LeapBuilder()._leap_script(
            chains, forcefield, tmp_path / "out.prmtop", tmp_path / "out.inpcrd"
        )

    def test_the_parameters_are_loaded_before_anything_else(self, forcefield, tmp_path):
        """One line per force field, and they come first."""
        lines = self.script(residue_chains(), forcefield, tmp_path).splitlines()
        assert lines[:2] == [
            f"source {forcefield.aptamer_source}",
            f"source {forcefield.ligand_source}",
        ]

    def test_each_chain_with_residues_gets_a_sequence_line(self, forcefield, tmp_path):
        """Two chains to build means two molecules to name."""
        lines = self.script(residue_chains(), forcefield, tmp_path).splitlines()
        assert sum(1 for line in lines if " = sequence {" in line) == 2

    def test_a_chain_with_no_residues_contributes_nothing_to_build(
        self, forcefield, tmp_path
    ):
        """A search starts from an empty strand, which has no molecule to make."""
        chains = (
            ResidueChain("aptamer", rna(), NucleotideSequence(())),
            residue_chains()[1],
        )
        lines = self.script(chains, forcefield, tmp_path).splitlines()
        assert sum(1 for line in lines if " = sequence {" in line) == 1

    def test_the_chains_are_joined_into_one_molecule(self, forcefield, tmp_path):
        """Both chains have to be in one structure for their contact to be
        scored."""
        script = self.script(residue_chains(), forcefield, tmp_path)
        assert "UNION = combine {CHAIN0 CHAIN1}" in script.splitlines()

    def test_the_result_is_written_out(self, forcefield, tmp_path):
        """Without this line the run would finish having saved nothing."""
        script = self.script(residue_chains(), forcefield, tmp_path)
        expected = (
            f"saveamberparm UNION {tmp_path / 'out.prmtop'} {tmp_path / 'out.inpcrd'}"
        )
        assert expected in script.splitlines()

    def test_the_script_ends_by_quitting(self, forcefield, tmp_path):
        """Without a closing instruction the program waits for input forever."""
        script = self.script(residue_chains(), forcefield, tmp_path)
        assert script.splitlines()[-1] == "quit"

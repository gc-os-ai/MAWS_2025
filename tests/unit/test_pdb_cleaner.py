"""
Tests for :mod:`maws.io.pdb_cleaner`.

Cleaning a structure file is where a run can go wrong most quietly. What comes
out is still a PDB file; LEaP still reads it, still builds something, and the
run still produces a strand. So most of what follows is not about what cleaning
keeps, but about what it must never do to what it keeps: reorder it, delete a
residue out of the middle of a chain, or join two chains together.

The files here are written out a line at a time, in the fixed-column layout the
format uses, so that every field a test depends on can be read off the string
above the assertion.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from maws.errors import ConfigurationError
from maws.io.pdb_cleaner import CleaningReport, clean_pdb_file, resolve_pdb_path


def atom(
    serial: int,
    name: str,
    residue: str,
    chain: str,
    seq: int,
    *,
    record: str = "ATOM",
    alt: str = " ",
    insertion: str = " ",
    occupancy: float = 1.0,
    element: str = "",
    x: float = 0.0,
) -> str:
    """Return one atom line of a PDB file.

    Parameters
    ----------
    serial : int
        The atom's number within the file.
    name : str
        The atom's name within its residue, e.g. ``"CA"``.
    residue : str
        The residue's name, e.g. ``"ALA"``.
    chain : str
        Which chain it belongs to.
    seq : int
        The residue's number within that chain.
    record : {"ATOM", "HETATM"}, default="ATOM"
        Which kind of line to write.
    alt : str, default=" "
        The alternate-position letter, or a space for none.
    insertion : str, default=" "
        The insertion code, or a space for none.
    occupancy : float, default=1.0
        What fraction of the time the atom is there.
    element : str, optional
        The chemical symbol. Left blank by default, as older files leave it.
    x : float, default=0.0
        The atom's position along x, in ångström. The other two are zero.

    Returns
    -------
    str
        An 80-character line, with a newline.

    Notes
    -----
    Atom names are placed from column 14 rather than 13, which is where a
    one-or-two-letter element's name goes. That is the convention the format
    uses to tell ``CA``, an alpha carbon, from ``CA``, a calcium ion.
    """
    return (
        f"{record:<6}{serial:>5} {name:<4}{alt}{residue:>3} {chain}{seq:>4}"
        f"{insertion}   {x:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{occupancy:>6.2f}"
        f"{0.0:>6.2f}          {element:>2}  \n"
    )


def write(tmp_path: Path, *lines: str) -> Path:
    """Write a PDB file out of the given lines and return where it is.

    Parameters
    ----------
    tmp_path : pathlib.Path
        A directory to write into.
    *lines
        The lines, each already ending in a newline.

    Returns
    -------
    pathlib.Path
        Where the file was written.
    """
    path = tmp_path / "target.pdb"
    path.write_text("".join(lines) + "END\n")
    return path


def clean(tmp_path: Path, source: Path, **options):
    """Clean a file into the same directory and return the result and report.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Where to write the tidied copy.
    source : pathlib.Path
        The file to clean.
    **options
        Passed through to :func:`~maws.io.pdb_cleaner.clean_pdb_file`.

    Returns
    -------
    lines : list of str
        The tidied file, split into lines.
    report : maws.io.pdb_cleaner.CleaningReport
        What was removed.
    """
    destination, report = clean_pdb_file(source, tmp_path / "clean.pdb", **options)
    return destination.read_text().splitlines(), report


def atom_lines(lines: list[str]) -> list[str]:
    """Return only the atom lines of a cleaned file.

    Parameters
    ----------
    lines : list of str
        Every line of the file.

    Returns
    -------
    list of str
        The ``ATOM`` and ``HETATM`` lines, in order.
    """
    return [line for line in lines if line.startswith(("ATOM", "HETATM"))]


def residue_order(lines: list[str]) -> list[str]:
    """Return the residue name and original position of each atom line, in order.

    Parameters
    ----------
    lines : list of str
        Every line of a cleaned file.

    Returns
    -------
    list of str
        One entry per atom, as ``"ALA/N"`` — the residue name and the atom
        name. Enough to see at a glance whether anything has been reordered.
    """
    return [
        f"{line[17:20].strip()}/{line[12:16].strip()}" for line in atom_lines(lines)
    ]


class TestNothingIsReordered:
    """The order of the lines is the molecule's chemistry.

    LEaP works out which atoms are bonded to which from the order they appear
    in. A cleaner that sorts the file — even as a step on the way to something
    else — hands back a different molecule, and one that still builds.
    """

    def test_a_file_with_alternate_positions_keeps_its_order(self, tmp_path):
        """Choosing between recorded positions must not disturb anything else.

        Four residues numbered 8 to 11, with residue 10 recorded in two
        positions. Resolving that is a decision about six atoms; the other
        residues, and the order of the atoms within each of them, are none of
        its business.
        """
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 8),
            atom(2, "CA", "ALA", "A", 8),
            atom(3, "C", "ALA", "A", 8),
            atom(4, "N", "GLY", "A", 9),
            atom(5, "CA", "GLY", "A", 9),
            atom(6, "N", "SER", "A", 10, alt="A", occupancy=0.4),
            atom(7, "N", "SER", "A", 10, alt="B", occupancy=0.6),
            atom(8, "N", "VAL", "A", 11),
            atom(9, "CA", "VAL", "A", 11),
        )
        lines, _ = clean(tmp_path, source)
        assert residue_order(lines) == [
            "ALA/N",
            "ALA/CA",
            "ALA/C",
            "GLY/N",
            "GLY/CA",
            "SER/N",
            "VAL/N",
            "VAL/CA",
        ]

    def test_residues_come_out_in_the_order_they_went_in(self, tmp_path):
        """Including when their numbers do not run in order in the original.

        A published file may number residues to match a reference sequence
        rather than the order they appear in. Sorting by that number would
        rearrange the chain.
        """
        source = write(
            tmp_path,
            atom(1, "N", "MET", "A", 12),
            atom(2, "N", "ALA", "A", 3),
            atom(3, "N", "GLY", "A", 7),
        )
        lines, _ = clean(tmp_path, source)
        assert residue_order(lines) == ["MET/N", "ALA/N", "GLY/N"]

    def test_the_end_of_a_chain_is_marked_after_its_last_atom(self, tmp_path):
        """A TER before the residue it closes would join two chains."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "N", "GLY", "B", 1),
        )
        lines, _ = clean(tmp_path, source)
        assert [line[:6].strip() for line in lines] == [
            "ATOM",
            "TER",
            "ATOM",
            "TER",
            "END",
        ]


class TestNoResidueIsDeletedFromTheMiddleOfAChain:
    """Deleting one joins its neighbours to each other, inventing a bond."""

    def test_a_residue_with_an_insertion_code_is_kept(self, tmp_path):
        """An insertion code says "this residue sits between two numbered ones".

        Antibody structures are numbered that way as a matter of course, and
        so is any structure with a loop inserted into a reference sequence. The
        residues carrying one are in the middle of the chain, not at its edges.
        """
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 52),
            atom(2, "N", "GLY", "A", 52, insertion="A"),
            atom(3, "N", "SER", "A", 53),
        )
        lines, _ = clean(tmp_path, source)
        assert residue_order(lines) == ["ALA/N", "GLY/N", "SER/N"]

    def test_insertion_codes_are_replaced_by_plain_numbering(self, tmp_path):
        """Renumbering is what makes it safe to keep them.

        Two residues sharing number 52 and told apart by a letter become
        residues 1 and 2. The letter is then not needed, and the column it sat
        in is cleared, so nothing downstream has to understand it.
        """
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 52),
            atom(2, "N", "GLY", "A", 52, insertion="A"),
        )
        lines, _ = clean(tmp_path, source)
        numbers = [line[22:27] for line in atom_lines(lines)]
        assert numbers == ["   1 ", "   2 "]

    def test_selenomethionine_is_renamed_rather_than_removed(self, tmp_path):
        """MSE is methionine with its sulfur swapped for selenium.

        It appears in a large share of published crystal structures, in the
        middle of chains, because the swap is what makes the structure
        solvable. No protein force field has an entry for it, so it has to
        become something — and the something is the residue it stands in for.
        """
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "N", "MSE", "A", 2),
            atom(3, "SE", "MSE", "A", 2, element="SE"),
            atom(4, "N", "GLY", "A", 3),
        )
        lines, report = clean(tmp_path, source)
        assert residue_order(lines) == ["ALA/N", "MET/N", "MET/SD", "GLY/N"]
        assert report.residues_renamed == 1

    def test_the_swapped_atom_becomes_the_one_it_stands_in_for(self, tmp_path):
        """A methionine has a sulfur where selenomethionine has a selenium.

        Leaving the atom named SE would leave the renamed residue with an atom
        methionine does not have, and LEaP would refuse it.
        """
        source = write(tmp_path, atom(1, "SE", "MSE", "A", 1, element="SE"))
        lines, _ = clean(tmp_path, source)
        line = atom_lines(lines)[0]
        assert line[12:16] == " SD "
        assert line[76:78] == " S"


class TestSeveralRecordedModels:
    """A structure measured by NMR is published as many complete copies."""

    def test_only_the_first_model_is_kept(self, tmp_path):
        """The models are alternative answers, not parts of one molecule.

        Read together they are the same molecule drawn twenty times over, all
        in the same place, with every residue number repeated.
        """
        path = tmp_path / "target.pdb"
        path.write_text(
            "MODEL        1\n"
            + atom(1, "N", "ALA", "A", 1, x=1.0)
            + "ENDMDL\n"
            + "MODEL        2\n"
            + atom(1, "N", "ALA", "A", 1, x=9.0)
            + "ENDMDL\nEND\n"
        )
        lines, report = clean(tmp_path, path)
        assert len(atom_lines(lines)) == 1
        assert report.models_dropped == 1
        assert "1.000" in atom_lines(lines)[0]

    def test_a_file_with_no_models_marked_is_left_alone(self, tmp_path):
        """An ordinary crystal structure has one model and does not say so."""
        source = write(tmp_path, atom(1, "N", "ALA", "A", 1))
        _, report = clean(tmp_path, source)
        assert report.models_dropped == 0


class TestAlternatePositions:
    """Choosing between the positions the experiment could not decide between."""

    def test_one_position_is_chosen_for_a_whole_residue(self, tmp_path):
        """Not atom by atom.

        Per-atom choices mix the two recorded positions, giving a shape that
        was never observed — usually one where the side chain runs into
        itself, because the two positions are alternatives precisely because
        they cannot both be true.
        """
        source = write(
            tmp_path,
            atom(1, "N", "SER", "A", 10, alt="A", occupancy=0.4),
            atom(2, "CA", "SER", "A", 10, alt="A", occupancy=0.9),
            atom(3, "N", "SER", "A", 10, alt="B", occupancy=0.6),
            atom(4, "CA", "SER", "A", 10, alt="B", occupancy=0.1),
        )
        lines, _ = clean(tmp_path, source)
        assert len(atom_lines(lines)) == 2

    def test_the_position_the_experiment_favoured_is_the_one_kept(self, tmp_path):
        """Occupancy is how much of the time the atoms are thought to be there.

        The two positions here are told apart by where they sit along x, so
        the one kept can be read off the coordinates.
        """
        source = write(
            tmp_path,
            atom(1, "N", "SER", "A", 10, alt="A", occupancy=0.3, x=1.0),
            atom(2, "N", "SER", "A", 10, alt="B", occupancy=0.7, x=9.0),
        )
        lines, _ = clean(tmp_path, source)
        assert "9.000" in atom_lines(lines)[0]

    def test_two_different_residues_at_one_position_do_not_get_mixed(self, tmp_path):
        """Sometimes the experiment cannot tell which amino acid is there.

        Choosing per atom takes the shared backbone from one and the side
        chain from the other, and one residue number comes out carrying two
        different residue names — which is not a molecule.
        """
        source = write(
            tmp_path,
            atom(1, "N", "SER", "A", 10, alt="A", occupancy=0.4),
            atom(2, "OG", "SER", "A", 10, alt="A", occupancy=0.4),
            atom(3, "N", "ALA", "A", 10, alt="B", occupancy=0.6),
            atom(4, "CB", "ALA", "A", 10, alt="B", occupancy=0.6),
        )
        lines, _ = clean(tmp_path, source)
        assert {line[17:20] for line in atom_lines(lines)} == {"ALA"}

    def test_the_letter_is_cleared_once_a_position_has_been_chosen(self, tmp_path):
        """A file still declaring position A says the choice was not made."""
        source = write(
            tmp_path,
            atom(1, "N", "SER", "A", 10, alt="A", occupancy=0.9),
            atom(2, "N", "SER", "A", 10, alt="B", occupancy=0.1),
        )
        lines, _ = clean(tmp_path, source)
        assert atom_lines(lines)[0][16] == " "

    def test_atoms_recorded_once_are_untouched(self, tmp_path):
        """Most of a structure has no alternate positions at all."""
        source = write(tmp_path, atom(1, "N", "ALA", "A", 1))
        _, report = clean(tmp_path, source)
        assert report.conformers_resolved == 0


class TestWhatIsKeptAndDropped:
    """Waters go; everything else stays unless it was asked to go."""

    def test_waters_are_dropped(self, tmp_path):
        """A crystal structure carries hundreds of them and they are not the
        target."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "O", "HOH", "A", 2, record="HETATM"),
        )
        lines, report = clean(tmp_path, source)
        assert report.waters_removed == 1
        assert len(atom_lines(lines)) == 1

    def test_a_bound_ligand_is_kept_by_default(self, tmp_path):
        """It is often the most interesting part of the structure.

        A co-crystallised inhibitor sits in the pocket an aptamer would have
        to compete for, so dropping it changes the problem.
        """
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "C1", "0G6", "A", 2, record="HETATM"),
        )
        lines, report = clean(tmp_path, source)
        assert len(atom_lines(lines)) == 2
        assert report.heteroatoms_removed == 0

    def test_a_metal_ion_is_kept_by_default(self, tmp_path):
        """One held inside a protein is part of what holds its shape."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "ZN", "ZN", "A", 2, record="HETATM", element="ZN"),
        )
        lines, _ = clean(tmp_path, source)
        assert len(atom_lines(lines)) == 2

    def test_heteroatoms_go_when_that_is_asked_for(self, tmp_path):
        """The option exists; it is simply not the default."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "C1", "0G6", "A", 2, record="HETATM"),
        )
        lines, report = clean(tmp_path, source, drop_heteroatoms=True)
        assert report.heteroatoms_removed == 1
        assert len(atom_lines(lines)) == 1

    def test_everything_after_the_last_chain_is_kept(self, tmp_path):
        """Ligands and ions are written after the protein, past its TER line.

        Anything that treats the last TER as the end of the file loses them
        all, whatever the settings say.
        """
        path = tmp_path / "target.pdb"
        path.write_text(
            atom(1, "N", "ALA", "A", 1)
            + "TER       2      ALA A   1\n"
            + atom(3, "C1", "0G6", "B", 2, record="HETATM")
            + "END\n"
        )
        lines, _ = clean(tmp_path, path)
        assert len(atom_lines(lines)) == 2

    def test_a_chain_with_no_end_marked_is_kept(self, tmp_path):
        """Files written by modelling tools often leave the last TER out."""
        path = tmp_path / "target.pdb"
        path.write_text(atom(1, "N", "ALA", "A", 1) + atom(2, "N", "GLY", "A", 2))
        lines, _ = clean(tmp_path, path)
        assert len(atom_lines(lines)) == 2


class TestRemovingHydrogens:
    """Dropping the experiment's hydrogens so the force field can place its own."""

    def test_hydrogens_go_when_that_is_asked_for(self, tmp_path):
        """Read from the element column, which modern files fill in."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1, element="N"),
            atom(2, "H", "ALA", "A", 1, element="H"),
        )
        lines, report = clean(tmp_path, source, remove_hydrogens=True)
        assert report.hydrogens_removed == 1
        assert len(atom_lines(lines)) == 1

    def test_a_file_with_no_element_column_still_loses_its_hydrogens(self, tmp_path):
        """Older files leave that column blank, and there are many of them.

        The name is read instead. ``HB2`` is a hydrogen; ``1HB`` is the same
        atom written the other way round, with the number first.
        """
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "HB2", "ALA", "A", 1),
            atom(3, "1HB", "ALA", "A", 1),
        )
        lines, report = clean(tmp_path, source, remove_hydrogens=True)
        assert report.hydrogens_removed == 2
        assert len(atom_lines(lines)) == 1

    def test_deuterium_counts_as_a_hydrogen(self, tmp_path):
        """It is a hydrogen with an extra neutron, used in neutron structures."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1, element="N"),
            atom(2, "D", "ALA", "A", 1, element="D"),
        )
        _, report = clean(tmp_path, source, remove_hydrogens=True)
        assert report.hydrogens_removed == 1

    def test_a_heavy_atom_whose_name_starts_with_h_is_kept(self, tmp_path):
        """Mercury is Hg, and its name does not make it a hydrogen."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1, element="N"),
            atom(2, "HG", "HG", "A", 2, record="HETATM", element="HG"),
        )
        lines, report = clean(tmp_path, source, remove_hydrogens=True)
        assert report.hydrogens_removed == 0
        assert len(atom_lines(lines)) == 2


class TestChoosingChains:
    """Keeping part of a structure that holds a whole complex."""

    def test_every_chain_is_kept_by_default(self, tmp_path):
        """Nothing is thrown away unless it was asked for."""
        source = write(
            tmp_path, atom(1, "N", "ALA", "A", 1), atom(2, "N", "GLY", "B", 1)
        )
        lines, report = clean(tmp_path, source)
        assert len(atom_lines(lines)) == 2
        assert report.chains_dropped == ()

    def test_named_chains_are_kept_and_the_rest_dropped(self, tmp_path):
        """The letters are read one at a time, so "AB" means two chains."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "N", "GLY", "B", 1),
            atom(3, "N", "SER", "C", 1),
        )
        lines, report = clean(tmp_path, source, keep_chains="AC")
        assert {line[21] for line in atom_lines(lines)} == {"A", "C"}
        assert report.chains_dropped == ("B",)

    def test_chain_letters_are_matched_exactly(self, tmp_path):
        """Upper and lower case are different chains, and large assemblies use
        both."""
        source = write(
            tmp_path, atom(1, "N", "ALA", "A", 1), atom(2, "N", "GLY", "a", 1)
        )
        lines, _ = clean(tmp_path, source, keep_chains="A")
        assert {line[21] for line in atom_lines(lines)} == {"A"}

    def test_keeping_one_chain_keeps_the_one_with_the_most_residues(self, tmp_path):
        """Counted in residues, not atoms.

        Counting atoms lets a short chain win by being made of larger
        residues, which is exactly what happens when one chain is a nucleic
        acid and another a protein.
        """
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "CA", "ALA", "A", 1),
            atom(3, "C", "ALA", "A", 1),
            atom(4, "O", "ALA", "A", 1),
            atom(5, "N", "GLY", "B", 1),
            atom(6, "N", "SER", "B", 2),
        )
        lines, _ = clean(tmp_path, source, keep_chains="one")
        assert {line[21] for line in atom_lines(lines)} == {"B"}

    def test_naming_a_chain_the_file_does_not_hold_is_refused(self, tmp_path):
        """Otherwise the answer is an empty file and a run against nothing."""
        source = write(tmp_path, atom(1, "N", "ALA", "A", 1))
        with pytest.raises(ConfigurationError, match="names no chain"):
            clean(tmp_path, source, keep_chains="Z")


class TestRefusingToHandOnSomethingBroken:
    """Cleaning either produces a usable file or says it could not."""

    def test_a_file_with_no_atoms_is_refused(self, tmp_path):
        """A header-only entry, or something that is not a PDB file."""
        path = tmp_path / "target.pdb"
        path.write_text("HEADER    NOTHING HERE\nEND\n")
        with pytest.raises(ConfigurationError, match="no atom records"):
            clean(tmp_path, path)

    def test_the_message_says_what_to_look_at(self, tmp_path):
        """A refusal is only useful if it points somewhere."""
        source = write(tmp_path, atom(1, "N", "ALA", "A", 1))
        with pytest.raises(ConfigurationError, match="case-sensitive"):
            clean(tmp_path, source, keep_chains="z")


class TestTheReport:
    """What was removed is counted and can be logged."""

    def test_a_file_that_needed_nothing_says_so(self, tmp_path):
        """Reading "kept all 2 atoms" is how a user knows cleaning was a no-op."""
        source = write(
            tmp_path, atom(1, "N", "ALA", "A", 1), atom(2, "CA", "ALA", "A", 1)
        )
        _, report = clean(tmp_path, source)
        assert report.summary() == "kept all 2 atoms"

    def test_what_was_removed_is_named_in_the_summary(self, tmp_path):
        """Rather than only a count, so a surprise is recognisable as one."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "O", "HOH", "A", 2, record="HETATM"),
        )
        _, report = clean(tmp_path, source)
        assert "1 water atoms" in report.summary()

    def test_the_atom_count_is_of_the_original_file(self, tmp_path):
        """So the two numbers in the summary can be compared."""
        source = write(
            tmp_path,
            atom(1, "N", "ALA", "A", 1),
            atom(2, "O", "HOH", "A", 2, record="HETATM"),
        )
        _, report = clean(tmp_path, source)
        assert (report.atoms_read, report.atoms_kept) == (2, 1)

    def test_a_report_cannot_be_changed_after_it_is_made(self):
        """It is a record of what happened, not working state."""
        report = CleaningReport(1, 1, 0, 0, 0, 0, (), 0, 0)
        with pytest.raises(AttributeError):
            report.atoms_kept = 5


class TestResolvingWhichFileToUse:
    """The wrapper the rest of MAWS calls."""

    def test_without_cleaning_the_file_is_used_as_given(self, tmp_path):
        """Both paths come back the same, so a caller can report either."""
        source = write(tmp_path, atom(1, "N", "ALA", "A", 1))
        assert resolve_pdb_path(source, "protein", clean_pdb=False) == (
            str(source),
            str(source),
        )

    def test_a_cleaned_copy_is_written_beside_the_original(self, tmp_path):
        """The original is never overwritten, so a run can be repeated."""
        source = write(tmp_path, atom(1, "N", "ALA", "A", 1))
        used, original = resolve_pdb_path(source, "protein", clean_pdb=True)
        assert used != original
        assert Path(used).exists()
        assert Path(original).exists()

    def test_a_target_that_is_not_a_protein_is_left_alone(self, tmp_path):
        """A small organic molecule is one residue with a name of its own.

        There is nothing here to do to it, and dropping its heteroatoms would
        drop the molecule.
        """
        source = write(tmp_path, atom(1, "C1", "LIG", "A", 1, record="HETATM"))
        used, _ = resolve_pdb_path(source, "organic", clean_pdb=True)
        assert used == str(source)

    def test_a_file_that_cannot_be_cleaned_stops_the_run(self, tmp_path):
        """Rather than carrying on with the untidied file.

        A run that believes its target was tidied, and is quietly working from
        the file that was not, produces an answer nobody can account for.
        """
        path = tmp_path / "target.pdb"
        path.write_text("HEADER    NOTHING HERE\nEND\n")
        with pytest.raises(ConfigurationError):
            resolve_pdb_path(path, "protein", clean_pdb=True)


class TestARealDownloadedStructure:
    """The same checks against a file as published, rather than one written here.

    ``1HAO`` is thrombin with a bound inhibitor: three chains, 149 waters, a
    30-atom inhibitor written after the last chain, and 309 atoms in residues
    carrying an insertion code, because antibody-style numbering is used
    throughout. Every one of those is a case a cleaner can quietly damage, and
    they are all present in one file.
    """

    @pytest.fixture
    def thrombin(self, data_dir):
        """The published structure, skipping the tests if it is not there."""
        path = data_dir / "1HAO.pdb"
        if not path.exists():
            pytest.skip(f"example structure not found: {path}")
        return path

    def test_only_the_waters_are_removed(self, thrombin, tmp_path):
        """2769 atoms in, 149 of them water, 2620 out. Nothing else goes."""
        _, report = clean(tmp_path, thrombin)
        assert report.atoms_read - report.waters_removed == report.atoms_kept

    def test_all_three_chains_are_still_separated(self, thrombin, tmp_path):
        """One TER per chain. Two would mean two chains had been joined."""
        lines, _ = clean(tmp_path, thrombin)
        assert sum(1 for line in lines if line.startswith("TER")) == 3

    def test_the_bound_inhibitor_survives(self, thrombin, tmp_path):
        """It is written after the last chain, and it is the interesting part.

        An aptamer designed against thrombin would be competing for the pocket
        this molecule sits in.
        """
        lines, _ = clean(tmp_path, thrombin)
        assert sum(1 for line in atom_lines(lines) if line[17:20] == "0G6") == 30

    def test_the_residues_with_insertion_codes_survive(self, thrombin, tmp_path):
        """309 of the file's atoms are in them, spread through the chains.

        They are not at the ends: this numbering puts them wherever the
        structure has more residues than the reference sequence it is numbered
        against. Dropping them would leave holes through the middle.
        """
        coded = sum(
            1
            for line in thrombin.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM")) and line[26] != " "
        )
        assert coded == 309
        _, report = clean(tmp_path, thrombin)
        assert report.atoms_kept == report.atoms_read - report.waters_removed

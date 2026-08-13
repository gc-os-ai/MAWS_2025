"""
maws.topology
=============

What to build, and what came out of building it.

MAWS designs an *aptamer*: a short strand of RNA or DNA that folds up against a
target molecule and sticks to it. Both molecules are described together as an
*assembly*, which is a list of chains — normally two, one for the aptamer and
one for the target. A chain is a run of *residues*, the repeating units a
molecule is built from; for an aptamer those are single nucleotides joined end
to end. The two ends of a nucleic acid strand are chemically different and are
called 5' and 3'; a sequence is written 5'-to-3'.

Two types split the work in half:

:class:`Assembly`
    A description of what to build. Frozen and cheap: making one opens no file
    and starts no program. A run can therefore be described, compared with
    another description, and looked up in a cache before anything is built.

:class:`BuiltSystem`
    The result of building one. Every field is filled in — coordinates,
    element symbols, masses — so there is no half-built state for anything to
    guard against.

A chain of an assembly is described in one of two ways, and its type says
which. :class:`ResidueChain` is a chain whose residues are known, so its atom
count is known too. :class:`PdbChain` is a chain that is still only a file on
disk: measured atom positions for a target molecule, with nothing counted and
no parameters worked out.

Getting from the second to the first means running AmberTools. ``antechamber``
and ``parmchk2`` work out the parameters for a molecule that the *force field*
— the collection of numbers saying how strongly each pair of atoms pulls on or
pushes away the other — does not already cover. ``tleap``, the program that
builds structures, then writes two files: a ``.prmtop`` holding the parameters
and the list of which atoms are bonded to which, and an ``.inpcrd`` holding the
starting coordinates.

See Also
--------
maws.build : Turns an :class:`Assembly` into a :class:`BuiltSystem`.
maws.pose : The coordinate and window types :class:`BuiltSystem` hands out.

Examples
--------
>>> from maws.libraries import rna
>>> assembly = (
...     Assembly().with_aptamer(rna(), sequence="G A").with_ligand_stub(n_atoms=40)
... )
>>> assembly.roles
('aptamer', 'ligand')
>>> assembly.n_atoms()
106
>>> assembly.with_sequence("aptamer", "G A U").n_atoms()  # a new Assembly
136
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

from maws.errors import ConfigurationError
from maws.forcefield import ForceField
from maws.pose import ChainView, Pose
from maws.values import AtomRange, NucleotideSequence, ResidueLibrary

if TYPE_CHECKING:  # pragma: no cover - only needed by type checkers
    from openmm.app import Topology

    from maws.energy import EnergyModel

__all__ = [
    "AmberArtifacts",
    "Assembly",
    "BuiltSystem",
    "ChainSpec",
    "PdbChain",
    "ResidueChain",
]

DEFAULT_STUB_ELEMENT = "C"
"""Element symbol given to every atom of a stub target.

A stub target is a stand-in with the right number of atoms and no chemistry at
all, used to exercise MAWS without AmberTools installed. Carbon is the choice
because it is unremarkable: nothing reading the element symbols has to make an
exception for it. See :meth:`Assembly.with_ligand_stub`.
"""

DEFAULT_STUB_MASS = 12.011
"""Mass in daltons given to every atom of a stub target, that of carbon."""


# ---------------------------------------------------------------------------
# Chain descriptions
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _ChainSpecBase:
    """Fields shared by every kind of chain description.

    Not part of the public interface. Use :class:`ResidueChain` or
    :class:`PdbChain`.

    Parameters
    ----------
    role : str
        The chain's name within its assembly, e.g. ``"aptamer"`` or
        ``"ligand"``. Chains are looked up by this name and never by position,
        so adding a third chain cannot quietly change which chain a piece of
        code is talking about.
    """

    role: str


@dataclass(frozen=True, slots=True)
class ResidueChain(_ChainSpecBase):
    """A chain whose residues are known, so its atom count is too.

    Because every residue's description says how many atoms it has, a chain
    like this can answer questions about its size without anything being built
    or any file being read.

    Parameters
    ----------
    role : str
        The chain's name within its assembly, e.g. ``"aptamer"``.
    library : maws.values.ResidueLibrary
        The descriptions the sequence's nucleotides are looked up in: how many
        atoms each has, and which of its bonds can be turned.
    sequence : maws.values.NucleotideSequence
        The chain's nucleotides as written, 5' end first.

    See Also
    --------
    PdbChain : A chain whose parameters have not been worked out yet.

    Examples
    --------
    >>> from maws.libraries import rna
    >>> chain = ResidueChain("aptamer", rna(), NucleotideSequence.parse("G A"))
    >>> chain.n_atoms
    66
    >>> chain.canonical
    ('G5', 'A3')
    """

    library: ResidueLibrary
    sequence: NucleotideSequence

    @property
    def n_atoms(self) -> int:
        """int : Total atoms across every residue of the sequence."""
        return self.library.atom_count(self.sequence.canonical(self.library))

    @property
    def canonical(self) -> tuple[str, ...]:
        """tuple of str : The name ``tleap`` knows each residue by.

        A written nucleotide such as ``"G"`` becomes a more specific name
        depending on where in the strand it sits, because a residue at an end
        of the strand carries an extra group that one in the middle does not.
        """
        return self.sequence.canonical(self.library)

    def with_sequence(self, sequence: NucleotideSequence) -> ResidueChain:
        """Return a copy of this chain carrying a different sequence.

        Parameters
        ----------
        sequence : maws.values.NucleotideSequence
            The nucleotides to use instead.

        Returns
        -------
        ResidueChain
            A new chain description. This one is unchanged.
        """
        return replace(self, sequence=sequence)


@dataclass(frozen=True, slots=True)
class PdbChain(_ChainSpecBase):
    """A chain that is still a PDB file, with no parameters worked out.

    A PDB file lists atoms and their measured positions. It does not say how
    those atoms behave, so before the molecule can be built or scored,
    ``antechamber`` and ``parmchk2`` have to work out its parameters and write
    them to a ``.lib`` file. That means running programs and writing to disk,
    which describing a run must not do — so a chain stays in this form until
    :class:`maws.build.LeapBuilder` is asked to do the work, and an
    :class:`Assembly` holding one can still be made and copied freely.

    How many atoms such a chain has is unknown until then, which is why
    :meth:`Assembly.n_atoms` refuses to answer while one is present.

    Parameters
    ----------
    role : str
        The chain's name within its assembly, e.g. ``"ligand"``.
    path : pathlib.Path
        Where the PDB file is.
    residue_name : str
        Name ``tleap`` will know this molecule by, and the stem of the ``.lib``
        file written for it. Two chains sharing a name would overwrite each
        other's parameters, so the name has to be unique per molecule.
    parameterized : bool
        Whether the target's force field already describes this molecule. When
        True, ``antechamber`` and ``parmchk2`` are skipped and the PDB file
        goes straight to ``tleap``.

    See Also
    --------
    ResidueChain : What this becomes once its parameters exist.
    maws.build.LeapBuilder.prepare : Does the work of converting one.
    """

    path: Path
    residue_name: str
    parameterized: bool


ChainSpec = ResidueChain | PdbChain
"""Either kind of chain description. Tell them apart with ``match`` or
``isinstance``."""


def default_residue_name(path: Path) -> str:
    """Return a residue name for a PDB file, derived from its contents.

    Parameters
    ----------
    path : pathlib.Path
        Where the PDB file is.

    Returns
    -------
    str
        An uppercase name of the form ``L<six hex digits>``, worked out from a
        digest of the file's bytes. If the file cannot be read, the path is
        digested instead, so a run can still be described for a file that has
        not been written yet.

    Notes
    -----
    The name decides what the molecule's parameter file is called, and two
    molecules sharing a name would overwrite each other's parameters. Taking
    the name from the contents makes that impossible: different files get
    different names, and the same file always gets the same name, so a
    parameter file already on disk is found again rather than rebuilt.

    Examples
    --------
    >>> name = default_residue_name(Path("nonexistent.pdb"))
    >>> len(name), name[0]
    (7, 'L')
    """
    import hashlib

    try:
        payload = path.read_bytes()
    except OSError:
        payload = str(path).encode()
    return "L" + hashlib.sha1(payload).hexdigest()[:6].upper()


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Assembly:
    """The chains to build, and nothing else.

    Every ``with_*`` method returns a new ``Assembly`` and leaves the one it
    was called on exactly as it was, in the same way that ``"abc".upper()``
    gives back a new string. A search relies on that: it tries several ways of
    extending the same aptamer and needs the description it started from still
    intact after each attempt. The same property makes an assembly safe to use
    as a cache key.

    Parameters
    ----------
    chains : tuple of ChainSpec, default=()
        The chains, in the order their atoms will be laid out. Atoms of the
        whole design are held in one array, and the first chain here takes the
        lowest indices in it.

    See Also
    --------
    BuiltSystem : What building one produces.

    Examples
    --------
    >>> from maws.libraries import rna
    >>> empty = Assembly()
    >>> one = empty.with_aptamer(rna(), sequence="G")
    >>> empty.roles, one.roles
    ((), ('aptamer',))
    """

    chains: tuple[ChainSpec, ...] = ()

    # -- inspection -------------------------------------------------------

    @property
    def roles(self) -> tuple[str, ...]:
        """tuple of str : The name of each chain, in the order they are laid out."""
        return tuple(chain.role for chain in self.chains)

    def __len__(self) -> int:
        return len(self.chains)

    def __iter__(self) -> Iterator[ChainSpec]:
        return iter(self.chains)

    def chain(self, role: str) -> ChainSpec:
        """Return the description of the chain with the given name.

        Parameters
        ----------
        role : str
            Name to look up, e.g. ``"aptamer"``.

        Returns
        -------
        ResidueChain or PdbChain
            The matching description.

        Raises
        ------
        maws.errors.ConfigurationError
            If no chain has that name. The message lists the names that exist.

        Examples
        --------
        >>> from maws.libraries import rna
        >>> assembly = Assembly().with_aptamer(rna(), sequence="G A")
        >>> assembly.chain("aptamer").n_atoms
        66
        """
        for chain in self.chains:
            if chain.role == role:
                return chain
        raise ConfigurationError(
            f"no chain named {role!r} in this assembly. Chains: "
            f"{', '.join(self.roles) or '(none)'}"
        )

    def n_atoms(self) -> int:
        """Return how many atoms the whole assembly comes to.

        Returns
        -------
        int
            The total across every chain.

        Raises
        ------
        maws.errors.ConfigurationError
            If any chain is still a :class:`PdbChain`. Nobody has counted the
            atoms in that file yet, so there is no total to give; build the
            assembly first, or use :meth:`with_ligand_stub` to stand in for the
            target.

        Examples
        --------
        >>> from maws.libraries import rna
        >>> Assembly().with_aptamer(rna(), sequence="G").n_atoms()
        33
        """
        total = 0
        for chain in self.chains:
            if isinstance(chain, PdbChain):
                raise ConfigurationError(
                    f"chain {chain.role!r} comes from {chain.path} and has not "
                    f"been parameterised, so its atom count is not known yet. "
                    f"Build the assembly first."
                )
            total += chain.n_atoms
        return total

    # -- construction -----------------------------------------------------

    def with_chain(self, spec: ChainSpec) -> Assembly:
        """Return a copy of this assembly with one more chain.

        Parameters
        ----------
        spec : ResidueChain or PdbChain
            The chain to add. It goes on the end, so it takes the highest atom
            indices of the design.

        Returns
        -------
        Assembly
            A new assembly. This one does not change.

        Raises
        ------
        maws.errors.ConfigurationError
            If a chain of that name is already there. Names are how chains are
            looked up, so two chains cannot share one.
        """
        if spec.role in self.roles:
            raise ConfigurationError(
                f"assembly already has a chain named {spec.role!r}"
            )
        return Assembly(self.chains + (spec,))

    def with_aptamer(
        self,
        library: ResidueLibrary,
        sequence: str | NucleotideSequence = "",
        *,
        role: str = "aptamer",
    ) -> Assembly:
        """Return a copy of this assembly with an aptamer chain added.

        Parameters
        ----------
        library : maws.values.ResidueLibrary
            The nucleotide descriptions to resolve `sequence` against, e.g.
            from :func:`maws.libraries.rna`. This is what decides whether the
            strand is RNA or DNA.
        sequence : str or maws.values.NucleotideSequence, default=""
            The nucleotides to start from, 5' end first, written as
            whitespace-separated tokens such as ``"G A U"``. Empty is the usual
            case: the search adds one nucleotide at a time.
        role : str, default="aptamer"
            Name to look this chain up by later.

        Returns
        -------
        Assembly
            A new assembly. This one does not change.

        Examples
        --------
        >>> from maws.libraries import rna
        >>> Assembly().with_aptamer(rna(), "G A").chain("aptamer").canonical
        ('G5', 'A3')
        """
        tokens = (
            NucleotideSequence.parse(sequence)
            if isinstance(sequence, str)
            else sequence
        )
        return self.with_chain(ResidueChain(role, library, tokens))

    def with_ligand(
        self,
        path: str | Path,
        forcefield: ForceField,
        *,
        role: str = "ligand",
        residue_name: str | None = None,
    ) -> Assembly:
        """Return a copy of this assembly with a target chain added from a PDB.

        Nothing is read and no program is run here. The chain is recorded as a
        :class:`PdbChain` and stays that way until it is built, so this is safe
        to call on a file that does not exist yet.

        Parameters
        ----------
        path : str or pathlib.Path
            Where the target's PDB file is.
        forcefield : maws.forcefield.ForceField
            The run's force fields. Only ``parameterized`` is read here, to
            record whether ``antechamber`` and ``parmchk2`` will have to work
            out this molecule's parameters before it can be built.
        role : str, default="ligand"
            Name to look this chain up by later.
        residue_name : str, optional
            Name ``tleap`` should know the molecule by. By default it is
            derived from the file's contents, so two different targets cannot
            end up writing to the same parameter file.

        Returns
        -------
        Assembly
            A new assembly. This one does not change.

        See Also
        --------
        with_ligand_stub : A stand-in target, for running without AmberTools.
        default_residue_name : How the default name is worked out.
        """
        resolved = Path(path)
        return self.with_chain(
            PdbChain(
                role=role,
                path=resolved,
                residue_name=residue_name or default_residue_name(resolved),
                parameterized=forcefield.parameterized,
            )
        )

    def with_ligand_stub(
        self,
        n_atoms: int,
        *,
        role: str = "ligand",
        residue_name: str = "LIG",
    ) -> Assembly:
        """Return a copy of this assembly with a stand-in target chain added.

        A stub target has the right number of atoms and no chemistry: none of
        its bonds can be turned and no parameters are worked out for it. It is
        there so that the search, the samplers and the bookkeeping can be run
        and tested on a machine with no AmberTools installed.

        Parameters
        ----------
        n_atoms : int
            How many atoms the stand-in should have. Pick roughly what the real
            target has, since it sets how much of the coordinate array the
            target occupies.
        role : str, default="ligand"
            Name to look this chain up by later.
        residue_name : str, default="LIG"
            Name for the stand-in's single residue. Nothing reads it from disk,
            so any name will do.

        Returns
        -------
        Assembly
            A new assembly. This one does not change.

        See Also
        --------
        with_ligand : Adds a real target, from a PDB file.
        maws.build.FakeBuilder : Builds an assembly of stand-ins into a system.

        Examples
        --------
        >>> Assembly().with_ligand_stub(n_atoms=40).n_atoms()
        40
        """
        library = ResidueLibrary.single(residue_name, n_atoms)
        sequence = NucleotideSequence((residue_name,))
        return self.with_chain(ResidueChain(role, library, sequence))

    def with_sequence(self, role: str, sequence: str | NucleotideSequence) -> Assembly:
        """Return a copy of this assembly with one chain's sequence replaced.

        This is how the search grows an aptamer: it makes one of these per
        candidate nucleotide, builds each, and keeps the best.

        Parameters
        ----------
        role : str
            Which chain to change, e.g. ``"aptamer"``.
        sequence : str or maws.values.NucleotideSequence
            The nucleotides to use instead, 5' end first.

        Returns
        -------
        Assembly
            A new assembly. This one does not change, which is what lets the
            search try several candidates against the same starting point.

        Raises
        ------
        maws.errors.ConfigurationError
            If no chain has that name, or the chain is a :class:`PdbChain`,
            which is a molecule read from a file and has no sequence to edit.

        Examples
        --------
        >>> from maws.libraries import rna
        >>> start = Assembly().with_aptamer(rna(), "G")
        >>> grown = start.with_sequence("aptamer", "G A")
        >>> start.n_atoms(), grown.n_atoms()
        (33, 66)
        """
        target = self.chain(role)
        if not isinstance(target, ResidueChain):
            raise ConfigurationError(
                f"chain {role!r} comes from a PDB file and has no editable sequence"
            )
        tokens = (
            NucleotideSequence.parse(sequence)
            if isinstance(sequence, str)
            else sequence
        )
        return Assembly(
            tuple(
                target.with_sequence(tokens) if chain.role == role else chain
                for chain in self.chains
            )
        )


def compute_spans(chains: Sequence[ResidueChain]) -> dict[str, AtomRange]:
    """Lay chains out end to end and return the atom range each one occupies.

    Every atom of a design lives in one array, so a chain is a run of rows in
    it. This puts the first chain at the start of the array and each later one
    immediately after the previous.

    Parameters
    ----------
    chains : sequence of ResidueChain
        The chains, in the order they are to be laid out. Each must know its
        own atom count, which is why a :class:`PdbChain` cannot appear here.

    Returns
    -------
    dict of str to maws.values.AtomRange
        For each chain name, the range of atom indices that chain covers.

    Notes
    -----
    The ranges depend on nothing but the chain lengths, and are worked out
    afresh every time. Changing one chain's sequence therefore cannot leave
    another chain's range stale — there is no stored offset to update.

    Examples
    --------
    >>> from maws.libraries import rna
    >>> compute_spans(
    ...     [
    ...         ResidueChain("aptamer", rna(), NucleotideSequence.parse("G")),
    ...         ResidueChain(
    ...             "ligand",
    ...             ResidueLibrary.single("LIG", 40),
    ...             NucleotideSequence(("LIG",)),
    ...         ),
    ...     ]
    ... )
    {'aptamer': AtomRange(start=0, stop=33), 'ligand': AtomRange(start=33, stop=73)}
    """
    spans: dict[str, AtomRange] = {}
    cursor = 0
    for chain in chains:
        length = chain.n_atoms
        spans[chain.role] = AtomRange(cursor, cursor + length)
        cursor += length
    return spans


# ---------------------------------------------------------------------------
# Built system
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AmberArtifacts:
    """The files ``tleap`` wrote when it built a system.

    A :class:`BuiltSystem` has one of these when it was built by ``tleap``, and
    none when it was made from stand-ins. That is exactly the difference
    between a system whose energy can be calculated and one whose cannot: the
    simulation package needs these files to know what the atoms are.

    Parameters
    ----------
    prmtop_path : pathlib.Path
        Where the ``.prmtop`` file is: the parameters for every atom, and the
        list of which atoms are bonded to which.
    inpcrd_path : pathlib.Path
        Where the ``.inpcrd`` file is: the coordinates ``tleap`` started the
        molecule off at.
    topology : openmm.app.Topology
        The same atom-and-bond listing, already read into the form the
        simulation package works with, so it does not have to be parsed again.
    """

    prmtop_path: Path
    inpcrd_path: Path
    topology: Topology


class BuiltSystem:
    """A built design: every chain laid out, with coordinates to match.

    All the chains' atoms are held together in one array, in the order the
    assembly listed them. Everything that follows — moving the aptamer,
    turning its bonds, calculating its energy — works on that array, and this
    object is what says which stretch of it belongs to which chain.

    Parameters
    ----------
    assembly : Assembly
        The description this was built from. Every chain in it must already
        know its residues.
    forcefield : maws.forcefield.ForceField
        The force fields the structure was built under. Kept so that whatever
        scores this system later uses the same physics it was built with.
    xyz : array_like
        Shape ``(N, 3)`` starting coordinates in ångström, one row per atom.
        An ångström is 10⁻¹⁰ m, roughly the width of one atom.
    elements : sequence of str
        Element symbol for each atom, e.g. ``"C"``. Samplers read these to
        judge how much room each atom takes up, and so how close another
        molecule may come.
    masses : sequence of float
        Mass of each atom in daltons, used to find the centre of mass of a
        molecule, which is heavier towards its heavier atoms rather than simply
        the middle of its atoms.
    amber : AmberArtifacts, optional
        The files ``tleap`` wrote, when ``tleap`` is what built this. Without
        them the system can be moved about but not scored.

    Attributes
    ----------
    pose : maws.pose.Pose
        The coordinates as built.
    spans : dict of str to maws.values.AtomRange
        For each chain name, the range of atom indices it covers.

    Raises
    ------
    maws.errors.ConfigurationError
        If the assembly still holds a :class:`PdbChain`, if `elements` or
        `masses` do not have exactly one entry per atom, or if the chains and
        the coordinates disagree about how many atoms there are.

    See Also
    --------
    maws.build.build : The usual way to obtain one.

    Examples
    --------
    >>> import numpy as np
    >>> from maws.forcefield import ForceField
    >>> from maws.libraries import rna
    >>> assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(10)
    >>> system = BuiltSystem(
    ...     assembly,
    ...     ForceField.for_target("RNA", "protein"),
    ...     np.zeros((43, 3)),
    ...     elements=["C"] * 43,
    ...     masses=[12.011] * 43,
    ... )
    >>> system.chain("ligand").span
    AtomRange(start=33, stop=43)
    """

    __slots__ = (
        "_amber",
        "_assembly",
        "_elements",
        "_forcefield",
        "_masses",
        "_pose",
        "_spans",
        "_views",
    )

    def __init__(
        self,
        assembly: Assembly,
        forcefield: ForceField,
        xyz: object,
        *,
        elements: Sequence[str],
        masses: Sequence[float],
        amber: AmberArtifacts | None = None,
    ) -> None:
        chains = []
        for chain in assembly.chains:
            if not isinstance(chain, ResidueChain):
                raise ConfigurationError(
                    f"chain {chain.role!r} is not parameterized yet, so it "
                    f"cannot be part of a built system"
                )
            chains.append(chain)

        self._assembly = assembly
        self._forcefield = forcefield
        self._spans = compute_spans(chains)
        self._amber = amber
        self._pose = Pose(xyz, self)
        self._elements = tuple(elements)
        self._masses = tuple(float(m) for m in masses)

        expected = self._pose.n_atoms
        for label, table in (("elements", self._elements), ("masses", self._masses)):
            if len(table) != expected:
                raise ConfigurationError(
                    f"{label} has {len(table)} entries but the system has "
                    f"{expected} atoms"
                )
        total = sum(len(span) for span in self._spans.values())
        if total != expected:
            raise ConfigurationError(
                f"chains account for {total} atoms but the coordinates hold "
                f"{expected}. The assembly and the build disagree."
            )

        self._views = {
            chain.role: ChainView.build(
                chain.role,
                self._spans[chain.role],
                chain.sequence,
                chain.library,
            )
            for chain in chains
        }

    # -- inspection -------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"<BuiltSystem {self.n_atoms} atoms, "
            f"{len(self._views)} chains: {', '.join(self._views)}>"
        )

    @property
    def assembly(self) -> Assembly:
        """Assembly : The description this system was built from."""
        return self._assembly

    @property
    def forcefield(self) -> ForceField:
        """maws.forcefield.ForceField : The physics this was built under."""
        return self._forcefield

    @property
    def pose(self) -> Pose:
        """maws.pose.Pose : The coordinates as built.

        Moving atoms about produces a new set of coordinates and never alters
        these, so this is the starting shape for as long as the system exists,
        however many shapes are derived from it.
        """
        return self._pose

    @property
    def spans(self) -> Mapping[str, AtomRange]:
        """Mapping of str to maws.values.AtomRange : Atom range per chain."""
        return dict(self._spans)

    @property
    def n_atoms(self) -> int:
        """int : How many atoms the whole system has."""
        return self._pose.n_atoms

    @property
    def elements(self) -> tuple[str, ...]:
        """tuple of str : Element symbol for each atom, e.g. ``"C"``."""
        return self._elements

    @property
    def masses(self) -> tuple[float, ...]:
        """tuple of float : Mass of each atom, in daltons."""
        return self._masses

    @property
    def amber(self) -> AmberArtifacts | None:
        """AmberArtifacts or None : The files ``tleap`` wrote, if it built this."""
        return self._amber

    def chain(self, role: str) -> ChainView:
        """Return a window onto one chain of this system.

        Parameters
        ----------
        role : str
            Which chain, e.g. ``"aptamer"``.

        Returns
        -------
        maws.pose.ChainView
            A window naming that chain's atoms and residues, which is what
            everything that turns a bond or moves a molecule is given.

        Raises
        ------
        maws.errors.ConfigurationError
            If no chain has that name. The message lists the ones that exist.
        """
        try:
            return self._views[role]
        except KeyError:
            raise ConfigurationError(
                f"no chain named {role!r} in this system. Chains: "
                f"{', '.join(self._views)}"
            ) from None

    def energy_model(
        self, *, platform: str | None = None, freeze: str | None = None
    ) -> EnergyModel:
        """Build an OpenMM energy model for this system.

        Parameters
        ----------
        platform : str, optional
            OpenMM platform name to force, e.g. ``"CPU"``. By default the
            fastest available is chosen.
        freeze : str, optional
            Name of a chain to hold in place while settling, normally the
            target. The search compares candidates by how each one sits
            against the same target, so a target that reshapes itself
            differently for each of them is not the same target. Nothing is
            held in place by default.

        Returns
        -------
        maws.energy.EnergyModel
            A model backed by an OpenMM context over this system's topology.

        Raises
        ------
        maws.errors.ConfigurationError
            If this system was not built by LEaP and so has no AMBER topology
            to hand to OpenMM.

        See Also
        --------
        maws.energy.StubEnergy : The analytic stand-in used in tests.
        """
        if self._amber is None:
            raise ConfigurationError(
                "this system has no AMBER topology, so OpenMM cannot score it. "
                "It was built without LEaP — use maws.energy.StubEnergy, or "
                "build with maws.build.LeapBuilder."
            )

        from maws.energy import OpenMMEnergy

        return OpenMMEnergy.from_prmtop(
            self._amber.prmtop_path,
            self._forcefield,
            platform=platform,
            frozen=None if freeze is None else self.chain(freeze).span,
        )

"""
maws.build
==========

Turning a description of a design into a structure with real coordinates.

An :class:`~maws.topology.Assembly` says which molecules a design contains: an
aptamer of so many *residues* — single nucleotides — and a target molecule to
design it against. It holds no coordinates, because working those out means
running AmberTools: a set of programs, installed separately from this package,
that know what a molecule made of given residues actually looks like. ``tleap``
— the command-line form of the modelling program LEaP — builds the structure,
while ``antechamber`` and ``parmchk2`` work out a *force field* for a molecule
nothing already describes. A force field is the collection of numbers saying
how strongly each pair of atoms pulls on or pushes away the other.

Building produces two files, and their extensions are used as their names
throughout: a *prmtop* holding the force field numbers for this particular
collection of atoms, and an *inpcrd* holding one position per atom. Both are
read back at once into a :class:`~maws.topology.BuiltSystem`.

Two builders are provided, and choosing between them is what separates a real
design run from a test of the search itself.

:class:`LeapBuilder`
    Runs AmberTools. This is the only place in MAWS that starts another program
    or writes a file, and the only way to get a structure that can be scored
    with real physics.

:class:`FakeBuilder`
    Puts atoms on a regular grid, one chain next to the other. The result is
    not a molecule and its energies mean nothing, but it has the right number
    of atoms in the right chains, so every part of MAWS that arranges and
    scores atoms can be exercised with nothing installed.

Both satisfy the :class:`Builder` protocol, so anything that builds takes one
as an argument rather than choosing for itself, and a design run can be tested
end to end by passing the other.

See Also
--------
maws.topology.Assembly : What goes in.
maws.topology.BuiltSystem : What comes out.

Examples
--------
>>> from maws.forcefield import ForceField
>>> from maws.libraries import rna
>>> from maws.topology import Assembly
>>> assembly = Assembly().with_aptamer(rna(), "G").with_ligand_stub(20)
>>> system = FakeBuilder().build(assembly, ForceField.for_target("RNA", "protein"))
>>> system.n_atoms
53
>>> system.chain("aptamer").span
AtomRange(start=0, stop=33)
"""

from __future__ import annotations

import hashlib
import json
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Protocol, runtime_checkable

import numpy as np

from maws.errors import BuildError, ConfigurationError
from maws.forcefield import ForceField
from maws.io.tools import run_tool
from maws.topology import (
    DEFAULT_STUB_ELEMENT,
    DEFAULT_STUB_MASS,
    AmberArtifacts,
    Assembly,
    BuiltSystem,
    PdbChain,
    ResidueChain,
)
from maws.values import NucleotideSequence, ResidueLibrary

__all__ = ["Builder", "FakeBuilder", "LeapBuilder", "build"]

DEFAULT_CACHE_DIR = Path(".maws_cache")
"""Directory built structures are cached in, relative to the working directory.

Building the same design twice gives the same structure, so the result is
stored under a name derived from the description and reused.
"""

GRID_SPACING = 1.5
"""Distance in ångström between neighbouring atoms of a :class:`FakeBuilder` grid.

Close to a real covalent bond length, so distances in a fake structure are at
least the right order of magnitude.
"""

CHAIN_SEPARATION = 25.0
"""Distance in ångström between the starting corners of two chains' grids.

This is the offset between grid origins, not the clearance between the chains:
each grid has a width of its own, so the nearest pair of atoms across two
chains sits closer than this. Far enough apart that the chains do not overlap,
near enough that a sampler looking for positions around one of them has
somewhere to put the other.
"""


@runtime_checkable
class Builder(Protocol):
    """Anything that can turn a design description into a structure.

    Two methods are needed. :meth:`prepare` takes a chain that is still just a
    filename and works out how many atoms it has, so the design knows how big
    it is. :meth:`build` takes a whole description and produces positions for
    every atom in it. Any object with those two methods can be passed wherever
    a builder is wanted; nothing needs to inherit from this class.

    Checked at run time, so ``isinstance(obj, Builder)`` works, though it only
    confirms the two methods exist and not what they do.

    See Also
    --------
    LeapBuilder : Runs AmberTools; the real one.
    FakeBuilder : Puts atoms on a grid; needs nothing installed.

    Examples
    --------
    >>> isinstance(FakeBuilder(), Builder)
    True
    >>> isinstance("not a builder", Builder)
    False
    """

    def prepare(self, chain: PdbChain, forcefield: ForceField) -> ResidueChain:
        """Work out how many atoms a molecule read from a file has.

        A design cannot say which stretch of its atom array belongs to the
        target until the target's size is known, and the only way to know that
        is to look at the file.

        Parameters
        ----------
        chain : maws.topology.PdbChain
            The chain to resolve, naming the file to read.
        forcefield : maws.forcefield.ForceField
            Which parameters to describe the molecule with.

        Returns
        -------
        maws.topology.ResidueChain
            The same chain, with its atom count known. Treated from here on as
            a single rigid residue, however large the molecule is.
        """
        ...

    def build(self, assembly: Assembly, forcefield: ForceField) -> BuiltSystem:
        """Turn a description into a structure with coordinates.

        Parameters
        ----------
        assembly : maws.topology.Assembly
            What to build. Any chain still described only by a filename is
            resolved with :meth:`prepare` first.
        forcefield : maws.forcefield.ForceField
            Which parameters to build it with.

        Returns
        -------
        maws.topology.BuiltSystem
            The built structure, with one position in ångström per atom.
        """
        ...


def build(
    assembly: Assembly,
    forcefield: ForceField,
    *,
    builder: Builder | None = None,
) -> BuiltSystem:
    """build(assembly, forcefield, *, builder=None) -> BuiltSystem

    Turn a design description into a structure with coordinates.

    A one-line convenience over calling a builder directly, so that code which
    has no opinion about how a structure is made does not have to name one.

    Parameters
    ----------
    assembly : maws.topology.Assembly
        What to build.
    forcefield : maws.forcefield.ForceField
        Which parameters to build it with.
    builder : Builder, optional
        Which builder to use. Defaults to a :class:`LeapBuilder`, which needs
        AmberTools installed. Pass a :class:`FakeBuilder` to build without it.

    Returns
    -------
    maws.topology.BuiltSystem
        The built structure.

    Raises
    ------
    maws.errors.ConfigurationError
        If the assembly has nothing in it to build.
    maws.errors.ToolchainError
        If the default builder is used and AmberTools is not installed.

    See Also
    --------
    LeapBuilder : What is used by default.
    FakeBuilder : What to pass when nothing is installed.

    Examples
    --------
    >>> from maws.forcefield import ForceField
    >>> from maws.libraries import rna
    >>> from maws.topology import Assembly
    >>> system = build(
    ...     Assembly().with_aptamer(rna(), "G A").with_ligand_stub(12),
    ...     ForceField.for_target("RNA", "protein"),
    ...     builder=FakeBuilder(),
    ... )
    >>> system.n_atoms
    78
    """
    chosen = builder if builder is not None else LeapBuilder()
    return chosen.build(assembly, forcefield)


def _resolved_chains(
    assembly: Assembly, forcefield: ForceField, builder: Builder
) -> tuple[ResidueChain, ...]:
    """Resolve every chain of an assembly to one with a known atom count.

    Parameters
    ----------
    assembly : maws.topology.Assembly
        The description whose chains are to be resolved.
    forcefield : maws.forcefield.ForceField
        Which parameters to describe unresolved molecules with.
    builder : Builder
        Whose :meth:`~Builder.prepare` does the resolving.

    Returns
    -------
    tuple of maws.topology.ResidueChain
        The chains, in the order they appear in the assembly.
    """
    resolved: list[ResidueChain] = []
    for chain in assembly.chains:
        if isinstance(chain, PdbChain):
            resolved.append(builder.prepare(chain, forcefield))
        else:
            resolved.append(chain)
    return tuple(resolved)


class FakeBuilder:
    """Builds a structure by putting atoms on a grid, with no chemistry at all.

    Each chain gets its own cubic grid of atoms, spaced :data:`GRID_SPACING`
    apart, and the grids are placed :data:`CHAIN_SEPARATION` apart along the
    x axis. Every atom is given carbon's symbol and mass.

    The result is not a molecule. Its energies are meaningless and it must
    never be used for a real design. What it is good for is everything else:
    the chain and residue windows are correct, the atom counts are correct, and
    the turnable bonds point at real atoms, so arranging, sampling, scoring
    with :class:`~maws.energy.StubEnergy`, and the whole growth search can be
    exercised in milliseconds with nothing installed.

    Parameters
    ----------
    spacing : float, default=1.5
        Distance in ångström between neighbouring grid atoms. Raising it
        spreads each chain out.
    separation : float, default=25.0
        Distance in ångström between one chain's grid and the next. Raising it
        moves the chains further apart, which lowers the score a repulsion-based
        stand-in energy gives.

    See Also
    --------
    LeapBuilder : The builder that produces real structures.
    maws.energy.StubEnergy : The scorer intended to go with this.

    Examples
    --------
    >>> from maws.forcefield import ForceField
    >>> from maws.libraries import rna
    >>> from maws.topology import Assembly
    >>> system = FakeBuilder().build(
    ...     Assembly().with_aptamer(rna(), "G").with_ligand_stub(5),
    ...     ForceField.for_target("RNA", "protein"),
    ... )
    >>> system.n_atoms
    38
    >>> system.pose.xyz[0]
    array([0., 0., 0.])
    """

    __slots__ = ("_separation", "_spacing")

    def __init__(
        self,
        *,
        spacing: float = GRID_SPACING,
        separation: float = CHAIN_SEPARATION,
    ) -> None:
        self._spacing = float(spacing)
        self._separation = float(separation)

    def __repr__(self) -> str:
        return f"<FakeBuilder spacing={self._spacing} sep={self._separation}>"

    def prepare(self, chain: PdbChain, forcefield: ForceField) -> ResidueChain:
        r"""Count the atoms in a PDB file, without running anything.

        Parameters
        ----------
        chain : maws.topology.PdbChain
            The chain to resolve. Its file is read and its atom lines counted.
            A PDB file is plain text with one line per atom; those lines start
            with ``ATOM`` for part of a protein or nucleic acid, and ``HETATM``
            for anything else, such as a small organic molecule. Both count.
        forcefield : maws.forcefield.ForceField
            Ignored. Accepted so this can stand in for :class:`LeapBuilder`.

        Returns
        -------
        maws.topology.ResidueChain
            The same chain as a single residue with that many atoms.

        Raises
        ------
        maws.errors.ConfigurationError
            If the file cannot be read, or holds no atom records.

        Examples
        --------
        >>> import tempfile
        >>> from pathlib import Path
        >>> from maws.forcefield import ForceField
        >>> from maws.topology import PdbChain
        >>> with tempfile.TemporaryDirectory() as tmp:
        ...     pdb = Path(tmp) / "target.pdb"
        ...     _ = pdb.write_text("ATOM      1  N\nHETATM    2  C\nEND\n")
        ...     chain = PdbChain(
        ...         role="ligand",
        ...         path=pdb,
        ...         residue_name="LIG",
        ...         parameterized=True,
        ...     )
        ...     resolved = FakeBuilder().prepare(
        ...         chain, ForceField.for_target("RNA", "protein")
        ...     )
        ...     resolved.n_atoms
        2
        """
        try:
            text = chain.path.read_text()
        except OSError as exc:
            raise ConfigurationError(
                f"cannot read {chain.path} to count its atoms: {exc}"
            ) from exc
        n_atoms = sum(
            1 for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))
        )
        if n_atoms == 0:
            raise ConfigurationError(
                f"{chain.path} holds no ATOM or HETATM records, so there is "
                f"nothing to build"
            )
        return ResidueChain(
            role=chain.role,
            library=ResidueLibrary.single(chain.residue_name, n_atoms),
            sequence=NucleotideSequence((chain.residue_name,)),
        )

    def build(self, assembly: Assembly, forcefield: ForceField) -> BuiltSystem:
        """Build a grid structure for a design description.

        Parameters
        ----------
        assembly : maws.topology.Assembly
            What to build. Chains still described by a file are resolved first
            by :meth:`prepare`.
        forcefield : maws.forcefield.ForceField
            Recorded on the result, but not otherwise used: a grid has no
            physics to apply parameters to.

        Returns
        -------
        maws.topology.BuiltSystem
            A structure with one grid per chain, in the order the chains appear
            in the assembly. Every atom is carbon, and no bonds are recorded.

        Raises
        ------
        maws.errors.ConfigurationError
            If the assembly has no chains, or a chain still described by a file
            cannot be read.

        Examples
        --------
        The atom counts and chain spans are the real ones, so anything that
        works in terms of those behaves as it would on a real structure:

        >>> from maws.forcefield import ForceField
        >>> from maws.libraries import rna
        >>> from maws.topology import Assembly
        >>> system = FakeBuilder().build(
        ...     Assembly().with_aptamer(rna(), "G A").with_ligand_stub(10),
        ...     ForceField.for_target("RNA", "protein"),
        ... )
        >>> system.chain("aptamer").span
        AtomRange(start=0, stop=66)
        >>> system.chain("ligand").span
        AtomRange(start=66, stop=76)

        The two chains start 25 Å apart along x:

        >>> system.pose.xyz[66][0]
        25.0
        """
        chains = _resolved_chains(assembly, forcefield, self)
        if not chains:
            raise ConfigurationError("cannot build an assembly with no chains")

        blocks = [
            _grid(chain.n_atoms, index * self._separation, self._spacing)
            for index, chain in enumerate(chains)
        ]
        xyz = np.concatenate(blocks, axis=0)
        return BuiltSystem(
            Assembly(chains),
            forcefield,
            xyz,
            elements=[DEFAULT_STUB_ELEMENT] * len(xyz),
            masses=[DEFAULT_STUB_MASS] * len(xyz),
        )


def _grid(n_atoms: int, offset: float, spacing: float) -> np.ndarray:
    """Lay out `n_atoms` points on a cubic grid starting at `offset` along x.

    Parameters
    ----------
    n_atoms : int
        How many points to produce.
    offset : float
        How far along the x axis the grid starts, in ångström.
    spacing : float
        Distance in ångström between neighbouring points.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_atoms, 3)`` positions in ångström. No two rows are equal, so
        any pair of them defines a usable rotation axis.

    Examples
    --------
    >>> _grid(3, 0.0, 1.0)
    array([[0., 0., 0.],
           [1., 0., 0.],
           [0., 1., 0.]])
    """
    side = max(1, int(np.ceil(n_atoms ** (1 / 3))))
    index = np.arange(n_atoms)
    points = np.stack(
        [index % side, (index // side) % side, index // (side * side)], axis=1
    ).astype(np.float64)
    points *= spacing
    points[:, 0] += offset
    return points


class LeapBuilder:
    """A builder that produces real structures by running AmberTools.

    Built structures are cached. The same description always produces the same
    structure, so a name is worked out from what is being built and the result
    is stored under it; a later build that finds files already under that name
    reads them back instead of running ``tleap`` again.

    Parameters
    ----------
    cache_dir : pathlib.Path, default=Path(".maws_cache")
        Where built structures are stored. Delete this directory to force
        everything to be rebuilt; nothing else has to be done, since a missing
        entry is simply built again.
    params_dir : pathlib.Path, optional
        Where parameter files worked out for target molecules are written.
        Defaults to a ``params`` directory inside `cache_dir`, which keeps them
        out of the directory holding the input files.

    Raises
    ------
    maws.errors.ToolchainError
        From :meth:`build` and :meth:`prepare`, if AmberTools is not installed
        or one of its programs fails.
    maws.errors.BuildError
        From :meth:`build`, if ``tleap`` finishes but writes no structure.

    See Also
    --------
    FakeBuilder : Builds without AmberTools, for tests.
    maws.topology.BuiltSystem : What a build produces.

    Notes
    -----
    Caching is what makes a design run practical. Building a fifteen-residue
    aptamer takes seconds, and the search that grows one rebuilds constantly:
    every step tries each of the four nucleotides at each of the two ends, so
    each step is eight more structures. Reading a stored pair of files back
    instead turns those seconds into milliseconds.

    The name a structure is stored under covers the force field names, every
    chain's residue names in order, and the *contents* of every parameter file
    that will be loaded — not those files' paths.

    .. warning::
        Hashing the paths would not be enough. A parameter file is named after
        the residue it describes, so two different target molecules given the
        same residue name produce files at the same path, and the second
        design would silently be handed the first one's structure. Hashing
        what is inside them means a different target always gets a different
        name.

    Cache entries are never removed, so the directory grows by one pair of
    files per distinct structure built, until it is deleted by hand.

    Examples
    --------
    >>> from maws.forcefield import ForceField  # doctest: +SKIP
    >>> from maws.libraries import rna  # doctest: +SKIP
    >>> from maws.topology import Assembly  # doctest: +SKIP
    >>> builder = LeapBuilder("/tmp/demo_cache")  # doctest: +SKIP
    >>> system = builder.build(  # doctest: +SKIP
    ...     Assembly().with_aptamer(rna(), "G A U"),
    ...     ForceField.for_target("RNA", "protein"),
    ... )
    >>> system.n_atoms  # doctest: +SKIP
    99

    Asking again reads the stored files instead of running ``tleap``:

    >>> again = builder.build(  # doctest: +SKIP
    ...     Assembly().with_aptamer(rna(), "G A U"),
    ...     ForceField.for_target("RNA", "protein"),
    ... )
    >>> again.n_atoms  # doctest: +SKIP
    99
    """

    __slots__ = ("_cache_dir", "_params_dir", "_resources")

    def __init__(
        self,
        cache_dir: Path | str = DEFAULT_CACHE_DIR,
        *,
        params_dir: Path | str | None = None,
    ) -> None:
        self._cache_dir = Path(cache_dir)
        self._params_dir = (
            Path(params_dir) if params_dir is not None else self._cache_dir / "params"
        )
        self._resources: dict[str, Path] = {}

    def __repr__(self) -> str:
        return f"<LeapBuilder cache={self._cache_dir}>"

    @property
    def cache_dir(self) -> Path:
        """pathlib.Path : Where built structures are stored, once built."""
        return self._cache_dir

    @property
    def params_dir(self) -> Path:
        """pathlib.Path : Where the ``.lib`` and ``.frcmod`` files are put."""
        return self._params_dir

    def prepare(self, chain: PdbChain, forcefield: ForceField) -> ResidueChain:
        """Work out parameters for a molecule read from a file.

        Runs ``antechamber`` and ``parmchk2`` to derive parameters when the
        force field does not already describe the molecule, then ``tleap`` to
        write a ``.lib`` file: a residue definition recording the molecule's
        atoms, their names and charges, and the bonds between them. That file
        stays in :attr:`params_dir` and is loaded again by every later build.

        However large the molecule is, it becomes one residue. MAWS never grows
        the target or turns any of its bonds, so there is nothing to gain from
        describing it in finer detail.

        Parameters
        ----------
        chain : maws.topology.PdbChain
            The chain to resolve, naming the file to read and the residue name
            to give it.
        forcefield : maws.forcefield.ForceField
            Which parameters to use. Its ``parameterized`` flag decides whether
            ``antechamber`` and ``parmchk2`` need to run at all.

        Returns
        -------
        maws.topology.ResidueChain
            The same chain as a single residue with a known atom count.

        Raises
        ------
        maws.errors.ToolchainError
            If AmberTools is not installed, or one of its programs fails.

        See Also
        --------
        maws.io.prepare.make_lib : Does the work of running the three programs.

        Examples
        --------
        >>> from maws.forcefield import ForceField  # doctest: +SKIP
        >>> from maws.topology import PdbChain  # doctest: +SKIP
        >>> chain = PdbChain(  # doctest: +SKIP
        ...     role="ligand",
        ...     path=Path("target.pdb"),
        ...     residue_name="LIG",
        ...     parameterized=True,
        ... )
        >>> resolved = LeapBuilder().prepare(  # doctest: +SKIP
        ...     chain, ForceField.for_target("RNA", "protein")
        ... )
        >>> resolved.n_atoms  # doctest: +SKIP
        1462
        """
        from maws.io.prepare import make_lib

        self._params_dir.mkdir(parents=True, exist_ok=True)
        n_atoms = make_lib(
            chain.path,
            chain.residue_name,
            output_dir=self._params_dir,
            force_field_aptamer=forcefield.aptamer_source,
            force_field_ligand=forcefield.ligand_source,
            parameterized=chain.parameterized,
            net_charge=chain.net_charge,
        )
        self._resources[chain.residue_name] = self._params_dir
        return ResidueChain(
            role=chain.role,
            library=ResidueLibrary.single(chain.residue_name, n_atoms),
            sequence=NucleotideSequence((chain.residue_name,)),
        )

    def build(self, assembly: Assembly, forcefield: ForceField) -> BuiltSystem:
        """Build a real structure, reusing a cached one where possible.

        Parameters
        ----------
        assembly : maws.topology.Assembly
            What to build. Chains still described by a file are resolved first
            by :meth:`prepare`.
        forcefield : maws.forcefield.ForceField
            Which parameters to build with. Two builds of the same chains under
            different force fields are cached separately.

        Returns
        -------
        maws.topology.BuiltSystem
            The built structure: positions in ångström, one element symbol and
            one mass in daltons per atom, all read back from the files
            ``tleap`` wrote.

        Raises
        ------
        maws.errors.ConfigurationError
            If the assembly has no chains with any residues in them.
        maws.errors.ToolchainError
            If ``tleap`` is not installed, or exits with an error.
        maws.errors.BuildError
            If ``tleap`` finishes but writes no usable structure.

        See Also
        --------
        FakeBuilder.build : The same call, without needing anything installed.

        Examples
        --------
        >>> from maws.forcefield import ForceField  # doctest: +SKIP
        >>> from maws.libraries import rna  # doctest: +SKIP
        >>> from maws.topology import Assembly  # doctest: +SKIP
        >>> system = LeapBuilder().build(  # doctest: +SKIP
        ...     Assembly().with_aptamer(rna(), "G"),
        ...     ForceField.for_target("RNA", "protein"),
        ... )
        >>> system.chain("aptamer").span  # doctest: +SKIP
        AtomRange(start=0, stop=33)
        """
        chains = _resolved_chains(assembly, forcefield, self)
        if not any(chain.sequence for chain in chains):
            raise ConfigurationError(
                "cannot build: no chain in this assembly has any residues"
            )

        key = self._cache_key(chains, forcefield)
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        prmtop = self._cache_dir / f"{key}.prmtop"
        inpcrd = self._cache_dir / f"{key}.inpcrd"
        if not (prmtop.exists() and inpcrd.exists()):
            self._run_leap(chains, forcefield, prmtop, inpcrd)

        return self._load(Assembly(chains), forcefield, prmtop, inpcrd)

    # -- internals --------------------------------------------------------

    def _leap_script(
        self,
        chains: Sequence[ResidueChain],
        forcefield: ForceField,
        prmtop: Path,
        inpcrd: Path,
    ) -> str:
        """Return the input script that tells ``tleap`` what to build.

        ``tleap`` is driven by a small command language rather than by
        arguments, so the whole build is expressed as a handful of lines: load
        the parameters, load any residue definitions worked out for target
        molecules, spell out each chain as a sequence of residue names, join
        them into one molecule, and save it.

        Parameters
        ----------
        chains : sequence of maws.topology.ResidueChain
            The chains to build, in order. A chain with no residues in it is
            skipped, since ``tleap`` has no way to write an empty sequence.
        forcefield : maws.forcefield.ForceField
            Supplies the two lines that load the parameters.
        prmtop, inpcrd : pathlib.Path
            Where ``tleap`` should write the force field numbers and the
            coordinates.

        Returns
        -------
        str
            The script, as newline-separated lines.
        """
        lines = list(forcefield.leap_preamble())
        for name, directory in sorted(self._resources.items()):
            lines.append(f"loadoff {directory / f'{name}.lib'}")
            frcmod = directory / f"{name}.frcmod"
            if frcmod.exists():
                lines.append(f"loadamberparams {frcmod}")

        built = []
        for index, chain in enumerate(chains):
            if not chain.sequence:
                continue
            label = f"CHAIN{index}"
            residues = " ".join(chain.canonical)
            lines.append(f"{label} = sequence {{{residues}}}")
            built.append(label)

        lines.append(f"UNION = combine {{{' '.join(built)}}}")
        lines.append(f"saveamberparm UNION {prmtop} {inpcrd}")
        lines.append("quit")
        return "\n".join(lines)

    def _run_leap(
        self,
        chains: Sequence[ResidueChain],
        forcefield: ForceField,
        prmtop: Path,
        inpcrd: Path,
    ) -> None:
        """Run ``tleap`` and check that it produced a structure.

        Parameters
        ----------
        chains : sequence of maws.topology.ResidueChain
            The chains to build.
        forcefield : maws.forcefield.ForceField
            Which parameters to build with.
        prmtop, inpcrd : pathlib.Path
            Where the results should end up.

        Raises
        ------
        maws.errors.ToolchainError
            If ``tleap`` is missing or exits with an error.
        maws.errors.BuildError
            If it exits cleanly but writes nothing usable. ``tleap`` reports
            many problems by printing a complaint and carrying on, so a clean
            exit is not on its own evidence that anything was built.
        """
        # tleap scatters intermediate files beside wherever it runs, so it runs
        # in a directory that is thrown away rather than in the cache.
        with tempfile.TemporaryDirectory() as workdir:
            script = Path(workdir) / "build.in"
            script.write_text(self._leap_script(chains, forcefield, prmtop, inpcrd))
            run_tool(["tleap", "-f", str(script)], cwd=workdir)

        if not (prmtop.exists() and inpcrd.exists()) or prmtop.stat().st_size == 0:
            raise BuildError(
                f"tleap exited cleanly but wrote no usable structure to "
                f"{prmtop}. Its own output above normally says which residue "
                f"or parameter it could not find."
            )

    def _load(
        self,
        assembly: Assembly,
        forcefield: ForceField,
        prmtop: Path,
        inpcrd: Path,
    ) -> BuiltSystem:
        """Read a built structure back off disk.

        Parameters
        ----------
        assembly : maws.topology.Assembly
            The description that was built.
        forcefield : maws.forcefield.ForceField
            Which parameters it was built with.
        prmtop, inpcrd : pathlib.Path
            The two files ``tleap`` wrote: force field numbers and positions.

        Returns
        -------
        maws.topology.BuiltSystem
            The structure, with positions in ångström and masses in daltons.
            Both files are kept on the result, so a later scoring step can hand
            them to OpenMM without rebuilding anything.
        """
        from openmm import app, unit

        parameters = app.AmberPrmtopFile(str(prmtop))
        coordinates = app.AmberInpcrdFile(str(inpcrd))
        atoms = list(parameters.topology.atoms())
        return BuiltSystem(
            assembly,
            forcefield,
            np.asarray(
                coordinates.positions.value_in_unit(unit.angstrom), dtype=np.float64
            ),
            elements=[atom.element.symbol for atom in atoms],
            masses=[atom.element.mass.value_in_unit(unit.dalton) for atom in atoms],
            amber=AmberArtifacts(
                prmtop_path=prmtop,
                inpcrd_path=inpcrd,
                topology=parameters.topology,
            ),
        )

    def _cache_key(self, chains: Sequence[ResidueChain], forcefield: ForceField) -> str:
        """Work out the name a built structure is stored under.

        Parameters
        ----------
        chains : sequence of maws.topology.ResidueChain
            The chains being built.
        forcefield : maws.forcefield.ForceField
            The parameters being built with.

        Returns
        -------
        str
            A 40-character hexadecimal digest, used as the filename stem for
            the pair of files this build produces.

        Notes
        -----
        The digest covers the force field names, each chain's residue names in
        order, and the *contents* of every parameter file that will be loaded.
        Anything that would change the structure changes the digest, and
        anything that would not, does not.

        .. warning::
            Hashing parameter file paths instead of their contents is not
            enough. Two different target molecules can produce files with the
            same name in the same place, and the two designs would then share
            one cache entry.
        """
        payload = {
            "aptamer": forcefield.aptamer_source,
            "ligand": forcefield.ligand_source,
            "chains": [list(chain.canonical) for chain in chains],
            "params": sorted(
                [
                    name,
                    _digest(directory / f"{name}.lib"),
                    _digest(directory / f"{name}.frcmod"),
                ]
                for name, directory in self._resources.items()
            ),
        }
        encoded = json.dumps(payload, sort_keys=True).encode()
        return hashlib.sha1(encoded).hexdigest()


def _digest(path: Path) -> str:
    """Return a digest of a file's contents, or ``""`` if it is not there.

    Parameters
    ----------
    path : pathlib.Path
        The file to read.

    Returns
    -------
    str
        A 40-character hexadecimal digest, or the empty string when the file
        does not exist. A ``.frcmod`` file is only written for a molecule whose
        parameters had to be worked out, so its absence is ordinary and is
        recorded rather than raised on.

    Examples
    --------
    >>> _digest(Path("no-such-file"))
    ''
    """
    if not path.exists():
        return ""
    return hashlib.sha1(path.read_bytes()).hexdigest()

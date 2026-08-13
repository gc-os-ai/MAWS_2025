"""
maws.values
===========

The plain data types the rest of MAWS is written in terms of.

MAWS designs an *aptamer*: a short strand of RNA or DNA that folds up against a
target molecule and sticks to it. A strand is a sequence of *residues* — single
nucleotides — joined end to end. The two ends of a strand are chemically
different and are conventionally called 5' and 3'; a sequence is written
5'-to-3', left to right.

Atoms in a design are held in one flat array covering every chain at once, so an
atom is identified by a single integer index into that array. The types here
describe pieces of that array: runs of it (:class:`AtomRange`), bonds within it
that can be turned (:class:`Torsion`), and the per-residue recipes that say
which bonds those are (:class:`ResidueTemplate`).

Nothing in this module reads a file, runs a program, or imports a simulation
package. It is all frozen dataclasses built from ``int``, ``float``, ``str`` and
``tuple``, which means anything here can be constructed by hand in a test.

Two coordinate frames
---------------------
A bare atom index means nothing until you know what it counts from. There are
two possibilities and a different type carries each:

- **Residue-local** — counted from the first atom of one residue. Residue
  templates are written this way, because a template describes a residue type
  and does not know where any particular copy of it will end up.
  :class:`BackboneTorsion` and :class:`LocalTorsion` hold residue-local indices.

- **Global** — counted from the first atom of the whole design.
  :class:`AtomRange` and :class:`Torsion` hold global indices.

:meth:`maws.pose.ResidueView.torsion` is the single place that converts from the
first to the second.

Examples
--------
>>> from maws.libraries import rna
>>> seq = NucleotideSequence.parse("G A U")
>>> seq.canonical(rna())
('G5', 'A', 'U3')
>>> str(seq.appended("C"))
'G A U C'
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

from maws.errors import ConfigurationError

__all__ = [
    "AliasSet",
    "AtomRange",
    "BackboneTorsion",
    "Connection",
    "Direction",
    "LocalTorsion",
    "NucleotideSequence",
    "ResidueLibrary",
    "ResidueTemplate",
    "Torsion",
    "TorsionTemplate",
]

Direction = Literal["3prime", "5prime"]
"""Which of a strand's two ends is meant. See the module docstring."""


# ---------------------------------------------------------------------------
# Atom spans
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AtomRange:
    """A run of consecutive atom indices, from `start` up to but not `stop`.

    The same half-open convention as a Python slice: ``AtomRange(3, 7)`` covers
    atoms 3, 4, 5 and 6, and has length 4.

    Parameters
    ----------
    start : int
        Index of the first atom in the run.
    stop : int
        Index one past the last atom in the run. Equal to `start` for an empty
        run.

    Raises
    ------
    maws.errors.ConfigurationError
        If `start` is negative, or `stop` comes before `start`.

    See Also
    --------
    Torsion : Pairs a range of moving atoms with the bond they turn about.

    Examples
    --------
    >>> span = AtomRange(10, 25)
    >>> len(span)
    15
    >>> 12 in span
    True
    >>> span.shifted(100)
    AtomRange(start=110, stop=125)
    """

    start: int
    stop: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ConfigurationError(
                f"atom indices count from zero and are never negative: "
                f"AtomRange({self.start}, {self.stop})"
            )
        if self.stop < self.start:
            raise ConfigurationError(
                f"a range cannot end before it begins: "
                f"AtomRange({self.start}, {self.stop})"
            )

    def __len__(self) -> int:
        return self.stop - self.start

    def __contains__(self, index: int) -> bool:
        return self.start <= index < self.stop

    def __iter__(self) -> Iterator[int]:
        return iter(range(self.start, self.stop))

    def shifted(self, offset: int) -> AtomRange:
        """Return the same run of atoms, moved `offset` places along the array.

        Parameters
        ----------
        offset : int
            How many indices to add to both ends. Typically the position at
            which a residue or chain begins.

        Returns
        -------
        AtomRange
            A new range. The one this was called on does not change.
        """
        return AtomRange(self.start + offset, self.stop + offset)

    def as_slice(self) -> slice:
        """Return this run as a :class:`slice`, ready to index a NumPy array.

        Returns
        -------
        slice
            ``slice(start, stop)``.

        Notes
        -----
        Indexing a NumPy array with a slice gives a *view* onto the original
        data and costs nothing. Indexing it with a list or a range copies the
        selected rows instead. Rotations run millions of times in one design,
        so this method exists to make sure the cheap form is the one that gets
        used.

        Examples
        --------
        >>> import numpy as np
        >>> xyz = np.zeros((10, 3))
        >>> xyz[AtomRange(2, 5).as_slice()].shape
        (3, 3)
        """
        return slice(self.start, self.stop)


# ---------------------------------------------------------------------------
# Torsions
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Torsion:
    """A single bond that can be turned, and the atoms that turn with it.

    A torsion is a rotatable bond. Two atoms define the bond, and therefore the
    line to rotate about. Turning it swings one part of the molecule around
    that line while the rest stays put; no bond length or bond angle changes,
    only the shape.

    All three fields use global atom indices, counted from the first atom of
    the whole design.

    Parameters
    ----------
    pivot : int
        One of the two atoms of the bond. The rotation is centred here, so this
        atom does not move even if it is listed in `moving`.
    bond : int
        The other atom of the bond. Together with `pivot` it fixes the line
        that everything turns about.
    moving : AtomRange
        The atoms that swing. Every atom outside this range keeps its position.

    See Also
    --------
    BackboneTorsion, LocalTorsion : Recipes that produce these.
    maws.pose.Pose.rotate : Applies one.

    Examples
    --------
    >>> t = Torsion(pivot=0, bond=3, moving=AtomRange(3, 34))
    >>> t.shifted(1000)
    Torsion(pivot=1000, bond=1003, moving=AtomRange(start=1003, stop=1034))
    """

    pivot: int
    bond: int
    moving: AtomRange

    def shifted(self, offset: int) -> Torsion:
        """Return the same torsion, moved `offset` places along the array.

        Parameters
        ----------
        offset : int
            How many indices to add to every atom index in the torsion.

        Returns
        -------
        Torsion
            A new torsion. The one this was called on does not change.
        """
        return Torsion(
            pivot=self.pivot + offset,
            bond=self.bond + offset,
            moving=self.moving.shifted(offset),
        )


@dataclass(frozen=True, slots=True)
class _TorsionTemplateBase:
    """Shared behaviour of the two kinds of torsion recipe.

    Not part of the public interface. Use :class:`BackboneTorsion` or
    :class:`LocalTorsion`, and see either for what the two indices mean.

    Parameters
    ----------
    pivot : int
        Residue-local index of one atom of the bond.
    bond : int
        Residue-local index of the other atom of the bond.
    """

    pivot: int
    bond: int

    def placed(self, offset: int, chain: AtomRange) -> Torsion:
        """Turn this recipe into a usable :class:`Torsion`, 3' side moving.

        Parameters
        ----------
        offset : int
            Global index of the first atom of the residue this recipe belongs
            to.
        chain : AtomRange
            Global span of the whole chain that residue sits in.

        Returns
        -------
        Torsion
            The same bond, in global atom indices.
        """
        raise NotImplementedError

    def reversed(self, offset: int, chain: AtomRange) -> Torsion:
        """Turn this recipe into a usable :class:`Torsion`, 5' side moving.

        Parameters
        ----------
        offset : int
            Global index of the first atom of the residue this recipe belongs
            to.
        chain : AtomRange
            Global span of the whole chain that residue sits in.

        Returns
        -------
        Torsion
            A torsion about the same bond, whose moving atoms are everything
            from the start of the chain up to that bond.

        Notes
        -----
        A bond has two sides, and turning either one changes the molecule's
        shape by the same amount — only which part of it stays still differs.
        So "rotate the 5' side" is a second torsion about the same bond, not a
        flag on the first one. Its axis runs the other way, so that the same
        angle produces the same change in shape either way round.

        .. warning::
            The two sides are complementary only when the bond separates the
            strand into exactly two pieces. That holds for a
            :class:`BackboneTorsion` and not for a :class:`LocalTorsion`, which
            overrides this.
        """
        bond = self.bond + offset
        return Torsion(
            pivot=bond,
            bond=self.pivot + offset,
            moving=AtomRange(chain.start, bond),
        )


@dataclass(frozen=True, slots=True)
class BackboneTorsion(_TorsionTemplateBase):
    """A recipe for a bond whose 3' side is the entire rest of the strand.

    These are the bonds along the backbone, the chain of atoms that runs
    through every residue and holds the strand together. Turning one of them
    swings every residue past it, so it changes the overall fold rather than
    the shape of a single residue.

    Both indices are residue-local: counted from the first atom of the residue
    this recipe describes.

    Parameters
    ----------
    pivot : int
        One atom of the bond.
    bond : int
        The other atom of the bond. Every atom from here to the end of the
        strand moves.

    See Also
    --------
    LocalTorsion : The other kind, bounded inside one residue.

    Examples
    --------
    A residue starting at atom 100, in a chain running to atom 300:

    >>> BackboneTorsion(pivot=0, bond=3).placed(100, AtomRange(100, 300))
    Torsion(pivot=100, bond=103, moving=AtomRange(start=103, stop=300))
    """

    def placed(self, offset: int, chain: AtomRange) -> Torsion:
        """Turn this recipe into a usable :class:`Torsion`, 3' side moving.

        Parameters
        ----------
        offset : int
            Global index of the first atom of the residue this recipe belongs
            to.
        chain : AtomRange
            Global span of the whole chain that residue sits in. The moving
            atoms run to the end of it.

        Returns
        -------
        Torsion
            The same bond, in global atom indices.
        """
        bond = self.bond + offset
        return Torsion(
            pivot=self.pivot + offset,
            bond=bond,
            moving=AtomRange(bond, chain.stop),
        )


@dataclass(frozen=True, slots=True)
class LocalTorsion(_TorsionTemplateBase):
    """A recipe for a bond that moves a fixed group of atoms in one residue.

    These are the bonds that reshape the sugar ring and swing the base of a
    single nucleotide. The group of atoms that moves is known in advance and
    lies entirely inside the residue, so turning one leaves the rest of the
    strand untouched.

    All three indices are residue-local: counted from the first atom of the
    residue this recipe describes.

    Parameters
    ----------
    pivot : int
        One atom of the bond.
    bond : int
        The other atom of the bond, and the first atom that moves.
    stop : int
        One past the last atom that moves.

    See Also
    --------
    BackboneTorsion : The other kind, which moves the rest of the strand.

    Examples
    --------
    >>> LocalTorsion(pivot=8, bond=10, stop=25).placed(100, AtomRange(100, 300))
    Torsion(pivot=108, bond=110, moving=AtomRange(start=110, stop=125))
    """

    stop: int

    def placed(self, offset: int, chain: AtomRange) -> Torsion:
        """Turn this recipe into a usable :class:`Torsion`, 3' side moving.

        Parameters
        ----------
        offset : int
            Global index of the first atom of the residue this recipe belongs
            to.
        chain : AtomRange
            Global span of the whole chain that residue sits in. Not needed
            here, since the moving atoms are bounded by `stop`; accepted so
            that both kinds of recipe can be called the same way.

        Returns
        -------
        Torsion
            The same bond, in global atom indices.
        """
        bond = self.bond + offset
        return Torsion(
            pivot=self.pivot + offset,
            bond=bond,
            moving=AtomRange(bond, self.stop + offset),
        )

    def reversed(self, offset: int, chain: AtomRange) -> Torsion:
        """Return the same torsion as :meth:`placed`.

        Parameters
        ----------
        offset : int
            Global index of the first atom of the residue this recipe belongs
            to.
        chain : AtomRange
            Global span of the whole chain that residue sits in.

        Returns
        -------
        Torsion
            Exactly what :meth:`placed` returns.

        Notes
        -----
        A bond only has two well-defined sides when cutting it would fall the
        molecule into two pieces. This kind of bond does not: its moving atoms
        stop at `stop`, and the atoms past `stop` are attached to the rest of
        the residue by a second path. Turning "everything before the bond"
        would move atoms on one side of that second path and not the other,
        which stretches a bond and destroys the molecule.

        So there is no 5' form to give, and the 3' form is returned instead.
        That is the right answer for what the 5' form is asked for: growing the
        strand at its 5' end, where the rest of the strand is already
        positioned against the target and only the new residue should move.
        This bond moves a group inside one residue whichever way it is read, so
        it already leaves the rest of the strand alone.
        """
        return self.placed(offset, chain)


TorsionTemplate = BackboneTorsion | LocalTorsion
"""Either kind of torsion recipe. Tell them apart with ``match`` or
``isinstance``."""


# ---------------------------------------------------------------------------
# Residues
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Connection:
    """Which atoms form the bond between one residue and its neighbour.

    Residues are joined into a strand by a single bond on each side. This
    records one of those two joins: one atom from this residue, one from the
    residue next to it, and how far apart the two should sit.

    Parameters
    ----------
    own_atom : int
        Index of the bonding atom within this residue, counted from its first
        atom.
    other_atom : int
        Index of the bonding atom within the neighbouring residue. May be
        negative, meaning "counted back from that residue's last atom" — useful
        because how many atoms the neighbour has is not known in advance.
    length : float
        How far apart the two atoms should sit, in ångström.
    """

    own_atom: int
    other_atom: int
    length: float


@dataclass(frozen=True, slots=True)
class AliasSet:
    """The four names one written nucleotide can have, by where it sits.

    A nucleotide is written as a single token, such as ``"G"``. The molecular
    modelling program that builds the structure needs a more specific name,
    because the same nucleotide is chemically different at the two ends of a
    strand than it is in the middle: an end residue carries an extra group that
    a middle one does not. So one written token maps to four possible names,
    and which one applies depends only on position.

    Parameters
    ----------
    alone : str
        Name to use when the strand is one residue long.
    start : str
        Name to use at the 5' end of a longer strand.
    middle : str
        Name to use anywhere that is not an end.
    end : str
        Name to use at the 3' end.

    See Also
    --------
    NucleotideSequence.canonical : Applies these rules to a whole sequence.

    Examples
    --------
    >>> from maws.libraries import rna
    >>> rna().alias("G")
    AliasSet(alone='GN', start='G5', middle='G', end='G3')
    """

    alone: str
    start: str
    middle: str
    end: str


@dataclass(frozen=True, slots=True)
class ResidueTemplate:
    """Everything about one type of residue: its size and its moving parts.

    A template describes a *kind* of residue, not any particular copy of one.
    So its atom indices are all counted from the residue's own first atom, and
    a copy placed somewhere in a strand is described by shifting them.

    Parameters
    ----------
    name : str
        The name the modelling program knows this residue by, e.g. ``"G5"``.
    n_atoms : int
        How many atoms the residue has.
    torsions : tuple of BackboneTorsion or LocalTorsion
        The bonds in this residue that can be turned, in a fixed order. Code
        that samples shapes refers to them by position in this tuple.
    head : Connection
        The bond joining this residue to its neighbour on the 5' side.
    tail : Connection
        The bond joining this residue to its neighbour on the 3' side.

    Raises
    ------
    maws.errors.ConfigurationError
        If `n_atoms` is not positive, or a torsion names an atom the residue
        does not have. Both are caught here, when the template is created,
        rather than later when something tries to use it.

    Examples
    --------
    >>> template = ResidueTemplate(
    ...     name="G",
    ...     n_atoms=34,
    ...     torsions=(BackboneTorsion(0, 3),),
    ...     head=Connection(0, -1, 1.6),
    ...     tail=Connection(-2, 0, 1.6),
    ... )
    >>> template.n_torsions
    1
    """

    name: str
    n_atoms: int
    torsions: tuple[TorsionTemplate, ...]
    head: Connection
    tail: Connection

    def __post_init__(self) -> None:
        if self.n_atoms <= 0:
            raise ConfigurationError(
                f"{self.name}: a residue must have at least one atom, "
                f"got n_atoms={self.n_atoms}"
            )
        for torsion in self.torsions:
            self._check_index(torsion.pivot, torsion)
            self._check_index(torsion.bond, torsion)
            if isinstance(torsion, LocalTorsion):
                self._check_index(torsion.stop, torsion, inclusive=True)
                if torsion.stop <= torsion.bond:
                    raise ConfigurationError(
                        f"{self.name}: torsion {torsion} would move no atoms, "
                        f"because stop ({torsion.stop}) is not past bond "
                        f"({torsion.bond})"
                    )

    def _check_index(
        self, index: int, torsion: TorsionTemplate, *, inclusive: bool = False
    ) -> None:
        """Raise if `index` names an atom this residue does not have.

        Parameters
        ----------
        index : int
            Residue-local atom index to check.
        torsion : BackboneTorsion or LocalTorsion
            The recipe the index came from. Quoted in the error message so the
            reader can find it.
        inclusive : bool, default=False
            Whether ``index == n_atoms`` is allowed. True for a `stop` bound,
            which points one past the last atom.

        Raises
        ------
        maws.errors.ConfigurationError
            If the index is out of range.
        """
        limit = self.n_atoms + 1 if inclusive else self.n_atoms
        if not 0 <= index < limit:
            raise ConfigurationError(
                f"{self.name}: torsion {torsion} names atom {index}, but the "
                f"residue only has atoms 0 to {self.n_atoms - 1}"
            )

    @property
    def n_torsions(self) -> int:
        """int : How many turnable bonds this residue has."""
        return len(self.torsions)


class ResidueLibrary(Mapping[str, ResidueTemplate]):
    """Every residue type MAWS can build with, plus how to name them.

    A library answers two questions. Given a modelling-program residue name
    such as ``"G5"``, it gives the :class:`ResidueTemplate` describing that
    residue — this is the mapping half, so ``lib["G5"]``, ``"G5" in lib``,
    ``len(lib)`` and iteration all work. And given a written token such as
    ``"G"``, :meth:`alias` gives the four names that token can take depending
    on where in a strand it sits.

    The two vocabularies are not the same. DNA sequences are written with the
    tokens ``G A T C``, while its residues are named ``DGN``, ``DA5``, ``DT``
    and so on.

    Parameters
    ----------
    templates : mapping of str to ResidueTemplate
        Residue descriptions, keyed by modelling-program name. Copied on the
        way in, so later changes to the argument have no effect.
    aliases : mapping of str to AliasSet, optional
        Written token to its four position-dependent names. Any residue name
        without an entry is taken to be its own token.

    See Also
    --------
    maws.libraries.rna, maws.libraries.dna : The two built-in libraries.

    Examples
    --------
    >>> lib = ResidueLibrary.single("LIG", n_atoms=40)
    >>> sorted(lib)
    ['LIG']
    >>> lib["LIG"].n_atoms
    40
    """

    __slots__ = ("_aliases", "_templates")

    def __init__(
        self,
        templates: Mapping[str, ResidueTemplate],
        aliases: Mapping[str, AliasSet] | None = None,
    ) -> None:
        self._templates = dict(templates)
        self._aliases = {
            name: AliasSet(name, name, name, name) for name in self._templates
        }
        if aliases:
            self._aliases.update(aliases)

    def __getitem__(self, name: str) -> ResidueTemplate:
        try:
            return self._templates[name]
        except KeyError:
            raise ConfigurationError(
                f"residue {name!r} is not in this library. It knows: "
                f"{', '.join(sorted(self._templates))}"
            ) from None

    def __iter__(self) -> Iterator[str]:
        return iter(self._templates)

    def __len__(self) -> int:
        return len(self._templates)

    def __repr__(self) -> str:
        return (
            f"<ResidueLibrary {len(self._templates)} residues, "
            f"{len(self._aliases)} tokens>"
        )

    def alias(self, token: str) -> AliasSet:
        """Return the four names a written token can take.

        Parameters
        ----------
        token : str
            A nucleotide as a user would write it, e.g. ``"G"``.

        Returns
        -------
        AliasSet
            The names for a lone residue, a 5' end, a middle position, and a
            3' end.

        Raises
        ------
        maws.errors.ConfigurationError
            If the library does not know this token.

        Examples
        --------
        >>> from maws.libraries import rna
        >>> rna().alias("G").start
        'G5'
        """
        try:
            return self._aliases[token]
        except KeyError:
            raise ConfigurationError(
                f"token {token!r} is not in this library. It knows: "
                f"{', '.join(sorted(self._aliases))}"
            ) from None

    @property
    def tokens(self) -> tuple[str, ...]:
        """tuple of str : Every token this library can resolve, sorted."""
        return tuple(sorted(self._aliases))

    def atom_count(self, names: Sequence[str]) -> int:
        """Return how many atoms a run of residues adds up to.

        Parameters
        ----------
        names : sequence of str
            Modelling-program residue names, in the order they appear in the
            strand.

        Returns
        -------
        int
            The total.

        Examples
        --------
        >>> lib = ResidueLibrary.single("LIG", n_atoms=40)
        >>> lib.atom_count(["LIG", "LIG"])
        80
        """
        return sum(self[name].n_atoms for name in names)

    @classmethod
    def single(
        cls,
        name: str,
        n_atoms: int,
        *,
        bond_length: float = 1.6,
    ) -> ResidueLibrary:
        """Build a library holding one residue and nothing else.

        This is the shape a target molecule takes. However large it is, MAWS
        treats it as a single rigid residue: it is never grown, and none of its
        bonds are turned.

        Parameters
        ----------
        name : str
            Residue name, matching the parameter file the modelling program
            will load for it.
        n_atoms : int
            How many atoms the molecule has.
        bond_length : float, default=1.6
            Nominal join distance in ångström. A one-residue strand has nothing
            to join to, so this only fills in the record.

        Returns
        -------
        ResidueLibrary
            A library containing exactly this residue, with no turnable bonds.
        """
        template = ResidueTemplate(
            name=name,
            n_atoms=n_atoms,
            torsions=(),
            head=Connection(0, -1, bond_length),
            tail=Connection(-2, 0, bond_length),
        )
        return cls({name: template})

    @classmethod
    def from_tables(
        cls,
        *,
        names: Sequence[str],
        lengths: Sequence[int],
        rotations: Sequence[tuple[str, int, int, int | None]],
        connect: Sequence[Sequence[object]],
        aliases: Sequence[Sequence[str]],
    ) -> ResidueLibrary:
        """Build a library from the flat tables in :mod:`maws.libraries`.

        The tables are written as parallel lists, some keyed by position and
        some by name. This method is where that ends: it checks they line up,
        resolves every negative index, and hands back one validated object. A
        table with a row missing fails here, with a message naming the problem,
        instead of much later with a complaint about a missing parameter.

        Parameters
        ----------
        names : sequence of str
            Residue names as the modelling program knows them.
        lengths : sequence of int
            Atom count per residue. Entry *i* describes ``names[i]``.
        rotations : sequence of tuple
            Rows of ``(residue, pivot, bond, stop_or_None)`` describing turnable
            bonds. A `stop_or_None` of ``None`` means the bond moves the rest of
            the strand and produces a :class:`BackboneTorsion`; a number bounds
            the moving atoms inside the residue and produces a
            :class:`LocalTorsion`. Negative indices count back from the
            residue's last atom and are resolved here.
        connect : sequence of sequence
            One row per residue, entry *i* describing ``names[i]``, shaped
            ``[[own_head, other_head], [own_tail, other_tail], head_length,
            tail_length]``.
        aliases : sequence of sequence of str
            Rows of ``[token, alone, start, middle, end]``. A token need not be
            a residue name.

        Returns
        -------
        ResidueLibrary
            A validated library.

        Raises
        ------
        maws.errors.ConfigurationError
            If the tables have different lengths, or a rotation row names a
            residue that has no atom count.

        Examples
        --------
        >>> lib = ResidueLibrary.from_tables(
        ...     names=["G", "A"],
        ...     lengths=[34, 33],
        ...     rotations=[("G", 0, 3, None), ("A", 0, 3, -7)],
        ...     connect=[[[0, -1], [-2, 0], 1.6, 1.6]] * 2,
        ...     aliases=[["G", "GN", "G5", "G", "G3"]],
        ... )
        >>> type(lib["G"].torsions[0]).__name__
        'BackboneTorsion'
        >>> type(lib["A"].torsions[0]).__name__
        'LocalTorsion'
        """
        if len(names) != len(lengths):
            raise ConfigurationError(
                f"residue tables do not line up: {len(names)} names but "
                f"{len(lengths)} atom counts. The shorter table stops at "
                f"index {min(len(names), len(lengths))}."
            )
        if len(names) != len(connect):
            raise ConfigurationError(
                f"residue tables do not line up: {len(names)} names but "
                f"{len(connect)} connectivity rows"
            )

        n_atoms_of = dict(zip(names, (int(n) for n in lengths), strict=True))
        torsions = _group_rotations(rotations, n_atoms_of)

        templates = {}
        for index, name in enumerate(names):
            head_pair, tail_pair, head_len, tail_len = connect[index]
            templates[name] = ResidueTemplate(
                name=name,
                n_atoms=n_atoms_of[name],
                torsions=torsions.get(name, ()),
                head=Connection(
                    own_atom=int(head_pair[0]),  # type: ignore[index]
                    other_atom=int(head_pair[1]),  # type: ignore[index]
                    length=float(head_len),  # type: ignore[arg-type]
                ),
                tail=Connection(
                    own_atom=int(tail_pair[0]),  # type: ignore[index]
                    other_atom=int(tail_pair[1]),  # type: ignore[index]
                    length=float(tail_len),  # type: ignore[arg-type]
                ),
            )
        return cls(templates, _parse_aliases(aliases))


def _resolve(index: int, n_atoms: int) -> int:
    """Turn an atom index that may count backwards into one that counts forwards.

    Parameters
    ----------
    index : int
        An index as written in a residue table. A negative value counts back
        from the residue's last atom, in the same way as Python list indexing.
    n_atoms : int
        How many atoms the residue has.

    Returns
    -------
    int
        The same atom, counted forwards from zero.

    Examples
    --------
    >>> _resolve(-7, 34)
    27
    >>> _resolve(3, 34)
    3
    """
    return index + n_atoms if index < 0 else index


def _group_rotations(
    rotations: Sequence[tuple[str, int, int, int | None]],
    n_atoms_of: Mapping[str, int],
) -> dict[str, tuple[TorsionTemplate, ...]]:
    """Sort flat rotation rows into one tuple of recipes per residue.

    Parameters
    ----------
    rotations : sequence of tuple
        Rows of ``(residue, pivot, bond, stop_or_None)``.
    n_atoms_of : mapping of str to int
        Atom count per residue name, needed to resolve negative indices.

    Returns
    -------
    dict of str to tuple of BackboneTorsion or LocalTorsion
        Recipes per residue, in the order the rows appeared. A residue with no
        rows is simply absent from the result, so there is no "this residue has
        no torsions" placeholder for anything else to check for.

    Raises
    ------
    maws.errors.ConfigurationError
        If a row names a residue that has no atom count.
    """
    grouped: dict[str, list[TorsionTemplate]] = {}
    for residue, start, bond, end in rotations:
        if residue not in n_atoms_of:
            raise ConfigurationError(
                f"the rotation table names residue {residue!r}, which does "
                f"not appear in the residue name and atom count tables"
            )
        n_atoms = n_atoms_of[residue]
        pivot = _resolve(int(start), n_atoms)
        axis = _resolve(int(bond), n_atoms)
        template: TorsionTemplate
        if end is None:
            template = BackboneTorsion(pivot=pivot, bond=axis)
        else:
            template = LocalTorsion(
                pivot=pivot, bond=axis, stop=_resolve(int(end), n_atoms)
            )
        grouped.setdefault(residue, []).append(template)
    return {name: tuple(items) for name, items in grouped.items()}


def _parse_aliases(aliases: Sequence[Sequence[str]]) -> dict[str, AliasSet]:
    """Turn ``[token, alone, start, middle, end]`` rows into an alias mapping.

    Parameters
    ----------
    aliases : sequence of sequence of str
        One row per written token.

    Returns
    -------
    dict of str to AliasSet
        One entry per row.

    Raises
    ------
    maws.errors.ConfigurationError
        If a row does not have exactly five entries. A short row means a typo
        in the table, which would otherwise show up as a structure built with
        the wrong residue at one end.
    """
    resolved: dict[str, AliasSet] = {}
    for row in aliases:
        if len(row) != 5:
            raise ConfigurationError(
                f"an alias row must be [token, alone, start, middle, end], "
                f"got {list(row)!r}"
            )
        token, alone, start, middle, end = row
        resolved[token] = AliasSet(alone, start, middle, end)
    return resolved


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class NucleotideSequence:
    """An aptamer strand, written as tokens, 5' end first.

    Growing a sequence returns a new one and leaves this one alone. The search
    relies on that: it tries several ways of extending the same strand and
    needs the original still intact after each attempt.

    Parameters
    ----------
    tokens : tuple of str
        Nucleotides as written, 5' to 3'. For RNA these are ``G``, ``A``, ``U``
        and ``C``; for DNA, ``G``, ``A``, ``T`` and ``C``.

    See Also
    --------
    ResidueLibrary.alias : What a single token maps to.

    Examples
    --------
    >>> seq = NucleotideSequence.parse("G A U")
    >>> str(seq)
    'G A U'
    >>> str(seq.appended("C"))
    'G A U C'
    >>> str(seq)
    'G A U'
    """

    tokens: tuple[str, ...]

    @classmethod
    def parse(cls, text: str) -> NucleotideSequence:
        """Read a sequence written as whitespace-separated tokens.

        Parameters
        ----------
        text : str
            Tokens separated by any run of whitespace. An empty or blank string
            gives an empty sequence, which is what a search starts from.

        Returns
        -------
        NucleotideSequence
            The parsed sequence.

        Examples
        --------
        >>> NucleotideSequence.parse("G  A").tokens
        ('G', 'A')
        >>> NucleotideSequence.parse("").tokens
        ()
        """
        return cls(tuple(text.split()))

    def appended(self, token: str) -> NucleotideSequence:
        """Return this sequence with one more nucleotide at the 3' end.

        Parameters
        ----------
        token : str
            The nucleotide to add.

        Returns
        -------
        NucleotideSequence
            A new sequence. This one does not change.
        """
        return NucleotideSequence(self.tokens + (token,))

    def prepended(self, token: str) -> NucleotideSequence:
        """Return this sequence with one more nucleotide at the 5' end.

        Parameters
        ----------
        token : str
            The nucleotide to add.

        Returns
        -------
        NucleotideSequence
            A new sequence. This one does not change.
        """
        return NucleotideSequence((token,) + self.tokens)

    def grown(self, token: str, direction: Direction) -> NucleotideSequence:
        """Return this sequence with one more nucleotide at the chosen end.

        Parameters
        ----------
        token : str
            The nucleotide to add.
        direction : {"3prime", "5prime"}
            Which end of the strand to add it to.

        Returns
        -------
        NucleotideSequence
            A new sequence. This one does not change.

        Examples
        --------
        >>> str(NucleotideSequence.parse("G").grown("A", "5prime"))
        'A G'
        """
        if direction == "3prime":
            return self.appended(token)
        return self.prepended(token)

    def canonical(self, library: ResidueLibrary) -> tuple[str, ...]:
        """Return the modelling-program name for each nucleotide in the strand.

        Each token has four possible names and the right one depends only on
        where the token sits: alone, at the 5' end, in the middle, or at the 3'
        end.

        Parameters
        ----------
        library : ResidueLibrary
            Library supplying the four names for each token.

        Returns
        -------
        tuple of str
            One name per token, in the same order.

        Raises
        ------
        maws.errors.ConfigurationError
            If the library does not know one of the tokens.

        Examples
        --------
        >>> from maws.libraries import rna
        >>> NucleotideSequence.parse("G").canonical(rna())
        ('GN',)
        >>> NucleotideSequence.parse("G A U").canonical(rna())
        ('G5', 'A', 'U3')
        """
        if not self.tokens:
            return ()
        if len(self.tokens) == 1:
            return (library.alias(self.tokens[0]).alone,)
        return (
            library.alias(self.tokens[0]).start,
            *(library.alias(token).middle for token in self.tokens[1:-1]),
            library.alias(self.tokens[-1]).end,
        )

    def __len__(self) -> int:
        return len(self.tokens)

    def __bool__(self) -> bool:
        return bool(self.tokens)

    def __iter__(self) -> Iterator[str]:
        return iter(self.tokens)

    def __str__(self) -> str:
        return " ".join(self.tokens)

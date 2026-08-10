"""
maws.pose
=========

Atom coordinates, and the named windows used to point at parts of them.

A :class:`Pose` is one set of atom positions for a whole design: an ``(N, 3)``
array of numbers in ångström, covering every chain at once. Moving atoms about
never changes a pose. Each method returns a *new* pose and leaves the one it
was called on exactly as it was, in the same way that ``"abc".upper()`` gives
back a new string.

That matters for the way MAWS works. It tries thousands of shapes for the same
strand, and after each attempt it needs the starting shape back. With a value
that cannot change, "get the starting shape back" is just using the variable
again — there is nothing to copy and nothing to restore.

A :class:`ChainView` names one chain within a pose, and a :class:`ResidueView`
names one residue within a chain. Neither holds coordinates. They exist so that
code can say "the first residue of the aptamer" instead of "atoms 0 through 32",
and so that turning a bond can be written without working out any array offsets
by hand.

The standard loop
-----------------
Almost everything MAWS does is the same four steps, repeated: take a candidate
position for the strand, give the strand a random shape, score it, and write
the score down. Written out, that is::

    energies = []
    for placement in islice(sampler, n_samples):
        pose = base.place(aptamer, placement)
        pose = pose.rotate_all(torsions, angles.sample())
        energies.append(model.evaluate(pose))

``base`` never changes, so every pass starts from the same place.

Units
-----
Positions are plain numbers in ångström everywhere in this module. The
simulation package MAWS scores with attaches units to its own arrays;
:meth:`Pose.to_openmm` and :meth:`Pose.from_openmm` are the only two places
that convert, so nothing else has to think about it.

Examples
--------
>>> import numpy as np
>>> from maws.values import AtomRange, Torsion
>>> pose = Pose(np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0]]))
>>> turned = pose.rotate(Torsion(0, 1, AtomRange(2, 3)), np.pi)
>>> np.round(turned.xyz[2], 6)
array([ 1., -1.,  0.])
>>> np.round(pose.xyz[2], 6)
array([1., 1., 0.])
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from maws.errors import ConfigurationError
from maws.values import (
    AtomRange,
    Direction,
    NucleotideSequence,
    ResidueLibrary,
    Torsion,
)

if TYPE_CHECKING:  # pragma: no cover - needed by type checkers only
    from openmm.unit import Quantity

    from maws.sampling import Placement

__all__ = ["ChainView", "Pose", "ResidueView", "rodrigues"]


def rodrigues(axis: np.ndarray, angle: float) -> np.ndarray:
    """Return the matrix that rotates by `angle` radians about `axis`.

    Parameters
    ----------
    axis : numpy.ndarray
        Shape ``(3,)``. The line to turn about. Its length is ignored, so it
        can be the difference between two atom positions as it comes.
    angle : float
        How far to turn, in radians.

    Returns
    -------
    numpy.ndarray
        Shape ``(3, 3)``. Positions are stored one atom per row, so apply it as
        ``xyz @ matrix.T``.

    Raises
    ------
    maws.errors.ConfigurationError
        If `axis` has zero length. That happens when the two atoms defining it
        sit at the same position, which leaves no line to turn about.

    Notes
    -----
    A positive angle turns counter-clockwise when you look down the axis from
    its far end towards its start — the right-hand rule.

    Examples
    --------
    Turning the x direction a quarter turn about z gives the y direction:

    >>> matrix = rodrigues(np.array([0.0, 0.0, 1.0]), np.pi / 2)
    >>> np.round(np.array([1.0, 0.0, 0.0]) @ matrix.T, 6)
    array([0., 1., 0.])
    """
    norm = float(np.linalg.norm(axis))
    if norm == 0.0:
        raise ConfigurationError(
            "cannot rotate about a zero-length axis: the two atoms defining "
            "it are at the same position"
        )
    x, y, z = np.asarray(axis, dtype=np.float64) / norm
    cross = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + np.sin(angle) * cross + (1.0 - np.cos(angle)) * (cross @ cross)


class Pose:
    """One set of atom positions for a whole design.

    Parameters
    ----------
    xyz : array_like
        Shape ``(N, 3)`` positions in ångström, one row per atom. Copied on the
        way in, so the array passed here can be reused freely afterwards.
    system : object, optional
        The :class:`~maws.topology.BuiltSystem` these positions belong to.
        Carried so that a pose can be written out or handed on without a second
        argument; it plays no part in any calculation here.

    Attributes
    ----------
    xyz : numpy.ndarray
        The positions. Read-only: attempting to write to it raises, which keeps
        the promise that a pose never changes.
    n_atoms : int
        How many atoms there are.

    Raises
    ------
    maws.errors.ConfigurationError
        If `xyz` is not shaped ``(N, 3)``.

    See Also
    --------
    ChainView : Names one chain's atoms within a pose.
    maws.energy.EnergyModel : Turns a pose into a number.

    Examples
    --------
    >>> import numpy as np
    >>> pose = Pose(np.zeros((4, 3)))
    >>> pose
    <Pose 4 atoms>
    >>> pose.xyz.flags.writeable
    False
    """

    __slots__ = ("_system", "_xyz")

    def __init__(self, xyz: Any, system: object | None = None) -> None:
        array = np.ascontiguousarray(xyz, dtype=np.float64)
        if array.ndim != 2 or array.shape[1] != 3:
            raise ConfigurationError(
                f"positions must be shaped (N, 3), one row per atom; "
                f"got shape {array.shape}"
            )
        array.flags.writeable = False
        self._xyz = array
        self._system = system

    # -- reading ----------------------------------------------------------

    @property
    def xyz(self) -> np.ndarray:
        """numpy.ndarray : Read-only ``(N, 3)`` positions in ångström."""
        return self._xyz

    @property
    def n_atoms(self) -> int:
        """int : How many atoms this pose holds positions for."""
        return int(self._xyz.shape[0])

    @property
    def system(self) -> object | None:
        """object or None : The built system these positions belong to."""
        return self._system

    def __len__(self) -> int:
        return self.n_atoms

    def __repr__(self) -> str:
        return f"<Pose {self.n_atoms} atoms>"

    def atoms(self, span: AtomRange | ChainView) -> np.ndarray:
        """Return the positions of one run of atoms.

        Parameters
        ----------
        span : maws.values.AtomRange or ChainView
            Which atoms to read. Passing a chain reads that whole chain.

        Returns
        -------
        numpy.ndarray
            Shape ``(len(span), 3)``, read-only. This looks into the pose's own
            array rather than copying it.

        Examples
        --------
        >>> import numpy as np
        >>> from maws.values import AtomRange
        >>> Pose(np.arange(12, dtype=float).reshape(4, 3)).atoms(AtomRange(1, 3))
        array([[3., 4., 5.],
               [6., 7., 8.]])
        """
        return self._xyz[_span_of(span).as_slice()]

    def centroid(self, span: AtomRange | ChainView | None = None) -> np.ndarray:
        """Return the average position of a run of atoms.

        Every atom counts the same, regardless of what element it is.

        Parameters
        ----------
        span : maws.values.AtomRange or ChainView, optional
            Which atoms to average. Defaults to all of them.

        Returns
        -------
        numpy.ndarray
            Shape ``(3,)``, in ångström.
        """
        chosen = self._xyz if span is None else self.atoms(span)
        return chosen.mean(axis=0)

    # -- moving atoms -----------------------------------------------------

    def _replace(self, xyz: np.ndarray) -> Pose:
        """Return a pose with these positions and the same owning system.

        Parameters
        ----------
        xyz : numpy.ndarray
            New ``(N, 3)`` positions.

        Returns
        -------
        Pose
            The new pose.
        """
        return Pose(xyz, self._system)

    def rotate(self, torsion: Torsion, angle: float) -> Pose:
        """Turn one bond, swinging the atoms attached to one side of it.

        Parameters
        ----------
        torsion : maws.values.Torsion
            The bond to turn, and which atoms swing with it.
        angle : float
            How far to turn, in radians.

        Returns
        -------
        Pose
            A new pose. This one does not change.

        See Also
        --------
        rotate_all : Turn several bonds in one call.

        Examples
        --------
        Three atoms in an L shape. Turning the first bond a quarter turn swings
        the third atom out of the plane:

        >>> import numpy as np
        >>> from maws.values import AtomRange, Torsion
        >>> pose = Pose(np.array([[0.0, 0, 0], [0, 0, 1], [1, 0, 1]]))
        >>> spun = pose.rotate(Torsion(0, 1, AtomRange(2, 3)), np.pi / 2)
        >>> np.round(spun.xyz[2], 6)
        array([0., 1., 1.])
        """
        origin = self._xyz[torsion.pivot]
        axis = self._xyz[torsion.bond] - origin
        moving = torsion.moving.as_slice()
        xyz = self._xyz.copy()
        xyz[moving] = (xyz[moving] - origin) @ rodrigues(axis, angle).T + origin
        return self._replace(xyz)

    def rotate_all(self, torsions: Sequence[Torsion], angles: Iterable[float]) -> Pose:
        """Turn several bonds in order, each by its own angle.

        Parameters
        ----------
        torsions : sequence of maws.values.Torsion
            The bonds to turn.
        angles : iterable of float
            How far to turn each one, in radians. Must supply exactly as many
            angles as there are bonds.

        Returns
        -------
        Pose
            A new pose with every turn applied. This one does not change.

        Raises
        ------
        maws.errors.ConfigurationError
            If the number of angles does not match the number of bonds.

        Notes
        -----
        The turns are applied one after another, so a bond that moves part of
        the strand also carries any bonds further along it. That is the same
        thing that happens in the molecule, so order matters and this method
        keeps it.

        Examples
        --------
        >>> import numpy as np
        >>> from maws.values import AtomRange, Torsion
        >>> pose = Pose(np.array([[0.0, 0, 0], [0, 0, 1], [1, 0, 1]]))
        >>> pose.rotate_all([Torsion(0, 1, AtomRange(2, 3))], [np.pi])
        <Pose 3 atoms>
        """
        angles = list(angles)
        if len(angles) != len(torsions):
            raise ConfigurationError(
                f"got {len(angles)} angles for {len(torsions)} bonds; "
                f"there must be exactly one angle per bond"
            )
        pose = self
        for torsion, angle in zip(torsions, angles, strict=True):
            pose = pose.rotate(torsion, angle)
        return pose

    def translate(self, span: AtomRange | ChainView, shift: np.ndarray) -> Pose:
        """Slide a run of atoms, without turning them.

        Parameters
        ----------
        span : maws.values.AtomRange or ChainView
            Which atoms to move. Passing a chain moves the whole chain.
        shift : numpy.ndarray
            Shape ``(3,)``. How far to move, in ångström, along each axis.

        Returns
        -------
        Pose
            A new pose. This one does not change.
        """
        xyz = self._xyz.copy()
        xyz[_span_of(span).as_slice()] += np.asarray(shift, dtype=np.float64)
        return self._replace(xyz)

    def rotate_about(
        self,
        span: AtomRange | ChainView,
        axis: np.ndarray,
        angle: float,
        origin: np.ndarray,
    ) -> Pose:
        """Turn a run of atoms about any line, not necessarily a bond.

        Parameters
        ----------
        span : maws.values.AtomRange or ChainView
            Which atoms to turn.
        axis : numpy.ndarray
            Shape ``(3,)``. The direction of the line to turn about; its length
            is ignored.
        angle : float
            How far to turn, in radians.
        origin : numpy.ndarray
            Shape ``(3,)``. A point the line passes through, in ångström.

        Returns
        -------
        Pose
            A new pose. This one does not change.

        See Also
        --------
        place : Puts a whole chain somewhere, using this.
        """
        moving = _span_of(span).as_slice()
        xyz = self._xyz.copy()
        xyz[moving] = (xyz[moving] - origin) @ rodrigues(axis, angle).T + origin
        return self._replace(xyz)

    def place(self, chain: AtomRange | ChainView, placement: Placement) -> Pose:
        """Move a whole chain to a new position and orientation.

        This is the first step of the standard loop: a sampler proposes
        somewhere for the strand to sit, and this puts it there.

        Parameters
        ----------
        chain : ChainView or maws.values.AtomRange
            The chain to move. Every other chain stays exactly where it is.
        placement : maws.sampling.Placement
            Where to move it and how to turn it, as produced by a sampler.

        Returns
        -------
        Pose
            A new pose. This one does not change.

        Notes
        -----
        The chain is first slid by ``placement.position``, then turned by
        ``placement.angle`` about ``placement.axis``. The turn is centred on
        the chain's first atom, so the slide decides roughly where the chain
        sits and the turn decides which way it faces.

        See Also
        --------
        maws.sampling.SurfaceSampler : Produces placements near a target.
        """
        span = _span_of(chain)
        moved = self.translate(span, placement.position)
        return moved.rotate_about(
            span,
            placement.axis,
            placement.angle,
            origin=moved.xyz[span.start],
        )

    def jittered(self, offsets: np.ndarray) -> Pose:
        """Nudge every atom by its own small displacement.

        Parameters
        ----------
        offsets : numpy.ndarray
            Shape ``(N, 3)``. How far to move each atom, in ångström.

        Returns
        -------
        Pose
            A new pose. This one does not change.

        Raises
        ------
        maws.errors.ConfigurationError
            If `offsets` does not have one row per atom.

        See Also
        --------
        maws.relax.perturb_and_minimize : Uses this to shake a structure loose
            from a poor local arrangement.
        """
        offsets = np.asarray(offsets, dtype=np.float64)
        if offsets.shape != self._xyz.shape:
            raise ConfigurationError(
                f"offsets must be shaped {self._xyz.shape}, one row per atom; "
                f"got {offsets.shape}"
            )
        return self._replace(self._xyz + offsets)

    def with_span(self, span: AtomRange | ChainView, xyz: np.ndarray) -> Pose:
        """Return a pose with one run of atoms given entirely new positions.

        Parameters
        ----------
        span : maws.values.AtomRange or ChainView
            Which atoms to overwrite.
        xyz : numpy.ndarray
            Shape ``(len(span), 3)`` replacement positions in ångström.

        Returns
        -------
        Pose
            A new pose. This one does not change.

        Raises
        ------
        maws.errors.ConfigurationError
            If `xyz` does not have one row per atom in `span`.
        """
        resolved = _span_of(span)
        xyz = np.asarray(xyz, dtype=np.float64)
        if xyz.shape != (len(resolved), 3):
            raise ConfigurationError(
                f"replacing {resolved} needs {len(resolved)} rows of "
                f"positions, got shape {xyz.shape}"
            )
        updated = self._xyz.copy()
        updated[resolved.as_slice()] = xyz
        return self._replace(updated)

    # -- talking to the simulation package --------------------------------

    def to_openmm(self) -> Quantity:
        """Return these positions with ångström units attached.

        Returns
        -------
        openmm.unit.Quantity
            The same numbers, tagged as ångström, which is the form OpenMM
            expects when it is given positions to score.

        See Also
        --------
        from_openmm : The other direction.
        """
        from openmm import unit

        return self._xyz * unit.angstrom

    @classmethod
    def from_openmm(cls, positions: Quantity, system: object | None = None) -> Pose:
        """Build a pose from OpenMM positions, converting to ångström.

        Parameters
        ----------
        positions : openmm.unit.Quantity
            Positions with length units attached, as OpenMM returns them.
        system : object, optional
            The built system these positions belong to.

        Returns
        -------
        Pose
            The same positions as plain numbers in ångström.
        """
        from openmm import unit

        return cls(np.asarray(positions.value_in_unit(unit.angstrom)), system)


def _span_of(target: AtomRange | ChainView) -> AtomRange:
    """Return the atom range a chain or range refers to.

    Parameters
    ----------
    target : maws.values.AtomRange or ChainView
        Either an atom range already, or a chain that covers one.

    Returns
    -------
    maws.values.AtomRange
        The range itself.
    """
    return target if isinstance(target, AtomRange) else target.span


@dataclass(frozen=True, slots=True)
class ChainView:
    """A named window onto one chain of a design.

    A design holds every atom of every chain in one array. A chain view records
    which stretch of that array belongs to one chain, what residues make it up,
    and where each residue starts. It holds no coordinates and refers to
    nothing that contains it, so one can be built by hand in a test from a
    sequence and a span.

    Parameters
    ----------
    role : str
        What this chain is called, e.g. ``"aptamer"`` or ``"ligand"``.
    span : maws.values.AtomRange
        Which atoms of the design belong to this chain.
    sequence : maws.values.NucleotideSequence
        The chain's nucleotides as written.
    library : maws.values.ResidueLibrary
        The residue descriptions the sequence is resolved against.
    canonical : tuple of str
        The modelling-program name of each residue, worked out from `sequence`.
    residue_offsets : tuple of int
        Where each residue's first atom sits, counted from this chain's first
        atom.

    See Also
    --------
    build : The usual way to make one.
    ResidueView : A window onto a single residue of a chain.

    Examples
    --------
    >>> from maws.libraries import rna
    >>> from maws.values import AtomRange, NucleotideSequence
    >>> chain = ChainView.build(
    ...     "aptamer", AtomRange(0, 66), NucleotideSequence.parse("G A"), rna()
    ... )
    >>> chain.n_residues
    2
    >>> chain.residue(1).span
    AtomRange(start=32, stop=66)
    """

    role: str
    span: AtomRange
    sequence: NucleotideSequence
    library: ResidueLibrary
    canonical: tuple[str, ...]
    residue_offsets: tuple[int, ...]

    @classmethod
    def build(
        cls,
        role: str,
        span: AtomRange,
        sequence: NucleotideSequence,
        library: ResidueLibrary,
    ) -> ChainView:
        """Work out a chain's residue names and offsets, and return the view.

        Parameters
        ----------
        role : str
            What this chain is called.
        span : maws.values.AtomRange
            Which atoms of the design belong to it.
        sequence : maws.values.NucleotideSequence
            The chain's nucleotides as written.
        library : maws.values.ResidueLibrary
            The residue descriptions to resolve the sequence against.

        Returns
        -------
        ChainView
            A view with residue names and offsets filled in.
        """
        canonical = sequence.canonical(library)
        offsets: list[int] = []
        tally = 0
        for name in canonical:
            offsets.append(tally)
            tally += library[name].n_atoms
        return cls(
            role=role,
            span=span,
            sequence=sequence,
            library=library,
            canonical=canonical,
            residue_offsets=tuple(offsets),
        )

    @property
    def n_residues(self) -> int:
        """int : How many residues this chain has."""
        return len(self.canonical)

    def __len__(self) -> int:
        return len(self.span)

    def __repr__(self) -> str:
        return (
            f"<ChainView {self.role!r} {len(self.span)} atoms, "
            f"{self.n_residues} residues: {self.sequence or '(empty)'}>"
        )

    def residue(self, index: int) -> ResidueView:
        """Return a window onto one residue of this chain.

        Parameters
        ----------
        index : int
            Which residue, counting from the 5' end. Negative values count back
            from the 3' end, as with a Python list, so ``-1`` is the last one.

        Returns
        -------
        ResidueView
            A window onto that residue.

        Raises
        ------
        maws.errors.ConfigurationError
            If there is no residue at that index.
        """
        resolved = index + self.n_residues if index < 0 else index
        if not 0 <= resolved < self.n_residues:
            raise ConfigurationError(
                f"chain {self.role!r} has {self.n_residues} residues, so "
                f"index {index} does not name one"
            )
        return ResidueView(chain=self, index=resolved)

    def residues(self) -> tuple[ResidueView, ...]:
        """Return a window onto every residue, 5' end first.

        Returns
        -------
        tuple of ResidueView
            One window per residue.
        """
        return tuple(self.residue(i) for i in range(self.n_residues))


@dataclass(frozen=True, slots=True)
class ResidueView:
    """A window onto one residue of a chain.

    Its job is to answer "which atoms is this residue" and "what are its
    turnable bonds, as positions in the whole design" without the caller doing
    any index arithmetic.

    Parameters
    ----------
    chain : ChainView
        The chain this residue is part of.
    index : int
        Which residue of that chain, counting from the 5' end. Already
        resolved, so never negative.

    See Also
    --------
    ChainView.residue : How one of these is normally obtained.
    """

    chain: ChainView
    index: int

    def __repr__(self) -> str:
        return f"<ResidueView {self.name} #{self.index} of {self.chain.role!r}>"

    @property
    def name(self) -> str:
        """str : The modelling-program name of this residue, e.g. ``"G5"``."""
        return self.chain.canonical[self.index]

    @property
    def span(self) -> AtomRange:
        """maws.values.AtomRange : Which atoms of the design this residue is."""
        start = self.chain.span.start + self.chain.residue_offsets[self.index]
        return AtomRange(start, start + self.chain.library[self.name].n_atoms)

    @property
    def n_torsions(self) -> int:
        """int : How many turnable bonds this residue has."""
        return self.chain.library[self.name].n_torsions

    def torsion(self, index: int, direction: Direction = "3prime") -> Torsion:
        """Return one of this residue's turnable bonds.

        The residue descriptions are written with atom numbers counted from the
        start of a residue. A pose numbers atoms from the start of the whole
        design. This method is where those two are reconciled, and it is the
        only place in MAWS that does so.

        Parameters
        ----------
        index : int
            Which bond, numbered as in the residue's description.
        direction : {"3prime", "5prime"}, default="3prime"
            Which side of the bond swings. ``"3prime"`` swings the part of the
            strand towards its 3' end; ``"5prime"`` swings the part towards its
            5' end. Either produces the same change in shape — they differ only
            in which part of the strand stays where it is, which matters when
            the rest of the strand is already positioned against a target.

        Returns
        -------
        maws.values.Torsion
            The bond, numbered from the start of the design.

        Raises
        ------
        maws.errors.ConfigurationError
            If there is no bond at that index, or `direction` is not one of the
            two allowed values.

        Examples
        --------
        >>> from maws.libraries import rna
        >>> from maws.values import AtomRange, NucleotideSequence
        >>> chain = ChainView.build(
        ...     "aptamer", AtomRange(0, 33), NucleotideSequence.parse("G"), rna()
        ... )
        >>> chain.residue(0).torsion(0)
        Torsion(pivot=0, bond=1, moving=AtomRange(start=1, stop=33))
        >>> chain.residue(0).torsion(0, "5prime")
        Torsion(pivot=1, bond=0, moving=AtomRange(start=0, stop=1))
        """
        template = self.chain.library[self.name]
        if not 0 <= index < template.n_torsions:
            raise ConfigurationError(
                f"residue {self.name!r} has {template.n_torsions} turnable "
                f"bonds, so index {index} does not name one"
            )
        if direction not in ("3prime", "5prime"):
            raise ConfigurationError(
                f"direction must be '3prime' or '5prime', got {direction!r}"
            )
        recipe = template.torsions[index]
        offset = self.chain.span.start + self.chain.residue_offsets[self.index]
        if direction == "3prime":
            return recipe.placed(offset, self.chain.span)
        return recipe.reversed(offset, self.chain.span)

    def torsions(
        self, direction: Direction = "3prime", *, limit: int | None = None
    ) -> tuple[Torsion, ...]:
        """Return several of this residue's turnable bonds at once.

        Handy for the standard loop, where the same set of bonds is turned to
        new random angles on every pass.

        Parameters
        ----------
        direction : {"3prime", "5prime"}, default="3prime"
            Which side of each bond swings. See :meth:`torsion`.
        limit : int, optional
            Return only the first this many bonds. Defaults to all of them.

        Returns
        -------
        tuple of maws.values.Torsion
            The bonds, numbered from the start of the design.

        Examples
        --------
        >>> from maws.libraries import rna
        >>> from maws.values import AtomRange, NucleotideSequence
        >>> chain = ChainView.build(
        ...     "aptamer", AtomRange(0, 33), NucleotideSequence.parse("G"), rna()
        ... )
        >>> len(chain.residue(0).torsions())
        4
        >>> len(chain.residue(0).torsions(limit=3))
        3
        """
        count = self.n_torsions if limit is None else min(limit, self.n_torsions)
        return tuple(self.torsion(i, direction) for i in range(count))

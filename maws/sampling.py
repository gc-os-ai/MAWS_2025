"""
maws.sampling
=============

Proposing places for the aptamer to sit, and shapes for it to take.

MAWS does not calculate where a strand should go. It guesses, many thousands of
times, and keeps track of how the guesses score. This module produces the
guesses.

There are two kinds. A :class:`Placement` says where the whole strand should sit
next to the target and which way it should face. A set of angles from
:class:`TorsionAngles` says how far to turn each of the strand's turnable bonds,
which decides its shape.

Samplers are iterable, in the same way a data loader is::

    for placement in islice(sampler, 5000):
        ...

so ``next(sampler)``, ``zip``, and a plain ``for`` loop all work without any
extra code.

Rejecting bad guesses
---------------------
A position drawn at random near the target is often *inside* it, which is
physically impossible. :class:`Excluder` throws those away. It uses the
solvent-accessible surface: imagine rolling a small ball over the target's
atoms, and keep only the positions where the ball's centre could sit without
overlapping any atom. The ball stands for a water molecule, which is why its
default radius is 1.4 ångström.

Repeatability
-------------
Every sampler takes a `rng` argument, a NumPy random generator. Pass one
created with a fixed seed and the whole run can be repeated exactly. Leave it
out and a fresh one is made, so runs differ.

Examples
--------
>>> import numpy as np
>>> angles = TorsionAngles(4, rng=np.random.default_rng(0))
>>> angles.sample().shape
(4,)
"""

from __future__ import annotations

import warnings
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, Protocol, runtime_checkable

import numpy as np
from scipy.spatial import KDTree

from maws.errors import ConfigurationError, SamplingError
from maws.geometry import centre_of_mass
from maws.values import AtomRange

if TYPE_CHECKING:  # pragma: no cover - needed by type checkers only
    from maws.topology import BuiltSystem

__all__ = [
    "Envelope",
    "Excluder",
    "FixedSampler",
    "Placement",
    "Sampler",
    "Sphere",
    "SurfaceFollowingSampler",
    "SurfaceSampler",
    "TorsionAngles",
    "make_sampler",
]

VDW_RADII: dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
    "Si": 2.10,
    "B": 1.92,
    "Se": 1.90,
    "Na": 2.27,
    "Mg": 1.73,
    "K": 2.75,
    "Ca": 2.31,
    "Fe": 2.00,
    "Zn": 1.39,
    "Cu": 1.40,
    "Mn": 2.00,
    "Ni": 1.63,
}
"""How much room each element takes up, in ångström, by chemical symbol.

These are van der Waals radii: roughly, how close another atom can get before
the two start pushing each other away. The values are Bondi's (1964), the same
set used by the common molecular viewers, so a surface drawn here matches one
drawn in those. They are written out rather than read from a chemistry package,
which keeps the surface a purely geometric idea and saves a dependency for
twenty-odd numbers.
"""

DEFAULT_VDW = 1.70
"""Radius in ångström assumed for an element not listed in :data:`VDW_RADII`.

The value for carbon, the most common element in these molecules.
"""


@dataclass(frozen=True)
class Placement:
    """Somewhere for the aptamer to sit, and which way it should face.

    Parameters
    ----------
    position : numpy.ndarray
        Shape ``(3,)``. How far to move the strand along each axis, in
        ångström.
    axis : numpy.ndarray
        Shape ``(3,)``, of length one. The line to turn the strand about after
        moving it.
    angle : float
        How far to turn it about that line, in radians.

    See Also
    --------
    maws.pose.Pose.place : Applies one of these.
    """

    position: np.ndarray
    axis: np.ndarray
    angle: float


@runtime_checkable
class Envelope(Protocol):
    """Anything that can propose a placement inside some region of space.

    Only used as the first, cheap stage of :class:`SurfaceSampler`, which then
    throws away proposals that land inside the target.
    """

    def sample(self) -> Placement:
        """Return one proposed placement.

        Returns
        -------
        Placement
            A position and orientation drawn from this region.
        """
        ...


@runtime_checkable
class Sampler(Protocol):
    """Anything that produces an unending stream of placements.

    Implementing :meth:`sample` and :meth:`accepts` is enough;
    :class:`SurfaceSampler` and the others also support ``next()`` and ``for``
    loops.
    """

    def sample(self) -> Placement:
        """Return one placement.

        Returns
        -------
        Placement
            The next placement in the stream.
        """
        ...

    def accepts(self, points: np.ndarray) -> bool:
        """Say whether a molecule sitting at these positions is acceptable.

        A placement says where the middle of the strand goes. The strand is
        not a point, so a middle that clears the target is not the same as a
        strand that clears it. This is asked once the placement has been
        applied and the atoms are known.

        Parameters
        ----------
        points : numpy.ndarray
            Shape ``(N, 3)``. Where the strand's atoms have ended up, in
            ångström.

        Returns
        -------
        bool
            True when the shape is acceptable and should be scored.
        """
        ...

    def __next__(self) -> Placement: ...

    def __iter__(self) -> Iterator[Placement]: ...


class _IterableSampler:
    """Gives a class with a ``sample`` method the behaviour of an iterator.

    Not part of the public interface. Mixing this in is what makes
    ``next(sampler)``, ``for placement in sampler`` and
    ``islice(sampler, 5000)`` work.
    """

    def sample(self) -> Placement:  # pragma: no cover - overridden everywhere
        """Return one placement.

        Returns
        -------
        Placement
            The next placement.
        """
        raise NotImplementedError

    def accepts(self, points: np.ndarray) -> bool:
        """Accept any arrangement of the atoms.

        The default for samplers that know nothing about a target, such as
        :class:`FixedSampler` and a bare :class:`Sphere`. The ones built around
        a target override this.

        Parameters
        ----------
        points : numpy.ndarray
            Shape ``(N, 3)``, in ångström. Ignored.

        Returns
        -------
        bool
            Always True.
        """
        return True

    def __next__(self) -> Placement:
        return self.sample()

    def __iter__(self) -> Iterator[Placement]:
        return self


# ---------------------------------------------------------------------------
# Angles
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TorsionAngles:
    """A source of random angles, one per turnable bond.

    Parameters
    ----------
    n : int
        How many angles to produce at a time. Match this to the number of bonds
        being turned.
    rng : numpy.random.Generator, optional
        Where the randomness comes from. Pass one built with a fixed seed to
        get the same angles every run.

    See Also
    --------
    maws.pose.Pose.rotate_all : Applies a set of angles to a set of bonds.

    Examples
    --------
    >>> import numpy as np
    >>> angles = TorsionAngles(4, rng=np.random.default_rng(0))
    >>> drawn = angles.sample()
    >>> drawn.shape, bool(((0 <= drawn) & (drawn < 2 * np.pi)).all())
    ((4,), True)

    A residue near the end of a strand can have fewer bonds worth turning than
    the setting asks for, so the count can be overridden per draw:

    >>> angles.sample(2).shape
    (2,)
    """

    n: int
    rng: np.random.Generator = field(default_factory=np.random.default_rng)

    def __post_init__(self) -> None:
        if self.n < 0:
            raise ConfigurationError(f"n must not be negative, got {self.n}")

    def sample(self, n: int | None = None) -> np.ndarray:
        """Return angles drawn evenly from a full turn.

        Parameters
        ----------
        n : int, optional
            How many angles to draw. Defaults to :attr:`n`. Pass the number of
            bonds actually being turned when that is smaller, which happens
            when some of a residue's bonds turn out to move nothing and are
            left out.

        Returns
        -------
        numpy.ndarray
            Shape ``(n,)``, each between 0 and 2π radians.

        Raises
        ------
        maws.errors.ConfigurationError
            If `n` is negative.
        """
        count = self.n if n is None else n
        if count < 0:
            raise ConfigurationError(f"n must not be negative, got {count}")
        return self.rng.uniform(0.0, 2.0 * np.pi, count)


# ---------------------------------------------------------------------------
# Regions of space
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Sphere:
    """A ball of space to draw positions from, evenly by volume.

    Parameters
    ----------
    radius : float
        How far the ball extends from its centre, in ångström.
    centre : numpy.ndarray
        Shape ``(3,)``. Where the ball is centred, in ångström.
    rng : numpy.random.Generator, optional
        Where the randomness comes from.

    Notes
    -----
    "Evenly by volume" means every cubic ångström inside the ball is equally
    likely, which takes two pieces of care. The distance from the centre is
    drawn as ``radius * u ** (1/3)`` rather than as a flat value, because a
    thin shell far out holds more volume than one near the middle. The
    direction is drawn by taking the height evenly between the poles rather
    than the polar angle, because otherwise proposals bunch up at the poles.

    Examples
    --------
    >>> import numpy as np
    >>> ball = Sphere(10.0, np.zeros(3), rng=np.random.default_rng(0))
    >>> bool(np.linalg.norm(ball.sample().position) <= 10.0)
    True
    """

    radius: float
    centre: np.ndarray
    rng: np.random.Generator = field(default_factory=np.random.default_rng)

    def sample(self) -> Placement:
        """Return one placement drawn evenly from inside the ball.

        Returns
        -------
        Placement
            A position inside the ball, with a direction and an angle drawn
            independently of it.
        """
        distance = self.radius * self.rng.uniform(0.0, 1.0) ** (1.0 / 3.0)
        return Placement(
            position=self.centre + distance * _random_direction(self.rng),
            axis=_random_direction(self.rng),
            angle=float(self.rng.uniform(0.0, 2.0 * np.pi)),
        )


def _random_direction(rng: np.random.Generator) -> np.ndarray:
    """Return a direction pointing evenly anywhere, of length one.

    Parameters
    ----------
    rng : numpy.random.Generator
        Where the randomness comes from.

    Returns
    -------
    numpy.ndarray
        Shape ``(3,)``, of length one.

    Notes
    -----
    Drawn by taking three numbers from a bell curve and scaling the result to
    length one. A bell curve in three dimensions is the same in every
    direction, so this covers directions evenly. Drawing three numbers from a
    flat range instead would favour the corners of a cube, because a cube
    reaches further along its diagonals than through its faces.
    """
    vector = rng.standard_normal(3)
    return vector / np.linalg.norm(vector)


class Excluder:
    """Decides whether a proposed position would sit inside the target.

    Works from the solvent-accessible surface: a position is acceptable when a
    ball of radius `probe` centred there would not overlap any of the target's
    atoms. Each atom is given the size its element normally takes up, so the
    test is ``distance to atom > that atom's radius + probe`` for every atom.

    Parameters
    ----------
    positions : numpy.ndarray
        Shape ``(N, 3)``. Where the target's atoms are, in ångström.
    elements : sequence of str
        Chemical symbol per atom, e.g. ``"C"``. Decides how much room each atom
        takes up. Symbols not listed in :data:`VDW_RADII` fall back to
        :data:`DEFAULT_VDW` and warn once.
    probe : float, default=1.4
        Radius in ångström of the ball rolled over the target. 1.4 Å is the
        size of a water molecule, the usual choice. A larger value keeps
        proposals further out; a smaller one lets them into tighter pockets.

    See Also
    --------
    SurfaceSampler : Combines this with a region to draw from.

    Examples
    --------
    One carbon atom at the origin. Carbon's radius is 1.70 Å, so with the
    default probe nothing within 3.1 Å of it is acceptable:

    >>> import numpy as np
    >>> excluder = Excluder(np.zeros((1, 3)), ["C"])
    >>> excluder.is_clear(np.array([2.0, 0.0, 0.0]))
    False
    >>> excluder.is_clear(np.array([5.0, 0.0, 0.0]))
    True
    """

    _warned: set[str] = set()

    __slots__ = ("_positions", "_reach", "_tree", "_widest")

    def __init__(
        self,
        positions: np.ndarray,
        elements: Sequence[str],
        *,
        probe: float = 1.4,
    ) -> None:
        positions = np.asarray(positions, dtype=np.float64)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ConfigurationError(
                f"positions must be shaped (N, 3), got {positions.shape}"
            )
        if len(elements) != positions.shape[0]:
            raise ConfigurationError(
                f"got {positions.shape[0]} positions but {len(elements)} "
                f"element symbols"
            )
        if probe < 0:
            raise ConfigurationError(f"probe must not be negative, got {probe}")

        radii = np.array([_radius_of(symbol) for symbol in elements], dtype=float)
        self._positions = positions
        self._reach = radii + probe
        self._widest = float(self._reach.max()) if len(self._reach) else 0.0
        self._tree = KDTree(positions)

    def __repr__(self) -> str:
        return f"<Excluder {len(self._positions)} atoms>"

    def is_clear(self, point: np.ndarray) -> bool:
        """Say whether a position is far enough from every target atom.

        Parameters
        ----------
        point : numpy.ndarray
            Shape ``(3,)``. The position to test, in ångström.

        Returns
        -------
        bool
            True when the position is acceptable.

        Notes
        -----
        Only atoms within the widest possible clearance are checked, found with
        a spatial index, so the cost does not grow with the size of the target.

        A position sitting exactly on the boundary counts as blocked. Since
        positions are drawn from a continuous range this never comes up in
        practice, but the rule has to be one or the other.
        """
        nearby = self._tree.query_ball_point(point, self._widest)
        if not nearby:
            return True
        gaps = self._positions[nearby] - np.asarray(point, dtype=np.float64)
        return bool(((gaps**2).sum(axis=1) > self._reach[nearby] ** 2).all())

    def all_clear(self, points: np.ndarray) -> bool:
        """Say whether every one of many positions clears the target.

        Parameters
        ----------
        points : numpy.ndarray
            Shape ``(N, 3)``. The positions to test, in ångström. An empty
            array is acceptable and passes.

        Returns
        -------
        bool
            True when no position is inside the target.

        See Also
        --------
        is_clear : The same test for one position.

        Notes
        -----
        This is what a whole molecule has to pass. Testing only where its
        middle went says nothing about the rest of it: a nucleotide is about
        ten ångström across, so it can have its middle in clear space and half
        its atoms buried in the target.

        Examples
        --------
        A carbon at the origin blocks anything within 3.1 Å of it, so a pair
        of positions is only clear when both of them are:

        >>> import numpy as np
        >>> excluder = Excluder(np.zeros((1, 3)), ["C"])
        >>> excluder.all_clear(np.array([[10.0, 0, 0], [20.0, 0, 0]]))
        True
        >>> excluder.all_clear(np.array([[10.0, 0, 0], [0.5, 0, 0]]))
        False
        """
        points = np.asarray(points, dtype=np.float64)
        if points.size == 0:
            return True
        for nearby, point in zip(
            self._tree.query_ball_point(points, self._widest), points, strict=True
        ):
            if not nearby:
                continue
            gaps = self._positions[nearby] - point
            if not ((gaps**2).sum(axis=1) > self._reach[nearby] ** 2).all():
                return False
        return True


def _radius_of(symbol: str) -> float:
    """Return how much room one element takes up, in ångström.

    Parameters
    ----------
    symbol : str
        Chemical symbol, e.g. ``"C"``. Case-sensitive, as chemical symbols are.

    Returns
    -------
    float
        The element's van der Waals radius, or :data:`DEFAULT_VDW` if it is not
        listed.

    Warns
    -----
    RuntimeWarning
        Once per unrecognised symbol per process, so a systematic problem is
        visible without one warning per atom.
    """
    radius = VDW_RADII.get(symbol)
    if radius is not None:
        return radius
    if symbol not in Excluder._warned:
        Excluder._warned.add(symbol)
        warnings.warn(
            f"unknown element {symbol!r}; assuming it takes up as much room "
            f"as carbon ({DEFAULT_VDW} Å)",
            RuntimeWarning,
            stacklevel=3,
        )
    return DEFAULT_VDW


# ---------------------------------------------------------------------------
# Samplers
# ---------------------------------------------------------------------------


@dataclass
class SurfaceSampler(_IterableSampler):
    """Proposes positions around a target, skipping any inside it.

    Each proposal goes through two tests. First the region: a position is drawn
    from `envelope`, normally a ball enclosing the target with some room to
    spare. Then the surface: `excluder` throws it away if it would sit inside
    the target's atoms. Drawing continues until one passes.

    Parameters
    ----------
    envelope : Envelope
        The region to draw positions from.
    excluder : Excluder
        The test that rejects positions inside the target.
    max_rejections : int, default=1000
        How many rejections in a row to accept before giving up. Stops a badly
        set-up run from looping forever; a region that lies entirely inside the
        target would otherwise never produce anything.

    See Also
    --------
    around : Build one for a target directly.
    SurfaceFollowingSampler : Concentrates proposals near the target's surface.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> sampler = SurfaceSampler(
    ...     Sphere(20.0, np.zeros(3), rng=rng),
    ...     Excluder(np.zeros((1, 3)), ["C"]),
    ... )
    >>> from itertools import islice
    >>> len(list(islice(sampler, 3)))
    3
    """

    envelope: Envelope
    excluder: Excluder
    max_rejections: int = 1000

    def sample(self) -> Placement:
        """Return one placement that clears the target.

        Returns
        -------
        Placement
            A position outside the target, with a random orientation.

        Raises
        ------
        maws.errors.SamplingError
            If `max_rejections` proposals in a row were all rejected.
        """
        for _ in range(self.max_rejections):
            proposal = self.envelope.sample()
            if self.excluder.is_clear(proposal.position):
                return proposal
        raise SamplingError(
            f"no position clear of the target after {self.max_rejections} "
            f"tries. The region being drawn from may lie inside the target: "
            f"try a larger reach, a smaller probe, or check the target's size."
        )

    def accepts(self, points: np.ndarray) -> bool:
        """Say whether the strand, once placed, clears the target.

        Parameters
        ----------
        points : numpy.ndarray
            Shape ``(N, 3)``. Where the strand's atoms ended up, in ångström.

        Returns
        -------
        bool
            True when no atom of the strand is inside the target.

        Notes
        -----
        :meth:`sample` already rejects positions inside the target, but the
        position it tests is where the middle of the strand goes. This is the
        same test applied to every atom, which is the question that actually
        matters and cannot be asked until the placement has been applied.
        """
        return self.excluder.all_clear(points)

    @classmethod
    def around(
        cls,
        system: BuiltSystem,
        role: str = "ligand",
        *,
        reach: float = 10.0,
        probe: float = 1.4,
        rng: np.random.Generator | None = None,
        max_rejections: int = 1000,
    ) -> SurfaceSampler:
        """Build a sampler that proposes positions around one chain.

        The region is a ball centred on the chain's balance point, big enough
        to contain every one of its atoms plus `reach` beyond the furthest.

        Parameters
        ----------
        system : maws.topology.BuiltSystem
            The built design. Its starting positions decide where the target
            is.
        role : str, default="ligand"
            Which chain to sample around.
        reach : float, default=10.0
            How far past the target's furthest atom the region extends, in
            ångström. Larger values consider positions further away, so the
            aptamer has more room but more proposals land in empty space.
        probe : float, default=1.4
            Radius of the ball used to find the target's surface, in ångström.
        rng : numpy.random.Generator, optional
            Where the randomness comes from.
        max_rejections : int, default=1000
            How many rejections in a row to accept before giving up.

        Returns
        -------
        SurfaceSampler
            A sampler ready to draw from.
        """
        rng = rng if rng is not None else np.random.default_rng()
        span, positions, elements, masses = _chain_atoms(system, role)
        centre = centre_of_mass(positions, masses)
        furthest = float(np.linalg.norm(positions - centre, axis=1).max())
        return cls(
            envelope=Sphere(furthest + reach, centre, rng=rng),
            excluder=Excluder(positions, elements, probe=probe),
            max_rejections=max_rejections,
        )


class SurfaceFollowingSampler(_IterableSampler):
    """Proposes positions in a thin layer wrapped around the target.

    Like :class:`SurfaceSampler`, but with a second rule: a proposal is also
    thrown away if it is further than `d_max` from the target's nearest atom.
    What survives is a shell of that thickness following the target's shape,
    including into its pockets. Proposals therefore concentrate where the
    aptamer could actually make contact, instead of spreading through the
    surrounding empty space.

    Parameters
    ----------
    positions : numpy.ndarray
        Shape ``(N, 3)``. Where the target's atoms are, in ångström.
    elements : sequence of str
        Chemical symbol per atom.
    masses : numpy.ndarray
        Shape ``(N,)``. Atomic mass per atom, in daltons, used to find the
        target's balance point.
    d_max : float, default=6.0
        How thick the shell is, in ångström. Smaller values hug the surface
        more closely but reject far more proposals to get there.
    probe : float, default=1.4
        Radius of the ball used to find the target's surface, in ångström.
    rng : numpy.random.Generator, optional
        Where the randomness comes from.
    max_rejections : int, default=50000
        How many rejections in a row to accept before giving up. Higher than
        :class:`SurfaceSampler`'s, because a thin shell is a small part of the
        ball proposals are drawn from.

    Raises
    ------
    maws.errors.ConfigurationError
        If `d_max` is not positive, which would leave no shell at all.
    """

    __slots__ = (
        "_bound",
        "_centre",
        "_d_max",
        "_excluder",
        "_max_rejections",
        "_rng",
        "_tree",
    )

    def __init__(
        self,
        positions: np.ndarray,
        elements: Sequence[str],
        masses: np.ndarray,
        *,
        d_max: float = 6.0,
        probe: float = 1.4,
        rng: np.random.Generator | None = None,
        max_rejections: int = 50_000,
    ) -> None:
        if d_max <= 0:
            raise ConfigurationError(f"d_max must be positive, got {d_max}")
        positions = np.asarray(positions, dtype=np.float64)
        self._centre = centre_of_mass(positions, masses)
        furthest = float(np.linalg.norm(positions - self._centre, axis=1).max())
        self._bound = furthest + d_max
        self._tree = KDTree(positions)
        self._excluder = Excluder(positions, elements, probe=probe)
        self._d_max = d_max
        self._max_rejections = max_rejections
        self._rng = rng if rng is not None else np.random.default_rng()

    def __repr__(self) -> str:
        return f"<SurfaceFollowingSampler shell {self._d_max} Å>"

    def sample(self) -> Placement:
        """Return one placement inside the shell around the target.

        Returns
        -------
        Placement
            A position within `d_max` of the target but outside its atoms.

        Raises
        ------
        maws.errors.SamplingError
            If `max_rejections` proposals in a row were all rejected.
        """
        for _ in range(self._max_rejections):
            distance = self._bound * self._rng.uniform(0.0, 1.0) ** (1.0 / 3.0)
            position = self._centre + distance * _random_direction(self._rng)
            nearest, _ = self._tree.query(position, k=1)
            if nearest > self._d_max:
                continue
            if not self._excluder.is_clear(position):
                continue
            return Placement(
                position=position,
                axis=_random_direction(self._rng),
                angle=float(self._rng.uniform(0.0, 2.0 * np.pi)),
            )
        raise SamplingError(
            f"no position in the shell around the target after "
            f"{self._max_rejections} tries. d_max may be too small, or the "
            f"target unusually densely packed."
        )

    def accepts(self, points: np.ndarray) -> bool:
        """Say whether the strand, once placed, clears the target.

        Parameters
        ----------
        points : numpy.ndarray
            Shape ``(N, 3)``. Where the strand's atoms ended up, in ångström.

        Returns
        -------
        bool
            True when no atom of the strand is inside the target.

        Notes
        -----
        Only the surface test is repeated here, not the shell test. Whether
        the middle of the strand sits within `d_max` of the target is what
        decides where to try putting it; once it is there, the far end of the
        strand reaching beyond the shell is not a reason to throw the shape
        away.
        """
        return self._excluder.all_clear(points)

    @classmethod
    def around(
        cls,
        system: BuiltSystem,
        role: str = "ligand",
        *,
        d_max: float = 6.0,
        probe: float = 1.4,
        rng: np.random.Generator | None = None,
        max_rejections: int = 50_000,
    ) -> SurfaceFollowingSampler:
        """Build a shell sampler around one chain of a built design.

        Parameters
        ----------
        system : maws.topology.BuiltSystem
            The built design.
        role : str, default="ligand"
            Which chain to sample around.
        d_max : float, default=6.0
            Shell thickness in ångström.
        probe : float, default=1.4
            Radius of the ball used to find the surface, in ångström.
        rng : numpy.random.Generator, optional
            Where the randomness comes from.
        max_rejections : int, default=50000
            How many rejections in a row to accept before giving up.

        Returns
        -------
        SurfaceFollowingSampler
            A sampler ready to draw from.
        """
        _, positions, elements, masses = _chain_atoms(system, role)
        return cls(
            positions,
            elements,
            masses,
            d_max=d_max,
            probe=probe,
            rng=rng,
            max_rejections=max_rejections,
        )


class FixedSampler(_IterableSampler):
    """Hands out placements from a list you wrote yourself.

    Useful in a test that needs to know exactly where the strand will be put,
    so the result can be worked out by hand and compared.

    Parameters
    ----------
    placements : sequence of Placement
        What to hand out, in order.
    cycle : bool, default=True
        What to do when the list runs out. True starts again from the
        beginning, so the sampler never stops. False raises instead, which
        catches a test that asked for more placements than it meant to.

    Raises
    ------
    maws.errors.ConfigurationError
        If `placements` is empty.

    Examples
    --------
    >>> import numpy as np
    >>> only = Placement(np.zeros(3), np.array([0.0, 0, 1]), 0.0)
    >>> sampler = FixedSampler([only])
    >>> next(sampler) is only and next(sampler) is only
    True
    """

    __slots__ = ("_cycle", "_index", "_placements")

    def __init__(self, placements: Sequence[Placement], *, cycle: bool = True) -> None:
        if not placements:
            raise ConfigurationError("a FixedSampler needs at least one placement")
        self._placements = list(placements)
        self._cycle = cycle
        self._index = 0

    def __repr__(self) -> str:
        return f"<FixedSampler {len(self._placements)} placements>"

    def sample(self) -> Placement:
        """Return the next placement from the list.

        Returns
        -------
        Placement
            The next one, starting again from the beginning if `cycle` is set.

        Raises
        ------
        maws.errors.SamplingError
            If the list has run out and `cycle` is False.
        """
        if self._index >= len(self._placements):
            if not self._cycle:
                raise SamplingError(
                    f"all {len(self._placements)} placements have been used "
                    f"and cycle is off"
                )
            self._index = 0
        placement = self._placements[self._index]
        self._index += 1
        return placement


def _chain_atoms(
    system: BuiltSystem, role: str
) -> tuple[AtomRange, np.ndarray, tuple[str, ...], np.ndarray]:
    """Pull out one chain's positions, element symbols and masses.

    Parameters
    ----------
    system : maws.topology.BuiltSystem
        The built design.
    role : str
        Which chain to pull out.

    Returns
    -------
    span : maws.values.AtomRange
        Which atoms of the design the chain is.
    positions : numpy.ndarray
        Shape ``(n, 3)``, in ångström.
    elements : tuple of str
        Chemical symbol per atom.
    masses : numpy.ndarray
        Shape ``(n,)``, in daltons.
    """
    span = system.chain(role).span
    positions = np.asarray(system.pose.atoms(span), dtype=np.float64)
    elements = system.elements[span.as_slice()]
    masses = np.asarray(system.masses[span.as_slice()], dtype=np.float64)
    return span, positions, elements, masses


def make_sampler(
    system: BuiltSystem,
    role: str = "ligand",
    *,
    mode: Literal["sphere", "surface-following"] = "sphere",
    reach: float = 10.0,
    d_max: float = 6.0,
    probe: float = 1.4,
    rng: np.random.Generator | None = None,
) -> Sampler:
    """Build the sampler named by `mode`, ready to draw from.

    Parameters
    ----------
    system : maws.topology.BuiltSystem
        The built design.
    role : str, default="ligand"
        Which chain to sample around.
    mode : {"sphere", "surface-following"}, default="sphere"
        Which sampler to build. ``"sphere"`` draws from a ball enclosing the
        target and skips positions inside it — simple and fast, and the right
        default. ``"surface-following"`` also skips positions far from the
        target, concentrating proposals near its surface at the cost of
        rejecting many more.
    reach : float, default=10.0
        Used by ``"sphere"``: how far past the target the ball extends, in
        ångström.
    d_max : float, default=6.0
        Used by ``"surface-following"``: shell thickness in ångström.
    probe : float, default=1.4
        Radius of the ball used to find the target's surface, in ångström.
        Used by both.
    rng : numpy.random.Generator, optional
        Where the randomness comes from.

    Returns
    -------
    Sampler
        Either a :class:`SurfaceSampler` or a :class:`SurfaceFollowingSampler`.

    Raises
    ------
    maws.errors.ConfigurationError
        If `mode` is not one of the two names, or `reach` or `probe` is
        negative.
    """
    # Every setting is checked whichever mode was chosen. Checking only the
    # ones the chosen mode reads would let a typo in an unused argument pass
    # silently, and the caller would never learn the value was ignored.
    if probe < 0:
        raise ConfigurationError(f"probe must not be negative, got {probe}")
    if reach < 0:
        raise ConfigurationError(f"reach must not be negative, got {reach}")
    if d_max <= 0:
        raise ConfigurationError(f"d_max must be positive, got {d_max}")
    if mode == "sphere":
        return SurfaceSampler.around(system, role, reach=reach, probe=probe, rng=rng)
    if mode == "surface-following":
        return SurfaceFollowingSampler.around(
            system, role, d_max=d_max, probe=probe, rng=rng
        )
    raise ConfigurationError(
        f"mode must be 'sphere' or 'surface-following', got {mode!r}"
    )

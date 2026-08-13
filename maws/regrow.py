r"""
maws.regrow
===========

Adding a nucleotide to a strand without losing the shape it already had.

The search grows an aptamer one residue at a time. Growing it means building a
new structure, because the strand now has more atoms than before and the
builder has to work out what the new ones are. But the shape found so far is
the whole result of the search up to that point, and a fresh build comes back
in whatever shape the builder happens to produce. Throwing the old shape away
and starting again would discard every step.

So the two are combined. The freshly built structure supplies correct internal
geometry for the new residue and its join to the rest of the strand. The old
structure supplies the positions of every atom that already existed. Fitting
one onto the other is a superposition: find the single rotation and shift that
best carries the new structure's copy of the old atoms onto where those atoms
actually were, apply it to the whole strand, then put the old atoms back
exactly.

The fit is Kabsch's algorithm, the standard way of superposing two sets of
matched points.

Examples
--------
>>> import numpy as np
>>> mobile = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
>>> target = mobile + np.array([5.0, 0, 0])
>>> rotation, shift = kabsch(mobile, target)
>>> np.allclose(mobile @ rotation.T + shift, target)
True
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from maws.errors import ConfigurationError
from maws.pose import ChainView, Pose
from maws.values import AtomRange, Direction

if TYPE_CHECKING:  # pragma: no cover - needed by type checkers only
    from maws.build import Builder
    from maws.topology import BuiltSystem

__all__ = ["Grown", "common_window", "grow_chain", "kabsch", "splice"]


class Grown(NamedTuple):
    """A structure that has just gained a residue, and its atom positions.

    Parameters
    ----------
    system : maws.topology.BuiltSystem
        The rebuilt structure. Its chains have new spans, because the strand
        now has more atoms.
    pose : maws.pose.Pose
        Positions for it: the old ones wherever an atom already existed, and
        freshly built ones for the residue that was added.
    """

    system: BuiltSystem
    pose: Pose


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r"""Find the rotation and shift that best carries `mobile` onto `target`.

    Parameters
    ----------
    mobile : numpy.ndarray
        Shape ``(N, 3)`` positions to be moved, in ångström.
    target : numpy.ndarray
        Shape ``(N, 3)`` positions to move them onto, in ångström. Row *i* must
        describe the same atom as row *i* of `mobile`.

    Returns
    -------
    rotation : numpy.ndarray
        Shape ``(3, 3)``. Apply to positions stored one atom per row as
        ``xyz @ rotation.T``.
    shift : numpy.ndarray
        Shape ``(3,)``, in ångström, applied after the rotation.

    Raises
    ------
    maws.errors.ConfigurationError
        If the two arrays are not the same shape, or hold fewer than three
        points. Fewer than three leaves the orientation undetermined: two
        points can be spun about the line joining them, and one about
        anything.

    Notes
    -----
    Kabsch's algorithm [1]_. The rotation minimising the sum of squared
    distances between matched points is found from the singular value
    decomposition of their covariance matrix,

    .. math::
        H = (\mathbf{m} - \bar{\mathbf{m}})^{T}(\mathbf{t} - \bar{\mathbf{t}})
        = U \Sigma V^{T},
        \qquad
        R = V \operatorname{diag}(1, 1, d) U^{T}

    where :math:`d = \operatorname{sign}(\det(V U^{T}))`. That last factor
    matters: without it the decomposition can return a reflection, which fits
    the points just as well numerically but turns the molecule into its mirror
    image.

    References
    ----------
    .. [1] Kabsch, W. (1976). "A solution for the best rotation to relate two
           sets of vectors". Acta Crystallographica A32, 922-923.

    Examples
    --------
    A quarter turn about z, recovered exactly:

    >>> import numpy as np
    >>> mobile = np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> target = np.array([[0.0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    >>> rotation, shift = kabsch(mobile, target)
    >>> np.allclose(mobile @ rotation.T + shift, target)
    True
    """
    mobile = np.asarray(mobile, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if mobile.shape != target.shape:
        raise ConfigurationError(
            f"the two point sets must have the same shape, got "
            f"{mobile.shape} and {target.shape}"
        )
    if mobile.shape[0] < 3:
        raise ConfigurationError(
            f"a superposition needs at least three points to pin down an "
            f"orientation, got {mobile.shape[0]}"
        )

    mobile_centre = mobile.mean(axis=0)
    target_centre = target.mean(axis=0)
    covariance = (mobile - mobile_centre).T @ (target - target_centre)
    left, _, right = np.linalg.svd(covariance)
    handedness = np.sign(np.linalg.det(right.T @ left.T))
    correction = np.diag([1.0, 1.0, handedness])
    rotation = right.T @ correction @ left.T
    return rotation, target_centre - mobile_centre @ rotation.T


def common_window(
    fresh: AtomRange, old: AtomRange, *, align: str
) -> tuple[AtomRange, AtomRange]:
    """Return the atoms a residue has in both structures, as matching spans.

    A residue at the end of a strand is not quite the same molecule as the same
    residue in the middle of one: an end residue carries an extra group. So
    when a strand grows, the residue that used to be outermost changes size,
    and its atoms in the old structure no longer line up one-for-one with its
    atoms in the new one.

    The atoms they do have in common are contiguous, because the group that
    appears or disappears sits at one end of the residue's atom list: a leading
    phosphate at the 5' end, a trailing hydrogen at the 3' end.

    Parameters
    ----------
    fresh : maws.values.AtomRange
        The residue in the new structure.
    old : maws.values.AtomRange
        The same residue in the old structure.
    align : {"start", "end"}
        Which end of the two spans to line up. Use ``"end"`` when the group
        that changed sits at the beginning of the residue's atoms, and
        ``"start"`` when it sits at the end.

    Returns
    -------
    fresh_window, old_window : maws.values.AtomRange
        Equal-length spans naming the same atoms in the two structures.

    Raises
    ------
    maws.errors.ConfigurationError
        If `align` is neither ``"start"`` nor ``"end"``.

    Examples
    --------
    A residue that lost its last atom, lined up from the start:

    >>> common_window(AtomRange(0, 30), AtomRange(0, 31), align="start")
    (AtomRange(start=0, stop=30), AtomRange(start=0, stop=30))

    One that gained three atoms at the front, lined up from the end:

    >>> common_window(AtomRange(0, 34), AtomRange(0, 31), align="end")
    (AtomRange(start=3, stop=34), AtomRange(start=0, stop=31))
    """
    if align not in ("start", "end"):
        raise ConfigurationError(f"align must be 'start' or 'end', got {align!r}")
    shared = min(len(fresh), len(old))
    if align == "start":
        return (
            AtomRange(fresh.start, fresh.start + shared),
            AtomRange(old.start, old.start + shared),
        )
    return (
        AtomRange(fresh.stop - shared, fresh.stop),
        AtomRange(old.stop - shared, old.stop),
    )


def splice(
    fresh: Pose,
    old: Pose,
    *,
    fresh_chain: AtomRange,
    matches: list[tuple[AtomRange, AtomRange]],
    anchor: int,
    others: list[tuple[AtomRange, AtomRange]],
) -> Pose:
    """Fit a freshly built strand onto the shape an older one already had.

    Parameters
    ----------
    fresh : maws.pose.Pose
        Positions from the new build. Correct internal geometry, arbitrary
        overall placement, and always fully extended.
    old : maws.pose.Pose
        Positions from before the strand grew. The shape worth keeping.
    fresh_chain : maws.values.AtomRange
        The whole grown chain, in the new structure.
    matches : list of (maws.values.AtomRange, maws.values.AtomRange)
        One entry per residue that existed before, giving the atoms it has in
        the new structure and the same atoms in the old one. Both spans of a
        pair must be the same length.
    anchor : int
        Which entry of `matches` to compute the fit from — the residue the new
        one is bonded to. A negative value means no fit is possible, which is
        the case when the strand was empty. Then, and when that entry names
        fewer than three atoms, `fresh` is returned as it came apart from the
        chains named in `others`.
    others : list of (maws.values.AtomRange, maws.values.AtomRange)
        For every chain that did not grow, where it sits in the new structure
        and where it sat in the old one, in that order. Their atoms are copied
        across unchanged.

    Returns
    -------
    maws.pose.Pose
        Positions for the new structure.

    Raises
    ------
    maws.errors.ConfigurationError
        If any pair of spans in `matches` describes different numbers of atoms.

    See Also
    --------
    common_window : Works out the spans for one residue.
    grow_chain : Assembles the arguments and calls this.

    Notes
    -----
    Every atom that existed before ends up at exactly the position it had, not
    merely close to it. The fit decides only where the new residue goes, and
    where the handful of atoms that have no counterpart go — the group that
    appears when a residue stops being at the end of the strand.

    .. warning::
        The fit is computed from the neighbouring residue alone, never from the
        whole strand. A fresh build is always fully extended, while the shape
        being fitted onto has been folded by many random turns. No single rigid
        motion carries one onto the other, so a fit over the whole strand
        spreads its error everywhere — and all of that error lands on the new
        residue, which is the one thing being positioned.
    """
    for index, (fresh_span, old_span) in enumerate(matches):
        if len(fresh_span) != len(old_span):
            raise ConfigurationError(
                f"residue {index} covers {len(fresh_span)} atoms in the new "
                f"structure and {len(old_span)} in the old one; the two spans "
                f"must name the same atoms"
            )

    xyz = fresh.xyz.copy()
    for new_span, old_span in others:
        xyz[new_span.as_slice()] = old.xyz[old_span.as_slice()]

    if not (0 <= anchor < len(matches) and len(matches[anchor][0]) >= 3):
        # Nothing to fit onto: either the strand was empty, or the residue the
        # new one joins has fewer than three atoms in common with its old self,
        # which leaves the orientation undetermined. Copying old positions in
        # anyway would drop a few atoms onto coordinates the rest of the fresh
        # build knows nothing about, tearing the residue apart. The fresh build
        # is used as it came instead.
        return Pose(xyz, fresh.system)

    fresh_anchor, old_anchor = matches[anchor]
    rotation, shift = kabsch(
        fresh.xyz[fresh_anchor.as_slice()], old.xyz[old_anchor.as_slice()]
    )
    chain = fresh_chain.as_slice()
    xyz[chain] = fresh.xyz[chain] @ rotation.T + shift

    for fresh_span, old_span in matches:
        xyz[fresh_span.as_slice()] = old.xyz[old_span.as_slice()]

    return Pose(xyz, fresh.system)


def grow_chain(
    system: BuiltSystem,
    pose: Pose,
    *,
    role: str,
    token: str,
    direction: Direction,
    builder: Builder,
) -> Grown:
    """Add one nucleotide to a strand, keeping the shape it already had.

    Parameters
    ----------
    system : maws.topology.BuiltSystem
        The structure before growing.
    pose : maws.pose.Pose
        The positions to keep. Normally the best shape found for `system` so
        far, which is not the same as the one it was built with.
    role : str
        Which chain to grow.
    token : str
        The nucleotide to add, as written, e.g. ``"G"``.
    direction : {"3prime", "5prime"}
        Which end of the strand to add it to.
    builder : maws.build.Builder
        What builds the new structure.

    Returns
    -------
    Grown
        The rebuilt structure and its positions.

    See Also
    --------
    splice : Does the fitting.
    maws.search.grow_aptamer : Calls this once per candidate.

    Notes
    -----
    Growing at the 5' end shifts every existing atom of the strand along by the
    size of the new residue, because atoms are numbered from the 5' end. The
    two ends therefore need different bookkeeping, which is why `direction` is
    needed here and not just when the sequence is edited.

    Examples
    --------
    >>> from maws.build import FakeBuilder
    >>> from maws.forcefield import ForceField
    >>> from maws.libraries import rna
    >>> from maws.topology import Assembly
    >>> builder = FakeBuilder()
    >>> forcefield = ForceField.for_target("RNA", "protein")
    >>> system = builder.build(
    ...     Assembly().with_aptamer(rna(), "G").with_ligand_stub(5), forcefield
    ... )
    >>> grown = grow_chain(
    ...     system,
    ...     system.pose,
    ...     role="aptamer",
    ...     token="A",
    ...     direction="3prime",
    ...     builder=builder,
    ... )
    >>> str(grown.system.chain("aptamer").sequence)
    'G A'
    """
    chain = system.chain(role)
    grown_sequence = chain.sequence.grown(token, direction)
    new_system = builder.build(
        system.assembly.with_sequence(role, grown_sequence), system.forcefield
    )

    old_span = chain.span
    new_chain = new_system.chain(role)
    new_span = new_chain.span
    n_added = len(new_span) - len(old_span)
    if n_added <= 0:
        raise ConfigurationError(
            f"growing chain {role!r} by {token!r} did not add any atoms; the "
            f"residue library may be missing that nucleotide"
        )

    matches, anchor = _match_residues(chain, new_chain, direction)
    others = [
        (new_system.chain(other).span, system.chain(other).span)
        for other in system.spans
        if other != role
    ]
    spliced = splice(
        new_system.pose,
        pose,
        fresh_chain=new_span,
        matches=matches,
        anchor=anchor,
        others=others,
    )
    return Grown(new_system, spliced)


def _match_residues(
    old_chain: ChainView, new_chain: ChainView, direction: Direction
) -> tuple[list[tuple[AtomRange, AtomRange]], int]:
    """Pair up each residue that already existed with its place in the rebuild.

    Parameters
    ----------
    old_chain : maws.pose.ChainView
        The strand before it grew.
    new_chain : maws.pose.ChainView
        The same strand after, one residue longer.
    direction : {"3prime", "5prime"}
        Which end the new residue was added at.

    Returns
    -------
    matches : list of (maws.values.AtomRange, maws.values.AtomRange)
        One entry per old residue, in old order, giving its atoms in the new
        structure and the same atoms in the old one.
    anchor : int
        Which entry to compute the fit from — the residue the new one is bonded
        to. ``-1`` when the strand was empty and there is nothing to fit onto.

    Notes
    -----
    Growing at the 3' end leaves the existing atoms where they are and adds
    after them, so old residue *i* is new residue *i*. Growing at the 5' end
    pushes everything along by one residue, because atoms are numbered from the
    5' end, so old residue *i* is new residue *i + 1*.

    One residue also changes size. The residue that was at the growing end
    stops being an end residue, and end residues carry an extra group: a
    trailing hydrogen at the 3' end, a leading phosphate at the 5' end. So the
    two spans for that residue are lined up from whichever end did not change.
    """
    shift = 0 if direction == "3prime" else 1
    changed = old_chain.n_residues - 1 if direction == "3prime" else 0

    matches = []
    for index in range(old_chain.n_residues):
        old_residue = old_chain.residue(index)
        new_residue = new_chain.residue(index + shift)
        # Only the residue at the growing end changes size, and the group it
        # gains sits at the end of its atom list going 3' and at the start
        # going 5'. Every other residue matches atom for atom, so the choice
        # of alignment makes no difference to them.
        align = "end" if index == changed and direction == "5prime" else "start"
        matches.append(common_window(new_residue.span, old_residue.span, align=align))

    anchor = changed if old_chain.n_residues else -1
    return matches, anchor

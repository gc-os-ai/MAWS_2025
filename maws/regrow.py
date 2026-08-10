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
from maws.pose import Pose
from maws.values import AtomRange, Direction

if TYPE_CHECKING:  # pragma: no cover - needed by type checkers only
    from maws.build import Builder
    from maws.topology import BuiltSystem

__all__ = ["Grown", "grow_chain", "kabsch", "splice"]


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


def splice(
    fresh: Pose,
    old: Pose,
    *,
    fresh_chain: AtomRange,
    fresh_kept: AtomRange,
    old_kept: AtomRange,
    others: list[tuple[AtomRange, AtomRange]],
) -> Pose:
    """Fit a freshly built strand onto the shape an older one already had.

    Parameters
    ----------
    fresh : maws.pose.Pose
        Positions from the new build. Correct internal geometry, arbitrary
        overall placement.
    old : maws.pose.Pose
        Positions from before the strand grew. The shape worth keeping.
    fresh_chain : maws.values.AtomRange
        The whole grown chain, in the new structure.
    fresh_kept : maws.values.AtomRange
        The part of that chain whose atoms already existed, in the new
        structure.
    old_kept : maws.values.AtomRange
        The same atoms, in the old structure. Must be the same length as
        `fresh_kept`.
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
        If `fresh_kept` and `old_kept` describe different numbers of atoms.

    Notes
    -----
    Every atom that already existed ends up at exactly the position it had
    before, not merely close to it. The superposition decides only where the
    newly added residue goes.

    When nothing was kept — the first residue of an empty strand — there is
    nothing to fit onto, and the freshly built positions are used as they are.
    """
    if len(fresh_kept) != len(old_kept):
        raise ConfigurationError(
            f"the kept part must be the same size in both structures, got "
            f"{len(fresh_kept)} and {len(old_kept)} atoms"
        )

    xyz = fresh.xyz.copy()
    for new_span, old_span in others:
        xyz[new_span.as_slice()] = old.xyz[old_span.as_slice()]

    if len(fresh_kept) >= 3:
        rotation, shift = kabsch(
            fresh.xyz[fresh_kept.as_slice()], old.xyz[old_kept.as_slice()]
        )
        chain = fresh_chain.as_slice()
        xyz[chain] = fresh.xyz[chain] @ rotation.T + shift
        xyz[fresh_kept.as_slice()] = old.xyz[old_kept.as_slice()]

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

    if direction == "3prime":
        fresh_kept = AtomRange(new_span.start, new_span.start + len(old_span))
        added = AtomRange(fresh_kept.stop, new_span.stop)
    else:
        added = AtomRange(new_span.start, new_span.start + n_added)
        fresh_kept = AtomRange(added.stop, new_span.stop)

    others = [
        (new_system.chain(other).span, system.chain(other).span)
        for other in system.spans
        if other != role
    ]
    spliced = splice(
        new_system.pose,
        pose,
        fresh_chain=new_span,
        fresh_kept=fresh_kept,
        old_kept=old_span,
        others=others,
    )
    return Grown(new_system, spliced)

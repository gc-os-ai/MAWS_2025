"""
maws.api
========

Designing an aptamer in one call.

This is the shortest way in. :func:`design` takes a target molecule's structure
file and gives back the strand it found, hiding the residue libraries, the
builder, the sampler and the scorer behind sensible defaults.

Everything it does is available separately if the defaults do not suit: see
:func:`maws.search.grow_aptamer` for the search itself, which reports its
progress step by step and can be stopped part-way.

Examples
--------
>>> result = design("target.pdb", length=15)  # doctest: +SKIP
>>> result.sequence  # doctest: +SKIP
'G A U C G A U C G A U C G A U'
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from maws.build import Builder, LeapBuilder
from maws.forcefield import AptamerType, ForceField, MoleculeType
from maws.io.pdb_cleaner import resolve_pdb_path
from maws.libraries import dna, rna
from maws.sampling import make_sampler
from maws.scoring import Scorer, entropy_score
from maws.search import (
    SearchFinished,
    StepEvent,
    grow_aptamer,
    openmm_energy,
)
from maws.topology import Assembly

if TYPE_CHECKING:  # pragma: no cover - needed by type checkers only
    from maws.pose import Pose
    from maws.topology import BuiltSystem

__all__ = ["MawsResult", "design"]

SamplingMode = Literal["sphere", "surface-following"]
"""Which region proposed positions for the aptamer are drawn from."""


@dataclass(frozen=True, slots=True)
class MawsResult:
    """The strand a design run settled on, and how it scored.

    Parameters
    ----------
    sequence : str
        The designed strand, written 5' to 3', e.g. ``"G A U C"``.
    energy : float
        Potential energy of the best arrangement found, in kJ/mol.
    entropy : float
        The score that chose this strand. Lower is better. See
        :func:`maws.scoring.entropy_score`.
    steps : int
        How many nucleotides were added.
    success : bool, default=True
        Whether the run finished as intended. False means `sequence` is
        whatever the run had reached when it stopped, which may be shorter than
        asked for.
    message : str, default=""
        What went wrong, when `success` is False. Empty otherwise.
    system : maws.topology.BuiltSystem, optional
        The final structure. Present after a real run; useful for writing the
        result out or scoring it again.
    pose : maws.pose.Pose, optional
        Atom positions of the best arrangement found.

    See Also
    --------
    design : Produces one of these.

    Notes
    -----
    `success` and `message` follow the same idea as
    :class:`scipy.optimize.OptimizeResult`: a run that found nothing usable
    should be distinguishable from one that worked, without inspecting the
    numbers.

    Examples
    --------
    >>> result = MawsResult("G A U", -1421.9, -0.072, steps=3)
    >>> result.length
    3
    >>> print(result)
    G A U  (3 nt, E=-1421.90 kJ/mol, S=-0.072000)
    """

    sequence: str
    energy: float
    entropy: float
    steps: int
    success: bool = True
    message: str = ""
    system: BuiltSystem | None = field(default=None, repr=False, compare=False)
    pose: Pose | None = field(default=None, repr=False, compare=False)

    @property
    def length(self) -> int:
        """int : How many nucleotides the designed strand has."""
        return len(self.sequence.split())

    def __str__(self) -> str:
        return (
            f"{self.sequence}  ({self.length} nt, "
            f"E={self.energy:.2f} kJ/mol, S={self.entropy:.6f})"
        )


def design(
    target: str | Path,
    *,
    length: int = 15,
    aptamer: AptamerType = "RNA",
    molecule: MoleculeType = "protein",
    samples: int = 5000,
    first_samples: int | None = None,
    beta: float = 0.01,
    salt_conc: float = 0.15,
    reach: float = 10.0,
    probe: float = 1.4,
    d_max: float = 6.0,
    sampling: SamplingMode = "sphere",
    relax_iterations: int = 50,
    seed: int | None = None,
    scorer: Scorer = entropy_score,
    builder: Builder | None = None,
    clean_pdb: bool = False,
    keep_chains: str = "all",
    remove_h: bool = False,
    drop_hetatm: bool = False,
    on_event: Callable[[StepEvent], None] | None = None,
    logger: logging.Logger | None = None,
) -> MawsResult:
    """design(target, *, length=15, aptamer="RNA", molecule="protein", ...)

    Design an aptamer against a target molecule.

    Parameters
    ----------
    target : str or pathlib.Path
        Structure file for the molecule the aptamer should bind to, in PDB
        format.
    length : int, default=15
        How many nucleotides the finished strand should have.
    aptamer : {"RNA", "DNA"}, default="RNA"
        Which nucleic acid to build the strand from.
    molecule : {"protein", "organic", "lipid"}, default="protein"
        What the target is. Decides which parameters describe it, and whether
        they have to be worked out first.
    samples : int, default=5000
        How many shapes to try per candidate at each growth step. This is the
        main cost: the total work is roughly ``8 * samples * length`` energy
        evaluations. Lower it to perhaps 50 for a quick check.
    first_samples : int, optional
        How many shapes to try on the first step, which also searches over
        where to put the strand. Defaults to `samples`.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ.
    salt_conc : float, default=0.15
        Concentration of dissolved salt in mol/L. Damps the pull between
        charged atoms; the default is roughly the saltiness of blood.
    reach : float, default=10.0
        How far past the target's furthest atom the strand may be placed, in
        ångström. Only used when `sampling` is ``"sphere"``.
    probe : float, default=1.4
        Radius in ångström of the ball rolled over the target to find its
        surface. 1.4 Å is the size of a water molecule, the usual choice.
    d_max : float, default=6.0
        Thickness in ångström of the shell positions are drawn from. Only used
        when `sampling` is ``"surface-following"``.
    sampling : {"sphere", "surface-following"}, default="sphere"
        Which region to draw positions from. ``"sphere"`` fills a ball around
        the target; ``"surface-following"`` keeps to a shell hugging its
        surface, which concentrates the search where contact is possible but
        rejects far more proposals.
    relax_iterations : int, default=50
        How many nudge-and-settle rounds to run after each nucleotide is joined
        on. Zero skips it and is much faster, at the cost of leaving strain at
        the join.
    seed : int, optional
        Fixes the randomness, so the same call gives the same answer. Left out,
        every run differs.
    scorer : maws.scoring.Scorer, default=entropy_score
        Reduces a candidate's energies to one number, lower being better.
    builder : maws.build.Builder, optional
        What builds the structures. Defaults to
        :class:`~maws.build.LeapBuilder`, which needs AmberTools installed.
    clean_pdb : bool, default=False
        Whether to tidy the target file before building: a downloaded
        structure often carries several copies of the molecule, water, and
        other records that would confuse the builder. Only applied when
        `molecule` is ``"protein"``.
    keep_chains : str, default="all"
        Which chains of the target to keep when cleaning: ``"all"``,
        ``"one"``, or a comma-separated list such as ``"A,B"``.
    remove_h : bool, default=False
        Whether to strip hydrogens while cleaning.
    drop_hetatm : bool, default=False
        Whether to drop records for anything that is not part of the molecule
        proper, such as bound water or ions.
    on_event : callable, optional
        Called with every :class:`~maws.search.StepEvent` as the search
        produces it. Use it to print progress or write files.
    logger : logging.Logger, optional
        Where to send progress messages. Defaults to this module's logger.

    Returns
    -------
    MawsResult
        The designed strand, its score, and the structure it was found in.

    Raises
    ------
    maws.errors.ConfigurationError
        If an argument is not one of the allowed values, or the target file
        cannot be read.
    maws.errors.ToolchainError
        If AmberTools is needed and is not installed.

    See Also
    --------
    maws.search.grow_aptamer : The search, reporting step by step.
    maws.AptamerDesigner : The same run configured as an object.

    Notes
    -----
    .. warning::
        A run at the default settings performs on the order of a million
        energy evaluations and takes hours. Try `samples` of 50 and `length`
        of 3 first, to confirm the target file and the installation are sound.

    Examples
    --------
    >>> result = design("target.pdb", length=4, samples=50, seed=0)
    ... # doctest: +SKIP
    >>> result.length  # doctest: +SKIP
    4
    """
    log = logger if logger is not None else logging.getLogger(__name__)
    rng = np.random.default_rng(seed)

    path, original = resolve_pdb_path(
        str(target),
        molecule,
        clean_pdb=clean_pdb,
        keep_chains=keep_chains,
        remove_h=remove_h,
        drop_hetatm=drop_hetatm,
        logger=log,
    )
    log.debug("target file: given %s, using %s", original, path)

    forcefield = ForceField.for_target(aptamer, molecule, salt_conc=salt_conc)
    library = rna() if aptamer == "RNA" else dna()
    assembly = Assembly().with_aptamer(library).with_ligand(path, forcefield)

    used_builder = builder if builder is not None else LeapBuilder()
    base = used_builder.build(assembly, forcefield)
    sampler = make_sampler(
        base,
        "ligand",
        mode=sampling,
        reach=reach,
        d_max=d_max,
        probe=probe,
        rng=rng,
    )

    events = grow_aptamer(
        base,
        energy=openmm_energy(),
        builder=used_builder,
        sampler=sampler,
        n_nucleotides=length,
        first_samples=first_samples if first_samples is not None else samples,
        samples=samples,
        beta=beta,
        scorer=scorer,
        relax_iterations=relax_iterations,
        rng=rng,
    )
    return collect(events, on_event=on_event)


def collect(
    events: Any,
    *,
    on_event: Callable[[StepEvent], None] | None = None,
) -> MawsResult:
    """Run a search to the end and return only its answer.

    Parameters
    ----------
    events : iterable of maws.search.StepEvent
        The stream :func:`maws.search.grow_aptamer` produces.
    on_event : callable, optional
        Called with each event as it arrives, before it is discarded.

    Returns
    -------
    MawsResult
        Built from the final :class:`~maws.search.SearchFinished` event. If the
        stream ends without one, the result has ``success=False`` and a
        `message` saying so.

    Examples
    --------
    >>> import numpy as np
    >>> from maws.build import FakeBuilder
    >>> from maws.forcefield import ForceField
    >>> from maws.libraries import rna
    >>> from maws.search import grow_aptamer, stub_energy
    >>> from maws.topology import Assembly
    >>> builder = FakeBuilder()
    >>> base = builder.build(
    ...     Assembly().with_aptamer(rna()).with_ligand_stub(20),
    ...     ForceField.for_target("RNA", "protein"),
    ... )
    >>> result = collect(
    ...     grow_aptamer(
    ...         base,
    ...         energy=stub_energy(),
    ...         builder=builder,
    ...         n_nucleotides=2,
    ...         first_samples=4,
    ...         samples=4,
    ...         rng=np.random.default_rng(0),
    ...     )
    ... )
    >>> result.success, result.length
    (True, 2)
    """
    finished: SearchFinished | None = None
    for event in events:
        if on_event is not None:
            on_event(event)
        if isinstance(event, SearchFinished):
            finished = event

    if finished is None:
        return MawsResult(
            sequence="",
            energy=float("nan"),
            entropy=float("nan"),
            steps=0,
            success=False,
            message="the search stopped before it finished",
        )

    winner = finished.winner
    return MawsResult(
        sequence=str(winner.sequence),
        energy=winner.energy,
        entropy=winner.entropy,
        steps=finished.steps,
        system=winner.system,
        pose=winner.pose,
    )

"""
maws.search
===========

Growing an aptamer one nucleotide at a time.

This is the algorithm MAWS is for. It starts with a target molecule and an
empty strand, and repeats a single step until the strand is as long as asked
for:

1. propose every way of extending the strand by one nucleotide — each of the
   four nucleotides, at either end, so eight candidates;
2. for each candidate, try thousands of shapes and score how promising it looks
   with :func:`~maws.scoring.entropy_score`;
3. keep the best candidate and throw the rest away.

The first step is different, because there is nothing to extend yet. Instead
each of the four nucleotides is tried on its own, and each shape tried also
puts the single residue somewhere new next to the target.

The result is reported as a stream of events rather than returned at the end.
That way one loop serves every purpose: the command-line program turns events
into progress lines and files, :func:`maws.api.design` keeps only the last one,
and a caller wanting to stop early simply stops reading.

Examples
--------
Running the whole search against stand-ins, with nothing installed:

>>> from maws.build import FakeBuilder
>>> from maws.forcefield import ForceField
>>> from maws.libraries import rna
>>> from maws.topology import Assembly
>>> import numpy as np
>>> builder = FakeBuilder()
>>> forcefield = ForceField.for_target("RNA", "protein")
>>> base = builder.build(
...     Assembly().with_aptamer(rna()).with_ligand_stub(20), forcefield
... )
>>> events = list(
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
>>> str(events[-1].winner.sequence).count(" ") + 1
2
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from itertools import islice

import numpy as np

from maws.build import Builder
from maws.energy import EnergyModel, StubEnergy
from maws.errors import ConfigurationError
from maws.pose import ChainView, Pose
from maws.regrow import grow_chain
from maws.relax import perturb_and_minimize
from maws.sampling import Sampler, SurfaceSampler, TorsionAngles
from maws.scoring import Scorer, entropy_score
from maws.topology import BuiltSystem
from maws.values import Direction, NucleotideSequence, Torsion

__all__ = [
    "Candidate",
    "CandidateScored",
    "EnergyFactory",
    "SearchFinished",
    "StepCompleted",
    "StepEvent",
    "StepStarted",
    "grow_aptamer",
    "openmm_energy",
    "stub_energy",
]

EnergyFactory = Callable[[BuiltSystem], EnergyModel]
"""Builds a scorer for a structure.

The strand gains atoms at every step, so each candidate is a different
structure and needs its own scorer. A factory is therefore passed to
:func:`grow_aptamer` rather than a scorer itself.

See Also
--------
openmm_energy : A factory that scores with real physics.
stub_energy : A factory that scores with a stand-in, for tests.
"""

DEFAULT_TORSIONS = 4
"""How many turnable bonds of a residue the search varies.

The four backbone bonds of a nucleotide. Varying these changes the strand's
fold; the remaining bonds move only the base, which matters much less to
whether the strand binds.
"""


def stub_energy(
    aptamer: str = "aptamer", ligand: str = "ligand", *, scale: float = 1.0
) -> EnergyFactory:
    """Return a factory producing stand-in scorers, for tests.

    Parameters
    ----------
    aptamer, ligand : str
        Names of the two chains whose separation is scored.
    scale : float, default=1.0
        Multiplies the score. See :class:`~maws.energy.StubEnergy`.

    Returns
    -------
    EnergyFactory
        A callable that builds a :class:`~maws.energy.StubEnergy` for whatever
        structure it is handed.

    Examples
    --------
    >>> factory = stub_energy()
    >>> callable(factory)
    True
    """

    def make(system: BuiltSystem) -> EnergyModel:
        return StubEnergy(
            system.chain(aptamer).span, system.chain(ligand).span, scale=scale
        )

    return make


def openmm_energy(*, platform: str | None = None) -> EnergyFactory:
    """Return a factory producing real scorers backed by OpenMM.

    Parameters
    ----------
    platform : str, optional
        Which OpenMM compute backend to use, e.g. ``"CPU"``. By default the
        fastest available is chosen.

    Returns
    -------
    EnergyFactory
        A callable that builds an :class:`~maws.energy.OpenMMEnergy` for
        whatever structure it is handed.

    Raises
    ------
    maws.errors.ConfigurationError
        When called with a structure that was not built by AmberTools, and so
        has no parameters for OpenMM to read.
    """

    def make(system: BuiltSystem) -> EnergyModel:
        return system.energy_model(platform=platform)

    return make


# ---------------------------------------------------------------------------
# What the search reports
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Candidate:
    """One way of extending the strand, and how well it scored.

    Parameters
    ----------
    sequence : maws.values.NucleotideSequence
        The strand this candidate would give.
    token : str
        The nucleotide that was added, as written.
    direction : {"3prime", "5prime"}
        Which end it was added at.
    entropy : float
        The score. Lower means more promising. See
        :func:`maws.scoring.entropy_score`.
    energy : float
        The lowest energy found among the shapes tried, in kJ/mol.
    system : maws.topology.BuiltSystem
        The structure this candidate was scored on.
    pose : maws.pose.Pose
        The positions that gave `energy`.
    """

    sequence: NucleotideSequence
    token: str
    direction: Direction
    entropy: float
    energy: float
    system: BuiltSystem
    pose: Pose


class StepEvent:
    """Base class for everything :func:`grow_aptamer` reports.

    Tell the four kinds apart with ``match`` or ``isinstance``.

    See Also
    --------
    StepStarted, CandidateScored, StepCompleted, SearchFinished
    """


@dataclass(frozen=True, slots=True)
class StepStarted(StepEvent):
    """A growth step is beginning.

    Parameters
    ----------
    step : int
        Which step, counting from 1.
    total : int
        How many steps the search will run in all.
    sequence : maws.values.NucleotideSequence
        The strand as it stands before this step. Empty on step 1.
    """

    step: int
    total: int
    sequence: NucleotideSequence


@dataclass(frozen=True, slots=True)
class CandidateScored(StepEvent):
    """One way of extending the strand has been tried and scored.

    Parameters
    ----------
    step : int
        Which step this candidate belongs to.
    candidate : Candidate
        The candidate and its score.
    """

    step: int
    candidate: Candidate


@dataclass(frozen=True, slots=True)
class StepCompleted(StepEvent):
    """A growth step has finished and its winner has been chosen.

    Parameters
    ----------
    step : int
        Which step finished.
    winner : Candidate
        The candidate with the lowest score, which the next step grows from.
    """

    step: int
    winner: Candidate


@dataclass(frozen=True, slots=True)
class SearchFinished(StepEvent):
    """The search has run to the requested length.

    Parameters
    ----------
    winner : Candidate
        The final strand, its structure and its positions.
    steps : int
        How many steps were run.
    """

    winner: Candidate
    steps: int


# ---------------------------------------------------------------------------
# The algorithm
# ---------------------------------------------------------------------------


def grow_aptamer(
    base: BuiltSystem,
    *,
    energy: EnergyFactory,
    builder: Builder,
    sampler: Sampler | None = None,
    n_nucleotides: int = 15,
    role: str = "aptamer",
    target_role: str = "ligand",
    alphabet: str | None = None,
    first_samples: int = 5000,
    samples: int = 5000,
    beta: float = 0.01,
    scorer: Scorer = entropy_score,
    n_torsions: int = DEFAULT_TORSIONS,
    relax_iterations: int = 0,
    relax_size: float = 0.5,
    rng: np.random.Generator | None = None,
) -> Iterator[StepEvent]:
    """grow_aptamer(base, *, energy, builder, n_nucleotides=15, ...)

    Grow an aptamer against a target, reporting progress as it goes.

    Parameters
    ----------
    base : maws.topology.BuiltSystem
        A built structure holding the target and an empty strand. Build one
        from an assembly whose aptamer chain has no sequence.
    energy : EnergyFactory
        Builds the scorer for each candidate structure. Use
        :func:`openmm_energy` for a real run and :func:`stub_energy` for a
        test.
    builder : maws.build.Builder
        What rebuilds the structure each time the strand grows.
    sampler : maws.sampling.Sampler, optional
        Proposes places to put the strand next to the target. Only used on step
        1, since later steps keep the position already found and vary only the
        shape. Defaults to a :class:`~maws.sampling.SurfaceSampler` around the
        target.
    n_nucleotides : int, default=15
        How long the finished strand should be.
    role : str, default="aptamer"
        Which chain of `base` to grow.
    target_role : str, default="ligand"
        Which chain is the target. Used to build the default sampler.
    alphabet : str, optional
        Which nucleotides may be added. Defaults to the four the force field
        describes: ``"GAUC"`` for RNA, ``"GATC"`` for DNA.
    first_samples : int, default=5000
        How many shapes to try per candidate on step 1. Each also puts the
        strand somewhere new, so step 1 explores more than later steps and
        usually wants a larger number.
    samples : int, default=5000
        How many shapes to try per candidate on every later step. This is the
        main cost of a run: the total work is roughly
        ``8 * samples * n_nucleotides`` energy evaluations.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ.
    scorer : maws.scoring.Scorer, default=entropy_score
        Reduces a candidate's energies to one number, lower being better.
    n_torsions : int, default=4
        How many of the new residue's turnable bonds to vary.
    relax_iterations : int, default=0
        How many nudge-and-settle rounds to run after joining a new residue.
        Zero skips it, which is right for a stand-in scorer that does not
        settle anything. Around 50 suits a real run.
    relax_size : float, default=0.5
        How far each atom may be nudged during that settling, in ångström.
    rng : numpy.random.Generator, optional
        Source of randomness. Pass one built with a fixed seed to make a run
        repeatable. Defaults to a fresh generator, so runs differ.

    Yields
    ------
    StepEvent
        In order: a :class:`StepStarted`, then one :class:`CandidateScored` per
        way of extending the strand, then a :class:`StepCompleted` — repeated
        once per step — and finally a single :class:`SearchFinished`.

    Raises
    ------
    maws.errors.ConfigurationError
        If `n_nucleotides` is less than 1, or `base` already has residues in
        the chain named by `role`.

    See Also
    --------
    maws.api.design : Runs this and returns only the final answer.

    Notes
    -----
    On step 1 the strand is a single residue, and every one of its turnable
    bonds is varied along with its position and orientation. On later steps the
    strand is already placed, so only the newly added residue is varied: its
    own bonds, plus the one bond of its neighbour that swings it. Which end the
    residue was added at decides which side of each bond moves, so that the
    part of the strand already sitting against the target stays where it is.

    Examples
    --------
    Stopping early once a candidate is good enough, which the event stream
    makes possible without any support from the search itself:

    >>> for event in grow_aptamer(base, ...):  # doctest: +SKIP
    ...     if isinstance(event, StepCompleted) and event.winner.entropy < -0.5:
    ...         break
    """
    if n_nucleotides < 1:
        raise ConfigurationError(
            f"n_nucleotides must be at least 1, got {n_nucleotides}"
        )
    start_chain = base.chain(role)
    if start_chain.n_residues:
        raise ConfigurationError(
            f"chain {role!r} already holds {start_chain.n_residues} residues; "
            f"grow_aptamer starts from an empty strand"
        )

    generator = rng if rng is not None else np.random.default_rng()
    letters = alphabet if alphabet is not None else base.forcefield.alphabet
    angles = TorsionAngles(n_torsions, rng=generator)
    placements = (
        sampler
        if sampler is not None
        else SurfaceSampler.around(base, target_role, rng=generator)
    )

    yield StepStarted(step=1, total=n_nucleotides, sequence=NucleotideSequence(()))
    scored = []
    for token in letters:
        candidate = _seed_candidate(
            base,
            token=token,
            role=role,
            builder=builder,
            energy=energy,
            sampler=placements,
            angles=angles,
            n_samples=first_samples,
            beta=beta,
            scorer=scorer,
            n_torsions=n_torsions,
        )
        scored.append(candidate)
        yield CandidateScored(step=1, candidate=candidate)
    best = _best_of(scored)
    yield StepCompleted(step=1, winner=best)

    for step in range(2, n_nucleotides + 1):
        yield StepStarted(step=step, total=n_nucleotides, sequence=best.sequence)
        scored = []
        for token in letters:
            for direction in ("3prime", "5prime"):
                candidate = _grow_candidate(
                    best,
                    token=token,
                    direction=direction,
                    role=role,
                    builder=builder,
                    energy=energy,
                    angles=angles,
                    n_samples=samples,
                    beta=beta,
                    scorer=scorer,
                    n_torsions=n_torsions,
                    relax_iterations=relax_iterations,
                    relax_size=relax_size,
                    rng=generator,
                )
                scored.append(candidate)
                yield CandidateScored(step=step, candidate=candidate)
        best = _best_of(scored)
        yield StepCompleted(step=step, winner=best)

    yield SearchFinished(winner=best, steps=n_nucleotides)


def _best_of(candidates: Sequence[Candidate]) -> Candidate:
    """Return the candidate with the lowest score.

    Parameters
    ----------
    candidates : sequence of Candidate
        The candidates tried in one step. Must not be empty.

    Returns
    -------
    Candidate
        The most promising one. Ties go to whichever was tried first.

    Raises
    ------
    maws.errors.ConfigurationError
        If `candidates` is empty, which means the alphabet was empty.
    """
    if not candidates:
        raise ConfigurationError(
            "no candidates were scored this step; the alphabet is empty"
        )
    return min(candidates, key=lambda candidate: candidate.entropy)


def _seed_candidate(
    base: BuiltSystem,
    *,
    token: str,
    role: str,
    builder: Builder,
    energy: EnergyFactory,
    sampler: Sampler,
    angles: TorsionAngles,
    n_samples: int,
    beta: float,
    scorer: Scorer,
    n_torsions: int,
) -> Candidate:
    """Score a single nucleotide placed against the target.

    Parameters
    ----------
    base : maws.topology.BuiltSystem
        The structure holding the target and an empty strand.
    token : str
        The nucleotide to try, as written.
    role : str
        Which chain is the strand.
    builder : maws.build.Builder
        What builds the one-residue structure.
    energy : EnergyFactory
        Builds the scorer for it.
    sampler : maws.sampling.Sampler
        Proposes where to put the residue.
    angles : maws.sampling.TorsionAngles
        Supplies the random shapes.
    n_samples : int
        How many shapes to try.
    beta : float
        How sharply lower energies are favoured.
    scorer : maws.scoring.Scorer
        Reduces the energies to one number.
    n_torsions : int
        How many turnable bonds to vary.

    Returns
    -------
    Candidate
        The nucleotide, its score, and the best positions found.
    """
    grown = grow_chain(
        base,
        base.pose,
        role=role,
        token=token,
        direction="3prime",
        builder=builder,
    )
    model = energy(grown.system)
    chain = grown.system.chain(role)
    torsions = chain.residue(0).torsions(limit=n_torsions)

    energies: list[float] = []
    best_energy = float("inf")
    best_pose = grown.pose
    for placement in islice(sampler, n_samples):
        pose = grown.pose.place(chain, placement)
        pose = pose.rotate_all(torsions, angles.sample(len(torsions)))
        value = model.evaluate(pose)
        energies.append(value)
        if value < best_energy:
            best_energy, best_pose = value, pose

    return Candidate(
        sequence=chain.sequence,
        token=token,
        direction="3prime",
        entropy=scorer(energies, beta=beta),
        energy=best_energy,
        system=grown.system,
        pose=best_pose,
    )


def _grow_candidate(
    previous: Candidate,
    *,
    token: str,
    direction: Direction,
    role: str,
    builder: Builder,
    energy: EnergyFactory,
    angles: TorsionAngles,
    n_samples: int,
    beta: float,
    scorer: Scorer,
    n_torsions: int,
    relax_iterations: int,
    relax_size: float,
    rng: np.random.Generator,
) -> Candidate:
    """Score one way of extending the strand by a nucleotide.

    Parameters
    ----------
    previous : Candidate
        The winner of the previous step, whose structure and positions this
        grows from.
    token : str
        The nucleotide to add, as written.
    direction : {"3prime", "5prime"}
        Which end to add it at.
    role : str
        Which chain is the strand.
    builder : maws.build.Builder
        What rebuilds the structure.
    energy : EnergyFactory
        Builds the scorer for the rebuilt structure.
    angles : maws.sampling.TorsionAngles
        Supplies the random shapes.
    n_samples : int
        How many shapes to try.
    beta : float
        How sharply lower energies are favoured.
    scorer : maws.scoring.Scorer
        Reduces the energies to one number.
    n_torsions : int
        How many turnable bonds to vary.
    relax_iterations : int
        How many nudge-and-settle rounds to run before sampling.
    relax_size : float
        How far each atom may be nudged during that settling, in ångström.
    rng : numpy.random.Generator
        Source of randomness for the settling.

    Returns
    -------
    Candidate
        The extended strand, its score, and the best positions found.
    """
    grown = grow_chain(
        previous.system,
        previous.pose,
        role=role,
        token=token,
        direction=direction,
        builder=builder,
    )
    model = energy(grown.system)

    start = grown.pose
    if relax_iterations:
        start = perturb_and_minimize(
            start,
            model,
            size=relax_size,
            iterations=relax_iterations,
            rng=rng,
        ).pose

    chain = grown.system.chain(role)
    torsions = _growth_torsions(chain, direction, n_torsions)

    energies: list[float] = []
    best_energy = float("inf")
    best_pose = start
    for _ in range(n_samples):
        pose = start.rotate_all(torsions, angles.sample(len(torsions)))
        value = model.evaluate(pose)
        energies.append(value)
        if value < best_energy:
            best_energy, best_pose = value, pose

    return Candidate(
        sequence=chain.sequence,
        token=token,
        direction=direction,
        entropy=scorer(energies, beta=beta),
        energy=best_energy,
        system=grown.system,
        pose=best_pose,
    )


def _growth_torsions(
    chain: ChainView, direction: Direction, n_torsions: int
) -> tuple[Torsion, ...]:
    """Choose which bonds to vary after a residue has been added.

    Parameters
    ----------
    chain : maws.pose.ChainView
        The strand, after growing.
    direction : {"3prime", "5prime"}
        Which end the new residue was added at.
    n_torsions : int
        How many bonds to return.

    Returns
    -------
    tuple of maws.values.Torsion
        The bonds to vary, in global atom indices. At most `n_torsions` of
        them, and fewer only when the residue does not have that many.

    Notes
    -----
    The rest of the strand is already positioned against the target, so the
    bonds are chosen and oriented to leave it where it is and move only the new
    residue.

    Both ends are treated the same way: at most ``n_torsions - 1`` bonds from
    inside the new residue, which change its own shape, plus the one bond that
    joins it to the residue it was added next to, which swings it as a whole.
    That last bond belongs to the neighbour, not to the new residue, and it is
    read in the direction that leaves the rest of the strand where it is.

    Keeping the two ends matched matters because the search decides which end
    to grow at by comparing the scores the two produce. An end offered fewer
    bonds would try a smaller set of shapes, find a worse best one, and lose
    that comparison for a reason that has nothing to do with the molecule.
    """
    if direction == "3prime":
        newest = chain.residue(-1)
        torsions = list(newest.torsions("3prime", limit=n_torsions - 1))
        if chain.n_residues > 1 and n_torsions > 1:
            neighbour = chain.residue(-2)
            index = min(n_torsions - 1, neighbour.n_torsions - 1)
            torsions.append(neighbour.torsion(index, "3prime"))
        return tuple(torsions)

    newest = chain.residue(0)
    torsions = list(newest.torsions("5prime", limit=n_torsions - 1))
    if chain.n_residues > 1 and n_torsions > 1:
        # The neighbour's first bond, read towards the 5' end, moves everything
        # in front of it — which is exactly the new residue. It is the mirror
        # image of the neighbour's last bond used above for the 3' end.
        torsions.append(chain.residue(1).torsion(0, "5prime"))
    return tuple(torsions)

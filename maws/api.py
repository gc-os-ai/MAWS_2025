"""
maws.api
========

Designing an aptamer in one call.

An *aptamer* is a short strand of RNA or DNA that folds up against another
molecule — the *target* — and sticks to it. :func:`design` searches for one.
Give it a PDB file for the target, the standard text format for a molecular
structure listing one atom per line with its position, and it hands back the
strand it settled on. The residue libraries, the program that assembles the
structures, the sampler that proposes shapes and the score that ranks them are
all chosen for you.

:class:`MawsResult` is what comes back. It carries a sequence whether or not
the run reached the requested length, so a run that gave up part-way is told
apart by its ``success`` flag rather than by inspecting the numbers.

Every stage is available on its own when the defaults do not suit.
:func:`maws.search.grow_aptamer` runs the same search while reporting each step
as it happens, and :func:`collect` reduces such a stream back to one result.

Examples
--------
>>> result = design("target.pdb", length=15)  # doctest: +SKIP
>>> print(result)  # doctest: +SKIP
G A U C G A U C G A U C G A U  (15 nt, E=-1421.90 kJ/mol, S=-0.072000)
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
from maws.scoring import Scorer, free_energy_score
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
"""Which region proposed positions for the aptamer are drawn from.

``"sphere"`` fills a ball centred on the target; ``"surface-following"`` keeps
to a thin shell hugging the target's surface.
"""


@dataclass(frozen=True, slots=True)
class MawsResult:
    """The strand a design run settled on, and how it scored.

    Parameters
    ----------
    sequence : str
        The designed strand, written 5' to 3' — the two ends of such a strand
        are chemically different and are named that way by convention. One
        whitespace-separated token per nucleotide, e.g. ``"G A U C"``.
    energy : float
        Potential energy of the best arrangement found, in kJ/mol. Comparable
        between runs against the same target only; the number includes the
        target's own internal energy, which differs from target to target.
    score : float
        The number that chose this strand, lower being better. With the
        default scorer it is a free energy in kJ/mol, covering all the shapes
        tried rather than only the best one. See
        :func:`maws.scoring.free_energy_score`.
    steps : int
        How many nucleotides were added. Equal to the requested length after a
        run that finished.
    success : bool, default=True
        Whether the run finished as intended. False means `sequence` is
        whatever the run had reached when it stopped, which may be shorter than
        asked for.
    message : str, default=""
        What went wrong, when `success` is False. Empty otherwise.
    system : maws.topology.BuiltSystem, optional
        The final structure, atoms and force-field parameters together. Present
        after a real run, and what to pass on to write the result out or score
        it again. Left out by runs that never built anything.
    pose : maws.pose.Pose, optional
        Atom positions of the best arrangement found, in ångström.

    See Also
    --------
    design : Produces one of these.
    maws.AptamerDesigner : Keeps one of these per target, in ``results_``.

    Notes
    -----
    Comparing two instances ignores `system` and `pose`, so two runs that
    arrived at the same sequence and scores compare equal even though their
    structures are separate objects.

    `success` and `message` follow the same idea as
    :class:`scipy.optimize.OptimizeResult`: a run that found nothing usable
    should be distinguishable from one that worked, without inspecting the
    numbers.

    Examples
    --------
    >>> result = MawsResult("G A U", -1421.9, -1215.4, steps=3)
    >>> result.length
    3
    >>> print(result)
    G A U  (3 nt, E=-1421.90 kJ/mol, score=-1215.40)

    A run that stopped early says so, and its sequence is short:

    >>> stopped = MawsResult(
    ...     "G A",
    ...     -900.0,
    ...     -780.0,
    ...     steps=2,
    ...     success=False,
    ...     message="the search stopped before it finished",
    ... )
    >>> stopped.success, stopped.length
    (False, 2)
    """

    sequence: str
    energy: float
    score: float
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
            f"E={self.energy:.2f} kJ/mol, score={self.score:.2f})"
        )


def design(
    target: str | Path,
    *,
    length: int = 15,
    aptamer: AptamerType = "RNA",
    molecule: MoleculeType = "protein",
    net_charge: int = 0,
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
    scorer: Scorer = free_energy_score,
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

    .. warning::
        A run at the default settings performs on the order of a million
        energy evaluations and takes hours. Try `samples` of 50 and `length`
        of 3 first, to confirm the target file and the installation are sound.

    Parameters
    ----------
    target : str or pathlib.Path
        PDB file for the molecule the aptamer should bind to. It is treated as
        a single unit: it is never grown and none of its bonds are turned. Its
        size sets the region the strand is placed in.
    length : int, default=15
        How many nucleotides the finished strand should have. Each one costs
        another full round of sampling, so the run time grows with it.
    aptamer : {"RNA", "DNA"}, default="RNA"
        Which nucleic acid to build the strand from. RNA strands are written
        with the letters ``G A U C`` and DNA with ``G A T C``, and the two are
        built from different residue libraries.
    molecule : {"protein", "organic", "lipid"}, default="protein"
        What the target is. This picks the force field used for it — the table
        of numbers that turns a set of atom positions into an energy. Proteins
        and lipids have ready-made parameters; an organic molecule does not, so
        its own are worked out at the start of the run, which takes extra
        minutes.
    net_charge : int, default=0
        The target's overall charge, in units of the electron charge: -2 for a
        doubly deprotonated acid, +1 for a protonated amine, and so on. Only
        used when `molecule` is ``"organic"``, where the parameters are worked
        out at the start of the run: the total is shared out across the atoms,
        so the wrong total puts the wrong charge on every one of them and
        every energy computed afterwards is wrong with it. Zero is right for a
        neutral molecule and wrong for everything else.
    samples : int, default=5000
        How many shapes to try per candidate at each growth step. This is the
        main cost: the total work is roughly ``8 * samples * length`` energy
        evaluations. Raising it searches each candidate more thoroughly at
        proportional cost. Lower it to perhaps 50 for a quick check.
    first_samples : int, optional
        How many shapes to try on the first step, which also searches over
        where next to the target to put the strand. Raising it beyond `samples`
        buys a better starting position, which every later step builds on.
        Defaults to `samples`.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ. Raising it makes
        the choice between candidates depend mostly on their few best shapes;
        at zero every shape weighs the same and no candidate can win.
    salt_conc : float, default=0.15
        Concentration of dissolved salt in mol/L. Dissolved ions gather around
        charged atoms and damp the pull between them, so raising this weakens
        every electrostatic interaction in the run. The default is roughly the
        saltiness of blood.
    reach : float, default=10.0
        How far past the target's furthest atom the sampling region extends, in
        ångström. Larger values consider positions further from the target, at
        the cost of more proposals landing in empty space. Only used when
        `sampling` is ``"sphere"``.
    probe : float, default=1.4
        Radius in ångström of the ball rolled over the target to find its
        surface. 1.4 Å is the size of a water molecule, the usual choice.
        Larger values smooth over narrow clefts, so the strand is not offered
        positions inside them.
    d_max : float, default=6.0
        Thickness in ångström of the shell positions are drawn from. Thinner
        shells hold the strand closer to the surface. Only used when `sampling`
        is ``"surface-following"``.
    sampling : {"sphere", "surface-following"}, default="sphere"
        Which region to draw positions from. ``"sphere"`` fills a ball around
        the target; ``"surface-following"`` keeps to a shell hugging its
        surface, which concentrates the search where contact is possible but
        rejects far more proposals.
    relax_iterations : int, default=50
        How many rounds of nudging the atoms and letting them settle to run
        after each nucleotide is joined on. Joining leaves the new residue
        strained against its neighbour, and each round displaces the atoms a
        little and minimises the energy again to work that out. Zero skips the
        whole thing and is much faster, at the cost of scoring strained shapes.
    seed : int, optional
        Fixes the randomness, so the same call gives the same answer. Left out,
        every run differs.
    scorer : maws.scoring.Scorer, default=free_energy_score
        Reduces one candidate's energies to the single number the candidates
        are ranked by, lower being better. Replace it to rank by something
        other than how tightly the energies cluster.
    builder : maws.build.Builder, optional
        What turns a sequence into a structure with atom positions. Defaults to
        :class:`~maws.build.LeapBuilder`, which drives ``tleap``, the
        structure-assembly program of the molecular modelling suite AmberTools,
        and so needs AmberTools installed. Pass
        :class:`~maws.build.FakeBuilder` to exercise the search on stand-in
        geometry without it.
    clean_pdb : bool, default=False
        Whether to tidy the target file before building: a downloaded
        structure often carries several copies of the molecule, water, and
        other records that would stop the structure being assembled. Only
        applied when `molecule` is ``"protein"``.
    keep_chains : str, default="all"
        Which chains of the target to keep when cleaning — a chain being one
        continuous molecule within the file, labelled with a letter. Either
        ``"all"``, ``"one"`` for the first, or a comma-separated list such as
        ``"A,B"``. Keeping fewer makes the run cheaper and changes what the
        strand can reach.
    remove_h : bool, default=False
        Whether to strip hydrogen atoms while cleaning. Useful when the file's
        hydrogens are named in a way the force field does not recognise, since
        fresh ones are added when the structure is assembled.
    drop_hetatm : bool, default=False
        Whether to drop records for anything that is not part of the molecule
        proper, such as bound water or ions. Keep them to design against a
        target that needs them, at the risk of one having no parameters.
    on_event : callable, optional
        Called with every :class:`~maws.search.StepEvent` as the search
        produces it. Use it to print progress or write files as the run goes,
        rather than waiting for the result.
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
        If AmberTools is needed and is not installed, or one of its programs
        fails.
    maws.errors.BuildError
        If the structure-assembly program runs but produces no structure, which
        usually means the target file holds something it cannot describe.

    See Also
    --------
    maws.AptamerDesigner : The same run configured as an object.
    maws.search.grow_aptamer : The same search, reporting step by step.
    MawsResult : What this returns.

    Notes
    -----
    Energies are in kJ/mol, lengths in ångström, and concentrations in mol/L.

    No result file is written. The structures built along the way are cached on
    disk by the builder, and `clean_pdb` leaves a tidied copy of the target
    next to the original, but everything the run found comes back in the return
    value or through `on_event`.

    Examples
    --------
    A short run, to check that a target file can be read and built:

    >>> result = design("target.pdb", length=4, samples=50, seed=0)  # doctest: +SKIP
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
    assembly = (
        Assembly()
        .with_aptamer(library)
        .with_ligand(path, forcefield, net_charge=net_charge)
    )

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
    """collect(events, *, on_event=None) -> MawsResult

    Run a search to the end and return only its answer.

    Parameters
    ----------
    events : iterable of maws.search.StepEvent
        The stream :func:`maws.search.grow_aptamer` produces. Reading it is
        what makes the search run, so this call is where the time goes.
    on_event : callable, optional
        Called with each event as it arrives, before it is discarded. The only
        chance to see the intermediate steps, since nothing but the last one
        is kept.

    Returns
    -------
    MawsResult
        Built from the final :class:`~maws.search.SearchFinished` event. If the
        stream ends without one, the result has ``success=False`` and a
        `message` saying so.

    See Also
    --------
    design : Sets up a search and calls this on it.
    maws.search.grow_aptamer : Produces the stream this consumes.

    Examples
    --------
    A whole search against stand-in geometry, which needs nothing installed:

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
            score=float("nan"),
            steps=0,
            success=False,
            message="the search stopped before it finished",
        )

    winner = finished.winner
    return MawsResult(
        sequence=str(winner.sequence),
        energy=winner.energy,
        score=winner.score,
        steps=finished.steps,
        system=winner.system,
        pose=winner.pose,
    )

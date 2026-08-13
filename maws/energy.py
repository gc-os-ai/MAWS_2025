"""
maws.energy
===========

Turning a set of atom positions into a single number.

MAWS chooses between candidate shapes by scoring each one, and the score starts
with the *potential energy* of the arrangement: a number saying how comfortable
the atoms are where they are. Atoms that are too close, or charges that repel,
push it up. Favourable contacts pull it down. Lower means a more settled
structure, and the units are kilojoules per mole throughout.

Anything with an :meth:`~EnergyModel.evaluate` and a
:meth:`~EnergyModel.minimize` method can act as the scorer. Two are provided:

:class:`OpenMMEnergy`
    The real one. It hands the positions to OpenMM, a molecular simulation
    package, and reads back the energy of the arrangement. This needs a
    structure that has actually been built by the AmberTools programs, so it
    only works when those are installed.

:class:`StubEnergy`
    A stand-in that computes a simple repulsion between two chains with a few
    lines of arithmetic. It has no chemistry in it, but it is a real function
    of the positions with a real minimum, so a search can be run against it and
    checked. It needs nothing installed and takes microseconds, which is what
    lets the whole search algorithm be tested.

Having both is the point of writing the scorer as an interface rather than as a
method on a molecule class: the choice of scorer is what separates "run this
design for real" from "check the search logic works".

Examples
--------
>>> import numpy as np
>>> from maws.pose import Pose
>>> from maws.values import AtomRange
>>> pose = Pose(np.array([[0.0, 0, 0], [5.0, 0, 0]]))
>>> model = StubEnergy(AtomRange(0, 1), AtomRange(1, 2))
>>> round(model.evaluate(pose), 4)
0.04
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple, Protocol, runtime_checkable

import numpy as np

from maws.errors import ConfigurationError
from maws.pose import Pose
from maws.values import AtomRange

if TYPE_CHECKING:  # pragma: no cover - needed by type checkers only
    from openmm.app import Topology

    from maws.forcefield import ForceField

__all__ = [
    "CompositeEnergy",
    "EnergyModel",
    "OpenMMEnergy",
    "Relaxed",
    "StubEnergy",
]


class Relaxed(NamedTuple):
    """The result of settling a structure into a nearby comfortable shape.

    Parameters
    ----------
    energy : float
        Potential energy of `pose`, in kilojoules per mole.
    pose : maws.pose.Pose
        The settled positions.

    Notes
    -----
    A named tuple, so it can be read either way round::

        energy, pose = model.minimize(pose)
        result = model.minimize(pose)
        result.energy
    """

    energy: float
    pose: Pose


@runtime_checkable
class EnergyModel(Protocol):
    """Anything that can score a set of atom positions.

    A class satisfies this by having the two methods below. There is no base
    class to inherit from and nothing to register.

    See Also
    --------
    StubEnergy : A fast stand-in with no chemistry, for tests.
    OpenMMEnergy : The real scorer.
    """

    def evaluate(self, pose: Pose) -> float:
        """Return the potential energy of an arrangement of atoms.

        Parameters
        ----------
        pose : maws.pose.Pose
            The positions to score.

        Returns
        -------
        float
            Potential energy in kilojoules per mole. Lower is better.
        """
        ...

    def minimize(self, pose: Pose, *, max_iterations: int = 100) -> Relaxed:
        """Settle an arrangement into the nearest comfortable shape.

        Parameters
        ----------
        pose : maws.pose.Pose
            The positions to start from.
        max_iterations : int, default=100
            How many adjustment steps to allow before stopping. More steps
            settle further but cost more time.

        Returns
        -------
        Relaxed
            The settled positions and their energy.
        """
        ...


class StubEnergy:
    """A scorer with no chemistry in it, for testing the search.

    The energy is the total inverse-square repulsion between the atoms of two
    named groups: every pair of atoms, one from each group, contributes
    ``scale / distance²``. So the two groups are pushed apart, arrangements
    that separate them score lower, and the score changes smoothly as atoms
    move.

    That is enough structure for a search to have something real to optimise,
    while costing microseconds and needing nothing installed.

    Parameters
    ----------
    first, second : maws.values.AtomRange
        The two groups of atoms whose separation is being scored. Usually the
        aptamer's atoms and the target's.
    scale : float, default=1.0
        Multiplies the whole score. Changes how strongly the two groups repel,
        and therefore how sharply the search prefers separated arrangements.

    See Also
    --------
    OpenMMEnergy : The real scorer, for actual designs.

    Examples
    --------
    Two single atoms 5 Å apart, then 10 Å apart:

    >>> import numpy as np
    >>> from maws.pose import Pose
    >>> from maws.values import AtomRange
    >>> model = StubEnergy(AtomRange(0, 1), AtomRange(1, 2))
    >>> round(model.evaluate(Pose(np.array([[0.0, 0, 0], [5.0, 0, 0]]))), 4)
    0.04
    >>> round(model.evaluate(Pose(np.array([[0.0, 0, 0], [10.0, 0, 0]]))), 4)
    0.01
    """

    __slots__ = ("_first", "_scale", "_second")

    def __init__(
        self, first: AtomRange, second: AtomRange, *, scale: float = 1.0
    ) -> None:
        self._first = first
        self._second = second
        self._scale = float(scale)

    def __repr__(self) -> str:
        return f"<StubEnergy {self._first} vs {self._second} scale={self._scale}>"

    def evaluate(self, pose: Pose) -> float:
        """Return the total repulsion between the two groups of atoms.

        Parameters
        ----------
        pose : maws.pose.Pose
            The positions to score.

        Returns
        -------
        float
            Sum of ``scale / distance²`` over every pair of atoms with one
            atom from each group. Squared distances are floored at a tiny value
            so that two atoms in exactly the same place give a large number
            rather than a division by zero.
        """
        first = pose.atoms(self._first)
        second = pose.atoms(self._second)
        gaps = first[:, None, :] - second[None, :, :]
        squared = (gaps**2).sum(axis=-1)
        return float(self._scale * np.reciprocal(np.maximum(squared, 1e-6)).sum())

    def minimize(self, pose: Pose, *, max_iterations: int = 100) -> Relaxed:
        """Return `pose` unchanged, with its score.

        Parameters
        ----------
        pose : maws.pose.Pose
            The positions to score.
        max_iterations : int, default=100
            Ignored. Accepted so this can stand in for a real scorer.

        Returns
        -------
        Relaxed
            The same positions and their energy.

        Notes
        -----
        Settling a structure is deliberately left as a no-op. A test that
        depends on this doing something is really a test of the real scorer,
        and should be marked as such.
        """
        return Relaxed(self.evaluate(pose), pose)


class OpenMMEnergy:
    """A scorer backed by OpenMM, the molecular simulation package.

    This is the only class in MAWS that builds an OpenMM simulation or holds
    one open, which is what keeps the rest of the package free of simulation
    state. Everything else works on plain arrays of numbers.

    The energy is computed with an *implicit solvent* model: water is treated
    as a smooth continuum surrounding the molecule rather than as individual
    molecules. Dissolved salt in that continuum damps the attraction and
    repulsion between charged atoms, which matters a great deal here because a
    nucleic acid backbone carries a charge on every residue.

    Parameters
    ----------
    prmtop : openmm.app.AmberPrmtopFile
        A parsed AMBER parameter file, describing which atoms exist, how they
        are bonded, and what parameters apply to them.
    forcefield : maws.forcefield.ForceField
        The run's force field settings. Only the salt concentration is used
        here; the force field names were already applied when the structure was
        built.
    platform : str, optional
        Which OpenMM compute backend to use, e.g. ``"CPU"`` or ``"CUDA"``. By
        default the fastest available is chosen.
    frozen : maws.values.AtomRange, optional
        Atoms that must not move when :meth:`minimize` settles the structure.
        Normally the target: the search compares candidates by how well each
        sits against the same target, and a target that reshapes itself
        differently for each one is not the same target. Nothing is frozen by
        default.

    See Also
    --------
    from_prmtop : Build one directly from a file path.
    StubEnergy : The stand-in used when AmberTools is not available.

    Notes
    -----
    Creating one of these builds an OpenMM simulation, which takes noticeably
    longer than a single evaluation. Build it once and reuse it.

    Freezing is done by setting an atom's mass to zero, which is how OpenMM is
    told an atom is held in place. Frozen atoms still contribute to the energy,
    exactly as they did before; they simply do not move. Their contribution is
    the same number for every candidate, so it shifts every energy by a
    constant and changes no comparison.
    """

    __slots__ = ("_simulation", "_system", "_topology")

    def __init__(
        self,
        prmtop: object,
        forcefield: ForceField,
        *,
        platform: str | None = None,
        frozen: AtomRange | None = None,
    ) -> None:
        import openmm as mm
        from openmm import app, unit

        self._system = prmtop.createSystem(  # type: ignore[attr-defined]
            nonbondedMethod=app.NoCutoff,
            nonbondedCutoff=5 * unit.angstrom,
            constraints=None,
            implicitSolvent=app.OBC1,
            implicitSolventSaltConc=forcefield.salt_conc * unit.molar,
        )
        if frozen is not None:
            n_atoms = self._system.getNumParticles()
            if frozen.stop > n_atoms:
                raise ConfigurationError(
                    f"cannot freeze atoms {frozen.start}-{frozen.stop}; the "
                    f"structure has only {n_atoms}"
                )
            for index in range(frozen.start, frozen.stop):
                self._system.setParticleMass(index, 0.0)
        self._topology: Topology = prmtop.topology  # type: ignore[attr-defined]
        integrator = mm.LangevinIntegrator(
            300.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
        )
        self._simulation = _make_simulation(
            self._topology, self._system, integrator, platform
        )

    def __repr__(self) -> str:
        return f"<OpenMMEnergy {self._system.getNumParticles()} atoms>"

    @classmethod
    def from_prmtop(
        cls,
        path: str | Path,
        forcefield: ForceField,
        *,
        platform: str | None = None,
        frozen: AtomRange | None = None,
    ) -> OpenMMEnergy:
        """Build a scorer by reading an AMBER parameter file from disk.

        Parameters
        ----------
        path : str or pathlib.Path
            Path to a ``.prmtop`` file, as written when the structure was
            built.
        forcefield : maws.forcefield.ForceField
            The run's force field settings.
        platform : str, optional
            Which OpenMM compute backend to use.
        frozen : maws.values.AtomRange, optional
            Atoms to hold in place when settling. See the class docstring.

        Returns
        -------
        OpenMMEnergy
            A scorer for that structure.
        """
        from openmm import app

        return cls(
            app.AmberPrmtopFile(str(path)),
            forcefield,
            platform=platform,
            frozen=frozen,
        )

    @property
    def topology(self) -> Topology:
        """openmm.app.Topology : Which atoms exist and how they are bonded."""
        return self._topology

    def evaluate(self, pose: Pose) -> float:
        """Return the potential energy of an arrangement of atoms.

        Parameters
        ----------
        pose : maws.pose.Pose
            The positions to score. Must hold one position per atom of the
            structure this scorer was built for.

        Returns
        -------
        float
            Potential energy in kilojoules per mole.
        """
        from openmm import unit

        self._simulation.context.setPositions(pose.to_openmm())
        state = self._simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    def minimize(self, pose: Pose, *, max_iterations: int = 100) -> Relaxed:
        """Settle an arrangement into the nearest comfortable shape.

        Parameters
        ----------
        pose : maws.pose.Pose
            The positions to start from.
        max_iterations : int, default=100
            How many adjustment steps OpenMM may take. More steps settle
            further but cost proportionally more time.

        Returns
        -------
        Relaxed
            The settled positions and their energy. Any atoms this scorer was
            built to freeze come back exactly where they were.
        """
        from openmm import unit

        self._simulation.context.setPositions(pose.to_openmm())
        self._simulation.minimizeEnergy(maxIterations=max_iterations)
        state = self._simulation.context.getState(getPositions=True, getEnergy=True)
        return Relaxed(
            energy=state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),
            pose=Pose.from_openmm(state.getPositions(), pose.system),
        )


class CompositeEnergy:
    """A scorer that adds up several other scorers, each with a weight.

    Useful when the potential energy is not the only thing worth optimising —
    for instance if a measure of how well the strand folds on its own is to be
    taken into account alongside how well it sits against the target.

    Parameters
    ----------
    *parts : tuple of (EnergyModel, float)
        Pairs of a scorer and the weight to multiply its result by.

    Raises
    ------
    maws.errors.ConfigurationError
        If no parts are given, which would make every arrangement score zero.

    Examples
    --------
    >>> import numpy as np
    >>> from maws.pose import Pose
    >>> from maws.values import AtomRange
    >>> plain = StubEnergy(AtomRange(0, 1), AtomRange(1, 2))
    >>> doubled = CompositeEnergy((plain, 2.0))
    >>> pose = Pose(np.array([[0.0, 0, 0], [5.0, 0, 0]]))
    >>> round(doubled.evaluate(pose) / plain.evaluate(pose), 6)
    2.0
    """

    __slots__ = ("_parts",)

    def __init__(self, *parts: tuple[EnergyModel, float]) -> None:
        if not parts:
            raise ConfigurationError(
                "a composite scorer needs at least one part, otherwise every "
                "arrangement scores zero"
            )
        self._parts = parts

    def __repr__(self) -> str:
        return f"<CompositeEnergy {len(self._parts)} parts>"

    def evaluate(self, pose: Pose) -> float:
        """Return the weighted sum of every part's score.

        Parameters
        ----------
        pose : maws.pose.Pose
            The positions to score.

        Returns
        -------
        float
            ``sum(weight * part.evaluate(pose))``.
        """
        return sum(weight * part.evaluate(pose) for part, weight in self._parts)

    def minimize(self, pose: Pose, *, max_iterations: int = 100) -> Relaxed:
        """Settle the arrangement using the first part, then score it fully.

        Parameters
        ----------
        pose : maws.pose.Pose
            The positions to start from.
        max_iterations : int, default=100
            How many adjustment steps the first part may take.

        Returns
        -------
        Relaxed
            The settled positions, scored by every part.

        Notes
        -----
        Only the first part does the settling, because settling means following
        the forces the first part defines. The returned energy is still the
        weighted sum over all parts, so it can be compared with
        :meth:`evaluate`.
        """
        first, _ = self._parts[0]
        settled = first.minimize(pose, max_iterations=max_iterations).pose
        return Relaxed(self.evaluate(settled), settled)


def _make_simulation(
    topology: Topology,
    system: object,
    integrator: object,
    platform: str | None = None,
) -> object:
    """Create an OpenMM simulation on the fastest backend available.

    Parameters
    ----------
    topology : openmm.app.Topology
        Which atoms exist and how they are bonded.
    system : openmm.System
        The forces acting on those atoms.
    integrator : openmm.Integrator
        How positions are advanced. Only used when settling a structure.
    platform : str, optional
        Force a particular backend by name. Overrides everything below.

    Returns
    -------
    openmm.app.Simulation
        A ready simulation.

    Notes
    -----
    Backends are tried in this order, and the first that works is used:

    1. the `platform` argument, if given;
    2. the ``MAWS_OPENMM_PLATFORM`` environment variable, if set;
    3. ``CUDA``, then ``OpenCL`` — graphics cards, much faster when present;
    4. whatever OpenMM picks by default, normally the CPU.

    A backend that is not installed raises when asked for, so each attempt is
    guarded and simply moves on to the next.
    """
    import openmm as mm
    from openmm import app

    mm.Platform.loadPluginsFromDirectory(mm.Platform.getDefaultPluginsDirectory())

    wanted = platform or os.getenv("MAWS_OPENMM_PLATFORM")
    if wanted:
        try:
            chosen = mm.Platform.getPlatformByName(wanted)
            return app.Simulation(topology, system, integrator, chosen)
        except Exception as exc:  # noqa: BLE001 - fall back to auto-detection
            warnings.warn(
                f"OpenMM backend {wanted!r} is unavailable ({exc}); "
                f"choosing one automatically instead",
                RuntimeWarning,
                stacklevel=2,
            )

    for name in ("CUDA", "OpenCL"):
        try:
            chosen = mm.Platform.getPlatformByName(name)
            return app.Simulation(topology, system, integrator, chosen)
        except Exception:  # noqa: BLE001,S112 - this backend is simply absent
            continue

    return app.Simulation(topology, system, integrator)

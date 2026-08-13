"""
maws
====

Designing aptamers against a target molecule, by simulation.

An *aptamer* is a short strand of RNA or DNA — a few dozen *nucleotides*, the
units such strands are built from, each written as a single letter — that folds
up against some other molecule and sticks to it. That other molecule is the
*target*. Aptamers are used much as antibodies are, to recognise one molecule in
a mixture, and the hard part of making one is deciding which sequence of
nucleotides to synthesise.

MAWS decides that on a computer. It starts from a PDB file for the target: the
standard text format for a molecular structure, listing one atom per line with
its position. From an empty strand it then repeats a single step. Every way of
adding one more nucleotide is tried — each of the four letters, at either end
of the strand, so eight candidates. Each candidate is given thousands of random
shapes, each shape is scored by the potential energy of the whole arrangement,
and the candidate whose energies cluster most tightly is kept. Nothing is
synthesised and nothing is measured: the answer is a sequence to go and test at
the bench.

Three ways in, in increasing order of control:

:func:`design`
    One call. Give it a target file and a length, get a sequence back.
:class:`AptamerDesigner`
    The same run as an object whose settings can be inspected, copied and
    changed one at a time, following scikit-learn's conventions.
:func:`grow_aptamer`
    The search itself, reporting each step as it happens so it can be watched,
    logged, or stopped part-way.

All three are assembled from parts that are usable on their own:
:class:`Assembly` describes what to build, :func:`build` builds it,
:class:`Pose` holds atom positions, and :class:`EnergyModel` scores them.

.. warning::
    A run at the default settings performs on the order of a million energy
    evaluations and takes hours. Try a short strand and a few dozen samples
    first, to confirm the target file and the installation are sound.

Examples
--------
>>> from maws import design
>>> result = design("target.pdb", length=15, seed=0)  # doctest: +SKIP
>>> result.sequence  # doctest: +SKIP
'G A U C G A U C G A U C G A U'

The same search run against stand-ins, which needs no molecular modelling
software installed and finishes in a moment:

>>> import numpy as np
>>> from maws import Assembly, FakeBuilder, ForceField, grow_aptamer, rna
>>> from maws.search import stub_energy
>>> builder = FakeBuilder()
>>> base = builder.build(
...     Assembly().with_aptamer(rna()).with_ligand_stub(20),
...     ForceField.for_target("RNA", "protein"),
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
>>> events[-1].steps
2
"""

from maws._designer import AptamerDesigner
from maws.api import MawsResult, design
from maws.build import Builder, FakeBuilder, LeapBuilder, build
from maws.energy import EnergyModel, OpenMMEnergy, Relaxed, StubEnergy
from maws.errors import (
    BuildError,
    ConfigurationError,
    MawsError,
    MissingLibrary,
    SamplingError,
    ToolchainError,
)
from maws.forcefield import ForceField
from maws.libraries import dna, rna
from maws.pose import ChainView, Pose, ResidueView
from maws.sampling import Placement, Sampler, SurfaceSampler, TorsionAngles
from maws.scoring import entropy_score
from maws.search import (
    Candidate,
    CandidateScored,
    SearchFinished,
    StepCompleted,
    StepEvent,
    StepStarted,
    grow_aptamer,
)
from maws.topology import Assembly, BuiltSystem
from maws.values import (
    AtomRange,
    NucleotideSequence,
    ResidueLibrary,
    ResidueTemplate,
    Torsion,
)

__version__ = "0.2.0.dev0"
"""Version of this package.

MAWS follows the usual scientific-Python deprecation policy: anything that is
going away warns for two minor releases first, and the warning names both the
replacement and the release it will be removed in.
"""

__all__ = [
    # one-call and object interfaces
    "AptamerDesigner",
    "MawsResult",
    "design",
    # the search, and what it reports
    "Candidate",
    "CandidateScored",
    "SearchFinished",
    "StepCompleted",
    "StepEvent",
    "StepStarted",
    "grow_aptamer",
    # describing and building a design
    "Assembly",
    "Builder",
    "BuiltSystem",
    "FakeBuilder",
    "ForceField",
    "LeapBuilder",
    "build",
    "dna",
    "rna",
    # positions and the windows onto them
    "ChainView",
    "Pose",
    "ResidueView",
    # scoring
    "EnergyModel",
    "OpenMMEnergy",
    "Relaxed",
    "StubEnergy",
    "entropy_score",
    # proposing positions and shapes
    "Placement",
    "Sampler",
    "SurfaceSampler",
    "TorsionAngles",
    # value types
    "AtomRange",
    "NucleotideSequence",
    "ResidueLibrary",
    "ResidueTemplate",
    "Torsion",
    # errors
    "BuildError",
    "ConfigurationError",
    "MawsError",
    "MissingLibrary",
    "SamplingError",
    "ToolchainError",
    "__version__",
]

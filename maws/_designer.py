"""
maws._designer
==============

A design run as a configurable object.

An *aptamer* is a short strand of RNA or DNA that folds up against a target
molecule and sticks to it; searching for one is what this package does.
:class:`AptamerDesigner` holds the settings for such a search, runs it in
:meth:`~AptamerDesigner.fit`, and leaves the answers on itself as attributes
whose names end in an underscore — ``sequences_``, ``energies_``. That trailing
underscore is a promise: the attribute does not exist until ``fit`` has run,
and asking for it before then raises :exc:`AttributeError`.

The estimator convention
------------------------
This layout is the one scikit-learn calls an *estimator*, and it is followed
here even though scikit-learn is not a dependency. Two rules define it:

- Every setting can be passed to ``__init__`` by name, and ``__init__`` stores
  each one unchanged under that same name and does nothing else. It validates
  nothing, reads no file, computes nothing.
- All the work happens in ``fit``, and the settings are checked there, at the
  moment they are used.

The first rule looks like a restriction and is really what makes the object
useful. Because each setting is still sitting there under the name it was
passed in by, the object can be taken apart into a plain dictionary of settings
with :meth:`~AptamerDesigner.get_params` and an equivalent one rebuilt from
that dictionary. Copying a configured run, changing one setting and trying
again is then two lines, and any tool that searches over settings — including
scikit-learn's own, where it is installed — can drive this class without
knowing anything about it. A constructor that quietly transformed or
rejected its arguments would break that round trip. The pleasant side effect is
that building a designer is instant and cannot fail, so a mistake in the
settings surfaces at ``fit``, next to the work it affects.

MAWS also fits the shape scikit-learn calls *transductive*: every target needs
its own search, and nothing is left over afterwards that could be applied to a
target never searched against. So there is :meth:`~AptamerDesigner.fit` and
:meth:`~AptamerDesigner.fit_predict`, and no ``predict``.

Examples
--------
>>> designer = AptamerDesigner(n_nucleotides=5, n_samples=50, random_state=0)
>>> designer.get_params()["n_nucleotides"]
5

Running it needs a real target file and the modelling software installed:

>>> designer.fit(["target.pdb"])  # doctest: +SKIP
>>> designer.sequences_[0]  # doctest: +SKIP
'G A U C G'
"""

from __future__ import annotations

import inspect
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from maws.api import MawsResult, SamplingMode, design
from maws.errors import ConfigurationError
from maws.forcefield import AptamerType, MoleculeType

__all__ = ["AptamerDesigner"]


class _BaseEstimator:
    """The part of the estimator convention every MAWS designer shares.

    Provides :meth:`get_params`, :meth:`set_params` and a readable ``repr``,
    all worked out by reading the arguments of the subclass's ``__init__``.
    Nothing has to be listed twice: adding a setting to ``__init__`` is enough
    for it to show up in all three.

    See Also
    --------
    AptamerDesigner : The one class built on this.

    Notes
    -----
    The rule that makes this work: ``__init__`` stores each argument under its
    own name and does nothing else — no validation, no computation, no reading
    files. Only then can an object be taken apart into its settings and rebuilt
    from them.

    Examples
    --------
    >>> designer = AptamerDesigner(n_nucleotides=3, n_samples=10)
    >>> designer.get_params()["n_samples"]
    10
    >>> AptamerDesigner(**designer.get_params()).n_nucleotides
    3
    """

    @classmethod
    def _param_names(cls) -> list[str]:
        """Return the names of this class's settings, sorted.

        Returns
        -------
        list of str
            Every argument of ``__init__`` except ``self``. Sorted, so that a
            ``repr`` and a settings dictionary are in a predictable order
            rather than in the order the arguments happen to be declared.
        """
        signature = inspect.signature(cls.__init__)
        return sorted(name for name in signature.parameters if name != "self")

    def get_params(self, deep: bool = True) -> dict[str, Any]:
        """Return this object's settings.

        Parameters
        ----------
        deep : bool, default=True
            Whether to include the settings of nested objects that have their
            own. No MAWS setting is such an object, so this makes no
            difference; it is accepted because tools that drive estimators pass
            it.

        Returns
        -------
        dict
            Setting name to value, in sorted order. Passing it back as
            ``**params`` to the same class builds an equivalent object.

        See Also
        --------
        set_params : Changes settings on an object that already exists.

        Examples
        --------
        >>> AptamerDesigner(n_nucleotides=3).get_params()["n_nucleotides"]
        3
        """
        return {name: getattr(self, name) for name in self._param_names()}

    def set_params(self, **params: Any) -> _BaseEstimator:
        """Change one or more settings in place.

        Parameters
        ----------
        **params
            Setting name to new value. Names not mentioned keep the value they
            have.

        Returns
        -------
        _BaseEstimator
            This object, changed in place, so calls can be chained. Not a copy.

        Raises
        ------
        maws.errors.ConfigurationError
            If a name is not one of this object's settings. The message lists
            the ones that are.

        See Also
        --------
        get_params : Reads the settings back out.

        Examples
        --------
        >>> designer = AptamerDesigner()
        >>> designer.set_params(n_nucleotides=7).n_nucleotides
        7

        Copy first to leave the original alone:

        >>> longer = AptamerDesigner(**designer.get_params()).set_params(beta=0.5)
        >>> designer.beta, longer.beta
        (0.01, 0.5)
        """
        allowed = self._param_names()
        for name, value in params.items():
            if name not in allowed:
                raise ConfigurationError(
                    f"{type(self).__name__} has no setting {name!r}. "
                    f"It has: {', '.join(allowed)}"
                )
            setattr(self, name, value)
        return self

    def __repr__(self) -> str:
        shown = ", ".join(
            f"{name}={value!r}" for name, value in sorted(self.get_params().items())
        )
        return f"{type(self).__name__}({shown})"


class AptamerDesigner(_BaseEstimator):
    """A configured aptamer design run, over one or more target molecules.

    .. warning::
        A run at the default settings performs on the order of a million energy
        evaluations per target and takes hours. Try ``n_nucleotides=3`` and
        ``n_samples=50`` first, to confirm the target files and the
        installation are sound.

    Parameters
    ----------
    n_nucleotides : int, default=15
        How many nucleotides each finished strand should have. Each one costs
        another full round of sampling, so the run time grows with it.
    aptamer : {"RNA", "DNA"}, default="RNA"
        Which nucleic acid to build the strands from. RNA strands are written
        with the letters ``G A U C`` and DNA with ``G A T C``.
    molecule : {"protein", "organic", "lipid"}, default="protein"
        What the targets are. This picks the force field used for them — the
        table of numbers that turns a set of atom positions into an energy.
        Proteins and lipids have ready-made parameters; an organic molecule has
        its own worked out at the start of each run, taking extra minutes.
    n_samples : int, default=5000
        How many shapes to try per candidate at each growth step. The main cost
        of a run, and roughly proportional to how long one takes. Raising it
        searches each candidate more thoroughly.
    n_first_samples : int, optional
        How many shapes to try on the first step, which also searches over
        where next to the target to put the strand. Raising it beyond
        `n_samples` buys a better starting position, which every later step
        builds on. Defaults to `n_samples`.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ. Raising it makes
        the choice between candidates depend mostly on their few best shapes;
        at zero every shape weighs the same and no candidate can win.
    salt_conc : float, default=0.15
        Concentration of dissolved salt in mol/L. Dissolved ions gather around
        charged atoms and damp the pull between them, so raising this weakens
        every electrostatic interaction. The default is roughly the saltiness
        of blood.
    reach : float, default=10.0
        How far past a target's furthest atom the sampling region extends, in
        ångström. Larger values consider positions further from the target, at
        the cost of more proposals landing in empty space. Only used when
        `sampling` is ``"sphere"``.
    probe : float, default=1.4
        Radius in ångström of the ball rolled over a target to find its
        surface. 1.4 Å is the size of a water molecule, the usual choice.
        Larger values smooth over narrow clefts, so no strand is offered a
        position inside one.
    sampling : {"sphere", "surface-following"}, default="sphere"
        Which region to draw positions from. ``"sphere"`` fills a ball around
        the target; ``"surface-following"`` keeps to a shell hugging its
        surface, which concentrates the search where contact is possible but
        rejects far more proposals.
    relax_iterations : int, default=50
        How many rounds of nudging the atoms and letting them settle to run
        after each nucleotide is joined on. Joining leaves the new residue
        strained against its neighbour, and each round works some of that out.
        Zero skips it and is much faster, at the cost of scoring strained
        shapes.
    random_state : int, optional
        Fixes the randomness, so the same settings give the same answer. Left
        out, every run differs.
    verbose : int, default=0
        How much to report while running. Zero is silent; anything else prints
        a line naming each target as it is started.

    Attributes
    ----------
    sequences_ : list of str
        The designed strand for each target, written 5' to 3', in the order the
        targets were given.
    energies_ : numpy.ndarray
        Shape ``(n_targets,)``. Best energy found for each, in kJ/mol. Energies
        for different targets are not comparable with each other.
    scores_ : numpy.ndarray
        Shape ``(n_targets,)``. The number that chose each strand, lower being
        better. With the default scorer it is a free energy in kJ/mol.
    results_ : list of maws.api.MawsResult
        The full result for each target, including its final structure and
        whether the run reached the requested length.
    n_steps_ : int
        The strand length the run aimed at, repeated here for convenience. A
        run that stopped early still leaves this at what was asked for; the
        matching entry of ``results_`` says what was actually reached.

    See Also
    --------
    maws.design : The same run as a single function call.
    maws.search.grow_aptamer : The same search, reporting step by step.
    maws.api.MawsResult : What each entry of ``results_`` is.

    Notes
    -----
    Every argument is stored exactly as given and nothing else happens until
    :meth:`fit` is called. Nothing is validated in the constructor, and no file
    is read, so building one of these is instant and cannot fail; a bad setting
    is reported by ``fit``.

    Targets are designed against one after another, in the order given. Each
    gets its own search, and `random_state` is used for every one of them, so
    two identical targets in the same list give identical answers.

    Examples
    --------
    >>> designer = AptamerDesigner(n_nucleotides=5, n_samples=50, random_state=0)
    >>> designer.n_nucleotides, designer.n_samples
    (5, 50)

    Copy it and change one setting, leaving the original alone:

    >>> longer = AptamerDesigner(**designer.get_params()).set_params(n_nucleotides=20)
    >>> designer.n_nucleotides, longer.n_nucleotides
    (5, 20)
    """

    def __init__(
        self,
        n_nucleotides: int = 15,
        *,
        aptamer: AptamerType = "RNA",
        molecule: MoleculeType = "protein",
        n_samples: int = 5000,
        n_first_samples: int | None = None,
        beta: float = 0.01,
        salt_conc: float = 0.15,
        reach: float = 10.0,
        probe: float = 1.4,
        sampling: SamplingMode = "sphere",
        relax_iterations: int = 50,
        random_state: int | None = None,
        verbose: int = 0,
    ) -> None:
        self.n_nucleotides = n_nucleotides
        self.aptamer = aptamer
        self.molecule = molecule
        self.n_samples = n_samples
        self.n_first_samples = n_first_samples
        self.beta = beta
        self.salt_conc = salt_conc
        self.reach = reach
        self.probe = probe
        self.sampling = sampling
        self.relax_iterations = relax_iterations
        self.random_state = random_state
        self.verbose = verbose

    def fit(self, X: Sequence[str | Path], y: object = None) -> AptamerDesigner:
        """Design an aptamer for every target in `X`.

        This is where the whole run happens; expect it to take hours at the
        default settings.

        Parameters
        ----------
        X : sequence of str or pathlib.Path
            Paths to the targets' structure files, in PDB format: the standard
            text format for a molecular structure, one atom per line with its
            position. One search is run per entry.
        y : ignored
            Accepted for consistency with the estimator convention, which
            passes labels here. A design run has none.

        Returns
        -------
        AptamerDesigner
            This object, changed in place, with its result attributes filled
            in. Returned so that ``designer.fit(targets).sequences_`` reads in
            one line.

        Raises
        ------
        maws.errors.ConfigurationError
            If a setting is not one of the allowed values, or `X` is empty.
        maws.errors.ToolchainError
            If AmberTools, the molecular modelling suite the structures are
            assembled with, is not installed.

        See Also
        --------
        fit_predict : Runs this and hands back the strands directly.

        Examples
        --------
        >>> designer = AptamerDesigner(n_nucleotides=5, n_samples=50)
        >>> designer.fit(["target.pdb"]).sequences_  # doctest: +SKIP
        ['G A U C G']

        Nothing to design against is an error, caught before any work starts:

        >>> designer.fit([])
        Traceback (most recent call last):
            ...
        maws.errors.ConfigurationError: fit needs at least one target file
        """
        targets = list(X)
        if not targets:
            raise ConfigurationError("fit needs at least one target file")
        self._check_params()

        results: list[MawsResult] = []
        for index, target in enumerate(targets):
            if self.verbose:
                print(f"[{index + 1}/{len(targets)}] designing against {target}")
            results.append(
                design(
                    target,
                    length=self.n_nucleotides,
                    aptamer=self.aptamer,
                    molecule=self.molecule,
                    samples=self.n_samples,
                    first_samples=self.n_first_samples,
                    beta=self.beta,
                    salt_conc=self.salt_conc,
                    reach=self.reach,
                    probe=self.probe,
                    sampling=self.sampling,
                    relax_iterations=self.relax_iterations,
                    seed=self.random_state,
                )
            )

        self.results_ = results
        self.sequences_ = [result.sequence for result in results]
        self.energies_ = np.array([result.energy for result in results])
        self.scores_ = np.array([result.score for result in results])
        self.n_steps_ = self.n_nucleotides
        return self

    def fit_predict(self, X: Sequence[str | Path], y: object = None) -> list[str]:
        """Design an aptamer for every target and return the strands.

        Parameters
        ----------
        X : sequence of str or pathlib.Path
            Paths to the targets' structure files, in PDB format.
        y : ignored
            Accepted for consistency with the estimator convention.

        Returns
        -------
        list of str
            One strand per target, in the order they were given. The same list
            :meth:`fit` leaves in ``sequences_``.

        See Also
        --------
        fit : Does the work, and keeps the energies and structures too.

        Notes
        -----
        There is no separate ``predict``, because a design run produces nothing
        that could be applied to a target it has not seen. Every target needs
        its own search.

        Examples
        --------
        >>> designer = AptamerDesigner(n_nucleotides=5, n_samples=50)
        >>> designer.fit_predict(["target.pdb"])  # doctest: +SKIP
        ['G A U C G']
        """
        return self.fit(X, y).sequences_

    def score(self, X: object = None, y: object = None) -> float:
        """Return how good the designs were, higher being better.

        Parameters
        ----------
        X : ignored
            Accepted for consistency with the estimator convention.
        y : ignored
            Accepted for consistency with the estimator convention. There is
            nothing to compare against: a design has no known right answer.

        Returns
        -------
        float
            Mean of the negated design scores across every target. Negated
            because a design score is lower-is-better, while a ``score``
            method is expected to be higher-is-better everywhere.

        Raises
        ------
        maws.errors.ConfigurationError
            If :meth:`fit` has not been called yet, so there is nothing to
            score.

        Examples
        --------
        Standing in for the attribute :meth:`fit` would have filled in, two
        targets scoring -800 and -400 average to 600 the other way up:

        >>> import numpy as np
        >>> designer = AptamerDesigner()
        >>> designer.scores_ = np.array([-800.0, -400.0])
        >>> designer.score()
        600.0
        """
        if not hasattr(self, "scores_"):
            raise ConfigurationError(
                "this designer has not been fitted yet; call fit first"
            )
        return float(-np.mean(self.scores_))

    def _check_params(self) -> None:
        """Check the settings, at the point they are about to be used.

        Raises
        ------
        maws.errors.ConfigurationError
            If a setting is out of range. Force field choices are left to
            :meth:`maws.forcefield.ForceField.for_target`, which already
            reports the allowed values.
        """
        if self.n_nucleotides < 1:
            raise ConfigurationError(
                f"n_nucleotides must be at least 1, got {self.n_nucleotides}"
            )
        if self.n_samples < 1:
            raise ConfigurationError(
                f"n_samples must be at least 1, got {self.n_samples}"
            )
        if self.beta < 0:
            raise ConfigurationError(f"beta must not be negative, got {self.beta}")

"""
maws._designer
==============

The design run as a configurable object, following scikit-learn's conventions.

:class:`AptamerDesigner` holds the settings for a run and does the work in
:meth:`~AptamerDesigner.fit`. Results appear as attributes whose names end in an
underscore, which is scikit-learn's mark for "this only exists once fit has
run".

The conventions are worth following even though scikit-learn is not a
dependency here. They give a familiar way to configure a run, copy it, change
one setting, and try again — and because the interface is the one scikit-learn
expects, its tools work on this class if it is installed.

MAWS fits the shape scikit-learn calls *transductive*: every target needs its
own search, and there is no model left over afterwards to apply to a new
target. So there is a :meth:`~AptamerDesigner.fit` and a
:meth:`~AptamerDesigner.fit_predict`, and no ``predict``.

Examples
--------
>>> designer = AptamerDesigner(n_nucleotides=5, n_samples=50, random_state=0)
>>> designer.get_params()["n_nucleotides"]
5
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
    """The part of scikit-learn's estimator interface MAWS needs.

    Provides :meth:`get_params`, :meth:`set_params` and a readable ``repr``,
    all derived from the arguments of ``__init__``. That is enough for
    scikit-learn's own ``clone`` and grid search to work on a subclass, without
    scikit-learn having to be installed to use one.

    Notes
    -----
    The rule that makes this work: ``__init__`` stores each argument under its
    own name and does nothing else — no validation, no computation, no reading
    files. Only then can an object be taken apart into its settings and rebuilt
    from them.
    """

    @classmethod
    def _param_names(cls) -> list[str]:
        """Return the names of this class's settings, sorted.

        Returns
        -------
        list of str
            Every argument of ``__init__`` except ``self``.
        """
        signature = inspect.signature(cls.__init__)
        return sorted(name for name in signature.parameters if name != "self")

    def get_params(self, deep: bool = True) -> dict[str, Any]:
        """Return this object's settings.

        Parameters
        ----------
        deep : bool, default=True
            Accepted for compatibility with scikit-learn. MAWS has no nested
            objects with their own settings, so it makes no difference.

        Returns
        -------
        dict
            Setting name to value.

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
            Setting name to new value.

        Returns
        -------
        _BaseEstimator
            This object, so calls can be chained.

        Raises
        ------
        maws.errors.ConfigurationError
            If a name is not one of this object's settings. The message lists
            the ones that are.

        Examples
        --------
        >>> designer = AptamerDesigner()
        >>> designer.set_params(n_nucleotides=7).n_nucleotides
        7
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
    """Design an aptamer against each of one or more target molecules.

    Parameters
    ----------
    n_nucleotides : int, default=15
        How many nucleotides each finished strand should have.
    aptamer : {"RNA", "DNA"}, default="RNA"
        Which nucleic acid to build the strand from.
    molecule : {"protein", "organic", "lipid"}, default="protein"
        What the targets are. Decides which parameters describe them.
    n_samples : int, default=5000
        How many shapes to try per candidate at each growth step. The main cost
        of a run.
    n_first_samples : int, optional
        How many shapes to try on the first step. Defaults to `n_samples`.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ.
    salt_conc : float, default=0.15
        Concentration of dissolved salt in mol/L.
    reach : float, default=10.0
        How far past a target's furthest atom the strand may be placed, in
        ångström.
    probe : float, default=1.4
        Radius in ångström of the ball rolled over a target to find its
        surface.
    sampling : {"sphere", "surface-following"}, default="sphere"
        Which region to draw positions from.
    relax_iterations : int, default=50
        How many nudge-and-settle rounds to run after each nucleotide is joined
        on.
    random_state : int, optional
        Fixes the randomness, so the same settings give the same answer.
    verbose : int, default=0
        How much to report while running. Zero is silent.

    Attributes
    ----------
    sequences_ : list of str
        The designed strand for each target, in the order the targets were
        given.
    energies_ : numpy.ndarray
        Shape ``(n_targets,)``. Best energy found for each, in kJ/mol.
    entropies_ : numpy.ndarray
        Shape ``(n_targets,)``. The score that chose each strand. Lower is
        better.
    results_ : list of maws.api.MawsResult
        The full result for each target, including its final structure.
    n_steps_ : int
        How many nucleotides were added per strand.

    See Also
    --------
    maws.design : The same run as a single function call.

    Notes
    -----
    Every argument is stored exactly as given and nothing else happens until
    :meth:`fit` is called. Nothing is validated in the constructor, and no file
    is read, so building one of these is instant and cannot fail.

    .. warning::
        A run at the default settings performs on the order of a million energy
        evaluations per target and takes hours.

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

        Parameters
        ----------
        X : sequence of str or pathlib.Path
            Paths to the targets' structure files, in PDB format.
        y : ignored
            Accepted for consistency with scikit-learn's interface, which
            passes labels here. MAWS has none.

        Returns
        -------
        AptamerDesigner
            This object, with its result attributes filled in.

        Raises
        ------
        maws.errors.ConfigurationError
            If a setting is not one of the allowed values, or `X` is empty.
        maws.errors.ToolchainError
            If AmberTools is not installed.
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
        self.entropies_ = np.array([result.entropy for result in results])
        self.n_steps_ = self.n_nucleotides
        return self

    def fit_predict(self, X: Sequence[str | Path], y: object = None) -> list[str]:
        """Design an aptamer for every target and return the strands.

        Parameters
        ----------
        X : sequence of str or pathlib.Path
            Paths to the targets' structure files.
        y : ignored
            Accepted for consistency with scikit-learn's interface.

        Returns
        -------
        list of str
            One strand per target, in the order they were given.

        Notes
        -----
        There is no separate ``predict``, because a design run produces nothing
        that could be applied to a target it has not seen. Every target needs
        its own search.
        """
        return self.fit(X, y).sequences_

    def score(self, X: object = None, y: object = None) -> float:
        """Return how good the designs were, higher being better.

        Parameters
        ----------
        X : ignored
            Accepted for consistency with scikit-learn's interface.
        y : ignored
            Accepted for consistency with scikit-learn's interface.

        Returns
        -------
        float
            Mean of the negated scores across every target. Negated because
            the design score is lower-is-better, while scikit-learn expects
            higher-is-better everywhere.

        Raises
        ------
        maws.errors.ConfigurationError
            If :meth:`fit` has not been called yet.
        """
        if not hasattr(self, "entropies_"):
            raise ConfigurationError(
                "this designer has not been fitted yet; call fit first"
            )
        return float(-np.mean(self.entropies_))

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

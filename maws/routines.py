"""
maws.routines
=============

Thermodynamic scoring function for the MAWS aptamer selection algorithm.

This module provides the entropy-based scoring function `S()` used to evaluate
aptamer-ligand binding configurations via Boltzmann-weighted sampling.

Functions
---------
S : Compute the entropy score from energy samples.

Notes
-----
Uses `mpmath` for arbitrary-precision arithmetic to avoid numerical underflow
when computing exponentials of large energy differences.

Examples
--------
>>> from maws.routines import S
>>> energies = [100.0, 150.0, 200.0, 175.0]
>>> score = S(energies, beta=0.01)
>>> hasattr(score, "__float__")  # S returns mpf (high-precision float)
True
"""

import numpy as np
from mpmath import mp as math

# Set high precision for numerical stability with large energy ranges
math.dps = 60


def _entropy(sample):
    """
    Compute Shannon entropy of a probability distribution.

    Internal helper - use :func:`S` as the public API.

    Parameters
    ----------
    sample : iterable
        Probability distribution P(x). Must sum to 1.

    Returns
    -------
    mpf
        Shannon entropy H = -sum(P * log(P * N)).
    """
    sample = list(sample)
    return -math.fsum(
        map(
            math.fmul,
            np.asarray(sample),
            map(
                math.log,
                map(math.fmul, np.asarray(sample), [len(sample)] * len(sample)),
            ),
        )
    )


def _zps(sample, beta):
    """
    Compute partition function, probabilities, and entropy from energy samples.

    Internal helper - use :func:`S` as the public API.

    Parameters
    ----------
    sample : array-like
        Energy values (kJ/mol) from conformational sampling.
    beta : float
        Inverse temperature parameter (1/kT in mol/kJ).

    Returns
    -------
    Z : mpf
        Partition function Z = sum(exp(-beta * E)).
    P : iterator
        Boltzmann probabilities P(i) = exp(-beta * E_i) / Z.
    entropy : mpf
        Shannon entropy of the Boltzmann distribution.
    """
    Z = math.fsum(
        map(math.exp, map(math.fmul, [-beta] * len(sample), np.asarray(sample)))
    )
    P = map(
        math.fdiv,
        map(math.exp, map(math.fmul, [-beta] * len(sample), np.asarray(sample))),
        np.asarray([Z] * len(sample)),
    )
    entropy = _entropy(P)
    return Z, P, entropy


def S(sample, beta=0.01):  # noqa: N802
    """
    Compute entropy score from energy samples.

    This is the primary scoring function used in MAWS. The entropy measures
    the "spread" of the Boltzmann distribution over sampled conformations.
    Lower (more negative) entropy indicates a more peaked distribution,
    suggesting stronger binding affinity.

    Parameters
    ----------
    sample : array-like
        Energy values (kJ/mol) from conformational sampling.
    beta : float, default=0.01
        Inverse temperature parameter. Higher beta = sharper distribution.

    Returns
    -------
    mpf
        Entropy of the Boltzmann distribution over sampled energies.

    Examples
    --------
    >>> energies = [100.0, 150.0, 200.0]
    >>> score = S(energies, beta=0.01)
    >>> hasattr(score, "__float__")  # Returns mpf (high-precision)
    True
    """
    return _zps(sample, beta)[2]

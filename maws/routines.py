"""
maws.routines
=============

Thermodynamic scoring function for the MAWS aptamer selection algorithm.

Functions
---------
entropy_score : Compute the entropy score from energy samples.

Notes
-----
Weights are evaluated in log-space via :func:`scipy.special.logsumexp`, so
energies spanning thousands of kJ/mol neither underflow nor overflow in double
precision.

Examples
--------
>>> from maws.routines import entropy_score
>>> round(entropy_score([100.0, 150.0, 200.0, 175.0], beta=0.01), 6)
-0.072433
"""

import numpy as np
from scipy.special import logsumexp
from scipy.stats import entropy as _relative_entropy


def _boltzmann(sample, beta):
    """
    Compute normalised Boltzmann probabilities from energy samples.

    Internal helper - use :func:`entropy_score` as the public API.

    Parameters
    ----------
    sample : array-like
        Energy values (kJ/mol) from conformational sampling.
    beta : float
        Inverse temperature parameter (1/kT in mol/kJ).

    Returns
    -------
    P : numpy.ndarray
        Boltzmann probabilities P(i) = exp(-beta * E_i) / Z, summing to 1.
    log_z : float
        Natural logarithm of the partition function Z = sum(exp(-beta * E)).
        Returned in log-space because Z itself overflows double precision for
        strongly negative energies.
    """
    log_weights = -beta * np.asarray(sample, dtype=float)
    log_z = float(logsumexp(log_weights))
    return np.exp(log_weights - log_z), log_z


def entropy_score(sample, beta=0.01):
    """
    Compute the entropy score of a set of sampled energies.

    This is the primary scoring function used in MAWS. It measures the spread
    of the Boltzmann distribution over sampled conformations: lower (more
    negative) scores indicate a more peaked distribution, suggesting stronger
    binding affinity.

    The score is ``-sum(P * log(P * N))``, the negative Kullback-Leibler
    divergence of the Boltzmann distribution from the uniform distribution. It
    is therefore at most 0, reaching 0 only when all sampled energies are equal.

    Parameters
    ----------
    sample : array-like
        Energy values (kJ/mol) from conformational sampling.
    beta : float, default=0.01
        Inverse temperature parameter. Higher beta = sharper distribution.

    Returns
    -------
    float
        Entropy of the Boltzmann distribution over sampled energies.

    Examples
    --------
    >>> round(entropy_score([100.0, 150.0, 200.0], beta=0.01), 6)
    -0.078421
    """
    probabilities, _ = _boltzmann(sample, beta)
    uniform = np.full(probabilities.size, 1.0 / probabilities.size)
    return -float(_relative_entropy(probabilities, uniform))

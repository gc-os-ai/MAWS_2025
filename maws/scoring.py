r"""
maws.scoring
============

The criterion MAWS ranks candidate nucleotides by.

A MAWS run grows an aptamer - a short strand of DNA or RNA meant to stick to
a chosen target molecule - one nucleotide at a time. At each step it tries
every nucleotide it could add, samples many shapes of each candidate strand
against the target, and records the potential energy of every shape. This
module turns one such list of energies into the single number the step is
decided on. The candidate scoring **lowest** is the one kept.

Where the criterion comes from
------------------------------

The score is not an invention of this package. It implements the Entropic
Fragment-Based Approach of Tseng et al. [1]_, which MAWS was built on [2]_,
and the formula here is the one those authors published.

The idea is that a nucleotide worth keeping is one whose shapes concentrate
into a narrow family rather than spreading over everything the sampler tried.
Concentration is measured as the distance of the Boltzmann distribution over
the sampled energies from a uniform distribution. A distribution that is
already uniform carries no information about where the strand prefers to sit,
and scores 0; a sharply peaked one scores far below 0.

Two properties of that choice are worth knowing before reading a score, and
both are described under :func:`entropy_score`: the score is unchanged by
shifting every energy by a constant, and a shape whose Boltzmann weight
underflows to zero leaves the sample entirely.

References
----------
.. [1] Tseng, C.-Y., Ashrafuzzaman, M., Mane, J. Y., Kapty, J., Mercer, J. R.,
       Tuszynski, J. A. (2011). "Entropic Fragment-Based Approach to Aptamer
       Design". Chemical Biology & Drug Design 78(1), 1-13. Also US patent
       8,484,010 B2.
.. [2] Kalinowski, M. et al. (2016). "MAWS - Making Aptamers Without SELEX".
       iGEM Heidelberg.

Examples
--------
>>> from maws.scoring import entropy_score
>>> round(entropy_score([100.0, 150.0, 200.0, 175.0], beta=0.01), 6)
-0.072433
"""

import numpy as np
from scipy.special import logsumexp
from scipy.stats import entropy as _relative_entropy


def _boltzmann(sample, beta):
    r"""
    Compute normalised Boltzmann probabilities from energy samples.

    Internal helper - use :func:`entropy_score` as the public API.

    Parameters
    ----------
    sample : array-like
        Energy values in kJ/mol, one per sampled shape.
    beta : float
        How sharply lower energies are favoured, in mol/kJ.

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
    r"""entropy_score(sample, beta=0.01) -> float

    Return how concentrated a candidate's sampled energies are.

    Turns the energies into a Boltzmann distribution and returns that
    distribution's distance from uniform, negated. The result is at most 0,
    reaching 0 when every sampled energy is equal, and falling towards
    ``-log N`` as the weight gathers onto a single shape. MAWS keeps the
    candidate scoring lowest.

    Parameters
    ----------
    sample : array-like
        Energy values in kJ/mol, one per sampled shape of one candidate
        strand. At least one is required.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ. Raising it makes
        the score depend mostly on the few lowest-energy shapes; at 0 every
        shape weighs the same and the score is 0.

    Returns
    -------
    float
        Zero when every shape is equally likely, and increasingly negative as
        the weight concentrates onto fewer shapes.

    Raises
    ------
    ValueError
        If `sample` is empty. There is no distribution to score, so no value
        would be meaningful.

    See Also
    --------
    maws.space.draw_clear_torsions : Keeps clashing shapes out of `sample`.
    maws.space.ClashFilter : Makes the accept/reject decision for a shape.

    Notes
    -----
    Each shape *i* is weighted by its Boltzmann factor, normalised to a
    probability, and the result is the negative relative entropy of that
    distribution against a uniform one over the same *N* shapes:

    .. math::
        p_i = \frac{e^{-\beta E_i}}{\sum_j e^{-\beta E_j}}
        \qquad
        S = -\sum_i p_i \ln(p_i N)

    The factor :math:`N` inside the logarithm puts the zero point at "every
    shape equally likely" whatever *N* is, so candidates sampled a different
    number of times stay comparable.

    Weights are evaluated in log-space via :func:`scipy.special.logsumexp`, so
    energies spanning thousands of kJ/mol neither underflow nor overflow.

    .. note::
        `beta` enters the source method as a Lagrange multiplier of the
        maximum-entropy derivation, not as a physical inverse temperature.
        The default of 0.01 is the value used there, where the ranking it
        produced held across every multiplier tried. Read as :math:`1/RT` it
        would correspond to about 12,000 K, which is why it does not match
        the 0.401 mol/kJ of a 300 K calculation.

    .. warning::
        The score reads the spread of `sample`, never the absolute energies.
        Adding the same constant to every energy leaves it unchanged, so a
        candidate whose shapes all sit at -5000 kJ/mol and one whose shapes
        all sit at +5000 kJ/mol score identically. It measures how tightly a
        candidate settles, not how strongly it binds.

    .. warning::
        A shape at around 1e8 kJ/mol, which is what atoms placed on top of
        each other cost, has a Boltzmann weight of exactly 0 in double
        precision. It drops out of the sum, so a sample holding such shapes
        scores as though it had been sampled fewer times, which lowers the
        score. Keep them out of `sample` rather than relying on the weighting
        to discount them.

    Examples
    --------
    Ten shapes, one far better than the rest, against ten much of a muchness.
    The first scores lower, meaning more promising.

    >>> one_clear_winner = entropy_score([0.0] + [1000.0] * 9)
    >>> nothing_to_choose = entropy_score([0.0] + [10.0] * 9)
    >>> one_clear_winner < nothing_to_choose
    True

    Equal energies score exactly zero.

    >>> abs(entropy_score([100.0, 100.0, 100.0])) < 1e-12
    True

    Shifting every energy by a constant changes nothing.

    >>> shifted = entropy_score([5100.0, 5150.0, 5200.0])
    >>> abs(entropy_score([100.0, 150.0, 200.0]) - shifted) < 1e-12
    True
    """
    energies = np.asarray(sample, dtype=float)
    if energies.size == 0:
        raise ValueError("sample must contain at least one energy value")

    probabilities, _ = _boltzmann(energies, beta)
    uniform = np.full(probabilities.size, 1.0 / probabilities.size)
    return -float(_relative_entropy(probabilities, uniform))

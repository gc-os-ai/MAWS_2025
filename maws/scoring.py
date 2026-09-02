r"""
maws.scoring
============

Score candidate nucleotides so a MAWS run can choose between them.

A MAWS run grows an aptamer one nucleotide at a time. An aptamer is a short
strand of DNA or RNA meant to stick to a chosen target molecule. At each step
the run tries every nucleotide it could add. It samples many shapes of each
candidate strand against the target and records the potential energy of every
shape. This module turns one such list of energies into the single number the
step is decided on. The candidate scoring **lowest** is the one kept.

Where the criterion comes from
------------------------------

The score implements the Entropic Fragment-Based Approach of Tseng et
al. [1]_, the method MAWS was built on [2]_. The formula here is the one
those authors published.

Sampling one candidate produces a list of energies, one for each shape
tried. The first step turns those energies into probabilities, weighting
each shape by :math:`e^{-\beta E}`. A low-energy shape comes out likely, a
high-energy one unlikely, and the probabilities across all the shapes add
up to 1.

Read those probabilities as a preference. Spread evenly over every shape,
the strand favours no particular way of sitting against the target. Piled
onto a few shapes, it favours those few. EFBA takes a strong preference as
its evidence of a good binder, on the reasoning that a strand which fits
the target settles into a small number of shapes.

The score reports how far the probabilities sit from evenly spread. Evenly
spread scores 0. The further the weight piles onto a few shapes, the
further below 0 the score falls, down to a floor of ``-log N`` for *N*
shapes. MAWS keeps the candidate scoring lowest.

Two things the score does not do are easy to assume it does. Read both
warnings under :func:`entropy_score` before you interpret one.

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

    Internal helper. :func:`entropy_score` is the public API.

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

    Return how tightly a candidate's sampled conformations cluster.

    A conformation is one shape the strand can take against the target.
    Sampling a candidate gives one energy per conformation. Those energies
    become a Boltzmann distribution, and the result is how far that
    distribution sits from uniform, with the sign flipped so that lower
    means more concentrated. MAWS keeps the candidate scoring lowest.

    Parameters
    ----------
    sample : array-like
        The energy of each sampled conformation of one candidate strand, in
        kJ/mol. At least one is required.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ. Raising it makes
        the score depend mostly on the few lowest-energy conformations. At 0
        every conformation weighs the same and the score is 0.

    Returns
    -------
    float
        A value between ``-log N`` and 0, for a sample of *N* conformations.
        It reaches 0 when every energy is equal, and approaches ``-log N``
        as the weight gathers onto a single conformation.

    Raises
    ------
    ValueError
        If `sample` is empty. There is no distribution to score, so no value
        would be meaningful.

    See Also
    --------
    maws.space.draw_clear_torsions : Keeps clashes out of `sample`.
    maws.run.MawsRunner : Runs the search this score decides.

    Notes
    -----
    Each conformation *i* carries its Boltzmann weight, normalised to a
    probability. The result is the negative relative entropy of that
    distribution against a uniform one over the same *N* conformations:

    .. math::
        p_i = \frac{e^{-\beta E_i}}{\sum_j e^{-\beta E_j}}
        \qquad
        S = -\sum_i p_i \ln(p_i N)

    Weights are evaluated in log-space via :func:`scipy.special.logsumexp`,
    so energies spanning thousands of kJ/mol neither underflow nor overflow.

    .. warning::
        Compare scores only between candidates sampled the same number of
        times. The factor :math:`N` fixes the zero point at "every
        conformation equally likely" for any *N*, but the far end of the
        range is ``-log N``. One dominant conformation therefore scores
        about -4.6 at *N* = 100 and about -6.4 at *N* = 1000, for the same
        situation. MAWS meets this by drawing a fixed number of
        conformations per step.

    .. warning::
        The score reads only the spread of `sample`. Adding the same
        constant to every energy leaves it unchanged. A candidate whose
        conformations all sit at -5000 kJ/mol and one whose conformations
        all sit at +5000 kJ/mol therefore score identically. The score
        measures how tightly a candidate settles. Binding strength needs a
        separate term.

    .. warning::
        Atoms placed on top of each other cost around 1e8 kJ/mol. In double
        precision that conformation's Boltzmann weight is exactly 0, so it
        leaves the distribution. A sample holding such conformations then
        scores as though it had been sampled fewer times, which lowers the
        score. Since MAWS keeps the lowest score, clashes make a candidate
        look better. Keep them out of `sample`. The weighting will not
        discount them for you.

    .. note::
        `beta` enters the source method as a Lagrange multiplier of the
        maximum-entropy derivation. The default of 0.01 is the value those
        authors used, and they report their nucleotide ranking held across
        every multiplier they tried. Read instead as :math:`1/RT`, 0.01
        corresponds to about 12,000 K. That is why it differs from the
        0.401 mol/kJ of a 300 K calculation.

    References
    ----------
    .. [1] Tseng, C.-Y., Ashrafuzzaman, M., Mane, J. Y., Kapty, J., Mercer,
           J. R., Tuszynski, J. A. (2011). "Entropic Fragment-Based Approach
           to Aptamer Design". Chemical Biology & Drug Design 78(1), 1-13.

    Examples
    --------
    Two samples of ten conformations each. In the first, one conformation
    sits 1000 kJ/mol below the other nine, so nearly all the weight lands
    on it. In the second the gap is 10 kJ/mol, so the weight stays spread
    across all ten. The concentrated sample scores lower, and MAWS keeps
    the lowest.

    >>> concentrated = entropy_score([0.0] + [1000.0] * 9)
    >>> spread = entropy_score([0.0] + [10.0] * 9)
    >>> concentrated < spread
    True

    Equal energies score exactly zero.

    >>> abs(entropy_score([100.0, 100.0, 100.0])) < 1e-12
    True

    Shifting every energy by a constant changes nothing.

    >>> shifted = entropy_score([5100.0, 5150.0, 5200.0])
    >>> abs(entropy_score([100.0, 150.0, 200.0]) - shifted) < 1e-12
    True

    The same situation, one dominant conformation, scores differently at
    different sample sizes. These two numbers are not comparable.

    >>> round(entropy_score([0.0] + [1000.0] * 99), 3)
    -4.556
    >>> round(entropy_score([0.0] + [1000.0] * 999), 3)
    -6.43
    """
    energies = np.asarray(sample, dtype=float)
    if energies.size == 0:
        raise ValueError("sample must contain at least one energy value")

    probabilities, _ = _boltzmann(energies, beta)
    uniform = np.full(probabilities.size, 1.0 / probabilities.size)
    return -float(_relative_entropy(probabilities, uniform))

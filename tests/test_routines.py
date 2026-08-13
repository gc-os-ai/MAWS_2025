"""
Tests for maws.routines module.

This module tests the entropy scoring function used in the MAWS algorithm:
- entropy_score(): Compute entropy from Boltzmann-weighted energy samples
- Lower entropy indicates stronger binding (sharper distribution)
- Evaluated in log-space so widely spread energies do not underflow
"""

import math

import numpy as np
import pytest

from maws.routines import _boltzmann, entropy_score


class TestEntropyScore:
    """Tests for the entropy_score() scoring function."""

    def test_returns_float(self):
        """entropy_score() returns a plain Python float."""
        score = entropy_score([100.0, 150.0, 200.0], beta=0.01)
        assert isinstance(score, float)

    def test_different_energies_finite_entropy(self):
        """A non-uniform Boltzmann distribution scores to a finite value."""
        score = entropy_score([100.0, 200.0, 300.0], beta=0.01)
        assert math.isfinite(score)

    def test_varying_beta(self):
        """Higher beta sharpens the distribution and lowers the score."""
        energies = [100.0, 150.0, 200.0]
        score_low_beta = entropy_score(energies, beta=0.001)
        score_high_beta = entropy_score(energies, beta=0.1)
        assert score_high_beta < score_low_beta

    def test_single_sample(self):
        """A single energy sample scores 0."""
        assert entropy_score([100.0], beta=0.01) == pytest.approx(0.0, abs=1e-10)

    def test_identical_energies(self):
        """Identical energies give a uniform distribution, which scores 0."""
        assert entropy_score([100.0] * 4, beta=0.01) == pytest.approx(0.0, abs=1e-10)

    @pytest.mark.parametrize("empty", [[], (), np.array([])])
    def test_rejects_empty_sample(self, empty):
        """An empty sample raises ValueError, not ZeroDivisionError."""
        with pytest.raises(ValueError, match="at least one energy value"):
            entropy_score(empty, beta=0.01)

    def test_is_non_positive(self):
        """The score is -KL(P || uniform), so it never exceeds 0."""
        rng = np.random.default_rng(0)
        for scale in (1.0, 1e3, 1e5):
            score = entropy_score(rng.normal(0.0, scale, 500).tolist(), beta=0.01)
            assert score <= 1e-12

    def test_matches_closed_form(self):
        """The score agrees with a direct sum(P * log(P * N)) evaluation."""
        energies = [-1500.0, -1490.0, -1450.0, -1200.0, -1199.5]
        weights = np.exp(-0.01 * np.asarray(energies))
        p = weights / weights.sum()
        expected = -np.sum(p * np.log(p * p.size))
        assert entropy_score(energies, beta=0.01) == pytest.approx(expected, rel=1e-12)

    def test_extreme_energies_do_not_underflow(self):
        """Energies large enough to underflow naive exp() still score finitely.

        At beta=0.01 these give exp(-1e4)-scale weights, which flush to zero in
        double precision unless the log-sum-exp shift is applied.
        """
        rng = np.random.default_rng(1)
        score = entropy_score(rng.normal(1e6, 5e5, 2000).tolist(), beta=0.01)
        assert math.isfinite(score)
        assert score < 0.0

    def test_is_shift_invariant(self):
        """Adding a constant to every energy leaves the score unchanged."""
        energies = [100.0, 150.0, 200.0, 175.0]
        shifted = [e + 5e5 for e in energies]
        assert entropy_score(shifted, beta=0.01) == pytest.approx(
            entropy_score(energies, beta=0.01), rel=1e-9
        )

    def test_accepts_numpy_array(self):
        """entropy_score() accepts an ndarray as well as a list."""
        energies = [100.0, 150.0, 200.0]
        assert entropy_score(np.asarray(energies), beta=0.01) == pytest.approx(
            entropy_score(energies, beta=0.01)
        )


class TestInternalFunctions:
    """Tests for internal helper functions."""

    def test_boltzmann_returns_normalised_probabilities(self):
        """_boltzmann() returns probabilities summing to 1 and a finite log Z."""
        p, log_z = _boltzmann([100.0, 150.0, 200.0], beta=0.01)
        assert p.shape == (3,)
        assert p.sum() == pytest.approx(1.0)
        assert math.isfinite(log_z)

    def test_boltzmann_lower_energy_higher_probability(self):
        """Lower energy states have higher Boltzmann probability."""
        p, _ = _boltzmann([100.0, 200.0], beta=0.1)
        assert p[0] > p[1]

    def test_boltzmann_probabilities_are_reusable(self):
        """P is a materialised array, not a one-shot iterator."""
        p, _ = _boltzmann([100.0, 150.0, 200.0], beta=0.01)
        assert p.sum() == pytest.approx(1.0)
        assert p.sum() == pytest.approx(1.0)

    def test_boltzmann_log_z_matches_direct_sum(self):
        """log_z equals log(sum(exp(-beta * E))) for moderate energies."""
        energies = [100.0, 150.0, 200.0]
        _, log_z = _boltzmann(energies, beta=0.01)
        expected = np.log(np.exp(-0.01 * np.asarray(energies)).sum())
        assert log_z == pytest.approx(expected, rel=1e-12)

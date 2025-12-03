"""
test_binaries
==========

This file contains tests for the binary population synthesis methods.
"""
import astropy.units as u
import numpy as np
from numpy.testing import assert_allclose

import arsenal_gear


class TestLogUniformPeriod:
    """Tests for the Log-Uniform Period distribution."""
    N = int(1e6)
    min_p = 1*u.d
    max_p = 1e4*u.d
    period = arsenal_gear.dist_funcs.binaries.LogUniformPeriod(min_p, max_p)

    def test_limits(self):
        """Check that the distribution limits are correct."""
        p = self.period.sample(self.N)
        assert p.max() <= self.max_p
        assert p.min() >= self.min_p

    def test_count(self):
        """Check that the number of samples is correct."""
        p = self.period.sample(self.N)
        assert len(p) == self.N

    def test_mean(self):
        """Check that the mean of the distribution is correct."""
        p = self.period.sample(self.N)
        assert_allclose(np.mean(p), u.d*self.period.mean(),
                        rtol=0, atol=5*self.period.var())

class TestLogUniformSemimajor:
    """Tests for the Log-Uniform Semimajor axis distribution."""
    N = int(1e6)
    min_a = 0.1*u.au
    max_a = 100*u.au
    sma = arsenal_gear.dist_funcs.binaries.LogUniformSemimajor(min_a, max_a)

    def test_limits(self):
        """Check that the distribution limits are correct."""
        a = self.sma.sample(self.N)
        assert a.max() <= self.max_a
        assert a.min() >= self.min_a

    def test_count(self):
        """Check that the number of samples is correct."""
        a = self.sma.sample(self.N)
        assert len(a) == self.N

    def test_mean(self):
        """Check that the mean of the distribution is correct."""
        a = self.sma.sample(self.N)
        assert_allclose(np.mean(a), u.au*self.sma.mean(),
                        rtol=0, atol=5*self.sma.var())

class TestUniformMassRatio:
    """Tests for the Uniform Mass Ratio distribution."""
    N = int(1e6)
    min_q = 0.1
    max_q = 1
    mratio = arsenal_gear.dist_funcs.binaries.UniformMassRatio(min_q, max_q)

    def test_limits(self):
        """Check that the distribution limits are correct."""
        q = self.mratio.sample(self.N)
        assert q.max() <= self.max_q
        assert q.min() >= self.min_q

    def test_count(self):
        """Check that the number of samples is correct."""
        q = self.mratio.sample(self.N)
        assert len(q) == self.N

    def test_mean(self):
        """Check that the mean of the distribution is correct."""
        q = self.mratio.sample(self.N)
        assert_allclose(np.mean(q), self.mratio.mean(),
                        rtol=0, atol=5*self.mratio.var())

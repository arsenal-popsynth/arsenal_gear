import astropy.units as u
import numpy as np
import pytest
from scipy.stats import sampling
from numpy.testing import assert_allclose

import arsenal_gear

class TestLogUniformSemimajor:
    N = int(1e6)
    min_a = 0.1*u.au
    max_a = 100*u.au
    sma = arsenal_gear.dist_funcs.binaries.LogUniformSemimajor(min_a, max_a)
    
    def test_limits(self):
        a = self.sma.sample(self.N)
        assert a.max() <= self.max_a
        assert a.min() >= self.min_a

    def test_count(self):
        a = self.sma.sample(self.N)
        assert len(a) == self.N

    def test_mean(self):
        a = self.sma.sample(self.N)
        assert_allclose(np.mean(a), u.au*self.sma.mean(), rtol=0, atol=5*self.sma.var())

class TestUniformMassRatio:
    N = int(1e6)
    min_q = 0.1
    max_q = 1
    mratio = arsenal_gear.dist_funcs.binaries.UniformMassRatio(min_q, max_q)
    
    def test_limits(self):
        q = self.mratio.sample(self.N)
        assert q.max() <= self.max_q
        assert q.min() >= self.min_q

    def test_count(self):
        q = self.mratio.sample(self.N)
        assert len(q) == self.N

    def test_mean(self):
        q = self.mratio.sample(self.N)
        assert_allclose(np.mean(q), self.mratio.mean(), rtol=0, atol=5*self.mratio.var())

class TestUniformEccentricity:
    N = int(1e6)
    min_e = 0
    max_e = 0.9
    ecc = arsenal_gear.dist_funcs.binaries.UniformEccentricity(min_e, max_e)
    
    def test_limits(self):
        e = self.ecc.sample(self.N)
        assert e.max() <= self.max_e
        assert e.min() >= self.min_e

    def test_count(self):
        e = self.ecc.sample(self.N)
        assert len(e) == self.N

    def test_mean(self):
        e = self.ecc.sample(self.N)
        assert_allclose(np.mean(e), self.ecc.mean(), rtol=0, atol=5*self.ecc.var())

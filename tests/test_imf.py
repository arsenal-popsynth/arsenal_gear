"""
test_imf
==========

This file contains unit tests for the Initial Mass Function (IMF) implementations.
"""
import astropy.units as u
import numpy as np
from numpy.testing import assert_allclose

import arsenal_gear


class TestSalpeter:
    """Tests for the Salpeter IMF implementation."""
    N = int(1e6)
    min_mass = 1*u.Msun
    max_mass = 100*u.Msun
    imf = arsenal_gear.dist_funcs.imf.Salpeter(min_mass=min_mass,
                                                    max_mass=max_mass, seed=1337)

    def test_limits(self):
        """Ensure sampled masses are within specified limits."""
        masses = self.imf.sample(self.N)
        print(masses)
        assert masses.max() <= self.max_mass
        assert masses.min() >= self.min_mass

    def test_count(self):
        """Ensure the correct number of stars are drawn."""
        masses = self.imf.sample(self.N)
        assert len(masses) == self.N

    def test_meanmass(self):
        """Ensure the distribution matches the expected mean stellar mass."""
        masses = self.imf.sample(self.N)
        assert_allclose(np.mean(masses), u.Msun*self.imf.mean(), rtol=0, atol=5*self.imf.var())

    def test_mass_sample(self):
        """Ensure the distribution matches the expected total stellar mass."""
        masses = self.imf.sample_mass(1e5*u.Msun)
        assert_allclose(np.sum(masses), 1e5*u.Msun, rtol=0.02) 

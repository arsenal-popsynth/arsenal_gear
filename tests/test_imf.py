"""
test_imf
==========

This file contains unit tests for the Initial Mass Function (IMF) implementations.
"""

import astropy.units as u
import numpy as np
import pytest

from arsenal_gear.dist_funcs.imf import Chabrier, Kroupa, MillerScalo, Salpeter


@pytest.mark.parametrize("imf_class", [Salpeter, Kroupa, MillerScalo, Chabrier])
@pytest.mark.parametrize(
    "limits",
    [
        (0.1, 100) * u.Msun,  # mass_min, mass_max -- Standard range
    ],
)
class TestIMFs:
    """All tests in this class will run for every IMF and every limit combo"""

    def test_normalization(self, imf_class, limits):
        """Run normalization check on every IMF implementation"""
        min_mass, max_mass = limits
        imf = imf_class(min_mass=min_mass, max_mass=max_mass)

        m_vals = np.logspace(np.log10(min_mass.value), np.log10(max_mass.value), 10000)
        pdf_vals = imf.pdf(m_vals)
        area = np.trapezoid(pdf_vals, m_vals)
        assert area == pytest.approx(1.0, abs=1e-3)

    def test_limits(self, imf_class, limits):
        """
        Ensure sampled masses are within specified limits
        for every IMF implementation
        """
        min_mass, max_mass = limits
        imf = imf_class(min_mass=min_mass, max_mass=max_mass)
        N = int(1e6)
        masses = imf.sample(N)

        expected_max = imf.max_mass
        expected_min = imf.min_mass
        assert masses.max().to(u.Msun).value <= expected_max * (1 + 1e-7)
        assert masses.min().to(u.Msun).value >= expected_min * (1 - 1e-7)

    def test_mean_mass(self, imf_class, limits):
        min_mass, max_mass = limits
        imf = imf_class(min_mass=min_mass, max_mass=max_mass)

        N = int(1e6)
        samples = imf.sample(N)

        actual_mean = np.mean(samples).to(u.Msun).value
        theoretical_mean = imf.mean().to(u.Msun).value

        assert actual_mean == pytest.approx(theoretical_mean, rel=2e-2)

    def test_count(self, imf_class, limits):
        min_mass, max_mass = limits
        imf = imf_class(min_mass=min_mass, max_mass=max_mass)
        N = int(1e3)
        masses = imf.sample(N)
        assert len(masses) == N

    def test_sample_mass(self, imf_class, limits):
        min_mass, max_mass = limits
        imf = imf_class(min_mass=min_mass, max_mass=max_mass)

        target_mtot = 1e5 * u.Msun
        masses = imf.sample_mass(target_mtot)
        actual_mtot = np.sum(masses).to_value(u.Msun)
        expected_mtot = target_mtot.to_value(u.Msun)
        assert actual_mtot == pytest.approx(expected_mtot, rel=0.1)

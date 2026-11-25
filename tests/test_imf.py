"""
test_imf
==========

This file contains unit tests for the Initial Mass Function (IMF) implementations.
"""

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.integrate import quad

from arsenal_gear.formation.dist_funcs.imf import (
    IMF,
    Chabrier,
    Kroupa1993,
    Kroupa2001,
    MillerScalo,
    PiecewisePowerLaw,
    Salpeter,
)

ALL_IMFS = [Salpeter, Kroupa1993, Kroupa2001, MillerScalo, Chabrier]


@pytest.mark.parametrize("imf_class", ALL_IMFS)
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
        masses = imf.rvs(size=N) * u.Msun

        expected_max = imf.max_mass
        expected_min = imf.min_mass

        assert masses.max().to(u.Msun).value <= expected_max
        assert masses.min().to(u.Msun).value >= expected_min

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


# --------------------------------------------------------------------------
# Mass ranges that cut through or skip IMF segments.
# The tests above only use (0.1, 100) Msun, and IMF.sample() clips to the
# limits, so sampling bugs elsewhere are hidden. These tests call the
# unclipped ppf directly.
# --------------------------------------------------------------------------
RANGES = [(0.1, 100), (0.08, 0.5), (0.6, 0.9), (2, 100), (20, 120)]


@pytest.mark.parametrize("imf_class", ALL_IMFS)
@pytest.mark.parametrize("lo, hi", RANGES)
class TestIMFRanges:
    """Sampling and normalization for mass ranges around the IMF breakpoints"""

    Q = np.linspace(0.0, 1.0, 2001)[1:-1]

    def test_normalization(self, imf_class, lo, hi):
        """The PDF integrates to one over [min_mass, max_mass]"""
        imf = imf_class(min_mass=lo * u.Msun, max_mass=hi * u.Msun)
        area, _ = quad(imf.pdf, lo, hi, points=[0.08, 0.5, 1.0, 10.0], limit=200)
        assert area == pytest.approx(1.0, abs=1e-3)

    def test_ppf_within_limits(self, imf_class, lo, hi):
        """The unclipped inverse CDF maps (0, 1) into [min_mass, max_mass]"""
        imf = imf_class(min_mass=lo * u.Msun, max_mass=hi * u.Msun)
        masses = imf.ppf(self.Q)
        assert np.all(np.isfinite(masses))
        assert np.all(masses >= lo * (1 - 1e-9))
        assert np.all(masses <= hi * (1 + 1e-9))

    def test_cdf_ppf_roundtrip(self, imf_class, lo, hi):
        """cdf(ppf(q)) == q, i.e. sampling matches the PDF"""
        imf = imf_class(min_mass=lo * u.Msun, max_mass=hi * u.Msun)
        q = self.Q[::50]
        assert_allclose(imf.cdf(imf.ppf(q)), q, atol=1e-3)


@pytest.mark.parametrize(
    "imf_class, breakpoints",
    [
        (Kroupa1993, [0.5, 1.0]),
        (Kroupa2001, [0.5]),
        (MillerScalo, [1.0, 10.0]),
        (Chabrier, [1.0]),
    ],
)
def test_pdf_continuous_at_breakpoints(imf_class, breakpoints):
    """Piecewise IMFs have no jumps at their segment boundaries"""
    imf = imf_class(min_mass=0.1 * u.Msun, max_mass=100 * u.Msun)
    eps = 1e-8
    for m in breakpoints:
        assert imf.pdf(m * (1 - eps)) == pytest.approx(imf.pdf(m * (1 + eps)), rel=1e-5)


# --------------------------------------------------------------------------
# Constructor validation and the IMF registry
# --------------------------------------------------------------------------
@pytest.mark.parametrize("imf_class", ALL_IMFS)
@pytest.mark.parametrize(
    "lo, hi, match",
    [
        (10, 1, "must be less than"),
        (5, 5, "must be less than"),
        (0, 100, "must both be positive"),
        (-1, 100, "must both be positive"),
    ],
)
def test_invalid_limits(imf_class, lo, hi, match):
    """Non-positive or inverted mass limits raise a ValueError"""
    with pytest.raises(ValueError, match=match):
        imf_class(min_mass=lo * u.Msun, max_mass=hi * u.Msun)


@pytest.mark.parametrize(
    "name, imf_class",
    [
        ("salpeter", Salpeter),
        ("SALPETER", Salpeter),
        ("kroupa1993", Kroupa1993),
        ("kroupa2001", Kroupa2001),
        ("millerscalo", MillerScalo),
        ("chabrier", Chabrier),
    ],
)
def test_get_imf(name, imf_class):
    """IMF.get_imf resolves case-insensitive names and forwards kwargs"""
    imf = IMF.get_imf(name, min_mass=0.5 * u.Msun, max_mass=50 * u.Msun)
    assert isinstance(imf, imf_class)
    assert imf.min_mass == pytest.approx(0.5)
    assert imf.max_mass == pytest.approx(50)


def test_get_imf_unknown():
    """Unknown IMF names raise a ValueError"""
    with pytest.raises(ValueError, match="Unknown IMF string"):
        IMF.get_imf("not_an_imf")


def test_get_imf_instance_passthrough():
    """Passing an IMF instance to get_imf returns it unchanged"""
    imf = Salpeter()
    assert IMF.get_imf(imf) is imf


@pytest.mark.parametrize("imf_class", ALL_IMFS)
def test_seed_reproducibility(imf_class):
    """The same seed gives the same sample"""
    a = imf_class(seed=11).sample(1000)
    b = imf_class(seed=11).sample(1000)
    assert_array_equal(a, b)


@pytest.mark.parametrize(
    "imf_class, expected",
    [
        (Salpeter, "Salpeter"),
        (Kroupa1993, "Kroupa1993"),
        (Kroupa2001, "Kroupa2001"),
        (MillerScalo, "MillerScalo"),
        (Chabrier, "Chabrier"),
    ],
)
def test_names(imf_class, expected):
    """Each IMF reports its own name"""
    assert imf_class().name == expected


# --------------------------------------------------------------------------
# Salpeter slope edge cases
# --------------------------------------------------------------------------
def test_salpeter_alpha_one():
    """alpha=1 (dN/dm ~ 1/m) is a valid power law and should be usable"""
    imf = Salpeter(alpha=1.0)
    assert np.isfinite(imf.mean().value)
    assert 0.08 < imf.ppf(0.5) < 100


def test_salpeter_alpha_two_mean():
    """alpha=2 has a finite, analytic mean: ln(m2/m1) / (1/m1 - 1/m2)"""
    imf = Salpeter(alpha=2.0)
    expected = np.log(100 / 0.08) / (1 / 0.08 - 1 / 100)
    assert imf.mean().value == pytest.approx(expected)


def test_salpeter_negative_alpha():
    """A negative slope should raise a ValueError (asserts vanish under -O)"""
    with pytest.raises(ValueError):
        Salpeter(alpha=-1.0)


@pytest.mark.parametrize(
    ("masses", "betas", "lo", "hi"),
    [
        ([1.0], [-1.0, -2.35], 0.1, 100.0),
        ([0.5, 2.0], [-1.3, -1.0, -2.7], 0.08, 50.0),
        ([], [-1.0], 0.2, 10.0),
    ],
)
def test_piecewise_power_law_beta_minus_one(masses, betas, lo, hi):
    """PiecewisePowerLaw handles beta == -1 segments cleanly across pdf, ppf, cdf, and mean"""
    imf = PiecewisePowerLaw(
        min_mass=lo * u.Msun,
        max_mass=hi * u.Msun,
        masses=masses,
        betas=betas,
    )

    area, _ = quad(imf.pdf, lo, hi, points=masses, limit=200)
    assert area == pytest.approx(1.0, abs=1e-3)

    q = np.linspace(0.01, 0.99, 200)
    masses_out = imf.ppf(q)
    assert np.all(np.isfinite(masses_out))
    assert np.all(masses_out >= lo * (1 - 1e-9))
    assert np.all(masses.out <= hi * (1 + 1e-9))

    assert_allclose(imf.cdf(imf.ppf(q)), q, rtol=1e-4, atol=1e-4)

    mean_mass = imf.mean()
    assert mean_mass.unit == u.Msun
    assert np.isfinite(mean_mass.value)
    assert lo < mean_mass.value < hi

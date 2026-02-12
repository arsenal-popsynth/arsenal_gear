"""
test_imf
==========

This file contains unit tests for the stellar population constructors.
"""

import astropy.units as u
import numpy as np
from numpy.testing import assert_array_equal

from arsenal_gear.formation.form_data_structures import SinglePop, BinaryPop


def test_single_pop():
    """Test the StarPopulation constructor for single stars."""
    N = 100
    mass = u.Msun * np.ones(N)
    Mtot = np.sum(mass)
    metallicity = 0.1 * u.dimensionless_unscaled
    star = SinglePop(Mtot=Mtot,
                     Nstar=N,
                     metallicity=metallicity,
                     imf=None,
                     mmin=0.1 * u.Msun,
                     mmax=100 * u.Msun,
                     discrete=True,
                     masses=mass
    )

    assert_array_equal(star.masses, mass)
    assert star.metallicity == metallicity


def test_binary_pop():
    """Test the BinaryPopulation constructor for binary stars."""
    N = 100
    masses = u.Msun * np.ones(N)
    Mtot = np.sum(masses) 
    metallicity = 0.1 * u.dimensionless_unscaled
    mmin = 0.1 * u.Msun
    mmax = 100 * u.Msun
    mrats = 0.5 * np.ones(N)
    periods = 10 * u.d * np.ones(N)
    semimajors = 100 * u.AU * np.ones(N)
    binary = BinaryPop(Mtot=Mtot,
                       Nstar=N,
                       metallicity=metallicity,
                       imf=None,
                       mmin=mmin,
                       mmax=mmax,
                       discrete=True,
                       masses=masses,
                       q_dist=None,
                       period_dist=None,
                       a_dist=None,
                       mrats=mrats,
                       periods=periods,
                       semimajors=semimajors
    )
    # Check the binary initilization
    assert_array_equal(binary.masses, masses)
    assert_array_equal(binary.mrats, mrats)
    assert_array_equal(binary.periods, periods)
    assert_array_equal(binary.semimajors, semimajors)
    assert binary.metallicity == metallicity

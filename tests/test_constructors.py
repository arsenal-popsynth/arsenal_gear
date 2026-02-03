"""
test_imf
==========

This file contains unit tests for the stellar population constructors.
"""

import astropy.units as u
import numpy as np
from numpy.testing import assert_array_equal

import arsenal_gear


def test_star():
    """Test the SSP constructor for single stars."""
    mass = u.Msun * np.ones(100)
    metals = 0.1
    tform = u.Myr
    star = arsenal_gear.population.SSP(mass=mass, metals=metals, tform=tform)

    assert_array_equal(star["mass"], mass)
    assert_array_equal(star["metals"], metals)
    assert_array_equal(star["tform"], tform)


def test_binary():
    """Test the BSP constructor for binary stars."""
    mass = u.Msun * np.ones(100)
    metals = 0.1
    tform = u.Myr
    star1 = arsenal_gear.population.SSP(mass=mass, metals=metals, tform=tform)
    star2 = arsenal_gear.population.SSP(mass=2 * mass, metals=metals, tform=tform)
    period = u.d * np.ones(100)
    eccentricity = 0.5 * np.ones(100)
    binary = arsenal_gear.population.BSP(
        primary=star1, secondary=star2, period=period, eccentricity=eccentricity
    )
    # Check the first star
    assert_array_equal(binary.primary["mass"], star1["mass"])
    assert_array_equal(binary.primary["metals"], star1["metals"])
    assert_array_equal(binary.primary["tform"], star1["tform"])

    # Check the second star
    assert_array_equal(binary.secondary["mass"], star2["mass"])
    assert_array_equal(binary.secondary["metals"], star2["metals"])
    assert_array_equal(binary.secondary["tform"], star2["tform"])
    # Check binary properties
    assert_array_equal(binary["period"], period)
    assert_array_equal(binary["eccentricity"], eccentricity)

"""
binary_fraction
==========

This submodule contains all the code required to sample from a mass-dependent binary fraction.
"""

from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.stats import rv_continuous, rv_histogram

from arsenal_gear.population import StarPopulation

class binary_fraction:
    """
    This class is the superclass of all mass-dependent binary fraction
    This assumes that the binary fraction is a step function

    :param fraction: Binary fraction of the mass bins
    :type fraction: float
    :param mass_bins: Limit of the mass bins
    :type mass_bins: astropy mass unit
    :param stars: Potential primaries
    :type stars: StarPopulation
    """
    def __init__(self, fraction: float, mass_bins: Quantity["mass"], stars: StarPopulation, name=""):
        self.fraction:  float = fraction
        self.mass_bins: float = mass_bins.to(u.Msun).value
        self.masses:    float = stars["mass"].to(u.Msun).value
        assert(np.min(self.fraction) >= 0)
        assert(np.max(self.fraction) <= 1)
        assert(np.min(self.mass_bins) >= 0)
        super().__init__(a=self.fraction, b=self.mass_bins, c=self.masses, name=name)


    def sample(self) -> bool:
        """
        Determine which stars are primaries

        :return: List of masses of stars
        :rtype: Quantity["mass"]
        """
        binary = np.heaviside(self.masses, self.mass_bins)
        print(binary)
        return binary
    



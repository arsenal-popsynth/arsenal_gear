"""
binaries
==========

This submodule contains all the code required to sample from mass-dependent a
binary fraction and distributions of orbital parameters.
"""

from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.stats import rv_continuous, rv_histogram

from arsenal_gear.population import StarPopulation

class Fraction:
    """
    This class is the superclass of all mass-dependent binary fraction
    This assumes that the binary fraction is a step function

    :param fraction: Binary fraction of the mass bins, of length k
    :type fraction: float
    :param mass_bins: Limit of the mass bins, of length k-1
    :type mass_bins: astropy mass unit
    :param stars: Potential primaries
    :type stars: StarPopulation
    :param name: Name of the binary fraction function
    :type name: str
    """
    def __init__(self, fraction: float, mass_bins: Quantity["mass"], stars: StarPopulation, name=""):
        self.fraction:  float = fraction
        self.mass_bins: float = mass_bins.to(u.Msun).value
        self.stars:     float = stars["mass"].to(u.Msun).value
        assert(np.min(self.fraction) >= 0)
        assert(np.max(self.fraction) <= 1)
        assert(np.min(self.mass_bins) >= 0)
        #assert(len(self.fraction)-len(self.mass_bins) == 1)
        #super().__init__(a=self.fraction, b=self.mass_bins, c=self.stars, name=name) # Why can't this be defined?

    
class StepFraction(Fraction):
    """
    A simple step function binary fraction, with a binary fraction
    of 0 below the changeover mass and 1 above the changeover mass

    :param mass: Changeover mass between binary fractions of 0 and 1
    :type mass: astropy mass unit
    """
    def __init__(self, fraction: float, mass_bins: Quantity["mass"], stars: StarPopulation):
        self.name = "Step"
        assert(len(fraction) == 2)
        super().__init__(fraction, mass_bins, stars, name=self.name)

    def binary_fraction(self) -> np.float64:
        prob = np.piecewise(self.stars, [self.stars < self.mass_bins, self.stars >= self.mass_bins], self.fraction)
        return prob
    
    def sample(self) -> bool:
        """
        Determine which stars are primaries

        :return: Boolean array
        :rtype: bool
        """
        _sample = np.random.rand(len(self.stars))
        _binary = np.zeros(len(self.stars), dtype=bool)
        _select = np.where(_sample <= self.binary_fraction())
        _binary[_select] = np.ones(len(_select), dtype=bool)
        return _binary
    
    
class SMA(rv_continuous):
    """
    This class is the superclass of all semi-major axis distributions

    :param stars: Primaries
    :type stars: StarPopulation
    :param name: Name for the scipy.stats rv_continuous instance
    :type name: str
    """
    def __init__(self, stars: StarPopulation, name=""):
        self.stars: float = stars["mass"].to(u.Msun).value
        super().__init__(a=self.stars, name=name)
        
        
class StepSMA(SMA):
    """
    A simple step function distribution of semi-major axes, with a log-uniform
    distribution above and below the changeover semi-major axis

    :param sma: Changeover semi-major axis
    :type sma: astropy length unit
    :param ratio: Ratio of probabilities for close and wide binaries
    :type ratio: float
    """

    #def sample_mass(self, mtot: Quantity["mass"]) -> Quantity["mass"]:
    #    """
    #    Draw a sample from the IMF with target total mass

    #    :param mtot: Targer total mass of the sample
    #    :type mtot: Quantity["mass"]
    #    :return: List of masses of stars
    #    :rtype: Quantity["mass"]
    #    """
    #    N_samp = round((mtot/self.mean()).value)
    #    return self.sample(N_samp)

    #def sample(self, N: int) -> Quantity["mass"]:
    #    """
    #    Draw a sample from the IMF with a specific number of stars
    #
    #    :param N: Number of stars to draw
    #    :type N: int
    #    :return: List of masses of stars
    #    :rtype: Quantity["mass"]
    #    """
    #    return self.rvs(size=N)*u.Msun

    



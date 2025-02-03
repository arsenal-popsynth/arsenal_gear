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
from scipy.stats import rv_continuous, rv_histogram, loguniform

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

    :param min_a: Minimum semi-major axis
    :type min_a: astropy length unit
    :param max_a: Maximum semi-major axis
    :type max_a: astropy length unit
    :param stars: Primaries
    :type stars: StarPopulation
    :param name: Name for the scipy.stats rv_continuous instance
    :type name: str
    """
    def __init__(self, min_a: Quantity["length"], max_a: Quantity["length"], name=""):
        self.min_a: float = min_a.to(u.au).value
        self.max_a: float = max_a.to(u.au).value
        assert self.min_a > 0
        assert self.max_a > self.min_a
        super().__init__(a=self.min_a, b=self.max_a, name=name)
        
    def sample(self, N: int) -> Quantity["length"]:
        """
        :param N: Number of stars to draw
        :type N: int
        :return: List of semi-major axes of stars
        :rtype: Quantity["length"]
        """
        return self.rvs(size=N)*u.au
        
        
class LogUniformSMA(SMA):
    """
    A simple loguniform distribution of semi-major axes
    """
    
    def __init__(self, min_a: Quantity["length"], max_a: Quantity["length"]):
        self.name = "loguniform"
        super().__init__(min_a, max_a, name=self.name)

    def _pdf(self, x: np.float64) -> np.float64:
        rv = loguniform(self.min_a, self.max_a)
        return rv.pdf(x)
   
    def _ppf(self, x: np.float64) -> np.float64:
        rv = loguniform(self.min_a, self.max_a)
        return rv.ppf(x)

    

    



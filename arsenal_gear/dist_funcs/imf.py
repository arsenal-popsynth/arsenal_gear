"""
imf
==========

This submodule contains all the code required to sample from an IMF.
"""

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.stats import rv_continuous

__all__ = ["IMF", "Salpeter"]

class IMF(rv_continuous):
    """
    This class is the superclass of all Initial Mass Functions (IMFs)

    :param min_mass: Least massive star in the IMF
    :type min_mass: astropy mass unit
    :param max_mass: Most massive star in the IMF
    :type max_mass: astropy mass unit
    :param name: Name for the scipy.stats rv_continuous instance
    :type name: str
    :param seed: Random seed for sampling
    :type seed: None, int, numpy.random.Generator, or numpy.random.RandomState
    """
    def __init__(self, min_mass: Quantity["mass"], max_mass: Quantity["mass"],
                 name="", seed=None):
        self.min_mass: float = min_mass.to(u.Msun).value
        self.max_mass: float = max_mass.to(u.Msun).value
        super().__init__(a=self.min_mass, b=self.max_mass, name=name, seed=seed)

    def sample_mass(self, mtot: Quantity["mass"]) -> Quantity["mass"]:
        """
        Draw a sample from the IMF with target total mass

        :param mtot: Targer total mass of the sample
        :type mtot: Quantity["mass"]
        :return: List of masses of stars
        :rtype: Quantity["mass"]
        """
        N_samp = round((mtot/self.mean()).value)
        return self.sample(N_samp)

    def sample(self, N: int) -> Quantity["mass"]:
        """
        Draw a sample from the IMF with a specific number of stars

        :param N: Number of stars to draw
        :type N: int
        :return: List of masses of stars
        :rtype: Quantity["mass"]
        """
        return self.rvs(size=N)*u.Msun

class Salpeter(IMF):
    """
    A simple, classic Salpeter 1955 (slope 2.35) IMF.

    :param min_mass: Least massive star in the IMF
    :type min_mass: astropy mass unit
    :param max_mass: Most massive star in the IMF
    :type max_mass: astropy mass unit
    :param alpha: The IMF slope (note that the slope applied is -alpha, so this
                  should be positive)
    :type alpha: float
    :param seed: Random seed for sampling
    :type seed: None, int, numpy.random.Generator, or numpy.random.RandomState
    """
    def __init__(self, min_mass: Quantity["mass"], max_mass: Quantity["mass"],
                 alpha:float=2.35, seed=None):
        self.alpha = alpha
        self.name = "Salpeter"
        assert alpha >= 0
        super().__init__(min_mass, max_mass, self.name, seed=seed)

    def _pdf(self, x: np.float64, *args) -> np.float64:
        upper = np.power(self.max_mass, 1-self.alpha)
        lower = np.power(self.min_mass, 1-self.alpha)
        norm = (upper-lower)/(1-self.alpha)
        return np.power(x, -self.alpha)/norm

    def _ppf(self, q: np.float64, *args) -> np.float64:
        upper = np.power(self.max_mass, 1-self.alpha)
        lower = np.power(self.min_mass, 1-self.alpha)
        return (q*(upper-lower)+lower)**(1./(1-self.alpha))

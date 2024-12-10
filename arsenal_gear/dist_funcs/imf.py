"""
imf
==========

This submodule contains all the code required to sample from an IMF.
"""

from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.stats import rv_continuous 


class IMF(rv_continuous):
    """
    This class is the superclass of all Initial Mass Functions (IMFs)

    :param min_mass: Least massive star in the IMF
    :type min_mass: astropy mass unit
    :param max_mass: Most massive star in the IMF
    :type max_mass: astropy mass unit
    :param name: Name for the scipy.stats rv_continuous instance
    :type name: str
    """
    def __init__(self, min_mass: Quantity["mass"], max_mass: Quantity["mass"],
                 name=""):
        self.min_mass: float = min_mass.to(u.Msun).value
        self.max_mass: float = max_mass.to(u.Msun).value
        super().__init__(a=self.min_mass, b=self.max_mass, name=name)

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
    """
    def __init__(self, min_mass: Quantity["mass"], max_mass: Quantity["mass"],
                 alpha:float=2.35):
        self.alpha = alpha
        self.name = "Salpeter"
        assert(alpha >= 0)
        super().__init__(min_mass, max_mass, self.name)

    def _pdf(self, x: np.float64) -> np.float64:
        upper = np.power(self.max_mass, 1-self.alpha)
        lower = np.power(self.min_mass, 1-self.alpha)
        norm = (upper-lower)/(1-self.alpha)
        return np.power(x, -self.alpha)/norm

"""
imf
==========

This submodule contains all the code required to sample from an IMF.
"""

from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from . import ProbDistFunc


class IMF(ProbDistFunc):
    """
    This class is the superclass of all Initial Mass Functions (IMFs)

    :param min_mass: Least massive star in the IMF
    :type max_mass: astropy mass unit
    :param pdf_max: Most massive star in the IMF
    :type pdf_max: astropy mass unit
    :param normalized: Should we return a normalized version of the
                       probability when called, or just a function proportional to it?
    :type normalized: bool
    """
    def __init__(self, min_mass: Quantity["mass"], max_mass: Quantity["mass"],
                 normalized:bool=False) -> None:
        self.min_mass: Quantity["mass"] = min_mass
        self.max_mass: Quantity["mass"] = max_mass
        self.min_mass_msun: float = self.min_mass.to(u.Msun).value
        self.max_mass_msun: float = self.max_mass.to(u.Msun).value
        super().__init__(min_mass.value, max_mass.value, normalized)

    # This is identical to base method, but we'd like type annotation here.
    def prob(self, masses: Quantity["mass"]) -> np.float64:
        return np.ones(masses.shape)/self.norm

    def __call__(self, x: Quantity["mass"]) -> np.float64:
        """
        Return the probability for value(s) x, normalized if the PDF is initialized
        with normalized = True

        :param x: The values to sample P(x) for.
        :type x: np.float64
        :return: The probability for x, normalized if desired.
        :rtype: np.float64
        """
        xmsun = x.to(u.Msun).value
        p = self.prob(x)
        (lb,hb) = (xmsun <= self.min_mass_msun, xmsun >= self.max_mass_msun)
        select_range = np.logical_and(lb, hb)
        p[select_range] = 0
        return p

class Salpeter(IMF):
    """
    A simple, classic Salpeter 1955 (slope 2.35) IMF.

    :param min_mass: Least massive star in the IMF
    :type max_mass: astropy mass unit
    :param pdf_max: Most massive star in the IMF
    :type pdf_max: astropy mass unit
    :param normalized: Should we return a normalized version of the probability
                       when called, or just a function proportional to it?
    :type normalized: bool
    :param alpha: The IMF slope (note that the slope applied is -alpha, so this
                  should be positive)
    :type alpha: float
    """
    def __init__(self, min_mass: Quantity["mass"], max_mass: Quantity["mass"],
                 normalized:bool=False, alpha:float=2.35) -> None:
        self.alpha = alpha
        super().__init__(min_mass, max_mass, normalized)

    def normalization(self) -> float:
        upper = np.power(self.max_mass_msun, 1-self.alpha)
        lower = np.power(self.min_mass_msun, 1-self.alpha)
        return (upper-lower)/(1-self.alpha)

    def prob(self, masses: Quantity["mass"]) -> np.float64:
        return np.power(masses.to(u.Msun).value, -self.alpha)/self.norm

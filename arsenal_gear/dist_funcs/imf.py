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
    :param normalized: Should we return a normalized version of the probability when called, or just a function proportional to it?
    :type normalized: bool
    """
    def __init__(self, min_mass: Quantity["mass"], max_mass: Quantity["mass"], normalized:bool=False) -> None:
        super().__init__(min_mass, max_mass, normalized)

    def prob(self, masses: Quantity["mass"]) -> np.float64: # Identical to base method, but we'd like type annotation here.
        return np.ones(masses.shape)

class Salpeter(IMF):
    """
    A simple, classic Salpeter 1955 (slope 2.35) IMF.

    :param min_mass: Least massive star in the IMF
    :type max_mass: astropy mass unit
    :param pdf_max: Most massive star in the IMF
    :type pdf_max: astropy mass unit
    :param normalized: Should we return a normalized version of the probability when called, or just a function proportional to it?
    :type normalized: bool
    :param alpha: The IMF slope (note that the slope applied is -alpha, so this should be positive)
    :type alpha: float
    """
    def __init__(self, min_mass: Quantity["mass"], max_mass: Quantity["mass"], normalized:bool=False, alpha:float=2.35) -> None:
        self.alpha = alpha
        super().__init__(min_mass, max_mass, normalized)

    def normalization(self) -> float:
        return (np.power(self.max/u.Msun, 1-self.alpha)-np.power(self.min/u.Msun, 1-self.alpha))/(1-self.alpha)

    def prob(self, masses: Quantity["mass"]) -> np.float64:
        return np.power(masses/u.Msun, -self.alpha)

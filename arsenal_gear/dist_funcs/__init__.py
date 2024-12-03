"""
dist_funcs
==========

Probability distribution functions for anything under the sun that arsenal
needs.
"""


from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity


class ProbDistFunc():
    """
    This class is the superclass of all PDFs that arsenal will use.

    :param pdf_min: Lower limit for the PDF.
                    Probabilities for values below this will be zero
    :type pdf_min: float
    :param pdf_max: Upper limit for the PDF.
                    Probabilities for values above this will be zero
    :type pdf_max: float
    :param normalized: Should we return a normalized version of the
                       probability when called, or just a function proportional to it?
    :type normalized: bool
    """
    def __init__(self, pdf_min:float, pdf_max:float, normalized:bool=False) -> None:
        self.min = pdf_min
        self.max = pdf_max
        self.norm = self.normalization() if normalized else 1

    def normalization(self) -> float:
        """
        Return the normalization for this PDF.
        The default base class returns a uniform distribution.

        :return: The normalization for the PDF (it's integral from self.min to self.max)
        :rtype: float
        """
        return self.max - self.min
    
    def pdf(self, x: np.float64) -> np.float64:
        """
        Return the normalized probability for value(s) x.

        :param x: The values to sample P(x) for.
        :type x: np.float64
        :return: The normalized probability for x
        :rtype: np.float64
        """
        p = np.ones(x.shape)
        p[np.logical_or(x < self.min, x > self.max)] = 0
        return p/self.norm

    def cdf(self, x: np.float64) -> np.float64:
        """
        Returns the value of the cumulative probability distribution function.

        :param x: The values to sample P(x) for.
        :type x: np.float64
        :return: The CDF value at x
        :rtype: np.float64
        """
        p = (x-self.min)/(self.max-self.min)
        p[np.logical_or(x < self.min, x > self.max)] = 0
        return p

    def __call__(self, x: np.float64) -> np.float64:
        """
        Simply calls the pdf method.
        """
        return self.prob(x)

from . import imf

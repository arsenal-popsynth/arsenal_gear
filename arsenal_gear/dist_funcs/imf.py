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
        self.mavg: Quantity["mass"] = 0.5*(self.min_mass + self.max_mass)
        super().__init__(min_mass.value, max_mass.value, normalized)

    def pdf(self, masses: Quantity["mass"]) -> np.float64:
        """
        Return the normalized probability for value(s) x.

        :masses: The values to sample P(x) for.
        """
        if masses.isscalar:
            masses = np.array([masses.value])*masses.unit
        p = np.ones(masses.shape)
        (lb,hb) = (xmsun <= self.min_mass_msun, xmsun >= self.max_mass_msun)
        select_range = np.logical_and(lb, hb)
        p[select_range] = 0
        if len(p) == 1:
            return p[0]/self.norm
        else:
            return p/self.norm

    def cdf(self, masses: Quantity["mass"]) -> np.float64:
        """
        Return the cumulative distribution function at x.

        :masses: The values to sample CDF(x) for.
        """
        if masses.isscalar:
            masses = np.array([masses.value])*masses.unit
        p = (masses - self.min_mass)/(self.max_mass-self.min_mass)
        p = p.to(" ").value
        xmsun = masses.to(u.Msun).value
        (lb,hb) = (xmsun <= self.min_mass_msun, xmsun >= self.max_mass_msun)
        select_range = np.logical_and(lb, hb)
        p[select_range] = 0
        if len(p) == 1:
            return p[0]/self.norm
        else:
            return np.ones(masses.shape)/self.norm

    def inv_cdf(self, c: float) -> Quantity["mass"]:
        """
        For a given value in (0,1) this returns the inverse
        of the CDF for that value, corresponding to a point in the measure space.
        TODO(ltlancas) : assure that c is in (0,1)

        :param c: float between 0 and 1
        :type c: np.float64
        :return: m such that CDF(X) = c
        :rtype: Quantity["mass"]
        """
        return c*self.max_mass + (1-c)*self.min_mass

    def sample(self, mtot: Quantity["mass"]) -> Quantity["mass"]:
        """
        Draw a sample from the IMF with target total mass

        :param mtot: Targer total mass of the sample
        :type mtot: Quantity["mass"]
        :return: List of masses of stars
        :rtype: Quantity["mass"]
        """
        N_samp = round((mtot/self.mavg).to(" ").value)
        c = np.random.uniform(0, 1, N_samp)
        return self.inv_cdf(c)

    def __call__(self, x: Quantity["mass"]) -> np.float64:
        """
        Simply calls the pdf method.
        """
        return self.pdf(x)

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
        self.mavg = self.get_mavg()

    def normalization(self) -> float:
        upper = np.power(self.max_mass_msun, 1-self.alpha)
        lower = np.power(self.min_mass_msun, 1-self.alpha)
        return (upper-lower)/(1-self.alpha)

    def get_mavg(self) -> Quantity["mass"]:
        upper1 = np.power(self.max_mass_msun, 1-self.alpha)
        lower1 = np.power(self.min_mass_msun, 1-self.alpha)
        upper2 = np.power(self.max_mass_msun, 2-self.alpha)
        lower2 = np.power(self.min_mass_msun, 2-self.alpha)
        t1 = (1-self.alpha)/(2-self.alpha)
        t2 = (upper2-lower2)/(upper1-lower1)
        return t1*t2*u.Msun

    def pdf(self, masses: Quantity["mass"]) -> np.float64:
        if masses.isscalar:
            masses = np.array([masses.value])*masses.unit
        xmsun = masses.to(u.Msun).value
        p = np.power(xmsun, -self.alpha)
        (lb,hb) = (xmsun <= self.min_mass_msun, xmsun >= self.max_mass_msun)
        select_range = np.logical_and(lb, hb)
        p[select_range] = 0
        if len(p) == 1:
            return p[0]/self.norm
        else:
            return p/self.norm

    def cdf(self, masses: Quantity["mass"]) -> np.float64:
        if masses.isscalar:
            masses = np.array([masses.value])*masses.unit
        upper = np.power(self.max_mass_msun, 1-self.alpha)
        lower = np.power(self.min_mass_msun, 1-self.alpha)
        xmsun = masses.to(u.Msun).value
        p = np.power(xmsun, 1-self.alpha)
        p = (p-lower)/(upper-lower)
        (lb,hb) = (xmsun <= self.min_mass_msun, xmsun >= self.max_mass_msun)
        (p[lb],p[hb]) = (0., 1.)
        if len(p) == 1:
            return p[0]
        else:
            return p

    def inv_cdf(self, c: float) -> Quantity["mass"]:
        upper = self.max_mass_msun**(1.-self.alpha)
        lower = self.min_mass_msun**(1.-self.alpha)
        m = (c*(upper-lower)+lower)**(1./(1-self.alpha))
        return m*u.Msun

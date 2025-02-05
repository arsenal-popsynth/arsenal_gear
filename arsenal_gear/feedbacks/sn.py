"""
sn
=========

Functions for supernovae
"""
import astropy.units as u
import numpy as np
from astropy.units import Quantity

from ..population import StarPopulation


def lifetimes_Raiteri(stars: StarPopulation) -> Quantity["time"]:
    """
    Stellar lifetimes calculated from Raiteri+ 1996
    (https://ui.adsabs.harvard.edu/abs/1996A%26A...315..105R/abstract)
    Equation 3.

    :param stars: Stellar Population
    :type stars: StarPopulation
    :return: lifetime of each star in the population
    :rtype: Quantity["time"]
    """
    a0 = 10.13 + 0.07547*np.log10(stars['metals']) - 0.008084*np.power(np.log10(stars['metals']), 2)
    a1 = -4.424 - 0.7939*np.log10(stars['metals']) - 0.1187*np.power(np.log10(stars['metals']), 2)
    a2 = 1.262 + 0.3385*np.log10(stars['metals']) + 0.05417*np.power(np.log10(stars['metals']), 2)
    return np.power(10, a0 + a1*np.log10(stars['mass'].to(u.Msun).value) +
                    a2*np.power(np.log10(stars['mass'].to(u.Msun).value), 2))*u.yr

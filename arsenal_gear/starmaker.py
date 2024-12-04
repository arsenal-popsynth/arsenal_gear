"""
starmaker
============

This submodule contains all the methods for generating populations of stars or binaries.
"""

from typing import Type, Optional

import numpy as np
import astropy.units as u
from scipy.stats import sampling
from astropy.units import Quantity

from .population import StarPopulation
from .dist_funcs import imf

__all__ = ['make_singles']

def make_singles(imf: imf.IMF, tform: Quantity["time"], metals: float, N: Optional[int] = None, 
                 M: Optional[Quantity["mass"]]= None):
    """
    Build a population of stars with either a fixed number N or a total mass M.

    :param imf:  Initial Mass Function to draw from
    :type imf:  Initial Mass Function to draw from
    :param tform:  Formation time for the stars
    :type tform:  Quantity["time"]
    :param metals:  Initial metallicity for the stars
    :type metals:  float
    :param N:  Number of stars to spawn
    :type N:  int or None
    :param M:  Total mass of stars to spawn
    :type M:  Quantity["mass"] or None
    """
    if (N is None) and (M is None):
        print("You must either request a number of stars N or a total mass M")
        return None
    if M is not None:
        print("Not Implemented")
        return None
    else:
        sampler = sampling.NumericalInverseHermite(imf)
        masses = sampler.rvs(N)
        return StarPopulation(mass=masses*u.Msun, metals=metals*np.ones(masses.shape), 
                              tform=tform*np.ones(masses.shape))

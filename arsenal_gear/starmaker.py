"""
starmaker
============

This submodule contains all the methods for generating populations of stars or binaries.
"""

from typing import Type

import numpy as np
import astropy.units as u
from scipy.stats import sampling
from astropy.units import Quantity

from .population import StarPopulation
from .dist_funcs import imf 

__all__ = ['make_singles']

def make_singles(imf, tform: Quantity["time"], metals: float, N=None, M=None):
    """
    Build a population of stars with either a fixed number N or a total mass M.
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

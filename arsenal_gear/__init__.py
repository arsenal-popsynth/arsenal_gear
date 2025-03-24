"""
arsenal-gear
============

A lightweight population synthesis code with an emphasis on the quantities
relevant for stellar feedback from massive stars.
"""
import time

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from . import dist_funcs, feedbacks, population, stellar_evolution

__version__ = '0.0.1'
__all__ = ['population', 'dist_funcs', 'feedbacks', 'stellar_evolution']

class StellarPopulation():
    """
    This class will act as the primary API for aresenal
    Ideally it will take an input yaml file or accept default values for
    certain parameters and return a population of stars with
    pre-calculated parameters for the feedback and radiaiton.
    """

    def __init__(self, verbose:bool=False, discrete:bool=True) -> None:
        self.Mtot = 1e6*u.Msun
        self.imf = dist_funcs.imf.Salpeter(0.08*u.Msun, 100*u.Msun, alpha=2.3)
        self.discrete = discrete
        # relative to solar
        self.metallicity = 1.0
        # generate masses
        if self.discrete:
            start_samp = time.time()
            self.masses = self.imf.sample_mass(self.Mtot)
            end_samp = time.time()
            if verbose:
                print("Time to sample masses: ", end_samp - start_samp)
        self.tmin = 0.0*u.Myr
        self.tmax = 40.0*u.Myr

        self.metallicity = 0.0
        self.iso = stellar_evolution.isochrone.MIST(self.metallicity,verbose=verbose)

    def get_NSN(self, t:Quantity["time"]) -> int:
        """
        Return the number of supernovae that have gone off by time t
        """
        Mmax = self.iso.get_Mmax(t)
        #if self.discrete:
        #    res = [len(np.where(self.masses.value >= max(8, mmi.value))[0]) for mmi in Mmax]
        #    return np.array(res)
        # fraction of stars that have exploded
        fexp_8 = 1 - self.imf.cdf(8*u.Msun)
        fexp = 1 - self.imf.cdf(Mmax)
        fexp = fexp*(fexp < fexp_8) + fexp_8*(fexp >= fexp_8)
        return (fexp*self.Mtot/self.imf.mean()).value

    def __call__(self, N:int) -> population.StarPopulation:
        """
        Return a
        """

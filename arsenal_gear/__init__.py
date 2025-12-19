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
__all__ = ['population', 'dist_funcs', 'feedbacks', 'stellar_evolution', 'StellarPopulation']

class StellarPopulation():
    """
    This class will act as the primary API for aresenal
    Ideally it will take an input yaml file or accept default values for
    certain parameters and return a population of stars with
    pre-calculated parameters for the feedback and radiaiton.
    """

    def __init__(self, **kwargs) -> None:
        # unpack kwarg parameters
        # total mass of the population
        self.Mtot = kwargs.get("Mtot",1e6*u.Msun)
        self.verbose = kwargs.get("verbose", False)
        self.discrete = kwargs.get("discrete", True)

        self.imf = dist_funcs.imf.Salpeter(0.08*u.Msun, 100*u.Msun, alpha=2.3)
        # expected number of stars
        self.Nstar = (self.Mtot/self.imf.mean()).value
        # log10(Z/Zsun)
        self.metallicity = kwargs.get("metallicity", 0.0)
        # generate masses
        if self.discrete:
            start_samp = time.time()
            self.masses = self.imf.sample_mass(self.Mtot)
            end_samp = time.time()
            if self.verbose:
                print("Time to sample masses: ", end_samp - start_samp)
        self.tmin = 0.0*u.Myr
        self.tmax = 40.0*u.Myr

        # initialize the isochrone
        self.iso = stellar_evolution.isochrone.MIST(**kwargs)

    def nsn(self, t:Quantity["time"]) -> int:
        """
        Return the number of supernovae that have gone off by time t
        """
        Mmax = self.iso.mmax(t)
        #if self.discrete:
        #    res = [len(np.where(self.masses.value >= max(8, mmi.value))[0]) for mmi in Mmax]
        #    return np.array(res)
        # fraction of stars that have exploded
        fexp_8 = 1 - self.imf.cdf(8*u.Msun)
        fexp = 1 - self.imf.cdf(Mmax)
        fexp = fexp*(fexp < fexp_8) + fexp_8*(fexp >= fexp_8)
        return fexp*self.Nstar

    def ndotsn(self, t:Quantity["time"]) -> int:
        """
        Return the rate of supernovae at time t: the derivative of nsn
        """
        Mmax = self.iso.mmax(t)
        Mmaxdot = self.iso.mmaxdot(t)
        return -self.imf.pdf(Mmax)*Mmaxdot*(Mmax.value>8)*self.Nstar

    def lbol(self, t:Quantity["time"]) -> Quantity["power"]:
        """
        Returns the bolometric luminosity of the population at time t
        """
        return np.sum(self.lbol_iso(t))

    def lbol_iso(self, t:Quantity["time"]) -> Quantity["power"]:
        """
        Returns the bolometric luminosity of each star in the population at time t
        """
        return self.iso.lbol(self.masses, t)

    def teff_iso(self, t:Quantity["time"]) -> Quantity["temperature"]:
        """
        Returns the effective temperature of each star in the population at time t
        """
        return self.iso.teff(self.masses, t)

    def __call__(self, N:int) -> population.StarPopulation:
        """
        Return a
        """

"""
arsenal-gear
============

A lightweight population synthesis code with an emphasis on the quantities
relevant for stellar feedback from massive stars.
"""

import time
import warnings

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from . import dist_funcs, element_yields, feedbacks, population, stellar_evolution

__version__ = "0.0.1"
__all__ = ["population", "dist_funcs", "feedbacks", "stellar_evolution"]


class StellarPopulation:
    """
    This class will act as the primary API for aresenal
    Ideally it will take an input yaml file or accept default values for
    certain parameters and return a population of stars with
    pre-calculated parameters for the feedback and radiaiton.
    """

    def __init__(self, **kwargs) -> None:
        # unpack kwarg parameters
        # total mass of the population
        self.Mtot = kwargs.get("Mtot", 1e6 * u.Msun)
        self.verbose = kwargs.get("verbose", False)
        self.discrete = kwargs.get("discrete", True)

        self.imf = dist_funcs.imf.Salpeter(0.08 * u.Msun, 100 * u.Msun, alpha=2.3)
        # expected number of stars
        self.Nstar = (self.Mtot / self.imf.mean()).value
        # log10(Z/Zsun)
        self.metallicity = kwargs.get("metallicity", 0.0) * u.dimensionless_unscaled
        self.rot = kwargs.get("rot", 0.0) * u.km * u.s**-1
        # generate masses
        if self.discrete:
            start_samp = time.time()
            self.masses = self.imf.sample_mass(self.Mtot)
            end_samp = time.time()
            if self.verbose:
                print("Time to sample masses: ", end_samp - start_samp)
        self.tmin = 0.0 * u.Myr
        self.tmax = 40.0 * u.Myr

        # initialize the isochrone
        self.iso = stellar_evolution.isochrone.MIST(**kwargs)

        # Yield model
        self.yield_model = kwargs.get("yield_model", None)

        if self.yield_model is None:
            from .element_yields import LimongiChieffi2018

            self.yield_model = LimongiChieffi2018(model="R")

    def __getitem__(self, key):
        """
        Allow dict-like access: pop['mass'], pop['metals'], etc.
        Introduced to support lifetimes_Raitieri function
        """
        if key == "mass":
            return self.masses
        if key == "metals":
            # Assuming uniform metallicity; return array matching mass shape
            return (
                np.full_like(self.masses.value, self.metallicity)
                * self.metallicity.unit
            )
        if key == "rot":
            # Assuming uniform rotation; return array matching mass shape
            return np.full_like(self.masses.value, 0.0) * u.km * u.s**-1
        if key in self.__dict__:
            return self.__dict__[key]

        raise KeyError(f"'{key}' not found in StellarPopulation")

    def nsn(self, t: Quantity["time"]) -> int:
        """
        Return the number of supernovae that have gone off by time t
        """
        Mmax = self.iso.mmax(t)
        # if self.discrete:
        #    res = [len(np.where(self.masses.value >= max(8, mmi.value))[0]) for mmi in Mmax]
        #    return np.array(res)
        # fraction of stars that have exploded
        fexp_8 = 1 - self.imf.cdf(8 * u.Msun)
        fexp = 1 - self.imf.cdf(Mmax)
        fexp = fexp * (fexp < fexp_8) + fexp_8 * (fexp >= fexp_8)
        return fexp * self.Nstar

    def ndotsn(self, t: Quantity["time"]) -> int:
        """
        Return the rate of supernovae at time t: the derivative of nsn
        """
        Mmax = self.iso.mmax(t)
        Mmaxdot = self.iso.mmaxdot(t)
        return -self.imf.pdf(Mmax) * Mmaxdot * (Mmax.value > 8) * self.Nstar

    def __call__(self, N: int) -> population.StarPopulation:
        """
        Return a
        """

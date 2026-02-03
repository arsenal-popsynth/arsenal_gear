"""
arsenal-gear
============

A lightweight population synthesis code with an emphasis on the quantities
relevant for stellar feedback from massive stars.
"""

import time
from abc import ABC, abstractmethod

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.integrate import trapezoid as trapz

from . import dist_funcs, element_yields, feedbacks, population, stellar_evolution
from .stellar_evolution.se_data_structures import Isochrone
from .utils import masked_power

__version__ = "0.0.1"
__all__ = [
    "population",
    "dist_funcs",
    "feedbacks",
    "stellar_evolution",
    "StellarPopulation",
    "DiscreteStellarPopulation",
]


class AbstractStellarPopulation(ABC):
    """
    This is the primary API for arsenal-gear.  Feed it what's needed to generate a population of stars,
    and it will allow you to evolve that population forward in time and query for various properties.

    :var t: Description
    :vartype t: the
    """

    def __init__(self, IMF, Mtot, metallicity, fbin) -> None:
        self.Mtot = Mtot
        self.imf = IMF
        self.metallicity = metallicity
        self.fbin = fbin

    @abstractmethod
    def sn_count(self, t0: Quantity["time"], t1: Quantity["time"]) -> int:
        """
        Return the number of supernovae that have gone off between time t0 and t1
        """


class DiscreteStellarPopulation(AbstractStellarPopulation):
    """
    A stellar population where individual stars are sampled from the IMF
    """

    def __init__(self, IMF, Mtot, metallicity, fbin) -> None:
        super().__init__(IMF, Mtot, metallicity, fbin)
        # sample stars until we reach Mtot
        self.masses = self.imf.sample_mass(self.Mtot)
        self.SingleStarPop = population.SSP(
            mass=self.masses, metals=self.metallicity, rot=0 * u.km / u.s
        )

    def sn_count(self, t0: Quantity["time"], t1: Quantity["time"]) -> int:
        explodability = feedbacks.sn.explodable_range(8 * u.Msun, 40 * u.Msun)
        return feedbacks.sn.get_sn_count(
            self.SingleStarPop, t0, t1, feedbacks.sn.lifetimes_Raiteri, explodability
        )


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
        self.metallicity = kwargs.get("metallicity", 0.0)
        # generate masses
        if self.discrete:
            start_samp = time.time()
            self.masses = self.imf.sample_mass(self.Mtot)
            end_samp = time.time()
            if self.verbose:
                print("Time to sample masses: ", end_samp - start_samp)
            self.Nstar = len(self.masses)
        self.tmin = 0.0 * u.Myr
        self.tmax = 40.0 * u.Myr

        # initialize the isochrone system
        self.iso = stellar_evolution.isochrone.IsochroneInterpolator(**kwargs)

    def _integrate_pop(self, iso: Isochrone, q: str) -> np.float64:
        """
        Integrate a given quantity over a population given an isochrone
        """
        return trapz(
            iso.qs[q] * self.imf.pdf(iso.qs[iso.mini_name]), iso.qs[iso.mini_name]
        )

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

    def lbol(self, t: Quantity["time"]) -> Quantity["power"]:
        """
        Returns the bolometric luminosity of the population at time t
        """
        if self.discrete:
            if np.isscalar(t):
                return np.sum(self.lbol_iso(t))
            else:
                return np.array([np.sum(self.lbol_iso(ti)).value for ti in t]) * u.Lsun
        else:
            if np.isscalar(t):
                iso = self.iso.construct_isochrone(t)
                iso.qs["L_bol"] = masked_power(10, iso.qs[self.iso.llbol_label])
                return self._integrate_pop(iso, "L_bol") * u.Lsun
            else:
                res = []
                for ti in t:
                    iso = self.iso.construct_isochrone(ti)
                    iso.qs["L_bol"] = masked_power(10, iso.qs[self.iso.llbol_label])
                    res.append(self._integrate_pop(iso, "L_bol"))
                return np.array(res) * u.Lsun

    def teff(self, t: Quantity["time"]) -> Quantity["power"]:
        """
        Returns the bolometric luminosity weighted
        effective temperature of the population at time t
        """
        if np.isscalar(t.value):
            teffs = self.teff_iso(t)
            lbols = self.lbol_iso(t)
            return np.sum(teffs * lbols) / np.sum(lbols)
        else:
            teff_arr = []
            for ti in t:
                teffs = self.teff_iso(ti)
                lbols = self.lbol_iso(ti)
                teff_arr.append((np.sum(teffs * lbols) / np.sum(lbols)).value)
            return np.array(teff_arr) * u.K

    def lbol_iso(self, t: Quantity["time"]) -> Quantity["power"]:
        """
        Returns the bolometric luminosity of each star in the population at time t
        """
        Lbols = self.iso.lbol(self.masses, t)
        return Lbols[np.logical_not(Lbols.mask)]

    def teff_iso(self, t: Quantity["time"]) -> Quantity["temperature"]:
        """
        Returns the effective temperature of each star in the population at time t
        """
        Teffs = self.iso.teff(self.masses, t)
        return Teffs[np.logical_not(Teffs.mask)]

    def __call__(self, N: int) -> population.SSP:
        """
        Return a
        """

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
from scipy.integrate import trapezoid as trapz

from . import formation
from .formation import SinglePop, BinaryPop
from . import element_yields, feedbacks, stellar_evolution


from .stellar_evolution.se_data_structures import Isochrone

__version__ = "0.0.1"
__all__ = [
    "element_yields",
    "feedbacks",
    "stellar_evolution",
    "SynthPop",
]


class SynthPop:
    """
    This class will act as the primary API for aresenal
    Ideally it will take an input yaml file or accept default values for
    certain parameters and return a population of stars with
    pre-calculated parameters for the feedback and radiaiton.
    """

    def __init__(self, **kwargs) -> None:
        # TO DO: unpack the kwargs into relevant sub-sets of arguments
        #       and change the below initialization so that it checks if users
        #       have already initialized the relevant submodules
        # total mass of the population
        self.form = formation.Formation(**kwargs)
        # initialize the isochrone interface/interpolator
        self.evol = stellar_evolution.Evolution(**kwargs)

    def _integrate_subpop(self, iso: Isochrone, pop: SinglePop, q: str) -> np.float64:
        """
        Integrate a given quantity over a single subpopulation given an isochrone
        """
        return pop.Nstar * trapz(iso.qs[q] * pop.imf.pdf(iso.mini), iso.mini)

    def nsn(self, t: Quantity["time"]) -> int:
        """
        Return the number of supernovae that have gone off by time t
        """
        Mmax = self.evol.se.mmax(t)
        NSN = np.zeros_like(Mmax.value)
        for pop in self.form.subpops:
            if pop.discrete:
                ms = pop.masses.value
                N_exp = [len(np.where((ms >= 8) & (ms >= m))[0]) for m in Mmax.value]
                NSN += np.array(N_exp)
            else:
                fexp_8 = 1 - pop.imf.cdf(8 * u.Msun)
                fexp = 1 - pop.imf.cdf(Mmax)
                fexp = fexp * (fexp < fexp_8) + fexp_8 * (fexp >= fexp_8)
                NSN += fexp * pop.Nstar
        if np.isscalar(t.value):
            return NSN[0]
        else:
            return NSN

    def ndotsn(self, t: Quantity["time"]) -> int:
        """
        Return: the rate of supernovae at time t, the derivative of nsn
                in Myr^-1
        """
        Mmax = self.evol.se.mmax(t)
        Mmaxdot = self.evol.se.mmaxdot(t)
        NSN_dot = 0
        for pop in self.form.subpops:
            NSN_dot += -pop.imf.pdf(Mmax)/u.Msun * Mmaxdot * (Mmax.value > 8) * pop.Nstar
        return NSN_dot

    def lbol(self, t: Quantity["time"]) -> Quantity["power"]:
        """
        Returns the bolometric luminosity of the population at time t
        """
        if np.isscalar(t.value):
            lbol_tot = 0 * u.Lsun
        else:
            lbol_tot = np.zeros_like(t.value) * u.Lsun

        for pop in self.form.subpops:
            if pop.discrete:
                if np.isscalar(t.value):
                    lbol_tot += np.sum(self.lbol_iso(t, pop))
                else:
                    lbol_tot += np.array([np.sum(self.lbol_iso(ti, pop)).value for ti in t]) * u.Lsun
            else:
                if np.isscalar(t.value):
                    iso = self.evol.se.construct_isochrone(t)
                    iso.qs["L_bol"] = iso.lbol
                    lbol_tot += self._integrate_subpop(iso, pop, "L_bol") * u.Lsun
                else:
                    res = []
                    for ti in t:
                        iso = self.evol.se.construct_isochrone(ti)
                        iso.qs["L_bol"] = iso.lbol
                        res.append(self._integrate_subpop(iso, pop, "L_bol"))
                    lbol_tot += np.array(res) * u.Lsun
        return lbol_tot

    def teff(self, t: Quantity["time"]) -> Quantity["power"]:
        """
        Returns the bolometric luminosity weighted
        effective temperature of the population at time t
        """
        if np.isscalar(t.value):
            teff_weight_tot = 0 * u.K * u.Lsun
        else:
            teff_weight_tot = np.zeros_like(t.value) * u.K * u.Lsun
        
        for pop in self.form.subpops:
            if pop.discrete:
                if np.isscalar(t.value):
                    teffs = self.teff_iso(t, pop)
                    lbols = self.lbol_iso(t, pop)
                    teff_weight_tot += np.sum(teffs * lbols)
                else:
                    teff_arr = []
                    for ti in t:
                        teffs = self.teff_iso(ti, pop)
                        lbols = self.lbol_iso(ti, pop)
                        teff_arr.append((np.sum(teffs * lbols)).value)
                    teff_weight_tot += np.array(teff_arr) * u.K * u.Lsun
            else:
                if np.isscalar(t.value):
                    iso = self.evol.se.construct_isochrone(t)
                    iso.qs["L_bol*Teff"] = iso.lbol.value*iso.teff.value
                    teff_weight_tot += self._integrate_subpop(iso, pop, "L_bol*Teff") * u.K * u.Lsun
                else:
                    res = []
                    for ti in t:
                        iso = self.evol.se.construct_isochrone(ti)
                        iso.qs["L_bol*Teff"] = iso.lbol.value*iso.teff.value
                        res.append(self._integrate_subpop(iso, pop, "L_bol*Teff"))
                    teff_weight_tot += np.array(res) * u.K * u.Lsun
        lbol = self.lbol(t)
        return (teff_weight_tot / lbol).to(u.K)

    def lbol_iso(self, t: Quantity["time"], pop:SinglePop) -> Quantity["power"]:
        """
        Returns the bolometric luminosity of each star in the population at time t
        """
        Lbols = self.evol.se.lbol(pop.masses, t)
        return Lbols[np.logical_not(Lbols.mask)]

    def teff_iso(self, t: Quantity["time"], pop:SinglePop) -> Quantity["temperature"]:
        """
        Returns the effective temperature of each star in the population at time t
        """
        Teffs = self.evol.se.teff(pop.masses, t)
        return Teffs[np.logical_not(Teffs.mask)]

    @property
    def masses(self) -> Quantity["mass"]:
        """
        Return the masses of all stars in the population as a 1D array
        If the population is not discrete, mass from these populations is not included
        """
        masses = []
        for pop in self.form.subpops:
            if pop.discrete:
                masses.append(pop.masses)
        return np.concatenate(masses) * u.Msun
    
    @property
    def Mtot(self) -> Quantity["mass"]:
        """
        Return the total mass of the population
        """
        return self.form.Mtot

    # ltlancas: commented this out for now, I can't think of what calling the SynthPop
    #           object should do by default right now...
    #def __call__(self, N: int) -> SinglePop:
    #    """
    #    Return a
    #    """

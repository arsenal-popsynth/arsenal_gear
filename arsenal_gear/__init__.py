"""
arsenal-gear
============

A lightweight population synthesis code with an emphasis on the quantities
relevant for stellar feedback from massive stars.
"""
import astropy.units as u
import numpy as np
from astropy.units import Quantity

from . import dist_funcs, feedbacks, population, stellar_evolution

__version__ = '0.0.1'
__all__ = ['population', 'dist_funcs', 'feedbacks', 'stellar_evolution']

class StarMaker():
    """
    This class will act as the primary API for aresenal
    Ideally it will take an input yaml file or accept default values for
    certain parameters and return a population of stars with
    pre-calculated parameters for the feedback and radiaiton.
    """

    def __init__(self, discrete:bool=True) -> None:
        self.Mtot = 1e6*u.Msun
        self.imf = dist_funcs.imf.Salpeter(0.08*u.Msun, 100*u.Msun, alpha=2.3)
        self.discrete = discrete
        # relative to solar
        self.metallicity = 1.0
        # generate masses
        if self.discrete:
            self.masses = self.imf.sample_mass(self.Mtot)
        self.tmin = 0.0*u.Myr
        self.tmax = 40.0*u.Myr

        mbase = "<path-to-mist>"
        isofname = mbase + "MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.0_full.iso"
        self.iso = stellar_evolution.isochrone.MIST(isofname)

    def get_NSN(self, t:Quantity["time"]) -> int:
        """
        Return the number of supernovae that have gone off by time t
        """
        Mmax = self.iso.get_Mmax(t)
        if self.discrete:
            return self.masses[self.masses >= max(8*u.Msun, Mmax)].size
        # fraction of stars that have exploded
        fexp_8 = 1 - self.imf.cdf(8*u.Msun)
        fexp = 1 - self.imf.cdf(Mmax)
        fexp = fexp*(fexp < fexp_8) + fexp_8*(fexp >= fexp_8)
        return (fexp*self.Mtot/self.imf.mean()).value

    def __call__(self, N:int) -> population.StarPopulation:
        """
        Return a
        """

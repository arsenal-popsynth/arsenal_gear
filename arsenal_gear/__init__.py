"""
arsenal-gear
============

A lightweight population synthesis code with an emphasis on the quantities
relevant for stellar feedback from massive stars.
"""

from . import dist_funcs, feedbacks, population, stellar_evolution

import astropy.units as u
from astropy.units import Quantity

__version__ = '0.0.1'
__all__ = ['population', 'dist_funcs', 'feedbacks']

class StellarPopulation():
    """
    This class will act as the primary API for aresenal
    Ideally it will take an input yaml file or accept default values for 
    certain parameters and return a population of stars with 
    pre-calculated parameters for the feedback and radiaiton.
    """

    def __init__(self) -> None:
        self.Mtot = 1e6*u.Msun
        self.imf = dist_funcs.imf.Salpeter(0.08*u.Msun, 100*u.Msun)
        # relative to solar
        self.metallicity = 1.0
        # generate masses
        self.masses = self.imf.sample(self.Mtot)
        self.tmin = 0.0*u.Myr
        self.tmax = 40.0*u.Myr

        mbase = "/Users/lachlanlancaster/Documents/Columbia/Research/pop_synth/wind_pop_synthesis/mist_iso_theory/MIST_v1.2_vvcrit0.0_full_isos/"
        isofname = mbase + "MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.0_full.iso"
        self.iso = stellar_evolution.isochrone.MIST(isofname)

    def __call__(self, N:int) -> population.StarPopulation:
        """
        Return a
        """
        pass
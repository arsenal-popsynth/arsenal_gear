"""
form_data_structures.py
=======================

Specifies containers for the basic properties of a stellar population
at the time of their formation. These classes are basically parsing
of input parameters and are used to construct the real stellar populations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from astropy.units import Quantity

if TYPE_CHECKING:
    from .dist_funcs.imf import IMF
    from .dist_funcs.binaries import MassRatio, Period, Semimajor

@dataclass
class SinglePop:
    """
    Data class for storing a population of single stars,
    BinaryPop extends this for binary star populations.
    """
    Mtot: Quantity["mass"]
    Nstar: int | float
    metallicity: Quantity["dimensionless"]
    imf: IMF
    mmin: Quantity["mass"]
    mmax: Quantity["mass"]
    discrete: bool
    # set to None if not discrete, otherwise a list of masses
    masses: Quantity["mass"]

@dataclass
class BinaryPop(SinglePop):
    """
    Data class for storing a binary star population.
    """
    q_dist: MassRatio
    period_dist: Period
    a_dist: Semimajor
    # set to None if not discrete
    # otherwise a list of mass ratios, periods, and semimajor axes
    mrats: Quantity["dimensionless"]
    periods: Quantity["time"]
    semimajors: Quantity["length"]
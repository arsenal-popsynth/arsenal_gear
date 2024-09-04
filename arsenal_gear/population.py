"""
population
==========

This submodule defines the classes for "populations", collections of single or binary stars.
"""

from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity

__all__ = ['population']

class StarPopulation(dict):
    """This class is used to instantiate populations of individual stars.

    :param mass: Initial mass of the stars
    :type mass: astropy mass unit
    :param metals: Metallicity (mass fraction) of the stars
    :type metals: float
    :param tform: Formation time of the stars, defaults to zero
    :type tform: astropy time unit
    """
    def __init__(self, mass: Quantity["mass"], metals: np.float64,
                 tform: Quantity["time"]=u.Quantity(0, u.Myr)) -> None:
        self['mass'] = mass
        self['metals'] = metals
        self['tform'] = tform

class BinaryPopulation(dict):
    """This class is used to instantiate a population of binary star systems.
    It is essentially a wrapper for a pair of
    :class:`arsenal_gear.population.StarPopulation` objects.

    :param primary: the primary stars in the binary pairs
    :type primary: class:arsenal_gear.population.StarPopulation
    :param secondary: the secondary stars in the binary pairs
    :type secondary: class:arsenal_gear.population.StarPopulation
    :param semimajor: the semimajor axis of the binary orbits
    :type semimajor: astropy length unit
    :param eccentricity: the eccentricity of the binary orbits
    :type eccentricity: float
    """
    def __init__(self, primary: StarPopulation, secondary: StarPopulation, semimajor:
                 Quantity["length"], eccentricity: np.float64) -> None:
        self.primary = primary
        self.secondary = secondary
        self['semimajor'] = semimajor
        self['eccentricity'] = eccentricity

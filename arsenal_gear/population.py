"""
population
==========

This submodule defines the classes for "populations", collections of single or binary stars.
"""

import numpy as np
from astropy.units import Quantity

__all__ = ["SSP", "BSP"]


class SSP(dict):
    """This class is used to instantiate populations of individual stars.

    :param mass: Initial mass of the stars
    :type mass: astropy mass unit
    :param metals: Metallicity (mass fraction) of the stars
    :type metals: float
    :param tform: Formation time of the stars, defaults to zero
    :type tform: astropy time unit
    """

    def __init__(
        self,
        mass: Quantity["mass"],
        metals: np.float64,
        tform: Quantity["time"] = 0,
    ) -> None:
        self["mass"] = mass
        self["metals"] = metals
        self["tform"] = tform


class BSP(dict):
    """This class is used to instantiate a population of binary star systems.
    It is essentially a wrapper for a pair of
    :class:`arsenal_gear.population.SSP` objects.

    :param primary: the primary stars in the binary pairs
    :type primary: class:arsenal_gear.population.SSP
    :param secondary: the secondary stars in the binary pairs
    :type secondary: class:arsenal_gear.population.SSP
    :param period: the orbital period of the binary orbits
    :type period: astropy length time
    :param eccentricity: the eccentricity of the binary orbits
    :type eccentricity: float
    """

    def __init__(
        self,
        primary: SSP,
        secondary: SSP,
        period: Quantity["time"],
        eccentricity: np.float64,
    ) -> None:
        self.primary = primary
        self.secondary = secondary
        self["period"] = period
        self["eccentricity"] = eccentricity

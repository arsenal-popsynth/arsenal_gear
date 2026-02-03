"""
population
==========

This submodule defines the classes for "populations", collections of single or binary stars.
"""

import astropy.units as u
import numpy as np
from astropy.units import Quantity

__all__ = ["SSP", "BSP"]


class SSP(Quantity):
    _unit = None
    metals = 0
    tform = 0 * u.yr
    """This class is used to instantiate populations of individual stars.  This
    "simple stellar population" has every star with equal metallicities and
    formation times.

    :param mass: the masses of the stars in the population
    :type mass: astropy mass unit
    :param metals: the metallicities of the stars in the population
    :type metals: float
    :param tform: the formation times of the stars in the population
    :type tform: astropy time unit
    """

    def __new__(
        subtype, mass: Quantity["mass"], metals=np.float64, tform=Quantity["time"]
    ):
        obj = Quantity(mass).view(subtype)
        obj._unit = mass.unit
        obj.metals = metals
        obj.tform = tform
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._unit = getattr(obj, "_unit", None)
        self.metals = getattr(obj, "metals", None)
        self.tform = getattr(obj, "tform", None)


class BSP(dict):
    """This class is used to instantiate a population of binary star systems.
    It is essentially a wrapper for a pair of
    :class:`arsenal_gear.population.SSP` objects.

    :param primary: the primary stars in the binary pairs
    :type primary: class:arsenal_gear.population.SSP
    :param secondary: the secondary stars in the binary pairs
    :type secondary: class:arsenal_gear.population.SSP
    :param period: the orbital period of the binary orbits
    :type period: astropy time unit
    :param eccentricity: the eccentricity of the binary orbits
    :type eccentricity: float
    """

    def __init__(
        self,
        primary: SSP,
        secondary: SSP,
        period: Quantity["time"],
        eccentricity: np.float64,
    ):
        self.primary = primary
        self.secondary = secondary
        self["period"] = period
        self["eccentricity"] = eccentricity

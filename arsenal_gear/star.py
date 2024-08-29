"""
star
====

This submodule defines the class for individual stars.
"""

from typing import Type
import astropy.units as u
from astropy.units import Quantity

__all__ = ['star']

class Star():
    """This class is used to instantiate individual stars.

    :param mass: Initial mass of the star
    :type mass: astropy mass unit
    :param metals: Metallicity (mass fraction) of the star
    :type metals: float
    :param age: Age of the star, defaults to zero
    :type age: astropy time unit
    """
    def __init__(self, mass: Quantity["mass"], metals: float, age: Quantity["time"]=u.Quantity(0, u.Myr)) -> None:
        self.mass = mass
        self.metals = metals
        self.age = age


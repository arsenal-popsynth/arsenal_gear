"""
binary
======

This submodule defines the class for binary pairs of stars.
"""

from .star import Star

from typing import Type
import astropy.units as u
from astropy.units import Quantity

__all__ = ['binary']

class Binary():
    """This class is used to instantiate individual binary star systems.  It is
    essentially a wrapper for a pair of :class:`arsenal_gear.star.Star` objects.

    :param star_1: the first star in the binary pair
    :type star_1: class:arsenal_gear.star.Star
    :param star_2: the second star in the binary pair
    :type star_2: class:arsenal_gear.star.Star
    :param radius: the separation of the pair 
    :type radius: astropy length unit
    :param eccentricity: the eccentricity of the binary orbit
    :type eccentricity: float
    """
    def __init__(self, star_1: Type[Star], star_2: Type[Star], radius:Quantity["length"], eccentricity:float) -> None:
        self.star_1 = star_1
        self.star_2 = star_2
        self.radius = radius
        self.eccentricity = eccentricity


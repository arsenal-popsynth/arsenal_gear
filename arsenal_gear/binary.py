"""
binary
====

This submodule defines the class for binary pairs of stars.
"""

from .star import Star

from typing import Type

__all__ = ['binary']

class Binary():
    def __init__(self, star_1: Type[Star], star_2: Type[Star], radius:float, eccentricity:float):
        self.star_1 = star_1
        self.star_2 = star_2
        self.radius = radius
        self.eccentricity = eccentricity


"""
arsenal-gear
============

A lightweight population synthesis code with an emphasis on the quantities
relevant for stellar feedback from massive stars.
"""

import time
from abc import ABC, abstractmethod

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.integrate import trapezoid as trapz

from . import dist_funcs, element_yields, feedbacks, population, stellar_evolution
from .stellar_evolution.se_data_structures import Isochrone
from .utils import masked_power

__version__ = "0.0.1"
__all__ = [
    "population",
    "dist_funcs",
    "feedbacks",
    "stellar_evolution",
    "FormationContext",
    "EvolutionContext",
]


class FormationContext:
    def __init__(self, **kwargs) -> None:
        self.imf = kwargs.get("imf", None)
        self.mass = kwargs.get("mass", None)
        self.metals = kwargs.get("metals", 0.012)  # Default to solar metallicity
        self.tform = kwargs.get(
            "tform", 0 * u.Myr
        )  # Default to formation time of 0 Myr
        self.N = kwargs.get("N", None)
        self.check_context()

    def check_context(self):
        if self.mass is None and self.N is None:
            raise ValueError(
                "Must provide either mass (total mass) or N (number of stars)"
            )

    def generate_population(self):
        return population.SSP(
            mass=self.imf.sample_mass(self.mass), metals=self.metals, tform=self.tform
        )


class EvolutionContext:
    def __init__(self, **kwargs) -> None:
        pass

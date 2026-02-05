"""
arsenal-gear
============

A lightweight population synthesis code with an emphasis on the quantities
relevant for stellar feedback from massive stars.
"""

import time
import warnings
from abc import ABC, abstractmethod
from typing import Type

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
    "element_yields",
    "stellar_evolution",
    "FormationContext",
    "EvolutionContext",
]


class FormationContext:
    """
    The FormationContext class encapsulates
    all the parameters, data, and models needed
    to generate a stellar population.  You must provide
    one of either a mass or a number of stars to generate.

    :param imf: An instance of an IMF class from dist_funcs.
    :type imf: dist_funcs.IMF
    :param mass: Total mass of the stellar population to generate.
    :type mass: astropy.units.Quantity, optional
    :param num: Number of stars in the stellar population to generate.
    :type num: int, optional
    """

    def __init__(self, imf: Type["IMF"], mass=0.0 * u.Msun, num=0, **kwargs) -> None:
        self.imf = imf
        self.mass = mass
        self.num = num
        self.metals = kwargs.get("metals", 0.012)  # Default to solar metallicity
        self.tform = kwargs.get(
            "tform", 0 * u.Myr
        )  # Default to formation time of 0 Myr
        self.check_context()

    def check_context(self):
        """
        Validate that the context has sufficient information to generate a population,
        and that there are no conflicting parameters.
        """
        if self.mass <= 0 and self.num <= 0:
            raise ValueError(
                "Must provide positive, non-zero mass (total mass) or num (number of stars)"
            )
        if self.mass > 0 and self.num > 0:
            raise ValueError("Cannot provide both mass and num; choose one to specify")

    def generate_population(self):
        """
        Actually generate a stellar population based on the context parameters.
        :return: Generated stellar population
        :rtype: population.SSP
        """
        return population.SSP(
            mass=self.imf.sample_mass(self.mass), metals=self.metals, tform=self.tform
        )


class EvolutionContext:
    """
    The EvolutionContext class encapsulates all the parameters, data, and models
    needed to evolve a stellar population and calculate feedback yields.
    """

    def __init__(self, mechanisms: list) -> None:
        """
        Initialize the EvolutionContext with a list of feedback mechanisms.
        :param mechanisms: List of feedback mechanism instances
        :type mechanisms: list of feedbacks.FBMechanism
        """
        self.mechanisms = mechanisms

    def mass(
        self, starpop: Type["SSP"], t0: Quantity["time"], t1: Quantity["time"]
    ) -> Quantity["mass"]:
        """
        Calculate the total mass returned by all feedback mechanisms
        between times t0 and t1 for the given stellar population.
        :param starpop: Stellar population to evolve
        :type starpop: population.SSP
        :param t0: Start time for feedback calculation
        :type t0: astropy.units.Quantity["time"]
        :param t1: End time for feedback calculation
        :type t1: astropy.units.Quantity["time"]
        :return: Total mass returned by feedback mechanisms
        :rtype: astropy.units.Quantity["mass"]
        """
        return sum(mech.mass(starpop, t0, t1) for mech in self.mechanisms)

    def metals_total(
        self, starpop: Type["SSP"], t0: Quantity["time"], t1: Quantity["time"]
    ) -> Quantity["mass"]:
        """
        Calculate the metal mass returned by all feedback mechanisms
        between times t0 and t1 for the given stellar population.
        :param starpop: Stellar population to evolve
        :type starpop: population.SSP
        :param t0: Start time for feedback calculation
        :type t0: astropy.units.Quantity["time"]
        :param t1: End time for feedback calculation
        :type t1: astropy.units.Quantity["time"]
        :return: Metal mass returned by feedback mechanisms
        :rtype: astropy.units.Quantity["mass"]
        """
        return sum(mech.metals_total(starpop, t0, t1) for mech in self.mechanisms)

    def metals_species(
        self,
        starpop: Type["SSP"],
        species: str,
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> Quantity["mass"]:
        """
        Calculate the metal mass of a given species returned by all feedback mechanisms
        between times t0 and t1 for the given stellar population.
        :param starpop: Stellar population to evolve
        :type starpop: population.SSP
        :param species: Chemical species to calculate (e.g., 'Fe', 'O')
        :type species: str
        :param t0: Start time for feedback calculation
        :type t0: astropy.units.Quantity["time"]
        :param t1: End time for feedback calculation
        :type t1: astropy.units.Quantity["time"]
        :return: Species mass returned by feedback mechanisms
        :rtype: astropy.units.Quantity["mass"]
        """
        return sum(
            mech.metals_species(starpop, species, t0, t1) for mech in self.mechanisms
        )

    def count(
        self, starpop: Type["SSP"], t0: Quantity["time"], t1: Quantity["time"]
    ) -> int:
        """
        Calculate the count of a discrete FB mechanism (SNII, SN Ia, stellar
        merger, etc) occuring
        between times t0 and t1 for the given stellar population.
        :param starpop: Stellar population to evolve
        :type starpop: population.SSP
        :param t0: Start time for feedback calculation
        :type t0: astropy.units.Quantity["time"]
        :param t1: End time for feedback calculation
        :type t1: astropy.units.Quantity["time"]
        :return: Total count of events from feedback mechanisms
        :rtype: int
        """
        return sum(mech.count(starpop, t0, t1) for mech in self.mechanisms)

    def energy(
        self, starpop: Type["SSP"], t0: Quantity["time"], t1: Quantity["time"]
    ) -> Quantity["energy"]:
        """
        Calculate the energy returned by all feedback mechanisms
        between times t0 and t1 for the given stellar population.
        :param starpop: Stellar population to evolve
        :type starpop: population.SSP
        :param t0: Start time for feedback calculation
        :type t0: astropy.units.Quantity["time"]
        :param t1: End time for feedback calculation
        :type t1: astropy.units.Quantity["time"]
        :return: Total energy returned by feedback mechanisms
        :rtype: astropy.units.Quantity["energy"]
        """
        return sum(mech.energy(starpop, t0, t1) for mech in self.mechanisms)

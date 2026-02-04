"""
feedbacks
=========

Functions for the returns from various feedback mechanisms (SN, winds, radiation, etc.)
"""

from abc import ABC, abstractmethod
from typing import Callable, Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from . import sn

__all__ = ["SNFeedbackMechanism", "sn"]


class MechanicalFBMechanism:
    """
    Abstract Base Class for Mechanical Feedback Mechanisms (SN, Winds, etc.)
    """

    @abstractmethod
    def mass(
        self,
        stars: Type["SSP"],
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> Quantity["mass"]:
        """
        Get the total mass yield from this mechanism between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: mass yield from this mechanism
        """

    @abstractmethod
    def metals_total(
        self,
        stars: Type["SSP"],
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> Quantity["mass"]:
        """
        Get the total metal yield from this mechanism between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: metal yield from this mechanism
        """

    @abstractmethod
    def metals_species(
        self,
        stars: Type["SSP"],
        species: str,
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> Quantity["mass"]:
        """
        Get the metal yield for a given species from this mechanism between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: species metal yield from this mechanism
        """


class SNFeedbackMechanism(MechanicalFBMechanism):
    """
    Supernova feedback mechanism
    """

    def __init__(
        self,
        energy_func: Callable[[Type["SSP"]], Quantity["energy"]],
        mass_func: Callable[[Type["SSP"]], Quantity["mass"]],
        metals_func: Callable[[Type["SSP"]], dict],
        explodability_func: Callable[
            [Type["SSP"], Quantity["time"], Quantity["time"]], Quantity["energy"]
        ],
    ) -> None:
        self.energy_ = energy_func
        self.mass_ = mass_func
        self.metals_ = metals_func
        self.explode = explodability_func

    def metals_species(self, stars, species, t0, t1):
        metals = self.metals_(stars[self.explode(stars, t0, t1)])
        if species not in metals:
            return 0.0 * u.Msun
        return np.sum(metals[species])

    def metals_total(self, stars, t0, t1):
        return np.sum(self.metals_(stars[self.explode(stars, t0, t1)])["metals"])

    def mass(self, stars, t0, t1):
        return np.sum(self.mass_(stars[self.explode(stars, t0, t1)]))

    def energy(self, stars, t0, t1):
        return np.sum(self.energy_(stars[self.explode(stars, t0, t1)]))

    def count(self, stars, t0, t1):
        return np.sum(self.explode(stars, t0, t1))

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

from . import analytic

__all__ = ["SNFeedbackMechanism", "analytic"]


class MechanicalFBMechanism(ABC):
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
        return 0 * u.Msun

    @abstractmethod
    def energy(
        self,
        stars: Type["SSP"],
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> Quantity["energy"]:
        """
        Get the total energy yield from this mechanism between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: energy yield from this mechanism
        """
        return 0 * u.erg

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
        return 0 * u.Msun

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
        return 0 * u.Msun

    @abstractmethod
    def count(
        self,
        stars: Type["SSP"],
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> int:
        """
        Get the number of events (for a discrete mechanism)
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: number of events from this mechanism
        :rtype: int
        """
        return 0


class SNFeedbackMechanism(MechanicalFBMechanism):
    """
    Supernova feedback mechanism.
    """

    def __init__(
        self,
        energy_func: Callable[[Type["SSP"]], Quantity["energy"]],
        mass_func: Callable[[Type["SSP"]], Quantity["mass"]],
        metals_func: Callable[[Type["SSP"]], dict],
        lifetime_func,
        explodability_func: Callable[[Type["SSP"]], bool],
    ) -> None:
        """
        Initialize the SN feedback mechanism.  The explodability function is used
        to select which stars from the SSP explode between times t0 and t1, and then the
        energy, mass, and metals functions are used to get the yields from those stars.

        :param energy_func: Function to get energy yield from an SSP
        :type energy_func: Callable[[SSP], Quantity["energy"]]
        :param mass_func: Function to get mass yield from an SSP
        :type mass_func: Callable[[SSP], Quantity["mass"]]
        :param metals_func: Function to get total and per-species metal yields from an SSP
        :type metals_func: Callable[[SSP], dict]
        :param lifetime_func: Function to get stellar lifetimes from an SSP
        :type lifetime_func: Callable[[SSP], Quantity["time"]] or IsochroneInterpolator
        :param explodability_func: Function to determine which stars explode
        :type explodability_func: Callable[[SSP], bool]
        """
        self.energy_ = energy_func
        if hasattr(mass_func, "ccsn_mass"):
            self.mass_ = mass_func.ccsn_mass
        else:
            self.mass_ = mass_func
        self.metals_ = metals_func
        if hasattr(lifetime_func, "mmax"):
            self.explode = lambda stars, t0, t1: np.bool(
                np.logical_and(
                    explodability_func(stars),
                    analytic.explodable_mass_range(
                        lifetime_func.mmax(t1), lifetime_func.mmax(t0)
                    )(stars),
                )
            )
        else:
            self.explode = lambda stars, t0, t1: np.bool(
                np.logical_and(
                    explodability_func(stars),
                    analytic.explodable_lifetime_range(t0, t1, lifetime_func)(stars),
                )
            )

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

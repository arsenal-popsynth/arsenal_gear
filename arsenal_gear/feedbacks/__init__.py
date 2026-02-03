"""
feedbacks
=========

Functions for the returns from various feedback mechanisms (SN, winds, radiation, etc.)
"""

from typing import Callable, Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from . import sn

__all__ = ["SNFeedbackMechanism", "sn"]


class SNFeedbackMechanism:
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
        self.energy_func = energy_func
        self.mass_func = mass_func
        self.metals_func = metals_func
        self.explodability_func = explodability_func

    def metals_species(
        self,
        stars: Type["SSP"],
        species: str,
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> int:
        """
        Get the total metal yield from supernovae between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param species: Chemical species to retrieve
        :type species: str
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: mass loss from supernovae
        """
        return np.sum(
            self.metals_func(stars)[species][self.explodability_func(stars, t0, t1)]
        )

    def metals_total(
        self,
        stars: Type["SSP"],
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> int:
        """
        Get the total metal yield from supernovae between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: mass loss from supernovae
        """
        return np.sum(
            self.metals_func(stars[self.explodability_func(stars, t0, t1)])["metals"]
        )

    def mass(
        self,
        stars: Type["SSP"],
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> int:
        """
        Get the number of supernovae that have gone off between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: mass loss from supernovae
        """
        return np.sum(self.mass_func(stars[self.explodability_func(stars, t0, t1)]))

    def energy(
        self,
        stars: Type["SSP"],
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> int:
        """
        Get the number of supernovae that have gone off between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: Energy from supernovae
        """
        return np.sum(self.energy_func(stars[self.explodability_func(stars, t0, t1)]))

    def count(
        self,
        stars: Type["SSP"],
        t0: Quantity["time"],
        t1: Quantity["time"],
    ) -> int:
        """
        Get the number of supernovae that have gone off between time t0 and t1.
        :param stars: Stellar Population
        :type stars: SSP
        :param t0: Start time
        :type t0: Quantity["time"]
        :param t1: End time
        :type t1: Quantity["time"]
        :return: Number of supernovae
        """
        return np.sum(self.explodability_func(stars, t0, t1))

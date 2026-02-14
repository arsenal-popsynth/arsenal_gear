"""
feedbacks
=========

Functions for the returns from various feedback mechanisms (SN, winds, radiation, etc.)
"""

from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity

__all__ = ["FBMechanism"]


class FBMechanism:
    """
    This class is inherited by a given feedback mechanism.
    """

    def __init__(self, start: Quantity["time"], end: Quantity["time"]) -> None:
        self.start = start
        self.end = end

    # What units should we put out here for each band? photons? energy?
    def radiation(self, bands: np.float64) -> np.float64:
        """
        Ionizing radiation output of the mechanism.
        """

    def ejecta_mass(self) -> Quantity["mass"]:
        """
        Mass ejected by the mechanism.
        """

    def ejecta_velocity(self) -> Quantity["velocity"]:
        """
        Velocity of the ejecta.
        """

    def metal_yields(self) -> Quantity["mass"]:
        """
        Metal fraction of the ejecta mass.
        """

    def species_yields(self, species: list[str]) -> Quantity["mass"]:
        """
        Individual species fractions of the ejecta mass.
        """

"""
feedbacks
=========

Functions for the returns from various feedback mechanisms (SN, winds, radiation, etc.)
"""


from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity


class FBMechanism():
    def __init__(self, start:Quantity["time"], end:Quantity["time"]) -> None:
        self.start = start
        self.end = end
        pass

    def radiation(self, bands:np.float64) -> np.float64: # What units should we put out here for each band? photons? energy?
        pass

    def ejecta_mass(self) -> Quantity["mass"]:
        pass

    def ejecta_velocity(self) -> Quantity["velocity"]:
        pass

    def metal_yields(self) -> Quantity["mass"]:
        pass

    def species_yields(self, species:list[str]) -> Quantity["mass"]:
        pass

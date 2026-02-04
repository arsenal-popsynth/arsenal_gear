"""
yieldtables
==========

This module contains class used consolidate the structure of different
yield tables available in the litterature.
"""

from pathlib import Path
from typing import List

from astropy.units import Quantity

from ..population import StarPopulation


def _not_implemented(channel: str, name: str):
    raise NotImplementedError(f"{channel} yields not implemented for {name}.")


class YieldTables:
    """Class for defining base structure of yield tables."""

    def __init__(self) -> None:
        self.name = "None"
        self.filedir = Path(__file__).parent.absolute()

        self.elements = None
        self.atomic_num = None

    def ccsn_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Placeholder function for core-collapse SNe yields."""
        _ = (elements, starPop, interpolate, extrapolate)
        return _not_implemented("CCSN", self.name)

    def snia_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Placeholder function for SNIa yields."""
        _ = (elements, starPop, interpolate, extrapolate)
        return _not_implemented("SNIa", self.name)

    def wind_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Placeholder function for wind yields (massive stars)."""
        _ = (elements, starPop, interpolate, extrapolate)
        return _not_implemented("Wind", self.name)

    def agb_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Placeholder function for AGB yields."""
        _ = (elements, starPop, interpolate, extrapolate)
        return _not_implemented("AGB", self.name)

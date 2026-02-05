"""
yields
==========

This module contains the abstract class used define a set
chemical yields that will interact with stellar feedback.
"""

from abc import ABC, abstractmethod
from typing import List

from astropy.units import Quantity

from ..population import SSP


class Yields(ABC):
    """Abstract base class for yield tables."""

    def __init__(self) -> None:

        self.ccsn_choice = None
        self.snia_choice = None
        self.wind_choice = None
        self.agb_choice = None

    @abstractmethod
    def ccsn_yields(
        self,
        elements: List[str],
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for core-collapse SNe yields.

        Args:
            elements: list of elements to get yields for
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Yields for the specified elements
        """

    @abstractmethod
    def snia_yields(
        self,
        elements: List[str],
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for SNe type Ia yields.

        Args:
            elements: list of elements to get yields for
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Yields for the specified elements
        """

    @abstractmethod
    def wind_yields(
        self,
        elements: List[str],
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for stellar wind yields.

        Args:
            elements: list of elements to get yields for
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Yields for the specified elements
        """

    @abstractmethod
    def agb_yields(
        self,
        elements: List[str],
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for AGB mass loss yields.

        Args:
            elements: list of elements to get yields for
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Yields for the specified elements
        """

    @abstractmethod
    def ccsn_mloss(
        self,
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for core-collapse SNe mloss.

        Args:
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Total mass loss from core-collapse SNe
        """

    @abstractmethod
    def snia_mloss(
        self,
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for SNe type Ia mloss.

        Args:
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Total mass loss from SNeIa
        """

    @abstractmethod
    def wind_mloss(
        self,
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for stellar wind mloss.

        Args:
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Total mass loss from stellar winds
        """

    @abstractmethod
    def agb_mloss(
        self,
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for AGB mass loss mloss.

        Args:
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Total mass loss from AGB
        """

    @abstractmethod
    def ccsn_relative(
        self,
        elements: List[str],
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for core-collapse SNe yields relative to
           total mass loss.

        Args:
            elements: list of elements to get yields for
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Fractional yields for the specified elements
        """

    @abstractmethod
    def snia_relative(
        self,
        elements: List[str],
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for SNe type Ia yields relative to
           total mass loss.


        Args:
            elements: list of elements to get yields for
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Fractional yields for the specified elements
        """

    @abstractmethod
    def wind_relative(
        self,
        elements: List[str],
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for stellar wind yields relative to
           total mass loss.


        Args:
            elements: list of elements to get yields for
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Fractional yields for the specified elements
        """

    @abstractmethod
    def agb_relative(
        self,
        elements: List[str],
        starPop: SSP,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Abstract method for AGB mass loss yields relative to
           total mass loss.


        Args:
            elements: list of elements to get yields for
            starPop: SSP object
            interpolate: interpolation method
            extrapolate: whether to extrapolate beyond bounds

        Returns:
            Fractional yields for the specified elements
        """

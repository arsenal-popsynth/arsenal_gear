"""
yieldtables
==========

This is basically a header file for the YieldTables class
it specifies the **maximal** implementation of a YieldTables
derived class, which will typically only implement a subset
of the available yield channels.
"""

from functools import wraps
from pathlib import Path
from typing import List

import numpy as np
from astropy import units as u

from ..formation import SinglePop


def _not_implemented(channel: str, name: str):
    raise NotImplementedError(f"{channel} yields not implemented for {name}.")


def requires_available(func):
    """Require the yield model to be available before executing a method."""

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if not self.available:
            raise RuntimeError(
                f"{self.__class__.__name__}: Yield model is unavailable."
            )

        return func(self, *args, **kwargs)

    return wrapper


class YieldTables:
    """Class for defining base structure of yield tables."""

    def __init__(self, interpolate: str = "nearest", extrapolate: bool = False) -> None:
        """
        Keyword Args:
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to limits if outside bound
        """
        self.interpolate = interpolate
        self.extrapolate = extrapolate

        self.name = "None"
        self.filedir = Path(__file__).parent.absolute()
        self.available = True

        self.elements = None
        self.atomic_num = None

        self.ccsn = None
        self.wind = None
        self.agb = None
        self.sagb = None
        self.snia = None

    def ccsn_args(self, starPop: SinglePop) -> List:
        """Placeholder function to define arguments going into ccsn yield,
        mass loss and relative. Also used to determine whether ccsn is
        available for the class."""
        _ = starPop
        return _not_implemented("CCSN", self.name)

    def wind_args(self, starPop: SinglePop) -> List:
        """Placeholder function to define arguments going into wind yield,
        mass loss and relative. Also used to determine whether wind is
        available for the class."""
        _ = starPop
        return _not_implemented("WIND", self.name)

    def agb_args(self, starPop: SinglePop) -> List:
        """Placeholder function to define arguments going into agb yield,
        mass loss and relative. Also used to determine whether agb is
        available for the class."""
        _ = starPop
        return _not_implemented("AGB", self.name)

    def sagb_args(self, starPop: SinglePop) -> List:
        """Placeholder function to define arguments going into sagb yield,
        mass loss and relative. Also used to determine whether sagb is
        available for the class."""
        _ = starPop
        return _not_implemented("SAGB", self.name)

    def snia_args(self, starPop: SinglePop) -> List:
        """Placeholder function to define arguments going into snia yield,
        mass loss and relative. Also used to determine whether snia is
        available for the class."""
        _ = starPop
        return _not_implemented("SNIA", self.name)

    @requires_available
    def ccsn_yields(
        self, elements: List[str], starPop: SinglePop
    ) -> u.Quantity["mass"]:
        """Interpolate yields from core-collapse supernovae for specified elements.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: SinglePop object
        Returns:
            Dictionary of yields matching provided element list

        """
        elements = np.atleast_1d(elements)
        args = self.ccsn_args(starPop)
        yld_array = self.ccsn.get_yld(
            elements, args, interpolate=self.interpolate, extrapolate=self.extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    @requires_available
    def wind_yields(
        self, elements: List[str], starPop: SinglePop
    ) -> u.Quantity["mass"]:
        """Interpolate yields from massive stars ejected as winds for specified elements.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: SinglePop object
        Returns:
            Dictionary of yields matching provided element list

        """
        elements = np.atleast_1d(elements)
        args = self.wind_args(starPop)
        yld_array = self.wind.get_yld(
            elements, args, interpolate=self.interpolate, extrapolate=self.extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    @requires_available
    def agb_yields(
        self,
        elements: List[str],
        starPop: SinglePop,
    ) -> u.Quantity["mass"]:
        """Interpolate yields from AGB stars ejected as winds for specified elements

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: SinglePop object
        Returns:
            Dictionary of yields matching provided element list

        """
        elements = np.atleast_1d(elements)
        args = self.agb_args(starPop)
        yld_array = self.agb.get_yld(
            elements, args, interpolate=self.interpolate, extrapolate=self.extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    @requires_available
    def sagb_yields(
        self,
        elements: List[str],
        starPop: SinglePop,
    ) -> u.Quantity["mass"]:
        """Interpolate yields from stars on the super asymptotic giant branch for specified elements.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: stellar population object
        Returns:
            Dictionary of yields matching provided element list

        """
        elements = np.atleast_1d(elements)
        args = self.sagb_args(starPop)
        yld_array = self.sagb.get_yld(
            elements, args, interpolate=self.interpolate, extrapolate=self.extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    @requires_available
    def snia_yields(
        self,
        elements: List[str],
        starPop: SinglePop,
    ) -> u.Quantity["mass"]:
        """Interpolate yields supernovae type Ia for specified elements.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: stellar population object
        Returns:
            Dictionary of yields matching provided element list

        """
        elements = np.atleast_1d(elements)
        args = self.snia_args(starPop)
        yld_array = self.snia.get_yld(
            elements, args, interpolate=self.interpolate, extrapolate=self.extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

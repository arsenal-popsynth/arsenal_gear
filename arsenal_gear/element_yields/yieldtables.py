"""
yieldtables
===========

Contains the YieldTables class which specifies the the **maximal**
implementation of a yield table derived class, which will typically
only implement a subset of the available yield channels.
"""

from functools import wraps
from pathlib import Path
from typing import List

import numpy as np
from astropy import units as u

from ..formation import SinglePop


def _not_implemented(channel: str, name: str):
    """Internal function to raise a NotImplementedError for missing yields."""
    raise NotImplementedError(f"{channel} yields not implemented for {name}.")


def requires_available(func):
    """Wrapper function to force available yield model before executing a method."""

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if not self.available:
            raise RuntimeError(
                f"{self.__class__.__name__}: Yield model is unavailable."
            )

        return func(self, *args, **kwargs)

    return wrapper


class YieldTables:
    """
    Class inherited when creating a new set of yield tables. The class implements
    the base structure for a set of yield tables available from some literature
    source and defines the methods to access the data. When inherited, functions
    *_args should be implemented for available nucleosynthetic sources, e.g.,
    ccsn, wind, and agb. If not implemented the yield class assumes that they are
    not part of the yield set and raises the appropriate error.

    Attributes:
        interpolate: interpolation method passed to scipy.interpolate.RegularGridInterpolator
        extrapolate: if False, then mass, metal, and rot are set to limits if outside bounds
        name: name of the yield table set
        filedir: path to tabulated data
        available: marks if yield set is available (used to mark errors during initialization)
        elements: list of chemical elements available
        atomic_num: list of atomic numbers corresponding to elements
        ccsn: Source class for core-collapse SNe
        wind: Source class for pre-SNe winds
        agb: Source class for asymptotic giant branch stars
        sagb: Source class for super-asymptotic giant branch stars
        snia: Source class for Type Ia SNe
    """

    def __init__(self, interpolate: str = "nearest", extrapolate: bool = False) -> None:
        """
        Initializes default values or placeholders for class attributes.

        Arguments:
            interpolate: interpolation method passed to scipy.interpolate.RegularGridInterpolator
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
        """
        Placeholder function to define arguments going into ccsn yield, mass loss and relative.
        Raises NotImplementedError if ccsn is not part of the yield set.

        Argument:
            starPop: SinglePop class used to access stellar data for interpolation
        """
        _ = starPop
        return _not_implemented("CCSN", self.name)

    def wind_args(self, starPop: SinglePop) -> List:
        """
        Placeholder function to define arguments going into wind yield, mass loss and relative.
        Raises NotImplementedError if wind is not part of the yield set.

        Argument:
            starPop: SinglePop class used to access stellar data for interpolation
        """
        _ = starPop
        return _not_implemented("WIND", self.name)

    def agb_args(self, starPop: SinglePop) -> List:
        """
        Placeholder function to define arguments going into agb yield, mass loss and relative.
        Raises NotImplementedError if agb is not part of the yield set.

        Argument:
            starPop: SinglePop class used to access stellar data for interpolation
        """
        _ = starPop
        return _not_implemented("AGB", self.name)

    def sagb_args(self, starPop: SinglePop) -> List:
        """
        Placeholder function to define arguments going into sagb yield, mass loss and relative.
        Raises NotImplementedError if sagb is not part of the yield set.

        Argument:
            starPop: SinglePop class used to access stellar data for interpolation
        """
        _ = starPop
        return _not_implemented("SAGB", self.name)

    def snia_args(self, starPop: SinglePop) -> List:
        """
        Placeholder function to define arguments going into snia yield, mass loss and relative.
        Raises NotImplementedError if snia is not part of the yield set.

        Argument:
            starPop: SinglePop class used to access stellar data for interpolation
        """
        _ = starPop
        return _not_implemented("SNIA", self.name)

    @requires_available
    def ccsn_yields(
        self, elements: List[str], starPop: SinglePop
    ) -> u.Quantity["mass"]:
        """
        Interpolate yields from core-collapse supernovae for specified elements.

        Arguments:
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
        """
        Interpolate yields from massive stars ejected as winds for specified elements.

        Arguments:
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
        """
        Interpolate yields from stars on the asymptotic giant branch for specified elements.

        Arguments:
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
        """
        Interpolate yields from stars on the super asymptotic giant branch for specified elements.

        Arguments:
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
        """Interpolate yields from supernovae Type Ia SNe for specified elements.

        Arguments:
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

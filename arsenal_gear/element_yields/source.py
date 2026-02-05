"""
source
==========

This module contains the class used for defining a source of yields
and interpolate tabulated data.
"""

import warnings

import numpy as np
from scipy.interpolate import RegularGridInterpolator


class Source:
    """
    Class used to store interpolators of yields from a given
    source (e.g., winds or core-collapse SNe). Makes use of
    scipy.interpolate.RegularGridInterpolator for flexibility
    and efficiency.

    Includes functions:
    get_yield - returns yields of given elements
    get_mloss - return sum of all yields in table.
    """

    def __init__(self, elements, params, yields) -> None:

        self.params = params
        self.yields = {
            element: RegularGridInterpolator(
                self.params, yields[i], fill_value=None, bounds_error=False
            )
            for i, element in enumerate(elements)
        }
        self.mloss = RegularGridInterpolator(
            self.params, np.sum(yields, axis=0), fill_value=None, bounds_error=False
        )

    def get_yld(self, elements, params, interpolate="nearest", extrapolate=False):
        """Interpolate yields from class.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            params: list of parameters of the table (e.g., mass, metallicity, rotation)
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then params are set to limits if outside bound.
        Returns:
            List of yields matching provided element list

        """
        elements = np.atleast_1d(elements)

        if len(params) != len(self.params):
            raise ValueError("Supplied parameters do not match yield set parameters.")

        points = self._convert2array(params)

        if extrapolate:
            warnings.warn(
                """Extrapolating yields might lead to problematic behaviour (e.g., negative yields).
                   Ensure that yields behave as expected."""
            )
        else:
            for ip in range(len(params)):
                points[:, ip][(points[:, ip] < np.min(self.params[ip]))] = np.min(
                    self.params[ip]
                )
                points[:, ip][(points[:, ip] > np.max(self.params[ip]))] = np.max(
                    self.params[ip]
                )

        if len(elements) == 1:
            try:
                return self.yields[elements[0]](points, method=interpolate)
            except KeyError:
                warnings.warn(f"Element {elements[0]} is not part of yield set.")
                return np.nan
        else:
            if all(element in list(self.yields.keys()) for element in elements):
                return np.array(
                    [
                        self.yields[element](points, method=interpolate)
                        for element in elements
                    ]
                )

            yld = []
            shape = self.yields["H"](points, method=interpolate).shape
            for element in elements:
                try:
                    yld.append(self.yields[element](points, method=interpolate))
                except KeyError:
                    warnings.warn(f"Element {element} is not part of yield set.")
                    yld.append(np.ones(shape) * np.nan)
            return yld

    def get_mloss(self, params, interpolate="nearest", extrapolate=False):
        """Interpolate sum of all yields (total mass) from class.

        Args:
            params: list of parameters of the table (e.g., mass, metallicity, rotation)
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then params are set to limits if outside bound.
        Returns:
            Total mass ejected by source (i.e., sum of all elements)

        """
        if len(params) != len(self.params):
            raise ValueError("Supplied parameters do not match yield set parameters.")

        points = self._convert2array(params)

        if extrapolate:
            warnings.warn(
                """Extrapolating yields might lead to problematic behaviour (e.g., negative yields).
                   Ensure that yields behave as expected."""
            )
        else:
            for ip in range(len(params)):
                points[:, ip][(points[:, ip] < np.min(self.params[ip]))] = np.min(
                    self.params[ip]
                )
                points[:, ip][(points[:, ip] > np.max(self.params[ip]))] = np.max(
                    self.params[ip]
                )

        return self.mloss(points, method=interpolate)

    @staticmethod
    def _convert2array(params):
        """Internal function to convert list of parameters to a format
        that RegularGridInterpolator can handle.

        Args:
            params: list of parameters of the table (e.g., mass, metallicity, rotation)
        Returns:
            Interpolation points for RegularGridInterpolator

        """
        plens = np.array([len(param) if isinstance(param, (list, np.ndarray))
                         else 1 for param in params])
        max_length = max(plens)

        # Ensure all arguments are lists or arrays of the same length
        if np.any((plens != 1) & (plens != max_length)):
            raise ValueError(
                "All lists of parameters must have the same length."
            )

        args = []
        for param in params:
            if isinstance(param, (list, np.ndarray)):
                args.append(np.array(param))
            else:
                # Fill with the same value
                args.append(np.full(max_length, param))

        # Combine arguments into points for interpolation
        return np.stack(args, axis=-1)

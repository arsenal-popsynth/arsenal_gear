import astropy.units as u
from astropy.units import Quantity
import numpy as np
from scipy.interpolate import RegularGridInterpolator

import os
import warnings
import re
from typing import Type
from typing import List

class Source:

    def __init__(self, elements, params, yields) -> None:

        self.params = params
        self.yields = {element: RegularGridInterpolator(self.params, yields[i],
                                                        fill_value=None, bounds_error=False)
                       for i, element in enumerate(elements)}
        self.mloss = RegularGridInterpolator(self.params, np.sum(yields, axis=0),
                                             fill_value=None, bounds_error=False)

    def get_yld(self, elements, params, interpolate='nearest', extrapolate=False):

        elements = np.atleast_1d(elements)

        if len(params) != len(self.params):
            raise ValueError(
                "Supplied parameters do not match yield set parameters.")

        points = self.convert2array(params)

        if extrapolate:
            warnings.warn("Extrapolating yields might lead to problematic behaviour (e.g., negative yields). Ensure that yields behave as expected.")
        else:
            for ip in range(len(params)):
                points[:, ip][(points[:, ip] < np.min(self.params[ip]))] = np.min(
                    self.params[ip])
                points[:, ip][(points[:, ip] > np.max(self.params[ip]))] = np.max(
                    self.params[ip])
        
        if len(elements) == 1:
            try:
                return self.yields[elements[0]](points, method=interpolate)
            except KeyError:
                warnings.warn("Element {elements[0]} is not part of yield set.")
                return np.nan
        else:
            if all(element in list(self.yields.keys()) for element in elements):
                return np.array([self.yields[element](points, method=interpolate) for element in elements])
            else:
                yld = []
                shape = self.yields['H'](points, method=interpolate).shape
                for element in elements:
                    try:
                        yld.append(self.yields[element](points, method=interpolate))
                    except KeyError:
                        warnings.warn("Element {elements[0]} is not part of yield set.")
                        yld.append(np.ones(shape) * np.nan)
                return yld

    def get_mloss(self, params, interpolate='nearest', extrapolate=False):

        if len(params) != len(self.params):
            raise ValueError(
                "Supplied parameters do not match yield set parameters.")

        points = self.convert2array(params)

        if extrapolate:
            warnings.warn("Extrapolating yields might lead to problematic behaviour (e.g., negative yields). Ensure that yields behave as expected.")
        else:
            for ip in range(len(params)):
                points[:, ip][(points[:, ip] < np.min(self.params[ip]))] = np.min(
                    self.params[ip])
                points[:, ip][(points[:, ip] > np.max(self.params[ip]))] = np.max(
                    self.params[ip])

        return self.mloss(points, method=interpolate)

    @staticmethod
    def convert2array(params):
        max_length = max([len(param) if isinstance(
            param, (list, np.ndarray)) else 1 for param in params])

        # Ensure all arguments are lists or arrays of the same length
        args = []
        for param in params:
            if isinstance(param, (list, np.ndarray)):
                if len(param) != max_length:
                    raise ValueError(
                        "All list of parameters must have the same length.")
                args.append(np.array(param))
            else:
                # Fill with the same value
                args.append(np.full(max_length, param))

        # Combine arguments into points for interpolation
        return np.stack(args, axis=-1)


class Yields:
    """ Header class for yield tables."""

    def __init__(self) -> None:

        warnings.warn(
            "This is a header class for inheritance and should not be used by itself")

        self.filedir = os.path.dirname(os.path.realpath(__file__))
        self.yield_tablefile = self.filedir + '<yield file name>'

        self.elements, self.atomic_num = None

    def ccsn_yields(self, 
                    elements: List[str], 
                    mass: Quantity["mass"], 
                    metal: Quantity["dimensionless"],
                    interpolate: str="nearest") -> Quantity["mass"]:
        

        raise ValueError("Core-collapse SNe is not part of this yield set.")

    def snia_yields(self,                     
                    elements: List[str], 
                    mass: Quantity["mass"], 
                    metal: Quantity["dimensionless"],
                    interpolate: str="nearest") -> Quantity["mass"]:
        raise ValueError("Type Ia SNe is not part of this yield set.")

    def wind_yields(self,
                    elements: List[str], 
                    mass: Quantity["mass"], 
                    metal: Quantity["dimensionless"],
                    interpolate: str="nearest") -> Quantity["mass"]:
        raise ValueError(
            "Stellar winds (main sequence) is not part of this yield set.")

    def agb_yields(self,
                    elements: List[str], 
                    mass: Quantity["mass"], 
                    metal: Quantity["dimensionless"],
                    interpolate: str="nearest") -> Quantity["mass"]:
        raise ValueError(
            "Stellar wind (asymptotic giant branch) is not part of this yield set.")

"""
formation
=================

Submodule to describe the properties of a stellar population
when it is first formed.
"""

import astropy.units as u
import numpy as np

# local imports
from . import dist_funcs
from .form_data_structures import SinglePop, BinaryPop

__all__ = ["Formation", "SinglePop", "BinaryPop"]


class Formation():
    """
    Abstract base class for describing the properties of a stellar population
    at the time of formation.
    """
    def __init__(self, **kwargs) -> None:
        # get population parameters from kwargs
        subpop_dicts = [kwargs.get(k) for k in kwargs if k.startswith("pop")]
        self.nsubpops = len(subpop_dicts)
        if (self.nsubpops == 0):
            default_pop = {
                "Mtot": 1e6 * u.Msun,
                "metallicity": 0.0 * u.dimensionless_unscaled,
                "imf": "salpeter",
                "mmin": 0.08,
                "mmax": 100.0,
                "discrete": True,
                "seed": kwargs.get("seed", None)}
            subpop_dicts = [default_pop]
            self.nsubpops = 1
        for i, pop_dict in enumerate(subpop_dicts):
            if not isinstance(pop_dict, dict):
                err_msg = f"Subpopulation {i} is not a dictionary.\
                           Please provide subpopulation parameters as dictionaries."
                raise ValueError(err_msg)
        
        subpops = []
        for pop_dict in subpop_dicts:
            pop_type = pop_dict.get("type", "single")
            if pop_type == "single":
                subpops.append(self._init_single_pop(pop_dict))
            elif pop_type == "binary":
                subpops.append(self._init_binary_pop(pop_dict))
            else:
                err_msg = f"Invalid population type: {pop_type}.\
                           Must be 'single' or 'binary'."
                raise ValueError(err_msg)
        self.subpops = subpops
        # Nstar can be a float if one of the populations is not discrete
        self.Nstar = np.sum([pop.Nstar for pop in self.subpops])
        self.Mtot = np.sum([pop.Mtot.value for pop in self.subpops]) * u.Msun
        self.metallicities = [pop.metallicity for pop in self.subpops]


    @staticmethod
    def _init_single_pop(pop_dict: dict) -> SinglePop:
        """
        Initialize a single population of stars from a dictionary of parameters.
        """
        Mtot = pop_dict.get("Mtot", 1e6 * u.Msun)
        metallicity = pop_dict.get("metallicity", 0.0) * u.dimensionless_unscaled
        imf = pop_dict.get("imf", "salpeter")
        if isinstance(imf, dist_funcs.imf.IMF):
            mmin = imf.min_mass * u.Msun
            mmax = imf.max_mass * u.Msun
        elif isinstance(imf, str):
            mmin = pop_dict.get("mmin", 0.08)
            mmax = pop_dict.get("mmax", 100.)
            if not isinstance(mmin, float) or not isinstance(mmax, float):
                raise ValueError("mmin and mmax must be floats specifying mass in Msun.")
            mmin *= u.Msun
            mmax *= u.Msun
            if imf == "salpeter":
                if "seed" in pop_dict and pop_dict["seed"] is not None:
                    imf = dist_funcs.imf.Salpeter(mmin,mmax, seed=pop_dict["seed"])
                else:
                    imf = dist_funcs.imf.Salpeter(mmin,mmax)
            else:
                err_msg = f"Invalid IMF string specified, {imf}. \
                            Currently only 'salpeter' is supported."
                raise ValueError(err_msg)
        else:
            err_msg = f"Invalid IMF specified, {imf}. \
                        Must be an IMF instance or string specifying IMF."
            raise ValueError(err_msg)

        Nstar = (Mtot / imf.mean()).value

        discrete = pop_dict.get("discrete", True)
        if not isinstance(discrete, bool):
            err_msg = f"Invalid discrete parameter: {discrete}. Must be a boolean."
            raise ValueError(err_msg)
        if discrete:
            masses = imf.sample_mass(Mtot)
            Nstar = len(masses)
        else:
            masses = None

        return SinglePop(Mtot,
                         Nstar,
                         metallicity,
                         imf,
                         mmin,
                         mmax,
                         discrete,
                         masses)
    
    @staticmethod
    def _init_binary_pop(pop_dict: dict) -> BinaryPop:
        """
        Initialize a binary population of stars from a dictionary of parameters.
        """
        raise NotImplementedError("Binary populations not yet implemented.")